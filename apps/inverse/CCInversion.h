#include <Potential.h>          // base class of the CC-inversion potentials

//------------------Reconstruction via the collective coordinates
/** \brief Collective-coordinate inversion potential: penalizes deviation of the structure factor from target values at a set of k-points. (CLAUDE) */
class FrankInvertCCPotential : public Potential
{
protected:
	const Configuration * pConfig;
public:
	/** \brief A single k-point constraint: target real/imaginary structure-factor amplitude with a rescaling weight. (CLAUDE) */
	struct KPointConstraint
	{
		GeometryVector k;//absolute k, not relative to basis vectors
		double V;//rescaling constant
		double targetRhoReal, targetRhoImag;
		/** \brief Construct a k-point constraint. (CLAUDE)
		 *  \param k              Absolute wavevector of the constraint.
		 *  \param V              Rescaling (weight) constant for this constraint.
		 *  \param targetRhoReal  Target real part of the collective coordinate at k.
		 *  \param targetRhoImag  Target imaginary part of the collective coordinate at k. */
		KPointConstraint(const GeometryVector & k, double V, double targetRhoReal, double targetRhoImag) : k(k), V(V), targetRhoReal(targetRhoReal), targetRhoImag(targetRhoImag)
		{
		}
	};
	//users or derived classes should fill this vector with k-constraints
	std::vector<KPointConstraint> constraints;

	//this value is added to the energy. It is initialized as 0, but user can change it.
	double Shift;

	//this class will manage these members, users can also read it.
	std::vector<double> RhoReal, RhoImag;
	bool RhoValid;
	/** \brief Compute and cache the real/imaginary parts of the collective coordinates at every constraint k-point. (CLAUDE) */
	virtual void GetRho(void)
	{
		assert(this->pConfig != nullptr);
		if (this->RhoValid)
		{
			return;
		}

		signed long end = this->constraints.size();
		this->RhoReal.resize(end);
		this->RhoImag.resize(end);

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for (signed long i = 0; i<end; i++)
		{
			double tr = 0, ti = 0;//real and imag parts of result
			for (size_t j = 0; j<pConfig->NumParticle(); j++)
			{
				GeometryVector c = pConfig->GetCartesianCoordinates(j);
				tr += std::cos(this->constraints[i].k.Dot(c));
				ti += std::sin(this->constraints[i].k.Dot(c));
			}
			this->RhoReal[i] = tr;
			this->RhoImag[i] = ti;
		}

		this->RhoValid = true;
	}

	/** \brief Construct a CC-inversion potential of the given dimension with an empty constraint set. (CLAUDE)
	 *  \param dimension  Spatial dimension. */
	FrankInvertCCPotential(DimensionType dimension) :Potential(dimension)
	{
		this->pConfig = nullptr;
		this->RhoValid = false;
		this->ParallelNumThread = 1;
		this->Shift = 0.0;
	}
	/** \brief Destructor. (CLAUDE) */
	virtual ~FrankInvertCCPotential()
	{
	}
	/** \brief Return a heap-allocated copy of this potential. (CLAUDE)
	 *  \return  A new FrankInvertCCPotential that copies this one. */
	FrankInvertCCPotential * clone()
	{
		return new FrankInvertCCPotential(*this);
	}
	/** \brief Return the total potential energy of the given configuration. */
	virtual double Energy()//return the total potential energy of the given configuration
	{
		if (this->RhoValid == false)
			this->GetRho();

		double result = 0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end = this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread = omp_get_thread_num();
#pragma omp for schedule(static)
			for (signed long i = 0; i < end; i++)
			{
				double temp1 = this->RhoReal[i] - this->constraints[i].targetRhoReal;
				double temp2 = this->RhoImag[i] - this->constraints[i].targetRhoImag;
				threadresult[numthread] += this->constraints[i].V*(temp1*temp1 + temp2*temp2);
			}
		}
		for (auto iter = threadresult.begin(); iter != threadresult.end(); iter++)
			result += *iter;

		return result + this->Shift;
	}
	/** \brief Calculate the force for the ith Particle.
	 *  \param[out] result  Force on particle i.
	 *  \param i            Particle index. */
	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		if (this->RhoValid == false)
			this->GetRho();

		result = GeometryVector(this->Dimension);

		signed long end = this->constraints.size();

		for (signed long j = 0; j<end; j++)
		{
			GeometryVector k = this->constraints[j].k;
			GeometryVector r = this->pConfig->GetCartesianCoordinates(i);
			double kr = k.Dot(r);
			double skr = std::sin(kr), ckr = std::cos(kr);
			double dRhoReal = RhoReal[j] - this->constraints[j].targetRhoReal;
			double dRhoImag = RhoImag[j] - this->constraints[j].targetRhoImag;
			result.AddFrom(k*(-2.0*(ckr*dRhoImag - skr*dRhoReal)));
		}
	}
	/** \brief Compute the forces on all particles. (CLAUDE)
	 *  \param[out] results  Resized to hold the force on every particle. */
	virtual void AllForce(std::vector<GeometryVector> & results)
	{
		if (this->RhoValid == false)
			this->GetRho();
		long end = this->pConfig->NumParticle();
		results.resize(end);
		unsigned long tempNumThread = 1;
		std::swap(tempNumThread, this->ParallelNumThread);
#pragma omp parallel for num_threads(tempNumThread) default(shared)
		for (long i = 0; i<end; i++)
			this->Force(results[i], i);
		std::swap(tempNumThread, this->ParallelNumThread);
	}
	/** \brief Set the configuration to evaluate and invalidate the cached structure factors. (CLAUDE)
	 *  \param Position  Configuration whose collective coordinates are constrained; its dimension must match. */
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->Dimension == Position.GetDimension());
		this->pConfig = &Position;
		this->RhoValid = false;
	}
};


//------------------Reconstruction via the collective coordinates 2 <- consider real and imaginary parts of collective coordinates separately...
/** \brief CC-inversion potential that constrains real and imaginary parts of the collective coordinates separately (handles an odd constraint count). (CLAUDE) */
class FrankInvertCCPotential_re : public FrankInvertCCPotential
{
protected:
	size_t numConstraints;
	size_t numkPoints;
	bool odd=true;
public:
	/** \brief Construct from a total number of (real+imaginary) constraints, deriving the number of k-points used. (CLAUDE)
	 *  \param d               Spatial dimension.
	 *  \param numConstraints  Total number of scalar constraints (real and imaginary parts counted separately). */
	FrankInvertCCPotential_re(DimensionType d, size_t numConstraints) : FrankInvertCCPotential(d)
	{
		this->numConstraints = numConstraints;
		this->numkPoints = (size_t)round(numConstraints/2.0 +0.1);
		if (numConstraints % 2 == 0)
			odd = false;
	}

	/** \brief Return the total potential energy of the given configuration. */
	virtual double Energy()//return the total potential energy of the given configuration
	{
		if (this->RhoValid == false)
			this->GetRho();

		double result = 0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end = this->numkPoints;
		//----for debug
		int count = 0;
#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread = omp_get_thread_num();
#pragma omp for schedule(static)
			for (signed long i = 0; i < end; i++)
			{
				double temp1 = this->RhoReal[i] - this->constraints[i].targetRhoReal;
				double temp2 = this->RhoImag[i] - this->constraints[i].targetRhoImag;
				count += 2;
				if (i == end - 1 && odd) { //if there are odd number of constraints.
					temp2 = 0;
					count -= 1;
				}
				threadresult[numthread] += this->constraints[i].V*(temp1*temp1 + temp2*temp2);
				

			}
		}
		for (auto iter = threadresult.begin(); iter != threadresult.end(); iter++)
			result += *iter;

		
		
		


		return result + this->Shift;
		std::cout << "num constraints=" << count <<std::endl;
	}
	
	/** \brief Calculate the force for the ith Particle.
	 *  \param[out] result  Force on particle i.
	 *  \param i            Particle index. */
	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		if (this->RhoValid == false)
			this->GetRho();

		result = GeometryVector(this->Dimension);

		signed long end = this->numkPoints;//this->constraints.size();

		for (signed long j = 0; j<end; j++)
		{
			GeometryVector k = this->constraints[j].k;
			GeometryVector r = this->pConfig->GetCartesianCoordinates(i);
			double kr = k.Dot(r);
			double skr = std::sin(kr), ckr = std::cos(kr);
			double dRhoReal = RhoReal[j] - this->constraints[j].targetRhoReal;
			double dRhoImag = RhoImag[j] - this->constraints[j].targetRhoImag;
			if (j == end - 1 && odd)//if there are odd number of constraints...
				dRhoImag = 0;
			result.AddFrom(k*(-2.0*(ckr*dRhoImag - skr*dRhoReal)));
		}
	}
};


//------------------Reconstruction via the collective coordinates 2 
//	consider real and imaginary parts of collective coordinates separately & possible to choose your constraints...
/** \brief CC-inversion potential with a selectable constraint type (real/imaginary/both) choosing which parts of the collective coordinates are constrained. (CLAUDE) */
class FrankInvertCCPotential2 : public FrankInvertCCPotential
{
protected:
	size_t numConstraints;
	size_t numkPoints;
	bool separateCount = true; //true: real and imaginary parts of a collective coordinate are considered seprately.
	bool realPartFirst = true; //true: if the number of constraints is odd, then consider real parts first.
	bool isNumConstraintsOdd = true;
public:
	/** \brief Construct from a dimension and number of constraints (the k-point count is finalized by SetConstraintType). (CLAUDE)
	 *  \param d               Spatial dimension.
	 *  \param numConstraints  Number of constraints requested. */
	FrankInvertCCPotential2(DimensionType d, size_t numConstraints) : FrankInvertCCPotential(d)
	{
		this->numConstraints = numConstraints;
	}

	/** \brief Select the constraint type and finalize the k-point count (separate real/imaginary counting only for d=1). (CLAUDE)
	 *  \param type  "real" (real parts first), "image" (imaginary parts first), or any other value for "both". */
	void SetConstraintType(const char * type) {
		if (strcmp(type, "real") == 0) {
			this->separateCount = true;
			this->realPartFirst = true;
		}
		else if (strcmp(type, "image") == 0) {
			this->separateCount = true;
			this->realPartFirst = false;
		}
		else {
			this->separateCount = false;
			this->realPartFirst = false;
		}
	
		//detailed set up...
		if (this->Dimension > 1 && this->separateCount ) {
			std::cout << "Seperate counting is not allowed when d>1 \n";
			this->separateCount = false;
		}
		else if (this->Dimension == 1 && this->separateCount) {
			this->isNumConstraintsOdd = (this->numConstraints % 2 == 1);
			this->numkPoints = (size_t)round((this->numConstraints + 0.8) / 2.0);
		}
		else {
			this->separateCount = false;
			this->numkPoints = numConstraints;
		}
	
	}

	/** \brief Return the number of k-points used by the current constraint type. (CLAUDE)
	 *  \return  Number of k-points. */
	size_t GetnumKPoints() {
		return this->numkPoints;
	}

	/** \brief Return the total potential energy of the given configuration. */
	virtual double Energy()//return the total potential energy of the given configuration
	{
		if (this->RhoValid == false)
			this->GetRho();

		double result = 0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end = this->numkPoints;
		//----for debug
		int count = 0;
#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread = omp_get_thread_num();
#pragma omp for schedule(static)
			for (signed long i = 0; i < end; i++)
			{
				double temp1 = this->RhoReal[i] - this->constraints[i].targetRhoReal;
				double temp2 = this->RhoImag[i] - this->constraints[i].targetRhoImag;
				count += 2;
				if (i == end - 1 && this->separateCount &&this->isNumConstraintsOdd ) { //if there are odd number of constraints => consider either image or real part first.
					if (this->realPartFirst)
						temp2 = 0;
					else
						temp1 = 0;
					count -= 1;
				}
				threadresult[numthread] += this->constraints[i].V*(temp1*temp1 + temp2*temp2);


			}
		}
		for (auto iter = threadresult.begin(); iter != threadresult.end(); iter++)
			result += *iter;


		return result + this->Shift;
		std::cout << "num constraints=" << count << std::endl;
	}

	/** \brief Calculate the force for the ith Particle.
	 *  \param[out] result  Force on particle i.
	 *  \param i            Particle index. */
	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		if (this->RhoValid == false)
			this->GetRho();

		result = GeometryVector(this->Dimension);

		signed long end = this->numkPoints;//this->constraints.size();

		for (signed long j = 0; j<end; j++)
		{
			GeometryVector k = this->constraints[j].k;
			GeometryVector r = this->pConfig->GetCartesianCoordinates(i);
			double kr = k.Dot(r);
			double skr = std::sin(kr), ckr = std::cos(kr);
			double dRhoReal = RhoReal[j] - this->constraints[j].targetRhoReal;
			double dRhoImag = RhoImag[j] - this->constraints[j].targetRhoImag;
			if (j == end - 1 && this->separateCount && this->isNumConstraintsOdd) {//if there are odd number of constraints...
				if (this->realPartFirst)
					dRhoImag = 0;
				else
					dRhoReal = 0;
			}
			result.AddFrom(k*(-2.0*(ckr*dRhoImag - skr*dRhoReal)));
		}
	}
};
