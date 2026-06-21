#include <Potential.h>          // base class of the CC-inversion potentials
#include "LatticeUtil.h"        // GetLattice_config (was in etc1.h, now retired)
//----------------------------------
/** \brief Harmonic potential of a lattice: nearest-neighbor springs acting on particle displacements from a reference. (CLAUDE) */
class SimpleHarmonic : public Potential //potential energy of the simple harmonic lattice
{
private:
	double K; //spring constant
	double a; //lattice constant
	Configuration ePositions; //original positions: presumably a lattice
	std::vector<std::vector<size_t > > nnList; // list of prt indices of the nearest neighbors.
											   //std::vector<GeometryVector> du; // list of displacements
	std::vector<GeometryVector> u; //list of displacements u

public:
	//const double a = 1; //lattice constant
	Configuration confi;

	/** \brief Construct a harmonic potential with the given dimension, spring constant, and lattice constant. (CLAUDE)
	 *  \param d                Spatial dimension.
	 *  \param springConstant   Spring constant K of the nearest-neighbor bonds.
	 *  \param latticeConstant  Lattice constant a (sets the nearest-neighbor search radius). */
	SimpleHarmonic(DimensionType d, double springConstant, double latticeConstant) : Potential(d), K(springConstant), a(latticeConstant)
	{}

	/** \brief Set the reference (equilibrium) configuration and build the nearest-neighbor list and zero displacements. (CLAUDE)
	 *  \param ref  Reference lattice configuration; must match this potential's dimension. */
	void SetReference(const Configuration & ref) {

		std::vector<std::vector<size_t > >& temp = nnList;

		if (this->Dimension == ref.GetDimension()) {
			this->ePositions = Configuration(ref);
			temp.reserve(this->ePositions.NumParticle());

			//make a list of nearest neighbors
			double r = a*1.1;
			for (int i = 0; i < this->ePositions.NumParticle(); i++) {
				std::vector<size_t> nn;
				this->ePositions.IterateThroughNeighbors(ePositions.GetRelativeCoordinates(i), r,
					[&nn, i](const GeometryVector &x, const GeometryVector &LatticeShift, const signed long * PeriodicShift, const size_t SourceParticle)->void
				{
					if (i != SourceParticle) {
						nn.push_back(SourceParticle);
					}
				}
				);
				temp.push_back(nn);
			}

			GeometryVector initializer(static_cast<int> (this->ePositions.GetDimension()));
			this->u.reserve(this->ePositions.NumParticle());
			for (int i = 0; i < this->ePositions.NumParticle(); i++) {
				this->u.push_back(initializer);
			}
		}
		else {
			std::cout << "SimpleHarmonic::SetReference::wrong dimensions\n";
		}
	}

	/** \brief Set the current configuration (must match the reference particle count and dimension). (CLAUDE)
	 *  \param Position  Configuration to evaluate the potential on. */
	virtual void SetConfiguration(const Configuration & Position) {
		if (this->ePositions.NumParticle() == Position.NumParticle() && this->Dimension == Position.GetDimension())
			this->confi = Configuration(Position);
		else
			std::cout << "SimpleHarmonic::setConfiguration::sth wrong\n";
	}

	/** \brief Update configuration based on the displacements. */
	void updateConfiguration() {//update configuration based on the displacements
		GeometryVector temp_position(static_cast<int>(this->confi.GetDimension()));
		for (int i = 0; i < this->confi.NumParticle(); i++) {
			temp_position = this->ePositions.GetCartesianCoordinates(i) + this->u[i];
			this->confi.MoveParticle(i, this->confi.CartesianCoord2RelativeCoord(temp_position));
		}
	}

	/** \brief Update un based from the configuration. */
	void updateUn() {
		GeometryVector U(static_cast<int>(this->confi.GetDimension()));
		for (int i = 0; i < this->confi.NumParticle(); i++) {
			U = this->confi.GetRelativeCoordinates(i) - this->ePositions.GetRelativeCoordinates(i);
			this->confi.RelativeCoordToMinimumImage(U);
			u[i] = this->confi.RelativeCoord2CartesianCoord(U);
		}
	}

	/** \brief Update a single displacement vector by the change du.
	 *  \param i   Index of the particle whose displacement is updated.
	 *  \param du  Increment added to the current displacement of particle i. */
	void updateUn(size_t i, GeometryVector du) {
		this->u[i] = this->u[i] + du;
	}

	/** \brief Return the current displacement vector of particle i. (CLAUDE)
	 *  \param i  Particle index.
	 *  \return   The displacement u[i] of particle i. */
	GeometryVector Get_u(size_t i) {
		return this->u[i];
	}
	/** \brief Compute the energy of ith particle.
	 *  \param i  Particle index.
	 *  \return   Harmonic energy contribution of particle i. */
	double Energy(size_t i) {
		double E = 0;
		GeometryVector u0 = this->u[i];

		for (int j = 0; j < this->nnList[i].size(); j++) {
			GeometryVector u1 = this->u[this->nnList[i].at(j)] - u0;
			E += u1.Modulus2();
		}

		return 0.5*(this->K)*E;
	}

	/** \brief Compute the change in the energy by the displacement dx from the current position.
	 *  \param i   Index of the particle being displaced.
	 *  \param dx  Displacement applied to particle i.
	 *  \return    Resulting change in the harmonic energy. */
	double Delta_energy(size_t i, GeometryVector & dx) {
		double dE = 0;
		GeometryVector u0 = this->u[i];

		GeometryVector sum_du(static_cast<int> (this->Dimension));
		for (int j = 0; j < this->nnList[i].size(); j++) {
			sum_du = sum_du + this->u[this->nnList[i].at(j)];
		}

		dE = dx.Dot((double)nnList[i].size() *(dx + 2.0* u0) - 2.0*sum_du);

		return 0.5*(this->K)*dE;
	}

	/** \brief Return the total harmonic energy of the current configuration. (CLAUDE)
	 *  \return  Sum of the per-particle harmonic energies. */
	virtual double Energy()
	{
		double E = 0;
#pragma omp parallel for schedule(guided) reduction (+:E)
		for (int i = 0; i < this->nnList.size(); i++) {
			E += this->Energy(i);
		}

		return E; //consider the double counting, and the half coefficient
	}

	/** \brief Compute the harmonic force on the ith particle. (CLAUDE)
	 *  \param[out] F  Force on particle i.
	 *  \param i       Particle index. */
	virtual void Force(GeometryVector &F, size_t i) {
		F = GeometryVector(static_cast<int> (this->Dimension));
		GeometryVector u0 = this->u[i];

		for (int j = 0; j < this->nnList[i].size(); j++) {
			GeometryVector u1 = this->u[this->nnList[i].at(j)];
			F = F + u1;
		}

		F = this->K * (F - (double)(this->nnList[i].size())*u0);
	}

	/** \brief Compute the harmonic forces on all particles. (CLAUDE)
	 *  \param[out] Fs  Resized to hold the force on every particle. */
	virtual void AllForce(std::vector<GeometryVector> & Fs) {
		GeometryVector F(static_cast<int>(this->Dimension));
		Fs.clear();
		Fs.resize(this->confi.NumParticle(), F);

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < confi.NumParticle(); i++) {
			GeometryVector temp_F;
			this->Force(temp_F, i);
			Fs.at(i) = temp_F;
		}
	}
};

/** \brief Harmonic lattice potential expressed in normal coordinates (displacements built from normal-mode amplitudes and phases). (CLAUDE) */
class SimpleHarmonicF : public Potential
{
private:
	double K; //spring constant
	double a; //lattice constant
	Configuration Ref; //original positions: presumably a lattice
	std::vector<std::vector<size_t > > nnList; // list of prt indices of the nearest neighbors.
											   //std::vector<GeometryVector> du; // list of displacements
public:
	//const double a = 1; //lattice constant
	Configuration confi;
	int numPrt;
	double preF;
	DimensionType d;
	std::vector<GeometryVector> q; //wavevectors of normal modes
	std::vector<GeometryVector> u_n; //displacements
	std::vector<GeometryVector> u_n_trial; //trial displacements by u_q
	std::vector<GeometryVector> u_q; //normal coordinates
	std::vector<GeometryVector> phi_q; //phi_q[i].x[j] = phase constant of j-th component of normal coordinates at wavevector q[i]

	/** \brief Construct a normal-coordinate harmonic potential with the given dimension, spring constant, and lattice constant. (CLAUDE)
	 *  \param d                Spatial dimension.
	 *  \param springConstant   Spring constant K of the nearest-neighbor bonds.
	 *  \param latticeConstant  Lattice constant a (sets the nearest-neighbor search radius). */
	SimpleHarmonicF(DimensionType d, double springConstant, double latticeConstant) : Potential(d), K(springConstant), a(latticeConstant)
	{}

	/** \brief Set the reference configuration and build the nearest-neighbor list and the first-Brillouin-zone q-vectors. (CLAUDE)
	 *  \param ref  Reference lattice configuration; must match this potential's dimension. */
	void SetReference(const Configuration & ref) {

		std::vector<std::vector<size_t > >& temp = nnList;

		if (this->Dimension == ref.GetDimension()) {
			this->Ref = Configuration(ref);
			temp.reserve(this->Ref.NumParticle());

			numPrt = ref.NumParticle();
			preF = 1.0 / sqrt(numPrt);
			d = ref.GetDimension();
			q.reserve(numPrt);
			u_n.reserve(numPrt);
			u_n_trial.reserve(numPrt);
			u_q.reserve(numPrt);
			phi_q.reserve(numPrt);
			confi = Configuration(ref);


			//-----------------------make a list of nearest neighbors
			double r = a*1.1;
			for (int i = 0; i < this->Ref.NumParticle(); i++) {
				std::vector<size_t> nn;
				this->Ref.IterateThroughNeighbors(Ref.GetRelativeCoordinates(i), r,
					[&nn, i](const GeometryVector &x, const GeometryVector &LatticeShift, const signed long * PeriodicShift, const size_t SourceParticle)->void
				{
					if (i != SourceParticle) {
						nn.push_back(SourceParticle);
					}
				}
				);
				temp.push_back(nn);
			}

			//-----------------------make a list of q-vectors
			int systSize = (int)round(pow(numPrt, 1.0 / (double)d));
			std::vector<GeometryVector> q_temp;
			for (int i = 0; i < d; i++) {
				q_temp.push_back(Ref.GetReciprocalBasisVector(i));
			}
			//make a list of q;
			Configuration temp_q_arrary = Configuration(GetLattice_config(q_temp, systSize));
			GeometryVector relativeCoord(static_cast<int>(this->d));
			for (int i = 0; i < temp_q_arrary.NumParticle(); i++) {
				//form the first Brillouin zone
				relativeCoord = temp_q_arrary.GetRelativeCoordinates(i);
				for (int j = 0; j < d; j++) {
					if (relativeCoord.x[j]>0.5) {
						relativeCoord.x[j] -= 1.0;
					}
				}
				q.push_back(temp_q_arrary.RelativeCoord2CartesianCoord(relativeCoord));
			}

			this->q.erase(this->q.begin());

			//initialize variable list
			GeometryVector opener(static_cast<int>(d));
			for (int i = 0; i < this->q.size(); i++) {
				u_n.push_back(opener);
				u_n_trial.push_back(opener);
				u_q.push_back(opener);
				phi_q.push_back(opener);
			}
			u_n.push_back(opener);
			u_n_trial.push_back(opener);
		}
		else {
			std::cout << "SimpleHarmonic::SetReference::wrong dimensions\n";
		}
	}

	/** \brief Compute the displacements changed by the trial u_q and phase constants.
	 *  \param i_q        Index of the normal mode (q-vector) being varied.
	 *  \param du_q       Trial increment to the normal-coordinate amplitude u_q[i_q].
	 *  \param phi_trial  Trial phase constants for mode i_q. */
	void update_u_trial(size_t i_q, GeometryVector & du_q, GeometryVector & phi_trial) {
		GeometryVector uq_new = this->u_q[i_q] + du_q;
#pragma omp parallel for schedule(guided)
		for (int n = 0; n < numPrt; n++) {
			double qR = q[i_q].Dot(this->Ref.GetCartesianCoordinates(n));
			//GeometryVector du_n(static_cast<int>(d));
#pragma omp critical
			{
				for (int j = 0; j < d; j++) {
					this->u_n_trial[n].x[j] = this->u_n[n].x[j] + (this->preF)*(uq_new.x[j] * sin(qR + phi_trial.x[j]) - (this->u_q[i_q].x[j])*sin(qR + phi_q[i_q].x[j]));
				}
			}

		}
	}

	/** \brief Accept a trial change of a normal mode: add du_q to its amplitude and set its phase. (CLAUDE)
	 *  \param i          Index of the normal mode (q-vector).
	 *  \param du_q       Increment added to the normal-coordinate amplitude u_q[i].
	 *  \param phi_trial  New phase constants for mode i. */
	void update_u_q(size_t i, GeometryVector & du_q, GeometryVector &phi_trial) {
		this->u_q[i] = this->u_q[i] + du_q;
		this->phi_q[i] = phi_trial;
	}

	/** \brief Synchronize u_n as u_n_trial. */
	void update_un() {
#pragma omp parallel for schedule(guided)
		for (int n = 0; n < numPrt; n++) {
			this->u_n[n] = this->u_n_trial[n];
		}
	}

	/** \brief Update the displacement vector of i-th particle.
	 *  \param n  Index of the particle whose real-space displacement u_n is recomputed from the normal coordinates. */
	void update_u(size_t n) {
		GeometryVector result(static_cast<int>(d));
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < q.size(); i++) {
#pragma omp critical
			{
				for (int j = 0; j < d; j++) {
					result.x[j] += u_q[i].x[j] * sin(q[i].Dot(this->Ref.GetCartesianCoordinates(n)) + phi_q[i].x[j]);
				}
			}

		}
		this->u_n[n] = (this->preF)*result;
	}
	/** \brief Update configuration based on the displacements; if option = 1 => recompute u_n for u_q.
	 *  \param option  If 1, recompute every real-space displacement u_n from the normal coordinates before moving particles. */
	void UpdateConfiguration(int option) {
		if (option == 1) {
			for (int i = 0; i < numPrt; i++) {
				update_u(i);
			}
		}

		//#pragma omp parallel for
		for (int i = 0; i < numPrt; i++) {
			this->confi.MoveParticle(i, this->Ref.CartesianCoord2RelativeCoord(this->Ref.GetCartesianCoordinates(i) + this->u_n[i]));
		}
	}


	/** \brief Set the current configuration (must match the reference particle count and dimension). (CLAUDE)
	 *  \param Position  Configuration to evaluate the potential on. */
	virtual void SetConfiguration(const Configuration & Position) {
		if (this->Ref.NumParticle() == Position.NumParticle() && this->Dimension == Position.GetDimension())
			this->confi = Configuration(Position);
		else
			std::cout << "SimpleHarmonic::setConfiguration::sth wrong\n";
	}

	/** \brief Compute the energy of ith particle.
	 *  \param i  Particle index.
	 *  \return   Harmonic energy contribution of particle i. */
	double Energy(size_t i) {
		double E = 0;
		GeometryVector u0 = this->u_n[i];

		for (int j = 0; j < this->nnList[i].size(); j++) {
			GeometryVector u1 = u_n[nnList[i].at(j)] - u0;
			E += u1.Modulus2();
		}
		return 0.5*(this->K)*E;
	}

	/** \brief Compute the energy of the system from u_n_trial.
	 *  \return  Total harmonic energy evaluated on the trial displacements u_n_trial. */
	double Energy_trial() {
		double E = 0;
#pragma omp parallel for schedule (guided) reduction (+:E)
		for (int n = 0; n < this->nnList.size(); n++) {
			double E_temp = 0;
			GeometryVector u0 = this->u_n_trial[n];
			for (int j = 0; j < this->nnList[n].size(); j++) {
				GeometryVector u1 = u_n_trial[nnList[n].at(j)] - u0;
				E_temp += u1.Modulus2();
			}
			E += E_temp;
		}
		return E*0.25*(this->K);
	}

	/** \brief Return the total harmonic energy of the current real-space displacements u_n. (CLAUDE)
	 *  \return  Total harmonic energy. */
	virtual double Energy()
	{
		double E = 0;
#pragma omp parallel for schedule (guided) reduction (+:E)
		for (int n = 0; n < this->nnList.size(); n++) {
			double E_temp = 0;
			GeometryVector u0 = this->u_n[n];
			for (int j = 0; j < this->nnList[n].size(); j++) {
				GeometryVector u1 = u_n[nnList[n].at(j)] - u0;
				E_temp += u1.Modulus2();
			}
			E += E_temp;
		}
		return E*0.25*(this->K);
	} //consider the double counting, and the half coefficient


	/** \brief Compute the harmonic force on the ith particle. (CLAUDE)
	 *  \param[out] F  Force on particle i.
	 *  \param i       Particle index. */
	virtual void Force(GeometryVector &F, size_t i) {
		F = GeometryVector(static_cast<int> (this->Dimension));
		GeometryVector u0 = this->u_n[i];

		for (int j = 0; j < this->nnList[i].size(); j++) {
			GeometryVector u1 = this->u_n[this->nnList[i].at(j)];
			F = F + u1;
		}

		F = this->K * (F - (double)(this->nnList[i].size())*u0);
	}

	/** \brief Compute the harmonic forces on all particles. (CLAUDE)
	 *  \param[out] Fs  Resized to hold the force on every particle. */
	virtual void AllForce(std::vector<GeometryVector> & Fs) {
		GeometryVector F(static_cast<int>(this->Dimension));
		Fs.clear();
		Fs.resize(this->confi.NumParticle(), F);

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < confi.NumParticle(); i++) {
			GeometryVector temp_F;
			this->Force(temp_F, i);
			Fs.at(i) = temp_F;
		}
	}
};

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
		KPointConstraint(const GeometryVector & k, double V, double targetRhoReal, double targetRhoImag) : k(k), V(V), targetRhoImag(targetRhoImag), targetRhoReal(targetRhoReal)
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
	void SetConstraintType(char * type) {
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



/** \brief Hard-sphere Metropolis Monte Carlo potential: energy is zero unless spheres overlap (then infinite/one). (CLAUDE) */
class HS_pot : public Potential
{
private:
	double maxR;
	double dx_max;
	bool monoDISPESE = true; //maximum displacement is identical for all particles
	RandomGenerator rng;
public:
	SpherePacking config;
	double Acceptance = 0.5;
	/** \brief Construct a 1D hard-sphere potential with zero radius and step size. (CLAUDE) */
	HS_pot() : Potential(1)
	{
		this->maxR = 0; this->dx_max = 0; 
	}
	/** \brief Construct a hard-sphere potential of dimension d with a time-seeded random generator. (CLAUDE)
	 *  \param d  Spatial dimension. */
	HS_pot(DimensionType d) : Potential(d) { this->rng.seed(std::time(NULL)); }
	/** \brief Construct a hard-sphere potential of dimension d with a given maximum sphere radius. (CLAUDE)
	 *  \param d     Spatial dimension.
	 *  \param rMax  Maximum sphere radius used for overlap and cell-size bookkeeping. */
	HS_pot(DimensionType d, double rMax) :Potential(d)
	{
		this->maxR = rMax;
		this->rng.seed(std::time(NULL));
	}
	/** \brief Construct a hard-sphere potential of dimension d with the random generator seeded from s. (CLAUDE)
	 *  \param d  Spatial dimension.
	 *  \param s  Extra seed offset combined with the current time. */
	HS_pot(DimensionType d, size_t s) :Potential(d)
	{
		this->rng.seed(s + std::time(NULL));
	}

	/** \brief Set the sphere packing to evolve and update its max radius, cell size, and maximum displacement. (CLAUDE)
	 *  \param syst  Sphere packing whose dimension must match this potential. */
	virtual void SetConfiguration(const SpherePacking & syst) {
		if (this->Dimension == syst.GetDimension()) {
			this->config = SpherePacking(syst);
			this->config.UpdateMaxRadius();
			this->maxR = this->config.GetMaxRadius();
			config.SetCellSize(2.0*this->maxR);
			this->dx_max = pow(this->config.PeriodicVolume() / (double)this->config.NumParticle(), 1.0 / this->config.GetDimension()) - this->maxR;
		}
		else {
			std::cout << "wrong dimension!!!" << std::endl;
		}
	}
	/** \brief Configuration overload rejected: HS_pot only accepts SpherePacking systems. (CLAUDE)
	 *  \param c  Ignored; only present to satisfy the base-class interface. */
	virtual void SetConfiguration(const Configuration & c) {
		std::cout << "wrong type of systems \n";
	}


	/** \brief Return the hard-sphere energy: 1 if any spheres overlap, 0 otherwise. (CLAUDE)
	 *  \return  1.0 when the packing has an overlap, 0 otherwise. */
	virtual double Energy() {
		if (config.CheckOverlap())
			return 1.0; //spheres are overlapped
		else
			return 0;
	}
	/** \brief Set whether the maximum displacement is identical for all particles (monodisperse mode). (CLAUDE)
	 *  \param yes  If true, use a single maximum displacement for all particles. */
	void SetDispersity(bool yes) { this->monoDISPESE = yes;	}
	/** \brief Set the maximum trial displacement per Monte Carlo move. (CLAUDE)
	 *  \param dx  Maximum displacement magnitude. */
	void SetMaxDisplacement(double dx) { this->dx_max = dx; }

	/** \brief Check whether the ith particle can move to x, x: relative coordinates.
	 *  \param i  Index of the trial particle.
	 *  \param x  Candidate relative coordinates of particle i.
	 *  \return   1 if the move would create an overlap, 0 otherwise. */
	double Energy_trial(size_t i, const GeometryVector & x) {//check whether the ith particle can move to x, x: relative coordinates
		bool overlap = false;
		std::vector<GeometryVector> data; //data[i].x[0]: index, data[i].x[1]:displace
		double r_i = config.GetCharacteristics(i);
		double D2 = this->maxR* this->maxR*4.0;
		double D = this->maxR * 2.0;
		//obtain indices of particles that can overlap with the test particle
		if (monoDISPESE) {
			config.IterateThroughNeighbors(x, D,
				[i, D2, &overlap](const GeometryVector &r, const GeometryVector &LatticeShift, const signed long * PeriodicShift, const size_t SourceParticle)
				->void {double dist2 = r.Modulus2();	if (SourceParticle != i && dist2 < D2) { overlap = true; return; }});

			if (overlap)
				return 1.0;
			else
				return 0;
		}

		else {
			config.IterateThroughNeighbors(x, this->maxR,
				[i, D2, &data](const GeometryVector &r, const GeometryVector &LatticeShift, const signed long * PeriodicShift, const size_t SourceParticle)
				->void {double dist2 = r.Modulus2();	if (SourceParticle != i && dist2 < D2) { data.push_back(GeometryVector(SourceParticle, sqrt(dist2))); }});

			for (int j = 0; j < data.size(); j++) {
				if (data[j].x[1] < r_i + config.GetCharacteristics(data[j].x[0]))
					return 1;
			}
			return 0;
		}
	}
	/** \brief Move particle i to relative coordinates x in the current packing. (CLAUDE)
	 *  \param i  Particle index.
	 *  \param x  New relative coordinates for particle i. */
	void UpdateConfig(size_t i, const GeometryVector & x) {
		this->config.MoveParticle(i, x);
	}

	/** \brief Perform one Metropolis sweep (one trial move per particle) and return the acceptance fraction. (CLAUDE)
	 *  \return  Fraction of accepted moves in the sweep. */
	double Evolve() {
		double count = 0;
		DimensionType d = this->Dimension;
		GeometryVector u(d);
		GeometryVector x0(d);
		for (int i = 0; i < this->config.NumParticle(); i++) {
			u = (this->dx_max)*RandomUnitVector(d, this->rng);
			x0 = this->config.CartesianCoord2RelativeCoord(u + this->config.GetCartesianCoordinates(i));
			if (Energy_trial(i, x0) < 0.5) {
				this->config.MoveParticle(i, x0);
				count++;
			}
		}
		return count / (double)this->config.NumParticle();
	}

	/** \brief Execute Evole() successively for a given numbers.
	 *  \param numOfCycles  Number of sweeps to perform.
	 *  \return             Mean acceptance fraction over the sweeps. */
	double Evolve(int numOfCycles) {//execute Evole() successively for a given numbers
		double result = 0;
		for (int i = 0; i < numOfCycles; i++)
			result += this->Evolve();

		return result / (double)numOfCycles;
	}

	/** \brief Tune the maximum displacement so the Metropolis acceptance fraction approaches the target Acceptance. (CLAUDE)
	 *  \param numtrials  Number of sweeps used per adjustment iteration (and 5x that for equilibration). */
	void AdjustDisplacement(size_t numtrials) {
		double acc = 0;
		//int numtrials = 10;
		//----equilibrate----
		for (int i = 0; i < 5 * numtrials; i++)
			acc += Evolve();

		acc /= 5.0*numtrials;
		double factor = 0.1;
		double dev = fabs(acc - this->Acceptance);

		//int l = 0;
		//char name[20] = {};
		while (dev>0.02) {
			if (dev > 0.1)
				factor = 0.1;
			else
				factor = 0.01;

			if (acc > this->Acceptance)
				this->dx_max *= 1 + factor;
			else if (acc < this->Acceptance)
				this->dx_max *= 1 - factor;

			acc = 0;
			for (int i = 0; i < numtrials; i++)
				acc += Evolve();

			acc /= (double)numtrials;
			//std::cout << "Acceptance = " << acc << "\t maximum displacements = " << this->dx_max << "\n";
			dev = abs(acc - this->Acceptance);
			//	l++;
			//	sprintf(name, "%d", l);
			//	PlotPacking2D(name, this->config);
		}
	}

	/** \brief Return the current sphere packing. (CLAUDE)
	 *  \return  A copy of the current configuration. */
	SpherePacking Output() { return this->config; }
	/** \brief Return the current maximum trial displacement. (CLAUDE)
	 *  \return  The maximum displacement per Monte Carlo move. */
	double GetMaxDisplacement() {		return this->dx_max;	}
	
};