/** Author: Ge Zhang
	email:	?
	Date:	?	
	Comments: Jaeuk Kim (June 2019). */

/** \file Potential.h
	\brief Header file to compute potential energy of point particle systems.*/

#ifndef POTENTIAL_INCLUDED
#define POTENTIAL_INCLUDED

/* In $(cores) */
#include <etc.h>
#include <GeometryVector.h>
#include <PeriodicCellList.h>

#include <vector>
#include <exception>
#include <cmath>
#include <cassert>
#include <memory>

/** \brief Exception thrown when EnergyDerivativeToBasisVectors is not implemented for a potential class. */
class EnergyDerivativeToBasisVectors_NotImplemented : public std::exception
{
public:
	/** \brief Construct the exception. */
	EnergyDerivativeToBasisVectors_NotImplemented() : exception()
	{}
	/** \brief Return the exception's descriptive message.
		\return A C-string describing the error. */
	virtual const char* what() const throw ()
	{
		return "EnergyDerivativeToBasisVectors is not implemented for this class!\n";
	}
};

/** \brief A virtual class to define potential energy. */
class Potential
{
public:
	DimensionType Dimension;	//< Space dimension.

	//some potentials might allow parallelization.
	//if so, this is the number of threads
	unsigned long ParallelNumThread;	//< The number of threads in openMP.

	/** \brief A default constructor.
		\param dimension  Spatial dimension of the system. */
	Potential(DimensionType dimension):Dimension(dimension)
	{
		this->ParallelNumThread=1;
	}

	/** \brief Destructor. */
	virtual ~Potential();

	/** \brief A virtual function to compute the total potential energy.
	 \return the total potential energy of the given configuration.*/
	virtual double Energy() =0;

	/** \brief A virtual function to compute force at the ith particle.
		\param[out] result	Force.
		\param[in]	i		Particle index.*/
	virtual void Force(GeometryVector & /*result*/, size_t /*i*/)
	{
		std::cerr<<"Error : virtual void Force is not implemented in this Potential!\n";
	}

	/** \brief resize the vector to the number of particles, and calculate the force for all Particles.
	 \param[out] results	A list of forces. */
	virtual void AllForce(std::vector<GeometryVector> & /*results*/)
	{
		//not implemented, because
		//1. Don't know how many particles are there
		//2. Implementations, if parallelized, should ensure thread safety
		std::cerr<<"Error : virtual void AllForce is not implemented in this Potential!\n";
	}
	//Energy derivative respect to basis vectors
	//grad[n*dim+m] is the energy derivative of mth (coordinate) component of nth basis vector
	/** \brief Energy derivative respect to basis vectors.
		\param[out] grad      Array of size dim*dim; grad[n*dim+m] is the energy derivative of the mth component of the nth basis vector.
		\param Pressure       Imposed pressure used to form the enthalpy-like derivative. */
	virtual void EnergyDerivativeToBasisVectors(double * /*grad*/, double /*Pressure*/)
	{
		throw EnergyDerivativeToBasisVectors_NotImplemented();
	}

	//set the potential to calculate energy and/or force for a specific configuration
	//the configuration must exist when calculating energy and force
	/** \brief Set the configuration for which energy and/or force is computed.
		\param Position  Configuration to bind; it must remain alive while energy and force are evaluated. */
	virtual void SetConfiguration(const Configuration & Position) =0;

	/** \brief Compute the total potential energy and the force on every particle of the bound configuration.
		Exactly Energy() followed by AllForce(forces). A potential may override it with one combined evaluation
		only if that returns bitwise the same energy and forces. Declared after every other virtual function so
		that their vtable slots are unchanged.
		\param[out] forces  A list of forces, as AllForce.
		\return The total potential energy, as Energy(). */
	virtual double EnergyAndAllForce(std::vector<GeometryVector> & forces)
	{
		const double e=this->Energy();
		this->AllForce(forces);
		return e;
	}

	/** \brief Print a number very close to 1 when Force() is consistent with Energy(), and c is not a local minimum of this potential.
		\param c    Configuration on which to test the force/energy consistency.
		\param Nbr  Index of the particle to perturb. */
	void Check(Configuration c, size_t Nbr = 0);

	/** \brief Print a number very close to 1 when EnergyDerivativeToBasisVectors() is consistent with Energy(), and c is not a local minimum of this potential.
		\param c         Configuration on which to test the consistency.
		\param Nbr       Index of the basis-vector component to perturb.
		\param Pressure  Imposed pressure used to form the enthalpy-like derivative. */
	void Check2(Configuration c, size_t Nbr = 0, double Pressure = 1.0);

};

//a potential that have good performance for Monte Carlo should inherent from this class and implement these two virtual functions
class forMCPotential
{
public:
	/** \brief A virtual member function to compute the change in potential energy by a trial move.
		The potential must have been set to the correct configuration before this function is called.
		\param prevCart	Current position in Cartesian coordinates.
		\param afterCart	A trial position in Cartesian coordinates.
		\return change in energy. */
	virtual double TryMove(GeometryVector & prevCart, GeometryVector & afterCart) = 0;

	/** \brief A virtual member function to accept the trial move and change the current state accordingly. */
	virtual void AcceptMove(void) = 0;
};

//add additional particles to the configuration in a potential
//useful in simulations with some particles fixed.
//this is a quick-and-dirty implementation
template<typename T> class AdditionalParticlesPotentialAdaptor : public T
{
private:
	std::vector<GeometryVector> AdditionalParticleRelatives;
	Configuration ConfigurationWithAdditionalParticles;
	//size_t NParticles;
public:
	/** \brief Construct an adaptor that appends fixed extra particles to a base potential.
		\param AdditionalParticleRelatives  Relative coordinates of the extra particles to add to every configuration.
		\param src                          Base potential of type T to wrap. */
	AdditionalParticlesPotentialAdaptor(const std::vector<GeometryVector> & AdditionalParticleRelatives, const T & src): T(src), AdditionalParticleRelatives(AdditionalParticleRelatives)
		//, NParticles(0)
	{
	}

	/** \brief Set the configuration, inserting the additional fixed particles before delegating to the base potential.
		\param Position  Configuration of the mobile particles. */
	virtual void SetConfiguration(const Configuration & Position)
	{
		//NParticles = Position.NumParticles();
		ConfigurationWithAdditionalParticles = Position;
		for (auto iter = AdditionalParticleRelatives.begin(); iter != AdditionalParticleRelatives.end(); ++iter)
			ConfigurationWithAdditionalParticles.Insert("Z", *iter);
		this->T::SetConfiguration(ConfigurationWithAdditionalParticles);
	}
	//virtual void AllForce() is not re-implemented because does not know how to parallelize
};

/** \brief A potential expressed as a sum of pairwise interactions. */
class PairPotential: public Potential
{
private:
	const Configuration * pConfig;
public:
	double Rcut2;//Rcut^2
	double Rcut;


	/** \brief Construct a pair potential.
		\param dimension  Spatial dimension of the system.
		\param Rcut       Cutoff distance beyond which the pair interaction vanishes. */
	PairPotential(DimensionType dimension, double Rcut);

	/** \brief Compute the pair force for a given separation offset.
		\param offset      Separation vector between the two particles.
		\param[out] result Resulting force vector (zero if beyond the cutoff). */
	void PairForce(const GeometryVector & offset, GeometryVector & result);

	/** \brief Compute the pair energy for a given separation offset.
		\param offset  Separation vector between the two particles.
		\return The pair interaction energy (zero if beyond the cutoff). */
	double PairEnergy(const GeometryVector & offset);

	/** \brief Destructor. */
	virtual ~PairPotential();

	/** \brief Set the configuration for which energy and force are computed.
		\param Position  Configuration to bind; it must remain alive while energy and force are evaluated. */
	virtual void SetConfiguration(const Configuration & Position);

	/** \brief Compute the total pair-potential energy of the bound configuration.
		\return The total potential energy. */
	virtual double Energy();

	/** \brief calculate the force for the ith Particle.
		\param[out] result Resulting force on particle i.
		\param i           Index of the particle. */
	virtual void Force(GeometryVector & result, size_t i);//calculate the force for the ith Particle

	//contrary to setting a configuration and then using Energy() or Force()
	//this version does not change the state of the object unless EnergyFormula() or ForceFormula() changes the state of the object.
	//can be handy in threads
	/** \brief return the total potential energy of the given configuration.
		\param conf  Configuration to evaluate without changing the bound state.
		\return The total potential energy of conf. */
	double Energy(const Configuration & conf);//return the total potential energy of the given configuration
	/** \brief calculate the force for the ith Particle.
		\param conf        Configuration to evaluate without changing the bound state.
		\param[out] result Resulting force on particle i.
		\param i           Index of the particle. */
	virtual void Force(const Configuration & conf, GeometryVector & result, size_t i);//calculate the force for the ith Particle
	/** \brief Compute the force on every particle of the bound configuration.
		\param[out] results  Resized list of per-particle forces. */
	virtual void AllForce(std::vector<GeometryVector> & results);

	/** \brief Pair-energy formula as a function of the separation offset.
		\param offset  Separation vector between the two particles.
		\return The pair interaction energy. */
	virtual double EnergyFormula(const GeometryVector &) =0;

	/** \brief Pair-force formula as a function of the separation offset.
		\param offset      Separation vector between the two particles.
		\param[out] result Resulting pair force vector. */
	virtual void ForceFormula(const GeometryVector & offset, GeometryVector & result) =0;

	/** \brief Numerically estimate a mixed second derivative of the pair energy.
		\param offset  Separation vector between the two particles.
		\param i       Coordinate component to perturb.
		\param j       Force component to sample.
		\return The mixed energy derivative with respect to offset components i and j, computed as minus the derivative of force component j with respect to component i. */
	virtual double SecondDerivative(const GeometryVector & offset, DimensionType i, DimensionType j);//(dFi)/(dxj), recommended to rewrite it

	//Energy derivative respect to basis vectors
	//grad[n*dim+m] is the energy derivative of mth (coordinate) component of nth basis vector
	/** \brief Energy derivative respect to basis vectors.
		\param[out] grad  Array of size dim*dim; grad[n*dim+m] is the energy derivative of the mth component of the nth basis vector.
		\param Pressure   Imposed pressure used to form the enthalpy-like derivative. */
	virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure);

	/** \brief Energy derivative respect to basis vectors for a given configuration.
		\param conf       Configuration to bind before computing the derivative.
		\param[out] grad  Array of size dim*dim; grad[n*dim+m] is the energy derivative of the mth component of the nth basis vector.
		\param Pressure   Imposed pressure used to form the enthalpy-like derivative. */
	void EnergyDerivativeToBasisVectors(const Configuration & conf, double * grad, double Pressure);
};

/** \brief A class derived from PairPotential class. It is useful when pair potential is isotropic. */
class IsotropicPairPotential: public PairPotential
{
public:
	double Ecut;	//< A shift in energy such that energy vanishes as r goes to infinity.
	bool EcutUnset; //< true = Ecut is not initialized.

	/** \brief A deafult constructor. This function do NOT initialize this->Ecut.
		\param dimension  Spatial dimension of the system.
		\param Rcut       Cutoff distance of the potential. */
	IsotropicPairPotential(DimensionType dimension, double Rcut);

	/** \brief Set this->Ecut and this->EcutUnset = false. */
	void SetEcut(void);

	/** \brief Destructor. */
	virtual ~IsotropicPairPotential();

	/** \brief Pair-energy formula as a function of the separation offset.
		\param offset  Separation vector between the two particles.
		\return The pair interaction energy. */
	virtual double EnergyFormula(const GeometryVector & offset);

	/** \brief Pair-force formula as a function of the separation offset.
		\param offset      Separation vector between the two particles.
		\param[out] result Resulting pair force vector. */
	virtual void ForceFormula(const GeometryVector & offset, GeometryVector & result);

	/** \brief Define pair potetial. Rcut and Ecut are not implemented.
		\param distance  Separation between the two particles.
		\return potential.	 */
	virtual double IsotropicEnergyFormula(double distance) =0;

	/** \brief Compute the force. By default, it performes numerical differentiation.
		It is strongly recommanded to redefine this function!.
		\param distance  Separation between the two particles.
		\return the magnitude of force. */
	virtual double IsotropicForceFormula(double distance);//\frac{- d \Psi(r)}{d r}, recommended to rewrite it

	//cutoff and shifted
	/** \brief Cutoff-and-shifted isotropic pair energy.
		\param distance  Separation between the two particles.
		\return The pair energy, shifted so it vanishes at the cutoff and zero beyond it. */
	inline double IsotropicPairEnergy(double distance)
	{
		if(distance>this->Rcut) return 0;
		if(this->EcutUnset)
			this->SetEcut();
		return this->IsotropicEnergyFormula(distance)-this->Ecut;
	}
};

/** \brief An isotropic pair potential that tabulates a source potential and interpolates it linearly. */
class LinearInterpolationIsotropicPairPotential : public IsotropicPairPotential
{
	std::vector<double> data;
	double Start, stepsize;
public:
	/** \brief Build the interpolation table by sampling a source potential.
		\param src     Source isotropic pair potential to tabulate.
		\param Start   Smallest separation at which the table starts.
		\param pieces  Number of interpolation intervals between Start and the cutoff. */
	LinearInterpolationIsotropicPairPotential(IsotropicPairPotential & src, double Start, size_t pieces) : IsotropicPairPotential(src.Dimension, src.Rcut)
	{
		this->Start=Start;
		this->data.resize(pieces+2);
		this->stepsize=(src.Rcut-Start)/pieces;
		for(size_t i=0; i<this->data.size(); i++)
			this->data[i]=src.IsotropicPairEnergy(Start+i*this->stepsize);
	}

	/** \brief Linearly interpolated pair energy at a given separation.
		\param distance  Separation between the two particles.
		\return The interpolated pair energy. */
	virtual double IsotropicEnergyFormula(double distance)
	{
		if(distance<this->Start)
			return this->data[0];
		double dindex=(distance-this->Start)/this->stepsize;
		size_t iindex=(size_t)(std::floor(dindex));
		double rindex=dindex-iindex;//remaining index
		return this->data[iindex]+(this->data[iindex+1]-this->data[iindex])*rindex;
	}
};

namespace
{
	const double AspectRatioSquared=4;
	const signed MaxI=6;
	const signed MaxJ=6;
	const double CorrectionCenter=1.0;
	const double CorrectionWidth=0.1;
	/** \brief Internal correction potential applied over a narrow window of separations. */
	class Correction : public IsotropicPairPotential
	{
		IsotropicPairPotential * orig;
		double V0;
	public:
		/** \brief Construct the correction by precomputing the reference offset V0 from lattice sums.
			\param original  Source isotropic pair potential to correct. */
		Correction(IsotropicPairPotential * original) : IsotropicPairPotential(original->Dimension, CorrectionCenter+2*CorrectionWidth)
		{
			this->orig=original;
			this->V0=0;
			for(signed i=(-1)*MaxI; i<=MaxI; i++)
			{
				for(signed j=(-1)*MaxJ; j<=MaxJ; j++)
				{
					double ii=i*i;
					double jj=j*j;
					if(i==0 && j==0)
						continue;
					V0-=orig->IsotropicPairEnergy(std::sqrt(ii+AspectRatioSquared*jj))*(ii-AspectRatioSquared*jj)/(ii+AspectRatioSquared*jj);
				}
			}
		}

		/** \brief Correction energy applied within the correction window.
			\param distance  Separation between the two particles.
			\return The correction energy (zero outside the correction window). */
		virtual double IsotropicEnergyFormula(double distance)
		{
			if(distance<CorrectionCenter-2*CorrectionWidth)
				return 0;
			if (distance>CorrectionCenter+2*CorrectionWidth)
				return 0;

			double V1=(-1)*this->V0;
			for(signed i=(-1)*MaxI; i<=MaxI; i++)
			{
				for(signed j=(-1)*MaxJ; j<=MaxJ; j++)
				{
					double ii=i*i;
					double jj=j*j;
					if(i==0 && j==0)
						continue;
					V1-=orig->IsotropicPairEnergy(distance*std::sqrt(ii+AspectRatioSquared*jj))*(ii-AspectRatioSquared*jj)/(ii+AspectRatioSquared*jj);
				}
			}
			V1/=2;
			double CorrectionCoeff;
			if(distance<CorrectionCenter-CorrectionWidth)
				CorrectionCoeff=( distance-(CorrectionCenter-2*CorrectionWidth) )/CorrectionWidth;
			else if(distance>CorrectionCenter+CorrectionWidth)
				CorrectionCoeff=( (CorrectionCenter+2*CorrectionWidth)-distance )/CorrectionWidth;
			else
				CorrectionCoeff=1;

			return V1*CorrectionCoeff;
		}
	};
};
/** \brief An isotropic pair potential equal to a source potential plus a localized correction term. */
class CorrectedIsotropicPairPotential : public IsotropicPairPotential
{
	std::unique_ptr<IsotropicPairPotential> orig;
	std::unique_ptr<LinearInterpolationIsotropicPairPotential> correction;

public:
	//the object (*original) will be freed when this class destruct
	/** \brief Construct from a source potential and build its tabulated correction.
		\param original  Source isotropic pair potential; it is owned and freed on destruction. */
	CorrectedIsotropicPairPotential(IsotropicPairPotential * original) : IsotropicPairPotential(original->Dimension, original->Rcut)
	{
		this->orig.reset(original);
		Correction corr(original);
		this->correction.reset(new LinearInterpolationIsotropicPairPotential(corr, CorrectionCenter-2*CorrectionWidth, 500));
	}

	/** \brief Corrected pair energy: source energy plus correction within the correction window.
		\param distance  Separation between the two particles.
		\return The corrected pair energy. */
	virtual double IsotropicEnergyFormula(double distance)
	{
		double result=orig->IsotropicEnergyFormula(distance);
		if(distance>CorrectionCenter-2*CorrectionWidth && distance<CorrectionCenter+2*CorrectionWidth)
			result+=this->correction->IsotropicEnergyFormula(distance);

		return result;
	}

	/** \brief Destructor; the unique_ptr members free the correction table and the owned source potential. */
	~CorrectedIsotropicPairPotential() = default;
};

/** \brief An isotropic pair potential equal to the sum of two isotropic pair potentials. */
class CombinedIsotropicPairPotential : public IsotropicPairPotential
{
private:
	std::unique_ptr<IsotropicPairPotential> p1;
	std::unique_ptr<IsotropicPairPotential> p2;
public:
	/** \brief Construct from two isotropic pair potentials of equal dimension.
		\param p1  First component potential; owned and freed on destruction.
		\param p2  Second component potential; owned and freed on destruction. */
	CombinedIsotropicPairPotential(IsotropicPairPotential * p1, IsotropicPairPotential * p2) : IsotropicPairPotential(p1->Dimension, p1->Rcut>p2->Rcut?p1->Rcut:p2->Rcut)
	{
		assert(p1->Dimension==p2->Dimension);
		this->p1.reset(p1);
		this->p2.reset(p2);
	}

	/** \brief Sum of the two component pair energies at a given separation.
		\param distance  Separation between the two particles.
		\return The combined pair energy. */
	virtual double IsotropicEnergyFormula(double distance)
	{
		return p1->IsotropicPairEnergy(distance)+p2->IsotropicPairEnergy(distance);
	}

	/** \brief Destructor; the unique_ptr members free both owned component potentials. */
	~CombinedIsotropicPairPotential() = default;
};
/** \brief A potential equal to a weighted sum of two arbitrary potentials. */
class CombinedPotential : public Potential
{
public:
	//p1/p2 are non-owning observers: this class never frees them (see the destructor).
	Potential * p1;
	Potential * p2;
	double Weight1, Weight2;
	/** \brief Construct from two potentials of equal dimension, each with unit weight.
		\param p1  First component potential.
		\param p2  Second component potential. */
	CombinedPotential(Potential * p1, Potential * p2): Potential(p1->Dimension)
	{
		assert(p1->Dimension==p2->Dimension);
		this->p1=p1;
		this->p2=p2;
		this->Weight1=1.0;
		this->Weight2=1.0;
	}

	/** \brief Destructor; does not free the component potentials. */
	virtual ~CombinedPotential()
	{
		//delete this->p1;
		//delete this->p2;
	}

	/** \brief return the total potential energy of the given configuration.
		\return The weighted sum of the two component energies. */
	virtual double Energy()//return the total potential energy of the given configuration
	{
		return this->p1->Energy()*this->Weight1 + this->p2->Energy()*this->Weight2;
	}

	/** \brief calculate the force for the ith Particle.
		\param[out] result  Weighted sum of the two component forces on particle i.
		\param i            Index of the particle. */
	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		GeometryVector f1, f2;
		this->p1->Force(f1, i);
		this->p2->Force(f2, i);
		result= f1*this->Weight1 + f2*this->Weight2;
	}

	/** \brief calculate the force for all Particles.
		\param[out] results  Weighted sum of the two component force lists. */
	virtual void AllForce(std::vector<GeometryVector> & results)//calculate the force for all Particles
	{
		std::vector<GeometryVector> temp;
		this->p1->AllForce(results);
		this->p2->AllForce(temp);
		for(size_t i=0; i<results.size(); i++)
			results[i]=results[i]*this->Weight1+temp[i]*this->Weight2;
	}
	//Energy derivative respect to basis vectors
	//grad[n*dim+m] is the energy derivative of mth (coordinate) component of nth basis vector
	/** \brief Energy derivative respect to basis vectors.
		\param[out] grad  Array of size dim*dim; grad[n*dim+m] is the energy derivative of the mth component of the nth basis vector.
		\param Pressure   Imposed pressure used to form the enthalpy-like derivative. */
	virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	{
		size_t dd=this->Dimension*this->Dimension;
		double * temp = new double[dd];
		this->p1->EnergyDerivativeToBasisVectors(grad, Pressure);
		this->p2->EnergyDerivativeToBasisVectors(temp, Pressure);
		for(size_t i=0; i<dd; i++)
			grad[i]=grad[i]*this->Weight1+temp[i]*this->Weight2;
		delete [] temp;
	}

	/** \brief Set the same configuration on both component potentials.
		\param Position  Configuration to bind. */
	virtual void SetConfiguration(const Configuration & Position)
	{
		this->p1->SetConfiguration(Position);
		this->p2->SetConfiguration(Position);
	}
};

/** \brief A cutoff-and-shifted isotropic pair potential whose unshifted energy is c * r^power. */
class RnPotential : public IsotropicPairPotential
{
public:
	double c, power;
	/** \brief Construct a power-law pair potential.
		\param dimension  Spatial dimension of the system.
		\param c          Prefactor of the power law.
		\param power      Exponent of the separation distance.
		\param Cutoff     Cutoff distance of the potential. */
	RnPotential(DimensionType dimension, double c, double power, double Cutoff) : IsotropicPairPotential(dimension, Cutoff), c(c), power(power)
	{}
	/** \brief Power-law pair energy c * Distance^power.
		\param Distance  Separation between the two particles.
		\return The pair energy. */
	virtual double IsotropicEnergyFormula(double Distance)
	{
		return c*std::pow(Distance, power);
	}

	/** \brief Signed radial power-law pair force -dV/dr.
		\param distance  Separation between the two particles.
		\return The signed radial force -c * power * distance^(power-1). */
	virtual double IsotropicForceFormula(double distance)
	{
		return (-1)*c*power*std::pow(distance, power-1.0);
	}
};





#endif
