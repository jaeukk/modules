/**
*	Author	: Ge Zhang
*	Email	:  - 	
*	Date	: 2015 */

/** \file CollectiveCoordinatePotential.h
 * \brief Header file to declare a class to compute collective-coordinate potential, including stealthy and nonstealthy variants (fixed simulation cell, plus the ForMC adapter).
 *  Split 2026-07-07: variable-cell machinery -> CCPotential_varBox.h; concrete kernels and named typedefs -> CCPotential_Kernels.h. */


#ifndef COLLECTIVECOORDINATEPOTENTIAL_INCLUDED
#define COLLECTIVECOORDINATEPOTENTIAL_INCLUDED

#include <omp.h>
#include <stdexcept>
#include <typeinfo>

#include "Potential.h"

//#define HAVETWOFIXEDPARTICLES


// wanted to write a class for Ewald summation. However, it seems to be hard to find a precise documentation
// that is more general than those online. Online descriptions usually only deal with 1/r^2 potential
// When doing it, it might be good to refer to Uche's 2004 collective coordinate paper
// Or refer to Sal's duality paper.

// \phi = N/(2V) \sum_{\mathbf k} u(\mathbf k) | \rho(\mathbf k) |^2 + this->Shift
/* References: 
	
	PHYSICAL REVIEW RESEARCH 5, 033190 (2023)
  */
// compared to definition in uche2004constraints, it is shifted and scaled

/** \brief Collective-coordinate potential: a weighted sum over k-space constraints of the structure factor. */
class ShiftedCCPotential : public Potential
{
protected:
	const Configuration * pConfig;
	double Volume;
public:
	/** \brief A single k-point constraint: a wavevector and its weight. */
	struct KPointConstraint
	{
		GeometryVector k;//absolute k, not relative to basis vectors
		double V;
		/** \brief Construct a k-point constraint.
			\param k  Absolute wavevector (not relative to basis vectors).
			\param V  Weight assigned to this k-point. */
		KPointConstraint(const GeometryVector & k, double V): k(k), V(V)
		{
		}
	};
	//users or derived classes should fill this vector with k-constraints
	std::vector<KPointConstraint> constraints;

	//this value is added to the energy. It is initialized as 0, but user can change it.
	double Shift;

	//a member function to define constraints
	/** \brief a member function to define constraints.
		\param k     Wavevector of the constraint to add.
		\param vals  Constraint values; vals[0] is the weight V(k). */
	virtual void AddConstraint(const GeometryVector & k, std::vector<double> vals){
		if(vals.empty())
			throw std::invalid_argument("ShiftedCCPotential::AddConstraint: vals must contain at least the weight V(k).");
		this->constraints.emplace_back(k, vals[0]);
	}

	/** \brief Print all k-point constraints and their weights.
		\param ofile  Output stream to print to (defaults to std::cout). */
	virtual void PrintConstraints(std::ostream & ofile = std::cout){
		for (size_t i =0; i<this->constraints.size(); i++){
			ofile << "@ k= " << this->constraints[i].k << ": v="<<this->constraints[i].V <<"\n";
		}
	}

	//this class will manage these members, users can also read it.
	std::vector<double> RhoReal, RhoImag;
	bool RhoValid;
	/** \brief Compute and cache the complex collective density rho(k) at every constraint. */
	virtual void GetRho(void)
	{
		assert(this->pConfig!=nullptr);
		if(this->RhoValid)
		{
			return;
		}
		this->Volume=this->pConfig->PeriodicVolume();

		signed long end=this->constraints.size();
		this->RhoReal.resize(end);
		this->RhoImag.resize(end);

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for(signed long i=0; i<end; i++)
		{
			double tr=0, ti=0;//real and imag parts of result
#ifdef HAVETWOFIXEDPARTICLES
			tr+=2;
#endif
			for(size_t j=0; j<pConfig->NumParticle(); j++)
			{
				GeometryVector c=pConfig->GetCartesianCoordinates(j);
				tr+=std::cos(this->constraints[i].k.Dot(c));
				ti+=std::sin(this->constraints[i].k.Dot(c));
			}
			this->RhoReal[i]=tr;
			this->RhoImag[i]=ti;
		}

		this->RhoValid=true;
	}

	/** \brief Construct an empty collective-coordinate potential.
		\param dimension  Spatial dimension of the system. */
	ShiftedCCPotential(DimensionType dimension):Potential(dimension)
	{
		this->pConfig=nullptr;
		this->RhoValid=false;
		this->ParallelNumThread=1;
		this->Shift=0.0;
	}

	/** \brief Destructor. */
	virtual ~ShiftedCCPotential()
	{
	}

	/** \brief Return a heap-allocated copy of this object's ShiftedCCPotential base subobject.
		\return A new ShiftedCCPotential copy. */
	virtual ShiftedCCPotential * clone()
	{
		//This copies only the ShiftedCCPotential base subobject, so a derived object would be sliced.
		//No derived class overrides clone() and there are no callers today; guard against silent misuse.
		assert(typeid(*this)==typeid(ShiftedCCPotential));
		return new ShiftedCCPotential(*this);
	}

	/** \brief return the total potential energy of the given configuration.
		\return The collective-coordinate energy plus Shift. */
	virtual double Energy()//return the total potential energy of the given configuration
	{
		if(this->RhoValid==false)
			this->GetRho();

		double result=0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end=this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long i=0; i<end; i++)
				threadresult[numthread]+=this->constraints[i].V*(this->RhoReal[i]*this->RhoReal[i] + this->RhoImag[i]*this->RhoImag[i]);
		}
		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result+=*iter;

		return result/(Volume)+this->Shift;
	}

	/** \brief calculate the force for the ith Particle.
		\param[out] result  Force on particle i.
		\param i            Index of the particle. */
	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		if(this->RhoValid==false)
			this->GetRho();

		result=GeometryVector(this->Dimension);

		signed long end=this->constraints.size();

			for(signed long j=0; j<end; j++)
			{
				GeometryVector k=this->constraints[j].k;
				GeometryVector r=this->pConfig->GetCartesianCoordinates(i);
				result.MinusFrom((2.0)*this->constraints[j].V*(this->RhoImag[j]*std::cos(k.Dot(r))-this->RhoReal[j]*std::sin(k.Dot(r)))*k);
			}

		result.MultiplyFrom(1.0/Volume);
	}

	/** \brief Compute the force on every particle.
		\param[out] results  Resized list of per-particle forces. */
	virtual void AllForce(std::vector<GeometryVector> & results)
	{
		if(this->RhoValid==false)
			this->GetRho();
		long end=this->pConfig->NumParticle();
		results.resize(end);
		unsigned long tempNumThread = 1;
		std::swap(tempNumThread, this->ParallelNumThread);
#pragma omp parallel for num_threads(tempNumThread) default(shared)
		for (long i = 0; i<end; i++)
			this->Force(results[i], i);
		std::swap(tempNumThread, this->ParallelNumThread);
	}
	//Energy derivative respect to basis vectors
	//grad[n*dim+m] is the energy derivative of mth (coordinate) component of nth basis vector
	//virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	//{
	//	std::cerr<<"Warning : CollectiveCoordinatePotential::EnergyDerivativeToBasisVectors is not implemented!\n";
	//}
	/** \brief Bind the configuration and invalidate the cached collective density.
		\param Position  Configuration to bind; must remain alive while energy and force are evaluated. */
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->Dimension==Position.GetDimension());
		this->pConfig=&Position;
		this->RhoValid=false;
	}
};

//same as ShiftedCCPotential, but allows user to specify a weight for each particle
//rho=\sum_j weight_j * exp (i*k*r_j)
/** \brief Collective-coordinate potential that weights each particle's contribution to rho(k). */
class ShiftedCCPotential_withWeight : public ShiftedCCPotential
{
protected:
	std::vector<double> weights;
public:
	/** \brief Construct from a base potential and per-particle weights.
		\param Pot      Base collective-coordinate potential to copy.
		\param weights  Weight for each particle in the collective density. */
	ShiftedCCPotential_withWeight(const ShiftedCCPotential & Pot, const std::vector<double> & weights) : ShiftedCCPotential(Pot), weights(weights)
	{
	}

	/** \brief Return a heap-allocated copy of this full derived object.
		\return A new ShiftedCCPotential_withWeight copy. */
	virtual ShiftedCCPotential_withWeight * clone()
	{
		return new ShiftedCCPotential_withWeight(*this);
	}

	/** \brief Compute and cache the weighted complex collective density rho(k). */
	virtual void GetRho(void)
	{
		assert(this->pConfig != nullptr);
		assert(this->pConfig->NumParticle() == weights.size());
		if (this->RhoValid)
		{
			return;
		}
		this->Volume = this->pConfig->PeriodicVolume();

		signed long end = this->constraints.size();
		this->RhoReal.resize(end);
		this->RhoImag.resize(end);

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for (signed long i = 0; i<end; i++)
		{
			double tr = 0, ti = 0;//real and imag parts of result
#ifdef HAVETWOFIXEDPARTICLES
			tr += 2;
#endif
			for (size_t j = 0; j<pConfig->NumParticle(); j++)
			{
				GeometryVector c = pConfig->GetCartesianCoordinates(j);
				tr += std::cos(this->constraints[i].k.Dot(c))*weights[j];
				ti += std::sin(this->constraints[i].k.Dot(c))*weights[j];
			}
			this->RhoReal[i] = tr;
			this->RhoImag[i] = ti;
		}
		this->RhoValid = true;
	}

	/** \brief calculate the force for the ith Particle (weighted).
		\param[out] result  Force on particle i.
		\param i            Index of the particle. */
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
			result.MinusFrom((2.0)*this->constraints[j].V*(this->RhoImag[j] * std::cos(k.Dot(r)) - this->RhoReal[j] * std::sin(k.Dot(r)))*k);
		}

		result.MultiplyFrom(weights[i] / Volume);
	}

	/** \brief Compute the force on every particle.
		\param[out] results  Resized list of per-particle forces. */
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

	/** \brief Bind the configuration and invalidate the cached collective density.
		\param Position  Configuration to bind. */
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->Dimension == Position.GetDimension());
		this->pConfig = &Position;
		this->RhoValid = false;
	}
};



// \phi = log ( 1/(2V) \sum_{\mathbf k} u(\mathbf k) | \rho(\mathbf k) |^2 )
/** \brief Collective-coordinate potential equal to the logarithm of the underlying ShiftedCCPotential energy. */
class LogCCPotential : public ShiftedCCPotential
{
public:
	/** \brief Construct from a base collective-coordinate potential.
		\param pot  Base potential whose energy is log-transformed. */
	LogCCPotential( const ShiftedCCPotential & pot ) : ShiftedCCPotential(pot)
	{
	}

	/** \brief Return a heap-allocated copy of this full derived object.
		\return A new LogCCPotential copy. */
	virtual LogCCPotential * clone()
	{
		return new LogCCPotential(*this);
	}

	/** \brief Return the logarithm of the base potential's energy.
		\return The log-transformed energy. */
	virtual double Energy()
	{
		return std::log( ShiftedCCPotential::Energy() );
	}

	/** \brief Compute the force of the log potential on particle i.
		\param[out] result  Force on particle i.
		\param i            Index of the particle. */
	virtual void Force(GeometryVector & result, size_t i)
	{
		ShiftedCCPotential::Force(result, i);
		result.MultiplyFrom( 1.0/ShiftedCCPotential::Energy() );
	}
};

// \phi = \sum_{\mathbf k} V(\mathbf k) ( | \rho(\mathbf k) |^2/N - S0(\mathbf k) )^2 + this->Shift
/** \brief Collective-coordinate potential penalizing the squared deviation of S(k) from a target S0(k). */
class ShiftedCCPotential_S2 : public ShiftedCCPotential
{
public:
	std::vector<double> S0;
	size_t N;
	/** \brief Construct an empty S2 collective-coordinate potential.
		\param d  Spatial dimension of the system. */
	ShiftedCCPotential_S2(DimensionType d) : ShiftedCCPotential(d)
	{
	}

	/** \brief Return a heap-allocated copy of this full derived object.
		\return A new ShiftedCCPotential_S2 copy. */
	virtual ShiftedCCPotential_S2 * clone()
	{
		return new ShiftedCCPotential_S2(*this);
	}

	//a member function to define constraints
	/** \brief a member function to define constraints.
		\param k     Wavevector of the constraint to add.
		\param vals  Constraint values; vals[0] is the weight V(k), vals[1] is the target S0(k). */
	virtual void AddConstraint(const GeometryVector & k, std::vector<double> vals){
		if(vals.size()<2)
			throw std::invalid_argument("ShiftedCCPotential_S2::AddConstraint: vals must contain the weight V(k) and the target S0(k).");
		this->constraints.emplace_back(k, vals[0]);
		this->S0.emplace_back(vals[1]);
	}

	/** \brief Print all k-point constraints with their weights and target structure factors.
		\param ofile  Output stream to print to (defaults to std::cout). */
	virtual void PrintConstraints(std::ostream & ofile = std::cout){
		for (size_t i =0; i<this->constraints.size(); i++){
			ofile << "@ k= " << this->constraints[i].k << ": v="<<this->constraints[i].V <<", S0(k)="<< S0[i]<<"\n";
		}
	}

	/** \brief return the total potential energy of the given configuration.
		\return The weighted sum of squared S(k) deviations plus Shift. */
	virtual double Energy() //return the total potential energy of the given configuration
	{
		if(this->RhoValid==false)
			this->GetRho();

		double result=0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end=this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long i=0; i<end; i++)
			{
				double s=(this->RhoReal[i]*this->RhoReal[i] + this->RhoImag[i]*this->RhoImag[i])/N;
				double deltaS=(s-this->S0[i]);
				threadresult[numthread]+=this->constraints[i].V*deltaS*deltaS;
			}
		}
		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result+=*iter;

		result+=this->Shift;

		return result;
	}

	/** \brief calculate the force for the ith Particle.
		\param[out] result  Force on particle i.
		\param i            Index of the particle. */
	virtual void Force(GeometryVector & result, size_t i) //calculate the force for the ith Particle
	{
		if(this->RhoValid==false)
			this->GetRho();

		result=GeometryVector(this->Dimension);
		std::vector<GeometryVector> threadresult(this->ParallelNumThread, result);

		signed long end=this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long j=0; j<end; j++)
			{
				GeometryVector k=this->constraints[j].k;
				GeometryVector r=this->pConfig->GetCartesianCoordinates(i);
				double s=(this->RhoReal[j]*this->RhoReal[j] + this->RhoImag[j]*this->RhoImag[j])/N;
				double deltaS=(s-this->S0[j]);
				threadresult[numthread].MinusFrom((4.0/this->N)*this->constraints[j].V*deltaS*(this->RhoImag[j]*std::cos(k.Dot(r))-this->RhoReal[j]*std::sin(k.Dot(r)))*k);
			}
		}

		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result.AddFrom(*iter);
	}

	/** \brief Bind the configuration and record the particle count N.
		\param Position  Configuration to bind. */
	virtual void SetConfiguration(const Configuration & Position)
	{
		this->N=Position.NumParticle();
		this->ShiftedCCPotential::SetConfiguration(Position);
	}
};

// \phi = \sum_{\mathbf k} V(\mathbf k) ( <S(\mathbf k)> - S0(\mathbf k) )^2 + this->Shift
// where <...> is average over sub-configurations.
//treat the configuration as NSubConfig sub-configurations.
//particles 1 to NumParticle/NSubConfig is in the first sub-configuration
//<S(\mathbf k)> means S(\mathbf k) avaraged over all sub-configurations 
/** \brief S2 collective-coordinate potential that targets S(k) averaged over sub-configurations. */
class ShiftedCCPotential_S2_MultiConfig : public ShiftedCCPotential
{
public:
	std::vector<double> S0;
	std::vector<double> Snow;//updated in GetRho
	size_t NSubConfig;//number of separate sub-configurations
	size_t Np;//number of particles per sub-configuration
	/** \brief Construct an empty multi-configuration S2 potential.
		\param d          Spatial dimension of the system.
		\param NSubConfig  Number of sub-configurations the system is split into. */
	ShiftedCCPotential_S2_MultiConfig(DimensionType d, size_t NSubConfig) : ShiftedCCPotential(d), NSubConfig(NSubConfig)
	{
	}

	/** \brief Return a heap-allocated copy of this full derived object.
		\return A new ShiftedCCPotential_S2_MultiConfig copy. */
	virtual ShiftedCCPotential_S2_MultiConfig * clone()
	{
		return new ShiftedCCPotential_S2_MultiConfig(*this);
	}

	/** \brief Compute rho(k) per sub-configuration and the sub-configuration-averaged S(k). */
	void GetRho_MultiConfig(void)
	{
		assert(this->pConfig!=nullptr);
		if(this->RhoValid)
		{
			return;
		}

		signed long end=this->constraints.size();
		this->RhoReal.resize(end*NSubConfig);
		this->RhoImag.resize(end*NSubConfig);
		this->Snow = std::vector<double>(end, 0.0);

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for(signed long nc=0; nc<(signed long)NSubConfig; nc++)
		{
			for(signed long i=0; i<end; i++)
			{
				double tr=0, ti=0;//real and imag parts of result
				for(size_t j=nc*Np; j<(nc+1)*Np; j++)
				{
					GeometryVector c=pConfig->GetCartesianCoordinates(j);
					tr+=std::cos(this->constraints[i].k.Dot(c));
					ti+=std::sin(this->constraints[i].k.Dot(c));
				}
				this->RhoReal[nc*end+i]=tr;
				this->RhoImag[nc*end+i]=ti;

#pragma omp atomic
				this->Snow[i] += (tr*tr+ti*ti);
			}
		}
		for (signed long i = 0; i < end; i++)
			this->Snow[i] /= (Np*NSubConfig);
		this->RhoValid = true;
	}

	/** \brief Compute and cache the collective density (delegates to GetRho_MultiConfig). */
	virtual void GetRho(void)
	{
		GetRho_MultiConfig();
	}

	/** \brief return the total potential energy of the given configuration.
		\return The weighted sum of squared deviations of the averaged S(k) from S0(k) plus Shift. */
	virtual double Energy() //return the total potential energy of the given configuration
	{
		if(this->RhoValid==false)
			this->GetRho_MultiConfig();

		double result=0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end=this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long i=0; i<end; i++)
			{
				double s = this->Snow[i];
				double deltaS=(s-this->S0[i]);
				threadresult[numthread]+=this->constraints[i].V*deltaS*deltaS;
			}
		}
		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result+=*iter;

		result+=this->Shift;

		return result;
	}

	/** \brief calculate the force for the ith Particle.
		\param[out] result  Force on particle i.
		\param i            Index of the particle. */
	virtual void Force(GeometryVector & result, size_t i) //calculate the force for the ith Particle
	{
		if(this->RhoValid==false)
			this->GetRho_MultiConfig();

		result=GeometryVector(this->Dimension);
		std::vector<GeometryVector> threadresult(this->ParallelNumThread, result);

		signed long end=this->constraints.size();

		size_t inc=i/this->Np;//the number of subconfig that contains particle i
#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long j=0; j<end; j++)
			{
				GeometryVector k=this->constraints[j].k;
				GeometryVector r=this->pConfig->GetCartesianCoordinates(i);
				double s = this->Snow[j];
				double deltaS=(s-this->S0[j]);
				threadresult[numthread].MinusFrom((4.0/this->Np)*this->constraints[j].V*deltaS*(this->RhoImag[j+inc*end]*std::cos(k.Dot(r))-this->RhoReal[j+inc*end]*std::sin(k.Dot(r)))*k);
			}
		}

		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result.AddFrom(*iter);
		result.MultiplyFrom(1.0/this->NSubConfig);
	}

	/** \brief Bind the configuration, validating that the particle count splits evenly into sub-configurations.
		\param Position  Configuration to bind. */
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->S0.size()==this->constraints.size());
		size_t N=Position.NumParticle();
		if(N%NSubConfig!=0)
		{
			std::cerr<<"Error in ShiftedCCPotential_S2_MultiConfig : the total number of particles is not multiple of the number of configurations!\n";
			std::cerr<<"N="<<N<<", NSubConfig="<<NSubConfig<<'\n';
			assert(false);
		}
		this->Np=N/this->NSubConfig;
		this->ShiftedCCPotential::SetConfiguration(Position);
		this->Volume=this->pConfig->PeriodicVolume();
	}
};


/** \brief A derived class to compute  a collective-coordinate potential to a form suitable for Monte Carlo simulations.
 *	This is derived from both ShiftedCCPotential and forMCPotential classes. */
class ShiftedCCPotential_ForMCMixin : public forMCPotential, public ShiftedCCPotential
{
protected:
	std::vector<double> RhoReal_AfterMove, RhoImag_AfterMove; //< The real and imaginary parts of collective coordinates when a trial move is accepted.					

public:
	/** \brief A constructor that converts a ShiftedCCPotential object. */
	ShiftedCCPotential_ForMCMixin(const ShiftedCCPotential & src) :ShiftedCCPotential(src)
	{
	}

	/** \brief Return a heap-allocated copy of this full derived object.
		\return A new ShiftedCCPotential_ForMCMixin copy. */
	virtual ShiftedCCPotential_ForMCMixin * clone()
	{
		return new ShiftedCCPotential_ForMCMixin(*this);
	}

	//parameters : Cartesian coordinates before and after move
	//the potential must have been set to the correct configuration before this function is called
	//return : delta energy
	/** \brief Compute the energy change of moving one particle, caching the trial collective density.
		\param prevCart   Cartesian coordinates of the particle before the move.
		\param afterCart  Cartesian coordinates of the particle after the move.
		\return The change in energy from the trial move. */
	virtual double TryMove(GeometryVector & prevCart, GeometryVector & afterCart)
	{
		this->RhoReal_AfterMove.resize(constraints.size());
		this->RhoImag_AfterMove.resize(constraints.size());

		double PrevE = this->Energy();
		signed long end=this->constraints.size();
#pragma omp parallel for num_threads(this->ParallelNumThread)
		for(signed long i=0; i<end; i++)
		{
			this->RhoReal_AfterMove[i]=this->RhoReal[i]-std::cos(this->constraints[i].k.Dot(prevCart))+std::cos(this->constraints[i].k.Dot(afterCart));
			this->RhoImag_AfterMove[i]=this->RhoImag[i]-std::sin(this->constraints[i].k.Dot(prevCart))+std::sin(this->constraints[i].k.Dot(afterCart));
		}
		std::swap(this->RhoImag, this->RhoImag_AfterMove);
		std::swap(this->RhoReal, this->RhoReal_AfterMove);
		double AfterE = this->Energy();
		std::swap(this->RhoImag, this->RhoImag_AfterMove);
		std::swap(this->RhoReal, this->RhoReal_AfterMove);
		return AfterE - PrevE;
	}

	/** \brief Commit the last trial move by adopting its cached collective density. */
	virtual void AcceptMove(void)
	{
		this->RhoImag=this->RhoImag_AfterMove;
		this->RhoReal=this->RhoReal_AfterMove;
	}
};

//typedef ShiftedCCPotential_ForMCMixin<ShiftedCCPotential> ShiftedCCPotential_ForMC;
typedef ShiftedCCPotential_ForMCMixin ShiftedCCPotential_ForMC;

#endif
