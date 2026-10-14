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
#include <cfloat>
#include <cstdint>
#include <limits>
#include <new>

#include "Potential.h"

// Optional FINUFFT Type 1 backend for ShiftedCCPotential::GetRho (off by default).
// Enable with -DCCPOTENTIAL_USE_FINUFFT, add the FINUFFT include directories, and link libfinufft;
// then select it at run time with ShiftedCCPotential::SetRhoBackend(RhoBackend::FinufftType1).
// Without the macro no FINUFFT header or library is needed and only the direct sum exists.
// The macro must be set identically in every translation unit of one program (ODR).
// See tests/finufft_type1/README.md for compile and run commands.
#ifdef CCPOTENTIAL_USE_FINUFFT
#include <complex>
#include <finufft.h>
#endif

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

	/** \brief Backend used by ShiftedCCPotential::GetRho to evaluate rho(k).
	 *  Scope: only ShiftedCCPotential::GetRho honors the selection, so classes that inherit it
	 *  (LogCCPotential, ShiftedCCPotential_S2, the _varBox variants, ShiftedCCPotential_ForMC) use it.
	 *  ShiftedCCPotential_withWeight and ShiftedCCPotential_S2_MultiConfig override GetRho and always
	 *  use their own direct sums; the inherited setting is ignored there. */
	enum class RhoBackend
	{
		Direct,		//< Direct O(N*Nk) summation (default, always available).
		FinufftType1	//< FINUFFT Type 1 on the reciprocal-lattice mode grid (needs -DCCPOTENTIAL_USE_FINUFFT).
	};

	/** \brief Value of GetLastFinufftStatus() when the last GetRho did not call FINUFFT. */
	static constexpr int FinufftNotCalled=-1;
	/** \brief FINUFFT is skipped when the mode grid exceeds this factor times (constraints + particles). */
	static constexpr double FinufftSparseGridFactor=64.0;

	/** \brief Report whether a backend was compiled into this build.
		\param backend  Backend to query.
		\return true for Direct always, and for FinufftType1 only when built with CCPOTENTIAL_USE_FINUFFT. */
	static bool IsRhoBackendAvailable(RhoBackend backend)
	{
		if(backend==RhoBackend::Direct)
			return true;
#ifdef CCPOTENTIAL_USE_FINUFFT
		return backend==RhoBackend::FinufftType1;
#else
		return false;
#endif
	}

	/** \brief Select the rho(k) backend and invalidate the cached rho.
		\param backend  Backend to use from the next GetRho call on.
		\throw std::invalid_argument if the backend is not compiled in; the object is then left unchanged. */
	void SetRhoBackend(RhoBackend backend)
	{
		if(!IsRhoBackendAvailable(backend))
			throw std::invalid_argument("ShiftedCCPotential::SetRhoBackend: backend not available; FinufftType1 requires compiling with -DCCPOTENTIAL_USE_FINUFFT and linking FINUFFT.");
		this->RhoBackendSelected=backend;
		this->RhoValid=false;
	}

	/** \brief Return the selected rho(k) backend.
		\return The backend chosen with SetRhoBackend (Direct by default). */
	RhoBackend GetRhoBackend(void) const
	{
		return this->RhoBackendSelected;
	}

	/** \brief Set the FINUFFT requested relative tolerance and invalidate the cached rho.
		\param eps  Tolerance passed to FINUFFT; must be positive and finite (default 1e-12).
		\throw std::invalid_argument for an invalid eps; the object is then left unchanged. */
	void SetFinufftTolerance(double eps)
	{
		if(!(eps>0.0) || !std::isfinite(eps))
			throw std::invalid_argument("ShiftedCCPotential::SetFinufftTolerance: eps must be positive and finite.");
		this->FinufftTolerance=eps;
		this->RhoValid=false;
	}

	/** \brief Return the FINUFFT requested tolerance.
		\return The tolerance set with SetFinufftTolerance. */
	double GetFinufftTolerance(void) const
	{
		return this->FinufftTolerance;
	}

	/** \brief Cap the number of FINUFFT output modes (memory guard) and invalidate the cached rho.
		The output array takes 16 bytes per mode and FINUFFT's fine grid roughly 2^d times more.
		\param maxModes  Largest allowed product of grid sizes (2K_d+1); must be >= 1 (default 2^22).
		\throw std::invalid_argument if maxModes < 1; the object is then left unchanged. */
	void SetFinufftMaxModes(std::int64_t maxModes)
	{
		if(maxModes<1)
			throw std::invalid_argument("ShiftedCCPotential::SetFinufftMaxModes: maxModes must be >= 1.");
		this->FinufftMaxModes=maxModes;
		this->RhoValid=false;
	}

	/** \brief Return the FINUFFT mode-count cap.
		\return The cap set with SetFinufftMaxModes. */
	std::int64_t GetFinufftMaxModes(void) const
	{
		return this->FinufftMaxModes;
	}

	/** \brief Set the number of threads FINUFFT uses (opts.nthreads) and invalidate the cached rho.
		The default is 1, independent of ParallelNumThread: only FINUFFT runs single-threaded, so that
		repeated runs can reproduce rho(k) bitwise. Earlier builds passed ParallelNumThread instead.
		A value > 1 is an optional performance choice; FINUFFT's multithreaded spreading may then change
		the summation order between runs, so no repeatability is promised. The setting affects only the
		FINUFFT call; the direct sums, energies, forces and other OpenMP regions still use ParallelNumThread.
		Repeatability means the same machine, build, settings and inputs; it does not mean identity between
		the FINUFFT and Direct backends or across platforms.
		Available in every build; without CCPOTENTIAL_USE_FINUFFT the value is stored but unused.
		\param nthreads  Number of FINUFFT threads; must be >= 1 (default 1).
		\throw std::invalid_argument if nthreads < 1; the object is then left unchanged. */
	void SetFinufftNumThreads(int nthreads)
	{
		if(nthreads<1)
			throw std::invalid_argument("ShiftedCCPotential::SetFinufftNumThreads: nthreads must be >= 1.");
		this->FinufftNumThreads=nthreads;
		this->RhoValid=false;
	}

	/** \brief Return the number of threads FINUFFT uses.
		\return The value set with SetFinufftNumThreads (1 by default). */
	int GetFinufftNumThreads(void) const
	{
		return this->FinufftNumThreads;
	}

	/** \brief Return the backend that actually produced the cached rho in the last ShiftedCCPotential::GetRho.
		Direct means either Direct was selected or the FINUFFT path fell back. Not updated by derived GetRho overrides.
		\return The backend used by the most recent base-class evaluation. */
	RhoBackend GetLastRhoBackendUsed(void) const
	{
		return this->RhoBackendLastUsed;
	}

	/** \brief Return the FINUFFT status of the last ShiftedCCPotential::GetRho.
		\return 0 on FINUFFT success, the nonzero FINUFFT code that caused a direct fallback, or FinufftNotCalled. */
	int GetLastFinufftStatus(void) const
	{
		return this->FinufftLastStatus;
	}

protected:
	RhoBackend RhoBackendSelected;	//< Backend requested by the user.
	RhoBackend RhoBackendLastUsed;	//< Backend that produced the cached rho.
	double FinufftTolerance;	//< FINUFFT requested tolerance.
	std::int64_t FinufftMaxModes;	//< Cap on the FINUFFT mode-grid size.
	int FinufftNumThreads;	//< FINUFFT opts.nthreads (default 1, independent of ParallelNumThread).
	int FinufftLastStatus;	//< Last FINUFFT return code, or FinufftNotCalled.

#ifdef CCPOTENTIAL_USE_FINUFFT
	/** \brief Try to fill RhoReal/RhoImag with one FINUFFT Type 1 transform.
		rho(k) = sum_j exp(+i k.r_j) with unit strengths and no normalization; since k = sum_d n_d b_d with
		n_d = k.a_d/(2 pi), k.r_j = sum_d n_d (2 pi s_jd) for fractional coordinates s_j. The mapping from
		constraints to integer modes is recomputed on every call because constraints is public and mutable.
		Snapping n_d to the nearest integer changes each phase by at most 2 pi dim max_d|n_d-round(n_d)|,
		so a mode is accepted only if the rounding residual |n_d-round(n_d)| is within
		max(budget, roundoff), where budget = min(1e-8, eps)/(8 pi dim) keeps the snapping error below
		min(1e-8, eps)/4 in phase and roundoff = 4*DBL_EPSILON*(1+sum_c|k_c a_dc|/(2 pi)) estimates the error
		of computing n_d itself. A residual at roundoff level means k is a reciprocal-lattice vector to double
		precision; snapping it is then no less accurate than the direct sum, whose phases k.r carry rounding
		of the same order, so large indices are accepted (earlier builds rejected every mode with
		roundoff > budget, i.e. indices above about 43/dim at eps=1e-12, silently falling back for large N).
		Every eps below about 2e-14*dim (4*DBL_EPSILON > budget) is still rejected before FINUFFT is called.
		Returns false, leaving the caller to use the direct sum, when: the dimension is not 1-3; there are
		no constraints or no particles; eps is below that floor; some mode fails the residual test above;
		|n_d| > 2^30; the mode grid exceeds
		FinufftMaxModes or FinufftSparseGridFactor*(Nk+N); allocation fails; or FINUFFT returns nonzero.
		All mapping checks run before the dense mode array is allocated.
		\return true if RhoReal/RhoImag hold the FINUFFT result, false if a direct fallback is required. */
	bool GetRho_FinufftType1(void)
	{
		const DimensionType dim=this->Dimension;
		const size_t nk=this->constraints.size();
		const size_t np=this->pConfig->NumParticle();
		if(dim<1 || dim>3 || nk==0 || np==0)
			return false;
		if(np>static_cast<size_t>(std::numeric_limits<std::int64_t>::max()))
			return false;

		const double twopi=2.0*::pi;
		const double maxIndex=static_cast<double>(std::int64_t(1)<<30);
		const double budget=std::min(1e-8, this->FinufftTolerance)/(4.0*twopi*static_cast<double>(dim));
		if(4.0*DBL_EPSILON>budget)
			return false;	//eps below what double-precision phases resolve
		try
		{
			GeometryVector basis[3];
			for(DimensionType d=0; d<dim; d++)
				basis[d]=this->pConfig->GetBasisVector(d);

			//integer reciprocal-lattice indices of every constraint
			std::vector<std::int64_t> modes(nk*dim);
			std::int64_t K[3]={0, 0, 0};
			for(size_t i=0; i<nk; i++)
			{
				const GeometryVector & k=this->constraints[i].k;
				for(DimensionType d=0; d<dim; d++)
				{
					double dot=0.0, scale=0.0;
					for(DimensionType c=0; c<dim; c++)
					{
						double t=k.x[c]*basis[d].x[c];
						dot+=t;
						scale+=std::abs(t);
					}
					double n=dot/twopi;
					if(!std::isfinite(n) || std::abs(n)>maxIndex)
						return false;
					double nr=std::round(n);
					double roundoff=4.0*DBL_EPSILON*(1.0+scale/twopi);
					if(std::abs(n-nr)>std::max(budget, roundoff))
						return false;
					std::int64_t m=static_cast<std::int64_t>(nr);
					modes[i*dim+d]=m;
					K[d]=std::max(K[d], m<0 ? -m : m);
				}
			}

			//grid sizes, guarded against overflow and unbounded allocation
			std::int64_t ms[3]={1, 1, 1};
			std::int64_t total=1;
			for(DimensionType d=0; d<dim; d++)
			{
				ms[d]=2*K[d]+1;
				if(ms[d]>this->FinufftMaxModes/total)
					return false;
				total*=ms[d];
			}
			if(static_cast<double>(total)>FinufftSparseGridFactor*(static_cast<double>(nk)+static_cast<double>(np)))
				return false;

			//fractional coordinates as angles in [0, 2 pi]
			std::vector<double> angles[3];
			for(DimensionType d=0; d<dim; d++)
				angles[d].resize(np);
			for(size_t j=0; j<np; j++)
			{
				const GeometryVector & s=this->pConfig->GetRelativeCoordinates(j);
				for(DimensionType d=0; d<dim; d++)
					angles[d][j]=twopi*(s.x[d]-std::floor(s.x[d]));
			}
			std::vector<std::complex<double> > strengths(np, std::complex<double>(1.0, 0.0));
			std::vector<std::complex<double> > f(static_cast<size_t>(total));

			finufft_opts opts;
			finufft_default_opts(&opts);
			opts.modeord=0;
			opts.nthreads=this->FinufftNumThreads;

			const std::int64_t nj=static_cast<std::int64_t>(np);
			int status;
			if(dim==1)
				status=finufft1d1(nj, angles[0].data(), strengths.data(), +1, this->FinufftTolerance, ms[0], f.data(), &opts);
			else if(dim==2)
				status=finufft2d1(nj, angles[0].data(), angles[1].data(), strengths.data(), +1, this->FinufftTolerance, ms[0], ms[1], f.data(), &opts);
			else
				status=finufft3d1(nj, angles[0].data(), angles[1].data(), angles[2].data(), strengths.data(), +1, this->FinufftTolerance, ms[0], ms[1], ms[2], f.data(), &opts);
			this->FinufftLastStatus=status;
			if(status!=0)
				return false;

			//modeord=0: mode n_d sits at n_d+K_d, first dimension fastest
			this->RhoReal.resize(nk);
			this->RhoImag.resize(nk);
			for(size_t i=0; i<nk; i++)
			{
				std::int64_t idx=0;
				for(int d=static_cast<int>(dim)-1; d>=0; d--)
					idx=idx*ms[d]+(modes[i*dim+d]+K[d]);
				double tr=f[static_cast<size_t>(idx)].real();
				double ti=f[static_cast<size_t>(idx)].imag();
#ifdef HAVETWOFIXEDPARTICLES
				tr+=2;
#endif
				this->RhoReal[i]=tr;
				this->RhoImag[i]=ti;
			}
		}
		catch(const std::bad_alloc &)
		{
			return false;
		}
		return true;
	}
#endif

public:
	//this class will manage these members, users can also read it.
	std::vector<double> RhoReal, RhoImag;
	bool RhoValid;
	/** \brief Compute and cache the complex collective density rho(k) at every constraint.
		Uses the backend chosen with SetRhoBackend. When FinufftType1 is selected but cannot be used
		(see GetRho_FinufftType1), this falls back to the direct sum, so rho is always valid on return. */
	virtual void GetRho(void)
	{
		assert(this->pConfig!=nullptr);
		if(this->RhoValid)
		{
			return;
		}
		this->Volume=this->pConfig->PeriodicVolume();
		this->FinufftLastStatus=FinufftNotCalled;

#ifdef CCPOTENTIAL_USE_FINUFFT
		if(this->RhoBackendSelected==RhoBackend::FinufftType1 && this->GetRho_FinufftType1())
		{
			this->RhoBackendLastUsed=RhoBackend::FinufftType1;
			this->RhoValid=true;
			return;
		}
#endif

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

		this->RhoBackendLastUsed=RhoBackend::Direct;
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
		this->RhoBackendSelected=RhoBackend::Direct;
		this->RhoBackendLastUsed=RhoBackend::Direct;
		this->FinufftTolerance=1e-12;
		this->FinufftMaxModes=std::int64_t(1)<<22;
		this->FinufftNumThreads=1;
		this->FinufftLastStatus=FinufftNotCalled;
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
