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
#include <algorithm>

#include "Potential.h"

// Optional FINUFFT Type 1 backend for ShiftedCCPotential::GetRho (off by default).
// Enable with -DCCPOTENTIAL_USE_FINUFFT, add the FINUFFT include directories, and link libfinufft;
// then select it at run time with ShiftedCCPotential::SetRhoBackend(RhoBackend::FinufftType1).
// Without the macro no FINUFFT header or library is needed and only the direct sum exists.
// The macro must be set identically in every translation unit of one program (ODR).
// See tests/finufft_type1/README.md for compile and run commands.
#ifdef CCPOTENTIAL_USE_FINUFFT
#include <complex>
#include <cstring>
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
	const Configuration * pViewConfig;	//< Full configuration of an indexed view (SetIndexedConfiguration), else nullptr.
	std::vector<size_t> ViewIndices;	//< Local particle j of the indexed view is particle ViewIndices[j] of *pViewConfig.

	/** \brief Configuration holding the bound particles and the cell.
		\return *pViewConfig while an indexed view is bound, otherwise *pConfig. */
	const Configuration & BoundConfig(void) const
	{
		return this->pViewConfig ? *this->pViewConfig : *this->pConfig;
	}

	/** \brief Number of bound particles.
		\return ViewIndices.size() while an indexed view is bound, otherwise pConfig->NumParticle(). */
	size_t NumBoundParticle(void) const
	{
		return this->pViewConfig ? this->ViewIndices.size() : this->pConfig->NumParticle();
	}

	/** \brief Index in BoundConfig() of a local particle.
		\param j  Local particle index.
		\return ViewIndices[j] while an indexed view is bound, otherwise j. */
	size_t BoundIndex(size_t j) const
	{
		return this->pViewConfig ? this->ViewIndices[j] : j;
	}
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
	 *  use their own direct sums; the inherited setting is ignored there.
	 *  FinufftType1 also makes ShiftedCCPotential::AllForce use a FINUFFT Type 2 transform (see
	 *  StartForces), but only on objects whose dynamic type is exactly ShiftedCCPotential and only
	 *  when a cost gate predicts a clear win over the direct loop (see FinufftForceGate). */
	enum class RhoBackend
	{
		Direct,		//< Direct O(N*Nk) summation (default, always available).
		FinufftType1	//< FINUFFT Type 1 on the reciprocal-lattice mode grid (needs -DCCPOTENTIAL_USE_FINUFFT).
	};

	/** \brief How ShiftedCCPotential::AllForce computed the collective-coordinate forces. */
	enum class ForceBackend
	{
		Direct,		//< Per-particle direct O(Nk) sums (Force).
		FinufftType2	//< One FINUFFT Type 2 execution with ntransf = dimension.
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

	/** \brief Return how the last ShiftedCCPotential::AllForce computed the forces.
		FinufftType2 requires the cost gate to accept (UseFinufftForces; for MultiStealthyPotential species the
		joint FinufftForceGate); Direct means the per-particle direct sum (selected, declined by the cost gate, or
		a fallback). Force(i) is always direct and does not update this value. In MultiStealthyPotential it
		describes each species' share of the last AllForce.
		\return The force backend used by the most recent base-class AllForce (Direct before any call). */
	ForceBackend GetLastForceBackendUsed(void) const
	{
		return this->ForceBackendLastUsed;
	}

	/** \brief Report whether ShiftedCCPotential::GetRho will run FINUFFT Type 1 for the bound configuration.
		Reads the eligibility of the persistent FINUFFT cache when it is current, and otherwise runs only the mode
		mapping of the FINUFFT path (MapFinufftModes), without a plan or transform, so the answer depends only
		on the constraints, the cell, the particle number and the FINUFFT settings (and a cached makeplan failure).
		A nonzero FINUFFT status at run time can still make GetRho fall back to the direct sum.
		Like the backend selection, it is not meaningful for derived classes that override GetRho.
		\return true if FinufftType1 is selected and compiled in, a configuration is bound, and every check passes. */
	bool FinufftEligible(void) const
	{
#ifdef CCPOTENTIAL_USE_FINUFFT
		return this->FinufftGridModes()>0;
#else
		return false;
#endif
	}

	/** \brief Report whether ShiftedCCPotential::AllForce takes the FINUFFT Type 2 path for the bound configuration.
		FinufftForceGate for this object alone with a team of ParallelNumThread: true only when built with
		CCPOTENTIAL_USE_FINUFFT, the object has exact type ShiftedCCPotential (derived classes may override Force
		or GetRho), FinufftEligible() holds, and the cost model predicts a clear win over the direct loop. The answer
		depends only on the configuration size, constraints, cell and settings, never on measured time. AllForce
		still falls back if an allocation or FINUFFT fails, or if rho fell back to the direct sum.
		MultiStealthyPotential decides its species jointly (FinufftForceGate over the concurrent species).
		\return true if AllForce will attempt the Type 2 transform. */
	bool UseFinufftForces(void) const
	{
		return FinufftForceGate({this}, this->ParallelNumThread)[0]!=0;
	}

	/** \brief Cost gate: choose which of the objects evaluated in one parallel region compute their forces with Type 2
		instead of the direct loop, and in how many particle blocks. An accepted object's Type 2 work runs as blocks, one
		task per block on one thread of a team of `threads`; the forces of the others are computed by the direct loop
		over the whole team (the model of MultiStealthyPotential::AllForce; a single object is
		ShiftedCCPotential::AllForce). With F the FINUFFT-eligible objects of exact type ShiftedCCPotential, each gets
		blocks = min(FinufftBlocks(2), max(1, floor(threads/F))): the input-only rule caps the repeated per-block work,
		and the share of the team keeps a busy team from paying for blocks it cannot run side by side. Every block
		interpolates the same fine grid point by point, so the block count never changes a force. Per object the cost
		model (FinufftForceCosts) predicts direct, its share of the direct loop over the team, and block, the time of
		its largest block on one thread, with work = blocks * block. With K the summed blocks of the F objects,
		candidates are those with direct >= FinufftForceMargin * work / min(K, threads). Taking candidates in order of
		increasing work, the longest prefix A is accepted whose summed savings are at least FinufftForceMargin times the
		wall time A adds to the region, estimated as min(sum_A work, ceil(K_A/threads) max_A block) with K_A the blocks
		of A: with K_A <= threads every block has its own thread, so A adds at most its longest block. For one unsplit
		object this is direct >= FinufftForceMargin * type2. Never measures time, so runs with the same inputs, settings
		and team size decide alike and stay bitwise repeatable.
		\param pots     Objects evaluated concurrently.
		\param threads  Team size of the region and of the direct loop.
		\return Per object, the number of Type 2 blocks if its forces should use Type 2, else 0 (all 0 without FINUFFT);
			the argument of PrepareFinufftForces. */
	static std::vector<size_t> FinufftForceGate(const std::vector<const ShiftedCCPotential *> & pots, unsigned long threads)
	{
		std::vector<size_t> accept(pots.size(), 0);
#ifdef CCPOTENTIAL_USE_FINUFFT
		const double team=static_cast<double>(std::max<unsigned long>(threads, 1));
		std::vector<double> direct(pots.size(), 0.0), block(pots.size(), 0.0), work(pots.size(), 0.0);
		std::vector<std::int64_t> grid(pots.size(), 0);
		std::vector<size_t> blocks(pots.size(), 0), order;
		for(size_t s=0; s<pots.size(); s++)
			grid[s]=typeid(*pots[s])==typeid(ShiftedCCPotential) ? pots[s]->FinufftGridModes() : 0;
		const double eligible=static_cast<double>(pots.size()-std::count(grid.begin(), grid.end(), 0));
		const size_t share=static_cast<size_t>(std::max(1.0, std::floor(team/std::max(eligible, 1.0))));
		double tasks=0.0;
		for(size_t s=0; s<pots.size(); s++)
		{
			if(grid[s]==0)
				continue;
			pots[s]->FinufftForceCosts(static_cast<double>(grid[s]), team, share, direct[s], block[s], blocks[s]);
			work[s]=static_cast<double>(blocks[s])*block[s];
			tasks+=static_cast<double>(blocks[s]);
		}
		for(size_t s=0; s<pots.size(); s++)
			if(blocks[s]>0 && direct[s]>=FinufftForceMargin*work[s]/std::min(tasks, team))
				order.push_back(s);
		std::stable_sort(order.begin(), order.end(), [&work](size_t a, size_t b) { return work[a]<work[b]; });
		double saved=0.0, added=0.0, taken=0.0, longest=0.0;
		size_t take=0;
		for(size_t j=0; j<order.size(); j++)
		{
			saved+=direct[order[j]];
			added+=work[order[j]];
			taken+=static_cast<double>(blocks[order[j]]);
			longest=std::max(longest, block[order[j]]);
			if(saved>=FinufftForceMargin*std::min(added, std::ceil(taken/team)*longest))
				take=j+1;
		}
		for(size_t j=0; j<take; j++)
			accept[order[j]]=blocks[order[j]];
#endif
		return accept;
	}

	/** \brief Refresh the persistent FINUFFT state serially and report whether GetRho will run FINUFFT.
		Rebuilds the mode map and the Type 1 plans if an input changed, exactly as GetRho would, so that a
		following GetRho only executes the plans. Call it before running several objects concurrently: plans are
		then never made inside a parallel region. The cached rho is not touched.
		\return FinufftEligible() after the refresh: true if GetRho will run FINUFFT unless it fails at run time. */
	bool PrepareFinufft(void)
	{
#ifdef CCPOTENTIAL_USE_FINUFFT
		try
		{
			size_t np=0;
			return this->RhoBackendSelected==RhoBackend::FinufftType1 && this->RefreshFinufftCache(np) && this->FinufftCache.eligible;
		}
		catch(const std::bad_alloc &)
		{
			return false;
		}
#else
		return false;
#endif
	}

	/** \brief Set up the force path of the next TryFinufftAllForce (or StartForces), serially, after PrepareFinufft:
		reset GetLastForceBackendUsed to Direct, then, for blocks > 0, keep or make one Type 2 plan per particle block
		(surplus plans are destroyed; a failed makeplan destroys them all and is not retried until an input changes),
		or destroy them for blocks = 0 (then a later request makes them again). blocks is the decision of
		FinufftForceGate (true, i.e. one block, is also accepted); it has no effect unless FINUFFT is eligible. The
		block count changes the work split, never a force.
		\param blocks  Number of Type 2 blocks (at most the particle count is used), or 0 for the direct sum.
		\return The number of Type 2 blocks TryFinufftAllForce will run (0: the direct sum), unless a block fails at
			run time or rho falls back. */
	size_t PrepareFinufftForces(size_t blocks)
	{
		this->ForceBackendLastUsed=ForceBackend::Direct;
#ifdef CCPOTENTIAL_USE_FINUFFT
		FinufftPlanCache & c=this->FinufftCache;
		if(blocks==0 || !c.built || !c.eligible)
		{
			if(!c.plan2.empty())
				c.status2=FinufftNotCalled;	//a deliberate destroy, not a failure: allow the plans again
			FinufftPlanCache::Keep(c.plan2, 0);
			return 0;
		}
		blocks=std::min(blocks, c.np);
		FinufftPlanCache::Keep(c.plan2, blocks);
		if(c.status2==FinufftNotCalled)
		{
			const int status=this->MakeFinufftPlans(2, blocks, c.plan2);
			if(status!=0)
				c.status2=status;	//cached until an input changes
		}
		return c.plan2.size();
#else
		return 0;
#endif
	}

	/** \brief Serial first step of a collective-density evaluation in Type 1 particle blocks. GetRho is StartRho(false),
		then RunRhoBlock(b) for every b below the returned count, then FinishRho; between the two serial steps the blocks
		of one object may run concurrently, one block per thread (MultiStealthyPotential), since each block has its own
		plan, particle range and output modes. Does nothing if rho is valid, except that with forces (objects of exact
		type ShiftedCCPotential only) a valid FINUFFT rho whose angles the cache no longer holds (a copied or assigned
		object, whose cache starts empty, or a rebuilt cache) is recomputed while Type 2 plans exist, so that Type 2
		sees the configuration of the cached rho; Type 1 is deterministic, so rho keeps its bits. Otherwise refreshes
		the persistent cache as PrepareFinufft does (which may make the Type 1 plans unless PrepareFinufft ran first).
		\param forces  The rho will feed Type 2 forces (StartForces).
		\return The number of Type 1 blocks to run; 0 if rho is valid or FINUFFT is not used (FinishRho then uses the
			direct sum). */
	size_t StartRho(bool forces)
	{
#ifdef CCPOTENTIAL_USE_FINUFFT
		FinufftPlanCache & c=this->FinufftCache;
		c.run1.clear();
		if(forces && this->RhoValid && this->RhoLacksAngles())
			this->RhoValid=false;
#endif
		if(this->RhoValid)
			return 0;
		this->Volume=this->BoundConfig().PeriodicVolume();
		this->FinufftLastStatus=FinufftNotCalled;
#ifdef CCPOTENTIAL_USE_FINUFFT
		try
		{
			size_t np=0;
			if(this->RhoBackendSelected==RhoBackend::FinufftType1 && this->RefreshFinufftCache(np))
			{
				if(!c.eligible)
					this->FinufftLastStatus=c.status;
				else
				{
					c.anglesOfRho=false;
					c.run1.assign(c.plan.size(), FinufftNotCalled);
					return c.plan.size();
				}
			}
		}
		catch(const std::bad_alloc &)
		{
			c.run1.clear();
		}
#endif
		return 0;
	}

	/** \brief Execute Type 1 particle block b of the evaluation begun by StartRho: fill the angles 2 pi s_jd of the
		block's particles (a contiguous range of the local order, FinufftPlanCache::BlockStart), pass them to the
		block's plan and execute it into the block's own mode array. Records a nonzero status (FinufftNotCalled for a
		non-finite fractional coordinate) that makes FinishRho fall back to the direct sum. Never makes a plan or
		resizes a buffer, so distinct blocks of one object may run concurrently.
		\param b  Block index, below the count StartRho returned. */
	void RunRhoBlock(size_t b)
	{
#ifdef CCPOTENTIAL_USE_FINUFFT
		FinufftPlanCache & c=this->FinufftCache;
		const DimensionType dim=this->Dimension;
		const size_t blocks=c.plan.size(), lo=c.BlockStart(b, blocks), hi=c.BlockStart(b+1, blocks);
		const size_t total=static_cast<size_t>(c.ms[0]*c.ms[1]*c.ms[2]);
		const double twopi=2.0*::pi;
		for(size_t j=lo; j<hi; j++)
		{
			const GeometryVector & s=this->BoundConfig().GetRelativeCoordinates(this->BoundIndex(j));
			for(DimensionType d=0; d<dim; d++)
			{
				if(!std::isfinite(s.x[d]))
					return;	//c.run1[b] stays FinufftNotCalled
				c.angles[d][j]=twopi*(s.x[d]-std::floor(s.x[d]));
			}
		}
		try
		{
			int status=finufft_setpts(c.plan[b], static_cast<std::int64_t>(hi-lo), c.Points(0, lo), c.Points(1, lo), c.Points(2, lo), 0, nullptr, nullptr, nullptr);
			if(status==0)
				status=finufft_execute(c.plan[b], c.strengths.data()+lo, c.f.data()+b*total);
			c.run1[b]=status;
		}
		catch(const std::bad_alloc &)
		{
		}
#endif
	}

	/** \brief Serial last step of a collective-density evaluation begun by StartRho. If every Type 1 block succeeded,
		rho(k) = sum_j exp(+i k.r_j) at each constraint is the sum of the blocks' modes in block order (unit strengths,
		no normalization; one block is the single transform of earlier builds, bitwise); otherwise, or without FINUFFT,
		the direct sum. GetLastFinufftStatus then reports the first nonzero block status in block order, or 0. Does
		nothing if rho is valid. */
	void FinishRho(void)
	{
		if(this->RhoValid)
			return;
#ifdef CCPOTENTIAL_USE_FINUFFT
		FinufftPlanCache & c=this->FinufftCache;
		if(!c.run1.empty())
		{
			int status=0;
			for(int s : c.run1)
				if(s!=0)
				{
					status=s;
					break;
				}
			c.run1.clear();
			this->FinufftLastStatus=status;
			if(status==0)
			{
				try
				{
					const size_t nk=this->constraints.size(), blocks=c.plan.size();
					const size_t total=static_cast<size_t>(c.ms[0]*c.ms[1]*c.ms[2]);
					this->RhoReal.resize(nk);
					this->RhoImag.resize(nk);
					for(size_t i=0; i<nk; i++)
					{
						const size_t slot=static_cast<size_t>(c.slot[i]);
						std::complex<double> rho=c.f[slot];
						for(size_t b=1; b<blocks; b++)
							rho+=c.f[b*total+slot];
						double tr=rho.real();
#ifdef HAVETWOFIXEDPARTICLES
						tr+=2;
#endif
						this->RhoReal[i]=tr;
						this->RhoImag[i]=rho.imag();
					}
					c.anglesOfRho=true;
					this->RhoBackendLastUsed=RhoBackend::FinufftType1;
					this->RhoValid=true;
					return;
				}
				catch(const std::bad_alloc &)
				{
				}
			}
		}
#endif
		this->GetRho_Direct();
	}

	/** \brief Serial first step of the Type 2 forces in particle blocks: TryFinufftAllForce is StartForces, then
		RunForceBlock(b, results) for every b below the returned count, then FinishForces; the blocks may run
		concurrently, one block per thread. Requires the Type 2 plans (PrepareFinufftForces) and a valid FINUFFT rho
		with its angles (StartRho(true) ... FinishRho). Builds the Type 2 input modes from rho (see RunForceBlock) and
		sizes results.
		\param[out] results  Resized to the particle count; RunForceBlock fills its blocks.
		\return The number of Type 2 blocks to run, 0 if Type 2 cannot be used (the caller then uses the direct sum). */
	size_t StartForces(std::vector<GeometryVector> & results)
	{
#ifdef CCPOTENTIAL_USE_FINUFFT
		FinufftPlanCache & c=this->FinufftCache;
		c.run2.clear();
		if(c.plan2.empty() || !this->RhoValid || this->RhoBackendLastUsed!=RhoBackend::FinufftType1 || !c.anglesOfRho)
			return 0;
		try
		{
			const DimensionType dim=this->Dimension;
			const size_t nk=this->constraints.size();
			const size_t total=static_cast<size_t>(c.ms[0]*c.ms[1]*c.ms[2]);
			c.coef.assign(dim*total, std::complex<double>(0.0, 0.0));
			for(size_t i=0; i<nk; i++)
			{
				const std::complex<double> w=this->constraints[i].V*std::complex<double>(this->RhoReal[i], -this->RhoImag[i]);
				for(DimensionType a=0; a<dim; a++)
					c.coef[a*total+static_cast<size_t>(c.slot[i])]+=this->constraints[i].k.x[a]*w;
			}
			c.out.resize(dim*c.np);
			results.resize(c.np);
			c.run2.assign(c.plan2.size(), FinufftNotCalled);
			return c.plan2.size();
		}
		catch(const std::bad_alloc &)
		{
			c.run2.clear();
		}
#endif
		return 0;
	}

	/** \brief Execute Type 2 particle block b of the forces begun by StartForces and write the block's forces.
		Force() sums F_j = -(2/V) sum_k V_k (Im rho_k cos(k.r_j) - Re rho_k sin(k.r_j)) k, which equals
		F_{j,a} = (2/V) Im[ sum_k V_k k_a conj(rho_k) exp(+i k.r_j) ]: one Type 2 transform per Cartesian component a
		(ntransf = dimension) with coefficient V_k k_a conj(rho_k) at the mode of constraint k (constraints sharing a mode
		add, in constraint order), isign=+1; every constraint enters once, as in Force(). The block passes the angles of
		its particles (filled by the Type 1 blocks of the cached rho) to its own plan; every block interpolates the same
		fine grid, point by point, so the forces do not depend on the block count. Never makes a plan or resizes a buffer.
		\param b             Block index, below the count StartForces returned.
		\param[out] results  Receives the forces of the block's particles if its execution succeeds. */
	void RunForceBlock(size_t b, std::vector<GeometryVector> & results)
	{
#ifdef CCPOTENTIAL_USE_FINUFFT
		FinufftPlanCache & c=this->FinufftCache;
		const DimensionType dim=this->Dimension;
		const size_t blocks=c.plan2.size(), lo=c.BlockStart(b, blocks), n=c.BlockStart(b+1, blocks)-lo;
		std::complex<double> * out=c.out.data()+dim*lo;
		try
		{
			int status=finufft_setpts(c.plan2[b], static_cast<std::int64_t>(n), c.Points(0, lo), c.Points(1, lo), c.Points(2, lo), 0, nullptr, nullptr, nullptr);
			if(status==0)
				status=finufft_execute(c.plan2[b], out, c.coef.data());
			c.run2[b]=status;
			if(status!=0)
				return;
		}
		catch(const std::bad_alloc &)
		{
			return;
		}
		const double scale=2.0/this->Volume;
		for(size_t j=0; j<n; j++)
		{
			results[lo+j]=GeometryVector(dim);
			for(DimensionType a=0; a<dim; a++)
				results[lo+j].x[a]=scale*out[a*n+j].imag();
		}
#endif
	}

	/** \brief Serial last step of the Type 2 forces begun by StartForces.
		\return true, recording FinufftType2 for GetLastForceBackendUsed, if every block succeeded; false (results then
			unspecified) if StartForces returned 0 or a block failed. */
	bool FinishForces(void)
	{
#ifdef CCPOTENTIAL_USE_FINUFFT
		FinufftPlanCache & c=this->FinufftCache;
		const bool ok=!c.run2.empty() && std::all_of(c.run2.begin(), c.run2.end(), [](int s) { return s==0; });
		c.run2.clear();
		if(ok)
			this->ForceBackendLastUsed=ForceBackend::FinufftType2;
		return ok;
#else
		return false;
#endif
	}

	/** \brief Predicted time of one particle block of the current Type 1 or Type 2 plans on one thread (cost model of
		FinufftBlocks, the largest block); used only to order tasks, never to choose a result.
		\param type  1 for the rho blocks, 2 for the force blocks.
		\return The time [us], or 0 without plans of that type. */
	double FinufftBlockTime(int type) const
	{
#ifdef CCPOTENTIAL_USE_FINUFFT
		const FinufftPlanCache & c=this->FinufftCache;
		const size_t blocks=type==1 ? c.plan.size() : c.plan2.size();
		if(blocks==0)
			return 0.0;
		double fixed=0.0, perParticle=0.0;
		this->FinufftExecutionCosts(type, static_cast<double>(c.ms[0]*c.ms[1]*c.ms[2]), fixed, perParticle);
		return fixed+perParticle*std::ceil(static_cast<double>(c.np)/static_cast<double>(blocks));
#else
		return 0.0;
#endif
	}

protected:
	RhoBackend RhoBackendSelected;	//< Backend requested by the user.
	RhoBackend RhoBackendLastUsed;	//< Backend that produced the cached rho.
	ForceBackend ForceBackendLastUsed;	//< Backend that produced the last AllForce result.
	double FinufftTolerance;	//< FINUFFT requested tolerance.
	std::int64_t FinufftMaxModes;	//< Cap on the FINUFFT mode-grid size.
	int FinufftNumThreads;	//< FINUFFT opts.nthreads (default 1, independent of ParallelNumThread).
	int FinufftLastStatus;	//< Last FINUFFT return code, or FinufftNotCalled.

#ifdef CCPOTENTIAL_USE_FINUFFT
	/** \brief Persistent FINUFFT state of one potential: the constraint-to-mode map, the guru plans and work buffers.
		RefreshFinufftCache rebuilds it lazily whenever an input it depends on changes (see FinufftCacheIsCurrent);
		an ineligible result is cached too, so a direct fallback costs only that check until an input changes.
		Each transform runs as particle blocks (contiguous ranges of the local particle order, BlockStart; the counts
		come from FinufftBlocks and never from the thread count), one guru plan per block. The Type 1 plans (rho) are
		made with the map; the Type 2 plans (forces, ntransf = dim) only by PrepareFinufftForces when the force gate
		accepts. Both read the same angles, which only the Type 1 blocks (RunRhoBlock) fill, so Type 2 always sees the
		configuration of the cached rho. A copy or assignment starts empty and is rebuilt on first use, so a plan
		handle is never shared or destroyed twice; the destructor destroys every plan. Invariant: every plan was made
		with the stored dim, ms, tol and nthreads, and plan is empty or holds one plan per Type 1 block. One object must
		not run two evaluations concurrently; only distinct blocks of one evaluation may run concurrently, between
		its serial start and finish steps (StartRho ... FinishRho, StartForces ... FinishForces). */
	struct FinufftPlanCache
	{
		bool built=false;	//< the members below describe the current inputs
		DimensionType dim=0;	//< fingerprint: Dimension
		size_t np=0;	//< fingerprint: particle count
		double tol=0.0;	//< fingerprint: FinufftTolerance
		std::int64_t maxModes=0;	//< fingerprint: FinufftMaxModes
		int nthreads=0;	//< fingerprint: FinufftNumThreads
		std::vector<double> key;	//< fingerprint: basis vectors, then constraint wavevectors (dim doubles each)
		bool eligible=false;	//< false: the direct sum is used until an input changes
		int status=FinufftNotCalled;	//< Type 1 finufft_makeplan code if that made the inputs ineligible
		int status2=FinufftNotCalled;	//< Type 2 finufft_makeplan code if it failed (not retried until a rebuild)
		bool anglesOfRho=false;	//< angles hold the configuration of the cached rho (set only by FinishRho)
		std::int64_t ms[3]={1, 1, 1};	//< mode-grid sizes 2K_d+1 (1 for unused dimensions)
		std::vector<std::int64_t> slot;	//< position of each constraint's mode in one block's modes
		std::vector<finufft_plan> plan;	//< owned guru Type 1 plans (rho), one per Type 1 block
		std::vector<finufft_plan> plan2;	//< owned guru Type 2 plans (forces), one per Type 2 block, while the gate accepts
		std::vector<int> run1, run2;	//< status of each block of the running evaluation (empty outside one)
		std::vector<double> angles[3];	//< nonuniform points 2 pi s_jd of the cached rho's configuration
		std::vector<std::complex<double> > strengths, f;	//< Type 1 unit strengths and output modes (one grid per block)
		std::vector<std::complex<double> > coef, out;	//< Type 2 input modes (dim grids) and output (dim x block size, per block)

		/** \brief Construct an empty cache. */
		FinufftPlanCache()=default;

		/** \brief Copy-construct as an empty cache; the source's plans are never shared. */
		FinufftPlanCache(const FinufftPlanCache &)
		{
		}

		/** \brief Assign as an empty cache: destroy this object's plans and never share the source's.
			\return *this. */
		FinufftPlanCache & operator=(const FinufftPlanCache &)
		{
			Keep(this->plan, 0);
			Keep(this->plan2, 0);
			this->run1.clear();
			this->run2.clear();
			this->built=false;
			this->anglesOfRho=false;
			return *this;
		}

		/** \brief Destroy the owned plans. */
		~FinufftPlanCache()
		{
			Keep(this->plan, 0);
			Keep(this->plan2, 0);
		}

		/** \brief Destroy owned plans from the back until at most count remain.
			\param plans  plan or plan2 of a cache.
			\param count  Plans to keep. */
		static void Keep(std::vector<finufft_plan> & plans, size_t count)
		{
			while(plans.size()>count)
			{
				if(plans.back()!=nullptr)
					finufft_destroy(plans.back());
				plans.pop_back();
			}
		}

		/** \brief First local particle of a block: block b of B covers [BlockStart(b, B), BlockStart(b+1, B)).
			\param b       Block index (b = B gives np).
			\param blocks  Block count B.
			\return b np / B. */
		size_t BlockStart(size_t b, size_t blocks) const
		{
			return b*this->np/blocks;
		}

		/** \brief Angles of one coordinate from a local particle on, for finufft_setpts.
			\param d   Coordinate (0-2).
			\param lo  First local particle.
			\return Pointer to angles[d][lo], or nullptr for a coordinate beyond the dimension. */
		double * Points(DimensionType d, size_t lo)
		{
			return this->angles[d].empty() ? nullptr : this->angles[d].data()+lo;
		}
	};
	FinufftPlanCache FinufftCache;	//< Persistent FINUFFT state; copies start empty (see FinufftPlanCache).

	/** \brief Read the inputs the FINUFFT path needs from the bound configuration.
		\param[out] basis  Receives the first Dimension cell basis vectors.
		\param[out] np     Receives the bound particle count.
		\return false, leaving the outputs unspecified, when nothing is bound, the dimension is not 1-3, or there
			are no constraints or no particles (FINUFFT is then not used). */
	bool FinufftInputs(GeometryVector basis[3], size_t & np) const
	{
		const DimensionType dim=this->Dimension;
		if((this->pConfig==nullptr && this->pViewConfig==nullptr) || dim<1 || dim>3 || this->constraints.empty())
			return false;
		np=this->NumBoundParticle();
		if(np==0 || np>static_cast<size_t>(std::numeric_limits<std::int64_t>::max()))
			return false;
		for(DimensionType d=0; d<dim; d++)
			basis[d]=this->BoundConfig().GetBasisVector(d);
		return true;
	}

	/** \brief Mode-grid size of the FINUFFT transforms for the bound configuration, without changing any state.
		Uses the persistent cache when it is current (including a cached makeplan failure) and otherwise runs
		only MapFinufftModes.
		\return The grid size G = prod_d (2K_d+1) if FinufftType1 is selected and FINUFFT is eligible, else 0. */
	std::int64_t FinufftGridModes(void) const
	{
		GeometryVector basis[3];
		size_t np=0;
		if(this->RhoBackendSelected!=RhoBackend::FinufftType1 || !this->FinufftInputs(basis, np))
			return 0;
		const FinufftPlanCache & c=this->FinufftCache;
		if(this->FinufftCacheIsCurrent(basis, np))
			return c.eligible ? c.ms[0]*c.ms[1]*c.ms[2] : 0;
		try
		{
			std::int64_t ms[3]={1, 1, 1};
			std::vector<std::int64_t> slot;
			return this->MapFinufftModes(basis, np, ms, slot) ? ms[0]*ms[1]*ms[2] : 0;
		}
		catch(const std::bad_alloc &)
		{
			return 0;
		}
	}

	/** \brief Report whether FinufftCache was built from the current inputs.
		Compares bit for bit everything the mode map, its eligibility and the plan depend on: Dimension, the
		particle count, FinufftTolerance, FinufftMaxModes, FinufftNumThreads, the cell basis vectors (SetConfiguration
		may bind a different cell) and every constraint wavevector (constraints is public and may be edited in
		place between calls, so it is compared on every call rather than trusted).
		\param basis  Basis vectors of the bound configuration (first Dimension entries are used).
		\param np     Particle count of the bound configuration.
		\return true if the cached map, eligibility and plan can be reused. */
	bool FinufftCacheIsCurrent(const GeometryVector * basis, size_t np) const
	{
		const FinufftPlanCache & c=this->FinufftCache;
		const DimensionType dim=this->Dimension;
		const size_t nk=this->constraints.size();
		const size_t bytes=dim*sizeof(double);
		if(!c.built || c.dim!=dim || c.np!=np || c.tol!=this->FinufftTolerance || c.maxModes!=this->FinufftMaxModes
			|| c.nthreads!=this->FinufftNumThreads || c.key.size()!=(dim+nk)*dim)
			return false;
		for(DimensionType d=0; d<dim; d++)
			if(std::memcmp(&c.key[d*dim], basis[d].x, bytes)!=0)
				return false;
		for(size_t i=0; i<nk; i++)
			if(std::memcmp(&c.key[(dim+i)*dim], this->constraints[i].k.x, bytes)!=0)
				return false;
		return true;
	}

	/** \brief Make guru plans on the cached mode grid (FinufftCache.ms) with the fixed options of both transforms:
		finufft_default_opts, then modeord=0, nthreads=FinufftNumThreads, upsampfac=2.0 and spread_sort=1, iflag +1
		and eps=FinufftTolerance. upsampfac and spread_sort are fixed instead of left to FINUFFT's density
		heuristics, so a transform never depends on point density. At the default eps=1e-12 they are also what
		those heuristics pick (FINUFFT 2.5.1; 1D skips sorting only beyond 1000 points per fine-grid cell), so there
		rho is bitwise identical to the simple interface finufft{1,2,3}d1 with default options, which earlier
		builds called.
		Appends plans until there are count; a failed makeplan destroys all of them, so plans is all or nothing.
		\param type         1 for rho (one transform), 2 for forces (Dimension transforms, one per component).
		\param count        Number of plans wanted.
		\param[in,out] plans  plan or plan2 of the cache (plans of this grid); count plans, or none on failure.
		\return 0, or the first nonzero finufft_makeplan status (a warning code may still return a plan; it is destroyed). */
	int MakeFinufftPlans(int type, size_t count, std::vector<finufft_plan> & plans)
	{
		finufft_opts opts;
		finufft_default_opts(&opts);
		opts.modeord=0;
		opts.nthreads=this->FinufftNumThreads;
		opts.upsampfac=2.0;
		opts.spread_sort=1;
		const int ntransf=type==1 ? 1 : static_cast<int>(this->Dimension);
		plans.reserve(count);	//push_back below cannot throw, so no made plan is lost
		while(plans.size()<count)
		{
			finufft_plan handle=nullptr;
			const int status=finufft_makeplan(type, static_cast<int>(this->Dimension), this->FinufftCache.ms, +1, ntransf, this->FinufftTolerance, &handle, &opts);
			if(handle!=nullptr)
				plans.push_back(handle);
			if(status!=0)
			{
				FinufftPlanCache::Keep(plans, 0);
				return status;
			}
		}
		return 0;
	}

	/** \brief Rebuild FinufftCache for the current inputs: mode map, eligibility, work buffers and Type 1 plans.
		All plans are kept when the grid sizes, tolerance and thread count are unchanged (for example after a
		particle-count change, or a rescaled cell whose constraints keep their integer modes) and are otherwise
		destroyed; then the Type 1 plans are trimmed or completed to FinufftBlocks(1) (the Type 2 plans only by
		PrepareFinufftForces). A nonzero Type 1 makeplan status makes the inputs ineligible and is kept in
		the cache; a cached Type 2 makeplan failure is cleared, so the next PrepareFinufftForces retries once.
		\param basis  Basis vectors of the bound configuration.
		\param np     Particle count of the bound configuration.
		\throw std::bad_alloc if a buffer cannot be allocated; the cache is then left unbuilt. */
	void RebuildFinufftCache(const GeometryVector * basis, size_t np)
	{
		FinufftPlanCache & c=this->FinufftCache;
		const DimensionType dim=this->Dimension;
		const size_t nk=this->constraints.size();
		c.built=false;
		std::int64_t ms[3]={1, 1, 1};
		const bool eligible=this->MapFinufftModes(basis, np, ms, c.slot);

		//keep the plans only if they were made for the same grid, tolerance and threads (preserves the invariant)
		if(!eligible || c.dim!=dim || c.tol!=this->FinufftTolerance || c.nthreads!=this->FinufftNumThreads
			|| c.ms[0]!=ms[0] || c.ms[1]!=ms[1] || c.ms[2]!=ms[2])
		{
			FinufftPlanCache::Keep(c.plan, 0);
			FinufftPlanCache::Keep(c.plan2, 0);
		}
		c.dim=dim;
		c.np=np;
		c.tol=this->FinufftTolerance;
		c.maxModes=this->FinufftMaxModes;
		c.nthreads=this->FinufftNumThreads;
		for(int d=0; d<3; d++)
			c.ms[d]=ms[d];
		c.key.resize((dim+nk)*dim);
		for(DimensionType d=0; d<dim; d++)
			std::memcpy(&c.key[d*dim], basis[d].x, dim*sizeof(double));
		for(size_t i=0; i<nk; i++)
			std::memcpy(&c.key[(dim+i)*dim], this->constraints[i].k.x, dim*sizeof(double));
		c.eligible=false;
		c.status=FinufftNotCalled;
		c.status2=FinufftNotCalled;
		c.anglesOfRho=false;

		if(eligible)
		{
			for(DimensionType d=0; d<3; d++)
				c.angles[d].resize(d<dim ? np : 0);
			c.strengths.assign(np, std::complex<double>(1.0, 0.0));
			const size_t total=static_cast<size_t>(ms[0]*ms[1]*ms[2]);
			const size_t blocks=this->FinufftBlocks(1, np, static_cast<double>(total));
			c.f.resize(blocks*total);
			FinufftPlanCache::Keep(c.plan, blocks);
			c.status=this->MakeFinufftPlans(1, blocks, c.plan);
			if(c.status!=0)
				FinufftPlanCache::Keep(c.plan2, 0);	//ineligible: no plan of these inputs stays
			else
				c.status=FinufftNotCalled;
			c.eligible=!c.plan.empty();
		}
		c.built=true;
	}

	/** \brief Map every constraint to its position in the FINUFFT mode grid and check eligibility.
		Since k = sum_d n_d b_d with n_d = k.a_d/(2 pi), k.r_j = sum_d n_d (2 pi s_jd) for fractional coordinates s_j.
		Snapping n_d to the nearest integer changes each phase by at most 2 pi dim max_d|n_d-round(n_d)|,
		so a mode is accepted only if the rounding residual |n_d-round(n_d)| is within
		max(budget, roundoff), where budget = min(1e-8, eps)/(8 pi dim) keeps the snapping error below
		min(1e-8, eps)/4 in phase and roundoff = 4*DBL_EPSILON*(1+sum_c|k_c a_dc|/(2 pi)) estimates the error
		of computing n_d itself. A residual at roundoff level means k is a reciprocal-lattice vector to double
		precision; snapping it is then no less accurate than the direct sum, whose phases k.r carry rounding
		of the same order, so large indices are accepted (earlier builds rejected every mode with
		roundoff > budget, i.e. indices above about 43/dim at eps=1e-12, silently falling back for large N).
		Every eps below about 2e-14*dim (4*DBL_EPSILON > budget) is still rejected before FINUFFT is called.
		Returns false when eps is below that floor, some mode fails the residual test above, |n_d| > 2^30,
		or the mode grid exceeds FinufftMaxModes or FinufftSparseGridFactor*(Nk+N).
		\param basis      Basis vectors of the bound configuration.
		\param np         Particle count of the bound configuration.
		\param[out] ms    Mode-grid sizes 2K_d+1 (left at 1 for unused dimensions).
		\param[out] slot  Position of each constraint's mode in the modeord=0 output (valid only on true).
		\return true if the constraints can use FINUFFT.
		\throw std::bad_alloc if the temporary mode list cannot be allocated. */
	bool MapFinufftModes(const GeometryVector * basis, size_t np, std::int64_t ms[3], std::vector<std::int64_t> & slot) const
	{
		const DimensionType dim=this->Dimension;
		const size_t nk=this->constraints.size();
		const double twopi=2.0*::pi;
		const double maxIndex=static_cast<double>(std::int64_t(1)<<30);
		const double budget=std::min(1e-8, this->FinufftTolerance)/(4.0*twopi*static_cast<double>(dim));
		if(4.0*DBL_EPSILON>budget)
			return false;	//eps below what double-precision phases resolve

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

		//modeord=0: mode n_d sits at n_d+K_d, first dimension fastest
		slot.resize(nk);
		for(size_t i=0; i<nk; i++)
		{
			std::int64_t idx=0;
			for(int d=static_cast<int>(dim)-1; d>=0; d--)
				idx=idx*ms[d]+(modes[i*dim+d]+K[d]);
			slot[i]=idx;
		}
		return true;
	}

	/** \brief Bring FinufftCache up to date with the bound configuration (FinufftInputs, then
		FinufftCacheIsCurrent, rebuilding with RebuildFinufftCache if an input changed).
		\param[out] np  Bound particle count.
		\return false, leaving the cache untouched, if FinufftInputs rejects the inputs; otherwise true, and
			FinufftCache.eligible tells whether FINUFFT can run.
		\throw std::bad_alloc if a rebuild cannot allocate; the cache is then left unbuilt. */
	bool RefreshFinufftCache(size_t & np)
	{
		GeometryVector basis[3];
		if(!this->FinufftInputs(basis, np))
			return false;
		if(!this->FinufftCacheIsCurrent(basis, np))
			this->RebuildFinufftCache(basis, np);
		return true;
	}

	/** \brief Report whether the cached rho came from FINUFFT but the cache no longer holds its angles (a copied or
		assigned object, whose cache starts empty, or a rebuilt cache) while Type 2 plans exist; the Type 2 forces
		then need rho recomputed first.
		\return true if Type 2 plans exist, rho came from FINUFFT and FinufftCache.anglesOfRho is false. */
	bool RhoLacksAngles(void) const
	{
		const FinufftPlanCache & c=this->FinufftCache;
		return !c.plan2.empty() && this->RhoBackendLastUsed==RhoBackend::FinufftType1 && !c.anglesOfRho;
	}

	// ponytail: machine-specific cost constants [microseconds], fitted on 2026-10-14 on the 8-thread (4-core SMT) WSL
	// development machine while another job loaded it: FINUFFT 2.5.1 persistent plans (setpts + execute plus the
	// coefficient and output copies of the Type 2 forces), one thread, eps=1e-12, g++ -O3. The Type 2 constants are
	// the tightest linear upper envelope of 142 cases (1-3D, N=100-10000, chi 0.02-0.8, orthogonal and skew cells);
	// the pair cost is the smallest single-thread rate of the long direct loops. They only choose the faster path,
	// never the result's accuracy. Ceiling: another CPU, compiler or FINUFFT version moves the crossover. Upgrade:
	// refit them on the production node with apps/CCOptimization/tests/finufft_gate (calibrate, fit, check; see its
	// README.md), never from timings taken during a run.
	static constexpr double FinufftCostCall[3]={1.8, 8.3, 129.5};	//< Type 2 fixed cost per execution, by dimension.
	static constexpr double FinufftCostParticle[3]={0.0150, 0.134, 1.95};	//< Type 2 cost per particle (all d transforms) at kernel width 13.
	static constexpr double FinufftCostMode[3]={0.0361, 0.0128, 0.0569};	//< Type 2 cost per mode-grid point.
	static constexpr double FinufftCostPair[3]={0.0126, 0.0189, 0.0212};	//< Direct Force() cost per particle-constraint pair.
	static constexpr double FinufftForceMargin=1.3;	//< Required predicted speedup before Type 2 replaces the direct loop.
	// ponytail: Type 1 cost constants [microseconds] and the block rule, fitted on 2026-10-14 on the same machine under
	// load (another job, load 4-7) with the same tools: tightest linear upper envelope of one Type 1 block (angles,
	// setpts and execute of RunRhoBlock) over the 142 calibration cases. FinufftBlockRatio and FinufftBlockMinWork were
	// chosen with a thread-CPU schedule model of whole evaluations (2D/3D, N=2000/10000, S=2/5, T=4/8, chi 0.4): a
	// Type 1 block must keep at least twice the work every block repeats, since its count cannot depend on the team
	// (it changes rho in the last bits); Type 2 blocks are also limited to the species' share of the team by the gate.
	// Ceiling and upgrade as above: refit with tests/finufft_gate on the production node.
	static constexpr double FinufftRhoCostCall[3]={3.04, 6.16, 84.4};	//< Type 1 fixed cost per execution, by dimension.
	static constexpr double FinufftRhoCostParticle[3]={0.0123, 0.046, 0.561};	//< Type 1 cost per particle at kernel width 13.
	static constexpr double FinufftRhoCostMode[3]={0.0375, 0.0282, 0.0477};	//< Type 1 cost per mode-grid point.
	static constexpr double FinufftBlockRatio[2]={2.0, 1.0};	//< Least particle work of a Type 1 / Type 2 block, in units of the work every block repeats.
	static constexpr double FinufftBlockMinWork=100.0;	//< Least particle work of a block [us].
	static constexpr size_t FinufftMaxBlocks=16;	//< Most particle blocks of one transform.

	/** \brief Cost model of one execution of a Type 1 or Type 2 plan on one thread with the persistent plan (no plan is
		made per call), in microseconds: fixed + perParticle * n for n points, with fixed = call + mode * G repeated by
		every particle block (plan execution, FFT and deconvolution on the mode grid) and perParticle = particle *
		(w/13)^d the spreading or interpolation work, w = ceil(-log10(eps/10)) in [2, 16] the estimated kernel width.
		\param type              1 (FinufftRhoCost* constants) or 2 (FinufftCost*, all Dimension transforms).
		\param gridModes         Mode-grid size G (> 0, so the dimension is 1-3).
		\param[out] fixed        Predicted fixed work of one execution.
		\param[out] perParticle  Predicted work per point. */
	void FinufftExecutionCosts(int type, double gridModes, double & fixed, double & perParticle) const
	{
		const int d=static_cast<int>(this->Dimension)-1;
		const double w=std::min(16.0, std::max(2.0, std::ceil(-std::log10(this->FinufftTolerance/10.0))));
		const double width=std::pow(w/13.0, d+1);
		if(type==1)
		{
			fixed=FinufftRhoCostCall[d]+FinufftRhoCostMode[d]*gridModes;
			perParticle=FinufftRhoCostParticle[d]*width;
		}
		else
		{
			fixed=FinufftCostCall[d]+FinufftCostMode[d]*gridModes;
			perParticle=FinufftCostParticle[d]*width;
		}
	}

	/** \brief Number of particle blocks of one transform: the largest B for which every block's predicted spreading or
		interpolation work, perParticle * N / B, is still at least FinufftBlockRatio[type-1] times the fixed work each block
		repeats (FinufftExecutionCosts) and at least FinufftBlockMinWork, clamped to [1, min(FinufftMaxBlocks, N)]. Splitting pays only on idle threads,
		and the repeated fixed work bounds its cost on a busy team. Depends only on the dimension, N, G and eps (and the
		constants), never on the thread count, so Type 1 rho, whose last bits depend on the blocks, is the same for
		every team; small species stay one block.
		\param type       1 for rho, 2 for forces.
		\param np         Particle count N.
		\param gridModes  Mode-grid size G.
		\return The block count. */
	size_t FinufftBlocks(int type, size_t np, double gridModes) const
	{
		double fixed=0.0, perParticle=0.0;
		this->FinufftExecutionCosts(type, gridModes, fixed, perParticle);
		const double fit=std::floor(perParticle*static_cast<double>(np)/std::max(FinufftBlockRatio[type-1]*fixed, FinufftBlockMinWork));
		return static_cast<size_t>(std::max(1.0, std::min({fit, static_cast<double>(FinufftMaxBlocks), static_cast<double>(np)})));
	}

	/** \brief Cost model of the two force paths of one object, in microseconds (see FinufftForceGate).
		direct = FinufftCostPair * N * Nk / threads, the object's share of the direct force loop with ideal parallel
		speedup (which favors the direct loop); the Type 2 path runs as blocks = min(FinufftBlocks(2), share) blocks, the
		largest of which takes block = fixed + perParticle * ceil(N / blocks) on one thread (FinufftExecutionCosts).
		\param gridModes    Mode-grid size G (> 0, so the dimension is 1-3).
		\param threads      Team size of the direct loop.
		\param share        Most blocks (the object's share of the team).
		\param[out] direct  Predicted direct-loop time.
		\param[out] block   Predicted time of one Type 2 block.
		\param[out] blocks  Number of Type 2 blocks. */
	void FinufftForceCosts(double gridModes, double threads, size_t share, double & direct, double & block, size_t & blocks) const
	{
		const int d=static_cast<int>(this->Dimension)-1;
		const size_t n=this->NumBoundParticle();
		const double np=static_cast<double>(n);
		direct=FinufftCostPair[d]*np*static_cast<double>(this->constraints.size())/threads;
		double fixed=0.0, perParticle=0.0;
		this->FinufftExecutionCosts(2, gridModes, fixed, perParticle);
		blocks=std::min(this->FinufftBlocks(2, n, gridModes), share);
		block=fixed+perParticle*std::ceil(np/static_cast<double>(blocks));
	}
#endif

	/** \brief Fill RhoReal/RhoImag with the direct sum rho(k) = sum_j exp(+i k.r_j), parallel over constraints with
		ParallelNumThread threads (each constraint sums its particles in order, so the team size does not change the
		bits), and mark rho valid with backend Direct. Volume must be set (StartRho). */
	void GetRho_Direct(void)
	{
		signed long end=this->constraints.size();
		this->RhoReal.resize(end);
		this->RhoImag.resize(end);
		//gather the bound coordinates once: the O(N*Nk) loop then reads them contiguously, as a sub-configuration would
		const Configuration & config=this->BoundConfig();
		std::vector<GeometryVector> coords(this->NumBoundParticle());
		for(size_t j=0; j<coords.size(); j++)
			coords[j]=config.GetCartesianCoordinates(this->BoundIndex(j));
		const size_t np=coords.size();

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for(signed long i=0; i<end; i++)
		{
			double tr=0, ti=0;//real and imag parts of result
#ifdef HAVETWOFIXEDPARTICLES
			tr+=2;
#endif
			for(size_t j=0; j<np; j++)
			{
				const GeometryVector & c=coords[j];
				tr+=std::cos(this->constraints[i].k.Dot(c));
				ti+=std::sin(this->constraints[i].k.Dot(c));
			}
			this->RhoReal[i]=tr;
			this->RhoImag[i]=ti;
		}

		this->RhoBackendLastUsed=RhoBackend::Direct;
		this->RhoValid=true;
	}

public:
	//this class will manage these members, users can also read it.
	std::vector<double> RhoReal, RhoImag;
	bool RhoValid;
	/** \brief Compute and cache the complex collective density rho(k) at every constraint.
		Uses the backend chosen with SetRhoBackend: StartRho(false), every Type 1 block (RunRhoBlock, in parallel over
		ParallelNumThread threads; the result does not depend on their order), then FinishRho.
		When FinufftType1 is selected but cannot be used, this falls back to the direct sum, so rho is always valid on
		return. FINUFFT is not used when the dimension is not 1-3, there are no constraints or no particles,
		MapFinufftModes rejects the constraints, or finufft_makeplan failed for these inputs (its code is reported again
		on every call until an input changes); a non-finite fractional coordinate, a failed allocation or a nonzero
		finufft_setpts or finufft_execute status falls back for that call only.
		Call it from one thread at a time per object (before any parallel force loop): it writes the cached rho
		and, for FinufftType1, executes this object's persistent FINUFFT plans. */
	virtual void GetRho(void)
	{
		assert(this->pConfig!=nullptr || this->pViewConfig!=nullptr);
		if(this->RhoValid)
		{
			return;
		}
		const long blocks=static_cast<long>(this->StartRho(false));
#pragma omp parallel for schedule(dynamic, 1) num_threads(this->ParallelNumThread) if(blocks>1)
		for(long b=0; b<blocks; b++)
			this->RunRhoBlock(b);
		this->FinishRho();
	}

	/** \brief Construct an empty collective-coordinate potential.
		\param dimension  Spatial dimension of the system. */
	ShiftedCCPotential(DimensionType dimension):Potential(dimension)
	{
		this->pConfig=nullptr;
		this->pViewConfig=nullptr;
		this->RhoValid=false;
		this->ParallelNumThread=1;
		this->Shift=0.0;
		this->RhoBackendSelected=RhoBackend::Direct;
		this->RhoBackendLastUsed=RhoBackend::Direct;
		this->ForceBackendLastUsed=ForceBackend::Direct;
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
		const GeometryVector & r=this->BoundConfig().GetCartesianCoordinates(this->BoundIndex(i));

			for(signed long j=0; j<end; j++)
			{
				GeometryVector k=this->constraints[j].k;
				result.MinusFrom((2.0)*this->constraints[j].V*(this->RhoImag[j]*std::cos(k.Dot(r))-this->RhoReal[j]*std::sin(k.Dot(r)))*k);
			}

		result.MultiplyFrom(1.0/Volume);
	}

	/** \brief Make rho valid, then compute every force with the persistent FINUFFT Type 2 plans if the last
		PrepareFinufftForces accepted them and rho came from FINUFFT Type 1: StartForces, every Type 2 block
		(RunForceBlock, in parallel over ParallelNumThread threads), FinishForces. A valid FINUFFT rho whose angles
		the cache no longer holds (a copied or assigned object, whose cache starts empty, or a rebuilt cache) is
		recomputed first; Type 1 is deterministic, so rho keeps its bits. Never makes a plan. Records the outcome for
		GetLastForceBackendUsed: FinufftType2 on success, otherwise Direct, in which case the caller must compute the
		forces with the direct sum of Force() (as AllForce does).
		\param[out] results  Holds the forces (resized to the particle count) when true is returned; unspecified otherwise.
		\return true if results holds the Type 2 forces, false if the direct sum is required. */
	bool TryFinufftAllForce(std::vector<GeometryVector> & results)
	{
		this->ForceBackendLastUsed=ForceBackend::Direct;
		if(this->RhoValid==false)
			this->GetRho();
#ifdef CCPOTENTIAL_USE_FINUFFT
		if(this->RhoLacksAngles())
		{
			this->RhoValid=false;
			this->GetRho();
		}
		const long blocks=static_cast<long>(this->StartForces(results));
#pragma omp parallel for schedule(dynamic, 1) num_threads(this->ParallelNumThread) if(blocks>1)
		for(long b=0; b<blocks; b++)
			this->RunForceBlock(b, results);
		return this->FinishForces();
#else
		return false;
#endif
	}

	/** \brief Compute the force on every particle.
		When UseFinufftForces() holds (FinufftType1 backend, exact type ShiftedCCPotential, and a predicted
		win of the cost gate), FINUFFT Type 2 gives all forces, in the particle blocks the gate chose, run in parallel
		(PrepareFinufft, PrepareFinufftForces, then TryFinufftAllForce); otherwise, or if it falls back, every particle
		uses the direct sum of Force().
		The exact-type test keeps derived classes that override Force or GetRho (S2, LogCC, MultiConfig, ...)
		on their own formulas. See GetLastForceBackendUsed.
		\param[out] results  Resized list of per-particle forces. */
	virtual void AllForce(std::vector<GeometryVector> & results)
	{
		this->PrepareFinufft();
		this->PrepareFinufftForces(FinufftForceGate({this}, this->ParallelNumThread)[0]);
		if(this->TryFinufftAllForce(results))
			return;
		long end=this->NumBoundParticle();
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
	/** \brief Bind the configuration, unbind any indexed view, and invalidate the cached collective density.
		\param Position  Configuration to bind; must remain alive while energy and force are evaluated. */
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->Dimension==Position.GetDimension());
		this->pConfig=&Position;
		this->pViewConfig=nullptr;
		this->RhoValid=false;
	}

	/** \brief Report whether SetIndexedConfiguration may bind this object.
		True only for ShiftedCCPotential itself: derived classes may read pConfig directly or cache state in their
		SetConfiguration override. A derived class may override this to return true only if it reads particles
		exclusively through the base-class GetRho, Force and AllForce and needs no SetConfiguration override.
		\return true if an indexed view is supported. */
	virtual bool SupportsIndexedView(void) const
	{
		return typeid(*this)==typeid(ShiftedCCPotential);
	}

	/** \brief Bind an indexed view and invalidate the cached collective density: local particle j is particle
		indices[j] of full. The base-class GetRho (both backends), Force and AllForce read particles through the
		view and use the cell of full, so they give the same bits as SetConfiguration on a sub-configuration of full
		holding those particles in that order. pConfig is nullptr while the view is bound, so code that reads it
		directly fails instead of reading the wrong particles. SetConfiguration unbinds the view.
		\param full     Configuration holding the particles; must stay alive at the same address while bound.
		\param indices  Local-to-global particle map (copied); every entry must be below full.NumParticle().
		\throw std::logic_error if SupportsIndexedView() is false; std::out_of_range for an invalid index.
			A failed check leaves the object unchanged. */
	void SetIndexedConfiguration(const Configuration & full, const std::vector<size_t> & indices)
	{
		if(!this->SupportsIndexedView())
			throw std::logic_error("ShiftedCCPotential::SetIndexedConfiguration: not supported by this derived class; use SetConfiguration.");
		assert(this->Dimension==full.GetDimension());
		for(size_t i : indices)
			if(i>=full.NumParticle())
				throw std::out_of_range("ShiftedCCPotential::SetIndexedConfiguration: particle index out of range.");
		this->ViewIndices=indices;
		this->pViewConfig=&full;
		this->pConfig=nullptr;
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
