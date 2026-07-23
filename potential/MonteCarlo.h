/**	Author	: Ge Zhang
 *	Email	:
 *	Date	:
 *  Annotated by Jaeuk Kim (2019). */

/** \file MonteCarlo.h
 *	\brief	Header file for performing the Monte Carlo simulations for general systems,
	including point configurations. Contains the sampling drivers (MonteCarlo,
	WangLandauMonteCarlo, ParallelTempering); the cooling schedules were split
	into CoolingSchedule.h (2026-07-07). The base MonteCarlo::Move / ::Anneal
	bodies are restored inline from the legacy c++/Potential/MonteCarlo.h
	(the split-out new_c++/Potential/MonteCarlo.cpp was ill-formed and never
	compiled, leaving them undefined). */
#ifndef MONTECARLO_INCLUDED
#define MONTECARLO_INCLUDED

//General-Purpose Monte Carlo, not necessary for particles

#include <etc.h>
#include <RandomGenerator.h>
#include "Move.h"
#include "CoolingSchedule.h"

#include <functional>
#include <cmath>

//General-Purpose Monte Carlo, not necessary for particles
template<typename System, typename AdditionalInfo> class MonteCarlo
{
protected:
	RandomGenerator gen;
public:

	double RecordEnergy;//this class keeps a log of the system energy, this is valid only if the potential is never changed
	size_t MoveCount;
	System * pSys;

	//Some trial moves, like particle displacement, has step size that needs to be adjusted.
	//The adjustment can be automatic. However, it must be turned off during sampling.
	//Set LockStepSize=true to turn off automatic adjustment of trial move parameters.
	bool LockStepSize;

	/** \brief Construct a Monte Carlo driver bound to a system.
		\param sys         System to simulate.
		\param RandomSeed  Seed for the random number generator.
		\param Energy      Initial recorded energy of the system. */
	MonteCarlo(System & sys, int RandomSeed, double Energy) : gen(RandomSeed)
	{
		this->MoveCount=0;
		this->pSys= & sys;
		this->RecordEnergy=Energy;
		this->LockStepSize=false;
	}

	/** \brief Destructor. */
	virtual ~MonteCarlo()
	{
	}


	/** \brief Monte Carlo move function at a given temperature.
		\param Temperature  Temperature at which the Metropolis acceptance is evaluated.
		\param Repeat The number of total trial moves.
		\param move	Definition of a trial move
		\param add  Additional information (e.g. the interaction potential) passed to the move. */
	virtual void Move(double Temperature, size_t Repeat, MCMove<System, AdditionalInfo> & move, AdditionalInfo add)
	{
		for(size_t i=0; i<Repeat; i++)
		{
			MoveCount++;
			double PreFactor=1.0;
			double dE=move.DeltaEnergy(*this->pSys, add, this->gen, PreFactor, LockStepSize, RecordEnergy);
			if(Temperature ==0)
			{
				if(dE>0)
					continue;
			}
			else if(gen.RandomDouble() > PreFactor*std::exp((-1)*dE/Temperature))
				continue;

			move.Accept(*this->pSys, add);
			this->RecordEnergy+=dE;
		}
	}


	/** \brief Perform simulated annealing with a cooling schedule.
		\param add  Additional information (e.g. the interaction potential) passed to the move.
		\param cool  Cooling schedule that controls the temperature sequence.
		\param move  Definition of a trial move.
		\param CallBack  Function called after each cooling stage with the energy, temperature, and system. */
	virtual void Anneal(AdditionalInfo add, CoolingSchedule & cool, MCMove<System, AdditionalInfo> & move, std::function<void(double & Energy, double Temperature, System & sys)> CallBack)
	{
		while(cool.Continue)
		{
			this->Move(cool.Temperature, cool.NumTry, move, add);
			CallBack(this->RecordEnergy, cool.Temperature, *pSys);

			cool.Report(this->RecordEnergy);
		}
	}
};
//class for Wang-Landau Monte Carlo
//class BinDecider should have the following member functions:
	//size_t NumBins()
	//long GetBin(double E)
	//double GetBinLowerBound(size_t NumBin)
	//double GetBinUpperBound(size_t NumBin)
//The user is responsible for adjusting SIncrease
/** \brief Wang-Landau Monte Carlo driver that estimates the density of states over energy bins. */
template<typename System, typename AdditionalInfo, typename BinDecider> class WangLandauMonteCarlo
{
private:
	std::vector<double> s;
	std::vector<size_t> Histogram;
	size_t MoveCount;
	RandomGenerator gen;
	BinDecider * pBinDecider;
public:
	bool LockStepSize;
	double RecordEnergy;
	System * pSys;
	double SIncrease;
	/** \brief Construct a Wang-Landau driver bound to a system and energy-bin decider.
		\param sys         System to simulate.
		\param RandomSeed  Seed for the random number generator.
		\param Energy      Initial energy of the system; must lie within the bin range.
		\param d           Bin decider that maps energies to histogram bins. */
	WangLandauMonteCarlo(System & sys, int RandomSeed, double Energy, BinDecider & d) :
		s(d.NumBins(), 0.0), Histogram(d.NumBins(), 0), gen(RandomSeed), RecordEnergy(Energy)
	{
		this->MoveCount=0;
		this->pSys= & sys;
		this->LockStepSize=false;
		this->pBinDecider=& d;
		long initBin = d.GetBin(Energy);//GetBin returns a signed bin index (may be -1 or NumBins)
		if(initBin >= static_cast<long>(d.NumBins()) || initBin<0 )
		{
			std::cerr<<"Error in WangLandauMonteCarlo : initial energy not within bounds! Bin="<<initBin<<"\n";
			assert(false);
		}
		this->SIncrease=1.0;
	}

	/** \brief Reset the visit histogram to zero. */
	void ClearHistogram(void)
	{
		for(auto iter=this->Histogram.begin(); iter!=this->Histogram.end(); iter++)
			*iter=0;
	}

	/** \brief Return the visit histogram as (bin-center energy, count) pairs.
		\return A vector of 2D points: bin-center energy and visit count. */
	std::vector<GeometryVector> GetHistogram(void)
	{
		std::vector<GeometryVector> result;
		for(size_t i=0; i<this->Histogram.size(); i++)
		{
			GeometryVector t(2);
			t.x[0]=(pBinDecider->GetBinLowerBound(i)+pBinDecider->GetBinUpperBound(i))*0.5;
			t.x[1]=this->Histogram[i];
			result.push_back(t);
		}
		return result;
	}

	/** \brief Return the estimated entropy (log density of states) as (bin-center energy, entropy) pairs.
		\return A vector of 2D points: bin-center energy and entropy. */
	std::vector<GeometryVector> GetEntropy(void)
	{
		std::vector<GeometryVector> result;
		for(size_t i=0; i<this->Histogram.size(); i++)
		{
			GeometryVector t(2);
			t.x[0]=(pBinDecider->GetBinLowerBound(i)+pBinDecider->GetBinUpperBound(i))*0.5;
			t.x[1]=this->s[i];
			result.push_back(t);
		}
		return result;
	}

	/** \brief Perform Wang-Landau trial moves, updating the entropy estimate and histogram.
		\param Repeat         Number of trial moves to perform.
		\param move           Definition of a trial move.
		\param add            Additional information (e.g. the interaction potential) passed to the move.
		\param GetEnergyFunc  Function returning the exact energy of a system, used to periodically re-sync RecordEnergy. */
	void Move(size_t Repeat, MCMove<System, AdditionalInfo> & move, AdditionalInfo add, std::function<double(const System & sys)> GetEnergyFunc)
	{
		long Nbins = static_cast<long>(pBinDecider->NumBins());//signed: GetBin() may return -1 or Nbins for out-of-range energies
		for(size_t i=0; i<Repeat; i++)
		{
			MoveCount++;
			double PreFactor=1.0;
			double dE=move.DeltaEnergy(*this->pSys, add, this->gen, PreFactor, LockStepSize, RecordEnergy);
			long PrevBin = pBinDecider->GetBin(RecordEnergy);
			long AfterBin = pBinDecider->GetBin(RecordEnergy+dE);


			if( (PrevBin+1 < 0) || (PrevBin > Nbins) )
			{
				std::cerr<<"Error in WangLandauMonteCarlo : PrevBin out of range!\n";
				std::cerr<<"PrevBin="<<PrevBin<<", RecordEnergy="<<RecordEnergy<<", MoveCount="<<MoveCount<<'\n';
				assert(false);
			}
			else if(PrevBin==-1)
				PrevBin=0;
			else if(PrevBin==Nbins)
				PrevBin=Nbins-1;
			if(AfterBin < 0 || AfterBin >= Nbins)
			{
				s[PrevBin]+=SIncrease;
				Histogram[PrevBin]++;
			}
			else if(gen.RandomDouble() < PreFactor*std::exp(s[PrevBin]-s[AfterBin]))
			{
				move.Accept(*this->pSys, add);
				this->RecordEnergy+=dE;
				s[AfterBin]+=SIncrease;
				Histogram[AfterBin]++;
			}
			else
			{
				s[PrevBin]+=SIncrease;
				Histogram[PrevBin]++;
			}
			if(MoveCount%1000==0)
			{
				double newE=GetEnergyFunc(*pSys);
				//scale-aware tolerance (absolute floor 1e-8 plus 1e-8 relative to the larger magnitude);
				//replaces the former newE<1e5 magic gate that silently trusted large energies
				double EScale=std::abs(newE)>std::abs(RecordEnergy)?std::abs(newE):std::abs(RecordEnergy);
				if(std::isfinite(newE) && std::abs(newE-RecordEnergy) > 1e-8*(1.0+EScale))
					std::cerr<<"Error in WangLandauMonteCarlo: Energy not consistent, newE="<<newE<<", record E="<<RecordEnergy<<'\n';

				this->RecordEnergy=newE;
			}
		}
	}
};

const size_t IndividualMoveTime=25;
//const double SwapMovePropability=0.2;
/** \brief Parallel-tempering Monte Carlo driver that runs replicas at rescaled temperatures and swaps them. */
template<typename System, typename AdditionalInfo> class ParallelTempering : public MonteCarlo<System, AdditionalInfo>
{
public:
	std::vector<double> vEnergyFactors;
	std::vector<MonteCarlo<System, AdditionalInfo> > vMonteCarlos;
	//one of the systems is given(others are allocated by this class), the pointer pSys points to it.
	size_t SwapSuccessCount, SwapFailCount;
	System * pOriginalSys;

	/** \brief Construct a parallel-tempering driver with several replicas.
		\param sys         System for the first replica; additional replicas are copied from it.
		\param RandomSeed  Base seed for the replicas' random number generators.
		\param NumCopy     Number of replicas (temperature copies).
		\param Energy      Initial energy assigned to every replica. */
	ParallelTempering(System & sys, int RandomSeed, size_t NumCopy, double Energy) : MonteCarlo<System, AdditionalInfo>(sys, RandomSeed, Energy)
	{
		if(NumCopy==0)
		{
			std::cerr<<"Error in Parallel Tempering : NumCopy is 0\n";
			assert(false);
		}
		{
			this->vEnergyFactors.clear();
			double nowRescale=1.0;
			for(size_t i=0; i<NumCopy; i++)
			{
				this->vEnergyFactors.push_back(nowRescale);
				nowRescale*=1.2;
			}
		}

		for(size_t i=0; i<this->vEnergyFactors.size(); i++)
		{
			if(i==0)
				this->vMonteCarlos.push_back(MonteCarlo<System, AdditionalInfo>(sys, (RandomSeed+i), Energy ));
			else
				this->vMonteCarlos.push_back(MonteCarlo<System, AdditionalInfo>(*(new System(sys)), (RandomSeed+i), Energy ));
			this->vMonteCarlos[i].RecordEnergy=Energy;
		}
		this->pSys = & sys;
		this->pOriginalSys= & sys;
		this->SwapFailCount=0;
		this->SwapSuccessCount=0;
	}

	/** \brief Destructor; frees the replica systems allocated by this class. */
	~ParallelTempering()
	{
		//debug temp
		//std::cerr<<"Parallel Tempering destructed, Fail times:"<<this->SwapFailCount<<", Success Count:"<<this->SwapSuccessCount<<'\n';

		for(size_t i=0; i<this->vMonteCarlos.size(); i++)
			if(this->vMonteCarlos[i].pSys!= this->pOriginalSys)
				delete this->vMonteCarlos[i].pSys;
	}

	//This driver owns the replica systems it allocated (freed above by pointer comparison),
	//so a shallow copy would double-delete them; copying a running driver is never needed.
	ParallelTempering(const ParallelTempering &) = delete;
	ParallelTempering & operator=(const ParallelTempering &) = delete;

	/** \brief Advance each replica independently at its rescaled temperature.
		\param Temperature  Base temperature; replica i runs at Temperature divided by its energy factor.
		\param Repeat       Number of trial moves per replica.
		\param moves        Per-replica trial move definitions.
		\param add          Additional information (e.g. the interaction potential) passed to the moves. */
	void IndividualMove(double Temperature, size_t Repeat, std::vector<MCMove<System, AdditionalInfo> * > & moves, AdditionalInfo add)
	{
		size_t stop = this->vEnergyFactors.size();
		assert(this->vMonteCarlos.size() == stop);
		assert(moves.size() == stop);

		//#pragma omp parallel for
		//(the former per-call NumThreads computation was dead: the pragma above is disabled)
		for(signed long i=0; i<(signed long)stop; i++)
			this->vMonteCarlos[i].Move(Temperature/this->vEnergyFactors[i], Repeat, *moves[i], add);
	}

	/** \brief Perform parallel-tempering moves, mixing replica swaps and individual moves.
		\param Temperature  Base temperature controlling swap acceptance and replica temperatures.
		\param Repeat       Number of moves to perform.
		\param moves        Per-replica trial move definitions.
		\param add          Additional information (e.g. the interaction potential) passed to the moves. */
	virtual void Move(double Temperature, size_t Repeat, std::vector<MCMove<System, AdditionalInfo> * > & moves, AdditionalInfo add)
	{
		//sync this->RecordEnergy
		for(auto iter = this->vMonteCarlos.begin(); iter!=this->vMonteCarlos.end(); iter++)
		{
			if(iter->pSys == this->pSys)
			{
				iter->RecordEnergy=this->RecordEnergy;
				break;
			}
		}

		for(size_t i=0; i<Repeat; i++)
		{
			if(this->gen.RandomDouble() < 0.9)
			{
				//swap moves (need at least two replicas; a single replica has nothing to swap,
				//and the index modulo (size-1) below would divide by zero when size==1)
				if(this->vMonteCarlos.size()>=2)
				for(size_t i=0; i<this->vMonteCarlos.size(); i++)
				{
					size_t m=static_cast<size_t>(this->gen.RandomDouble()*100000000)%(this->vMonteCarlos.size()-1);
					size_t n=m+1;
					if(n>=this->vMonteCarlos.size())
					{
						std::cerr<<"error in ParallelTempering : index out of range\n";
						continue;
					}

					double DeltaE=this->vMonteCarlos[n].RecordEnergy*this->vEnergyFactors[m]+this->vMonteCarlos[m].RecordEnergy*this->vEnergyFactors[n]-this->vMonteCarlos[n].RecordEnergy*this->vEnergyFactors[n]-this->vMonteCarlos[m].RecordEnergy*this->vEnergyFactors[m];

					if(Temperature ==0)
					{
						if(DeltaE>0)
						{
							this->SwapFailCount++;
							continue;
						}
					}
					else if(this->gen.RandomDouble() > std::exp((-1)*DeltaE/Temperature))
					{
						this->SwapFailCount++;
						continue;
					}

					std::swap(this->vMonteCarlos[m].pSys, this->vMonteCarlos[n].pSys);
					std::swap(this->vMonteCarlos[m].RecordEnergy, this->vMonteCarlos[n].RecordEnergy);
					this->SwapSuccessCount++;

					//debug temp
					if(this->pSys!=this->vMonteCarlos.back().pSys)
					{
						this->pSys=this->vMonteCarlos.back().pSys;
						//std::cerr<<"pSys changed!";
					}

				}

				//modify energy factors if necessary
				size_t count=this->SwapFailCount+this->SwapSuccessCount;
				if(count>1000 && this->vEnergyFactors.size()>1)
				{
					//debug temp
					//std::cout<<"Swap success:"<<this->SwapSuccessCount<<", fail:"<<this->SwapFailCount<<'\n';
					if(static_cast<double>(this->SwapSuccessCount)/count<0.1)
					{
						double factor=std::pow(this->vEnergyFactors[1], 0.7);
						for(size_t i=0; i<this->vEnergyFactors.size(); i++)
							this->vEnergyFactors[i]=std::pow(factor, static_cast<double>(i));

						//debug temp
						std::cout<<"Energy factor has been changed to:"<<factor<<'\n';
					}
					this->SwapFailCount=0;
					this->SwapSuccessCount=0;
				}
			}
			else
			{
				this->IndividualMove(Temperature, IndividualMoveTime, moves, add);
				i+=IndividualMoveTime;
			}
		}
		this->MoveCount += Repeat;

		//sync this->RecordEnergy
		for(auto iter = this->vMonteCarlos.begin(); iter!=this->vMonteCarlos.end(); iter++)
		{
			if(iter->pSys == this->pSys)
			{
				this->RecordEnergy=iter->RecordEnergy;
				break;
			}
		}
	}

	/** \brief Perform parallel-tempering moves using a single move template cloned per replica.
		\param Temperature  Base temperature controlling swap acceptance and replica temperatures.
		\param Repeat       Number of moves to perform.
		\param move         Trial move template cloned once per replica.
		\param add          Additional information (e.g. the interaction potential) passed to the moves. */
	virtual void Move(double Temperature, size_t Repeat, MCMove<System, AdditionalInfo>  & move, AdditionalInfo add)
	{
		std::vector<MCMove<System, AdditionalInfo> * > moves;
		for(size_t i=0; i<this->vMonteCarlos.size(); i++)
			moves.push_back(move.clone());

		this->Move(Temperature, Repeat, moves, add);
		//delete stuff
		for(auto iter=moves.begin(); iter!=moves.end(); iter++)
			delete *iter;
	}

	/** \brief Run parallel-tempering simulated annealing with a cooling schedule.
		\param add       Additional information (e.g. the interaction potential) passed to the moves.
		\param cool      Cooling schedule that controls the temperature sequence.
		\param move      Trial move template cloned once per replica.
		\param CallBack  Function called after each cooling stage with the energy, temperature, and system. */
	virtual void Anneal(AdditionalInfo add, CoolingSchedule & cool, MCMove<System, AdditionalInfo> & move, std::function<void(double & Energy, double Temperature, System & sys)> CallBack)
	{
		std::vector<MCMove<System, AdditionalInfo> * > moves;
		for(size_t i=0; i<this->vMonteCarlos.size(); i++)
			moves.push_back(move.clone());

		while(cool.Continue)
		{
			this->Move(cool.Temperature, cool.NumTry, moves, add);
			CallBack(this->RecordEnergy, cool.Temperature, *this->pSys);

			//report the average energy of all configurations to the cooling schedule
			double AverageEnergy=0;
			for(size_t i=0; i<this->vMonteCarlos.size(); i++)
				AverageEnergy+=this->vEnergyFactors[i]*this->vMonteCarlos[i].RecordEnergy;
			AverageEnergy/=this->vMonteCarlos.size();
			cool.Report(AverageEnergy);

			//debug temp
			//for(size_t i=0; i<cool.NumTry; i++)
			//{
			//	this->Move(cool.Temperature, 1, move, add);
			//	CallBack(this->RecordEnergy, cool.Temperature, *pSys);
			//}
			//cool.Report(this->RecordEnergy);
		}
		//delete stuff
		for(auto iter=moves.begin(); iter!=moves.end(); iter++)
			delete *iter;

		//the one with lowest temperature is returned
		if(this->pOriginalSys!=this->pSys)
			(*this->pOriginalSys)=(*this->vMonteCarlos.back().pSys);
	}
};

#endif
