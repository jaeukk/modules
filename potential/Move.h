/** Author: Ge Zhang
	email:	?
	Date:	?
	Comments: Jaeuk Kim (June 2019). */

/** \file Move.h
	\brief Header file for the MCMove class and its derivative classes.
	Those classes are designed to perform the Monte Carlo Simulations for point configurations as well as any type of systems.*/


#ifndef MCMOVE_INCLUDED
#define MCMOVE_INCLUDED

#include <etc.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>

#include "Potential.h"

#include <iostream>
#include <cassert>
#include <exception>
#include <memory>
const size_t AcceptanceSampleSize=1000;
const double MaxAcceptance=0.7;
const double MinAcceptance=0.3;


const double MaxAspectRatioOfCell=1.9;
const double MaxCellSideLengthSquared=10000.0;

/** \brief A virtual class to support general MC simulations. 
 *	System: a class to define systems of interest.
	AdditionalInfo:	a class to define interactions. 
	To understand a usage, please see the AtomMove class. */
template <typename System, typename AdditionalInfo>
class MCMove
{
public:
	/** \brief Generate a random move, return Energy difference. the PreFactor is also set if it is not 1.
	 	\param sys	Define current state of a system.
		\param add	Define interaction
		\param gen	Random number generator.
		\param PreFactor	Acceptance prefactor; set when it differs from 1.
		\param LockStepSize	If true, do not auto-adjust the trial-move step size.
		\param CurrentEnergy	Energy of the current state.
		\return Change in energy from the trial move. */
	virtual double DeltaEnergy(System & sys, AdditionalInfo add, RandomGenerator & gen,
		double & PreFactor, bool LockStepSize, double CurrentEnergy) =0;

	/** \brief A pure virtual function to determine whether a random trial is accepted or not.
		\param sys	Define the current system state.
		\param add	Define interaction. */
	virtual void Accept(System & sys, AdditionalInfo add) =0;

	/** \brief Return a heap-allocated copy of this move.
		\return A new clone of this move. */
	virtual MCMove<System, AdditionalInfo> * clone() =0;

	/** \brief Virtual destructor so derived moves are destroyed correctly through a base pointer. */
	virtual ~MCMove() {}
};

/** \brief A derivative of MCMove specialized in point configurations in canonical ensemble. */
class AtomMove : public MCMove<Configuration, Potential *>
{
private:
	size_t TrialCount,			//< The number of total number of trial moves.
		AcceptCount;			//< The number of accepted moves.
	size_t moved;				//< Index of the particle with a trial move. 
	GeometryVector newcoord;	//< A trial move in Relative coordinates.

public:
	double AcceptRatio;			//< the most recently calculated acceptance ratio, updated by the class
	double sigma;				//< 
	double MyMinAcceptance;		//<

	double MinDistance;			//< abandon this move immediately if particles are within this distance
	
	/** \brief A default constructor. */
	AtomMove()
	{
		this->sigma=0.1;
		this->TrialCount=0;
		this->AcceptCount=0;
		this->moved=0;
		this->AcceptRatio=0;
		this->MyMinAcceptance= ::MinAcceptance;
		this->MinDistance=0;
	}

	//calculate the delta E for a given move
	//move is : move the particle (this->moved) so that its new relative coordinate is newcoord
	/** \brief Calculate the change in energy when the newcoord is accepted.
		\param list	Define the current state.
		\param pot	Define interaction.	This object can be one of three classes: PairPotential, forMCPotential, and CombinedPotential.
		\param newcoord	A trial move in relative coordiantes.
		\param CurrentEnergy	energy of the current state.
		\return Change in energy. */
	double CalculateDeltaEnergy(Configuration & list, Potential * pot, const GeometryVector & newcoord, double CurrentEnergy)
	{
		//Determine which class pot belongs to. 
		PairPotential * pPair = dynamic_cast<PairPotential *>(pot);
		forMCPotential * pForMCPotential = dynamic_cast<forMCPotential *>(pot);
		CombinedPotential * pComb = dynamic_cast<CombinedPotential *>(pot);
		if(pPair!=nullptr)
		{
			double deltaEnergy=0;
			list.IterateThroughNeighbors(newcoord, pPair->Rcut, [&](const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * /*PeriodicShift*/, const size_t SourceAtom) -> void 
			{
				if(SourceAtom!=moved)
					deltaEnergy+=pPair->EnergyFormula(shift);
			});
			list.IterateThroughNeighbors(list.GetRelativeCoordinates(moved), pPair->Rcut, [&](const GeometryVector & shift, const GeometryVector & /*LatticeShift*/, const signed long * /*PeriodicShift*/, const size_t SourceAtom) -> void 
			{
				if(SourceAtom!=moved)
					deltaEnergy-=pPair->EnergyFormula(shift);
			});
			return deltaEnergy;
		}
		else if (pForMCPotential!=nullptr)
		{
			GeometryVector prevCart=list.GetCartesianCoordinates(this->moved);
			GeometryVector afterCart=list.RelativeCoord2CartesianCoord(this->newcoord);
			return pForMCPotential->TryMove(prevCart, afterCart);
		}
		else if (pComb!=nullptr)
		{
			//recursion
			return this->CalculateDeltaEnergy(list, pComb->p1, newcoord, CurrentEnergy)+this->CalculateDeltaEnergy(list, pComb->p2, newcoord, CurrentEnergy);
		}
		else
		{	//When pot belongs to none of the three classes.
			GeometryVector backup=list.GetRelativeCoordinates(this->moved);
			pot->SetConfiguration(list);
			double E1=pot->Energy();

			list.MoveParticle(this->moved, this->newcoord); 
			pot->SetConfiguration(list);
			double E2=pot->Energy();
			list.MoveParticle(this->moved, backup); 

			return E2-E1;
		}
	}

	/** \brief Generate a random single-particle displacement and return its energy difference.
		\param list Current point configuration.
		\param pot Interaction potential.
		\param gen Random number generator.
		\param PreFactor Acceptance prefactor; left unchanged by this move.
		\param LockStepSize If true, do not auto-adjust the displacement step size.
		\param CurrentEnergy Energy of the current state.
		\return Change in energy from the trial displacement, or MaxEnergy if a too-near neighbor is found. */
	virtual double DeltaEnergy(Configuration & list, Potential * pot, RandomGenerator & gen,
		double & /*PreFactor*/, bool LockStepSize, double CurrentEnergy)//Generate a random move, return Energy difference
	{
		if(this->TrialCount>=::AcceptanceSampleSize )
		{
			//TODO : THIS MAY BREAK DETAILED BALANCE. FIND A WAY TO FIX IT.
			AcceptRatio=static_cast<double>(this->AcceptCount)/this->TrialCount;//acceptance rate
			if(LockStepSize==false)
			{
				if(AcceptRatio> ::MaxAcceptance && this->sigma<0.5)
				{
					this->sigma*=1.1;
					//std::cout<<"In AtomMove, step size changed to"<<this->sigma<<"\n";
				}
				else if(AcceptRatio < this->MyMinAcceptance)
				{
					this->sigma*=0.9;
					//std::cout<<"In AtomMove, step size changed to"<<this->sigma<<"\n";
				}
				//else
					//std::cout<<"In AtomMove, keep step size"<<this->sigma<<"\n";
			}
			//debug temp
			//std::cout<<"Sigma is:"<<this->sigma<<", Accept Count is:"<<this->AcceptCount<<'\n';
			this->TrialCount=0;
			this->AcceptCount=0;
		}
		this->TrialCount++;

		bool TooNearNeighborFound=true;

		this->moved=list.GetRandomParticle(gen);
		this->newcoord=list.GetRelativeCoordinates(this->moved);
		DimensionType dim=list.GetDimension();
		for(int i=0; i<dim; i++)
			this->newcoord.x[i]+=this->sigma*(gen.RandomDouble()-0.5);

		TooNearNeighborFound=false;
		auto moved=this->moved;
		list.IterateThroughNeighbors(newcoord, this->MinDistance, [&TooNearNeighborFound, &moved](const GeometryVector & /*shift*/, const GeometryVector & /*LatticeShift*/, const signed long * /*PeriodicShift*/, const size_t SourceAtom) -> void 
		{
			if(SourceAtom!=moved)
				TooNearNeighborFound=true;
		}, &TooNearNeighborFound);

		//new_idea
		if(TooNearNeighborFound)
			return ::MaxEnergy;
		//restore?ReTrialCount++;

		return CalculateDeltaEnergy(list, pot, newcoord, CurrentEnergy);
	}

	//check if the potential needs specific operation for accepting the move
	/** \brief Notify the potential (recursing into combined potentials) that the move was accepted.
		\param pot  Interaction potential to update on acceptance. */
	void AcceptCheck(Potential * pot)
	{
		forMCPotential * pForMCPotential = dynamic_cast<forMCPotential *>(pot);
		CombinedPotential * pComb= dynamic_cast<CombinedPotential *>(pot);
		if(pForMCPotential!=nullptr)
			pForMCPotential->AcceptMove();
		else if(pComb!=nullptr)
		{
			this->AcceptCheck(pComb->p1);
			this->AcceptCheck(pComb->p2);
		}
	}

	/** \brief Accept the last trial displacement, applying it to the configuration.
		\param list  Configuration to update with the accepted move.
		\param pot   Interaction potential to notify of the acceptance. */
	virtual void Accept(Configuration & list, Potential * pot)
	{
		this->AcceptCount++;

		list.MoveParticle(this->moved, this->newcoord); 

		this->AcceptCheck(pot);
	}

	/** \brief Return a heap-allocated copy of this move.
		\return A new clone of this AtomMove. */
	virtual MCMove<Configuration, Potential *> * clone()
	{
		return new AtomMove(*this);
	}
};

/** \brief A MCMove that proposes changes to the simulation-cell basis vectors (isobaric ensemble). */
class CellMove : public MCMove<Configuration, Potential *>
{
private:
	size_t TrialCount, AcceptCount;
	double Pressure;

	GeometryVector newvec[::MaxDimension];
public:
	double sigma;
	double MinDistance;//abandon this move immediately if particles are within this distance
	/** \brief A default constructor (zero pressure). */
	CellMove()
	{
		this->sigma=3e-2;
		this->TrialCount=0;
		this->AcceptCount=0;
		this->Pressure=0;
		this->MinDistance=0;
	}

	/** \brief Construct a cell move at a given pressure.
		\param pressure  Pressure used to form the enthalpy in the acceptance criterion. */
	CellMove(double pressure)
	{
		this->sigma=3e-2;
		this->TrialCount=0;
		this->AcceptCount=0;
		this->Pressure=pressure;
		this->MinDistance=0;
	}

	/** \brief Generate a random change of the cell basis vectors and return the enthalpy difference.
		\param list	Current point configuration.
		\param pot	Interaction potential.
		\param gen	Random number generator.
		\param[out] PreFactor	Acceptance prefactor accounting for the volume change.
		\param LockStepSize	If true, do not auto-adjust the cell step size.
		\param CurrentEnergy	Current enthalpy of the state.
		\return Change in enthalpy from the trial cell move (MaxEnergy if the move is rejected outright). */
	virtual double DeltaEnergy(Configuration & list, Potential * pot, RandomGenerator & gen, double & PreFactor, bool LockStepSize, double CurrentEnergy)//Generate a random move, return Energy difference
	{
		if(this->TrialCount>=::AcceptanceSampleSize && LockStepSize==false)
		{
			double arate=static_cast<double>(this->AcceptCount)/this->TrialCount;//acceptance rate
			//std::cout<<"arate is:"<<arate<<'\n';
			if(arate> ::MaxAcceptance)
			{
				//std::cout<<"Warning in CellMove : step size changed!\n";
				this->sigma*=1.1;
				if(this->sigma>0.3)
					this->sigma=0.3;
				//std::cout<<"Sigma is:"<<this->sigma<<'\n';
			}
			else if(arate < ::MinAcceptance)
			{
				//std::cout<<"Warning in CellMove : step size changed!\n";
				this->sigma*=0.9;
				//std::cout<<"Sigma is:"<<this->sigma<<'\n';
			}
			this->TrialCount=0;
			this->AcceptCount=0;
			//std::cout<<"Sigma is:"<<this->sigma<<'\n';
		}
		if(this->TrialCount==::AcceptanceSampleSize/10)
		{
			double arate=static_cast<double>(this->AcceptCount)/this->TrialCount;//acceptance rate
			//std::cout<<"arate is:"<<arate<<'\n';
			if(arate < ::MinAcceptance/2)
			{
				this->sigma*=0.5;
				this->TrialCount=0;
				this->AcceptCount=0;
				//std::cout<<"Sigma is:"<<this->sigma<<'\n';
			}
		}
		this->TrialCount++;

		auto dim=list.GetDimension();
		bool CellWeired;
		bool TooNearNeighborFound=true;
		Configuration newlist(list);
		size_t ReTrialTime =0;//don't try the following loop too much, it's not terrible even if we don't have a perfect move
		//pot->SetConfiguration(list);
		//double OldH = pot->Energy() + list.PeriodicVolume()*this->Pressure;
		double OldH = CurrentEnergy;
		do
		{
			CellWeired=false;
			double L2Max=0, L2Min=1000000000;
			for(DimensionType i=0; i<dim; i++)
				newvec[i]=list.GetBasisVector(i);
			for(DimensionType i=0; i<dim; i++)
			{
				for(int j=0; j<list.GetDimension(); j++)
				{
					double length=std::sqrt(newvec[i].Modulus2());
                    newvec[i].x[j]+=this->sigma*(gen.RandomDouble()-0.5)*(length);
				}
				double dtemp=newvec[i].Modulus2();
				if(dtemp>L2Max)
					L2Max=dtemp;
				if(dtemp<L2Min)
					L2Min=dtemp;
			}
			double OrigVol=list.PeriodicVolume();
			double NewVol=::Volume(newvec, dim);
			if(std::abs((NewVol-OrigVol)/OrigVol) > this->sigma )
				CellWeired=true;
			if( (L2Max/L2Min)> (MaxAspectRatioOfCell)*(MaxAspectRatioOfCell) )
				CellWeired=true;
			if(L2Max > MaxCellSideLengthSquared)
				CellWeired=true;

			//search for too near neighbors
			newlist.ChangeBasisVector(newvec);
			TooNearNeighborFound=false;
			for(size_t i=0; i<newlist.NumParticle(); i++)
			{
				//Configuration::particle * pa=newlist.GetParticle(i);
				newlist.IterateThroughNeighbors(newlist.GetRelativeCoordinates(i), this->MinDistance, [&TooNearNeighborFound, &i](const GeometryVector & /*shift*/, const GeometryVector & LatticeShift, const signed long * /*PeriodicShift*/, const size_t SourceAtom) -> void 
				{
					if(SourceAtom!=i || LatticeShift.Modulus2()!=0)
						TooNearNeighborFound=true;
				}, &TooNearNeighborFound);
			}
			//new_idea
			ReTrialTime++;
			if( (TooNearNeighborFound||CellWeired) && ReTrialTime>100 )
				return ::MaxEnergy;
		}
		while( (CellWeired || TooNearNeighborFound) );

		double NewH;
		// newlist is a local trial configuration that dies on return. Once the potential is bound
		// to it, EVERY exit from here (normal return, the OOM abandon path, or any other exception
		// propagating out) must rebind it to the caller's live configuration `list`. An RAII guard
		// does this on destruction; it is never dismissed because the move is not applied yet, so
		// DeltaEnergy always leaves the potential bound to `list` (Accept rebinds after applying it).
		struct RebindGuard
		{
			Potential * pot;
			const Configuration * live;
			bool dismissed;
			~RebindGuard() { if(!dismissed) pot->SetConfiguration(*live); }
		} rebindGuard{pot, &list, false};
		try
		{
			pot->SetConfiguration(newlist);
			NewH = pot->Energy() + newlist.PeriodicVolume()*this->Pressure;
		}
		catch(std::bad_alloc & a)
		{
			std::cerr<<"Warning : not enough memory in cell resize, abandon move!";
			return ::MaxEnergy;//rebindGuard restores the binding
		}

		PreFactor=std::pow(newlist.PeriodicVolume()/list.PeriodicVolume(), (double)(list.NumParticle()));

		double DeltaE= NewH-OldH;
		//std::cout <<" \t"<< PreFactor << " \t" << DeltaE << '\n';
		return DeltaE;//rebindGuard restores the binding
	}

	/** \brief Accept the last trial cell move, applying the new basis vectors.
		\param list  Configuration to update with the accepted basis vectors.
		\param pot   Interaction potential; rebound to the accepted configuration. */
	virtual void Accept(Configuration & list, Potential * pot)
	{
		list.ChangeBasisVector(newvec);
		list.TryRefineBasisVectors();
		pot->SetConfiguration(list);//the basis just changed; rebind so the potential tracks the accepted configuration
		this->AcceptCount++;
	}

	/** \brief Return a heap-allocated copy of this move.
		\return A new clone of this CellMove. */
	virtual MCMove<Configuration, Potential *> * clone()
	{
		return new CellMove(*this);
	}
};


/** \brief A MCMove that randomly selects between two sub-moves with a fixed probability. */
template <typename System, typename AdditionalInfo>
class CombinedMove : public MCMove<System, AdditionalInfo>
{
	std::unique_ptr<MCMove<System, AdditionalInfo> > m1;
	std::unique_ptr<MCMove<System, AdditionalInfo> > m2;
	double Move1Probability;
	RandomGenerator gen;
	bool IsMoving1;

public:
	/** \brief Construct from two moves and the probability of choosing the first.
		\param m1               First sub-move; cloned and owned by this object.
		\param m2               Second sub-move; cloned and owned by this object.
		\param Move1Probability  Probability in [0,1] of attempting the first move.
		\param RandomSeed       Seed for the internal random number generator. */
	CombinedMove(MCMove<System, AdditionalInfo> & m1, MCMove<System, AdditionalInfo> & m2, double Move1Probability, int RandomSeed) : m1(m1.clone()), m2(m2.clone()), Move1Probability(Move1Probability), gen(RandomSeed)
	{
		assert(Move1Probability>=0);
		assert(Move1Probability<=1);
	}

	/** \brief Destructor; the owned sub-moves are freed automatically by their unique_ptrs. */
	~CombinedMove() = default;

	/** \brief Randomly select one sub-move, perform it, and return the energy difference.
		\param sys	Current system state.
		\param add	Additional information (e.g. the interaction potential).
		\param gen	Random number generator; also selects which sub-move is attempted (the member generator is only used to seed clones).
		\param[out] PreFactor	Acceptance prefactor set by the selected sub-move.
		\param LockStepSize	If true, do not auto-adjust step sizes.
		\param CurrentEnergy	Energy of the current state.
		\return Change in energy from the selected sub-move. */
	virtual double DeltaEnergy(System & sys, AdditionalInfo add, RandomGenerator & gen, double & PreFactor, bool LockStepSize, double CurrentEnergy)//Generate a random move, return Energy difference
	{
		if(gen.RandomDouble()<this->Move1Probability)
		{
			this->IsMoving1=true;
			return m1->DeltaEnergy(sys, add, gen, PreFactor, LockStepSize, CurrentEnergy);
		}
		else
		{
			this->IsMoving1=false;
			return m2->DeltaEnergy(sys, add, gen, PreFactor, LockStepSize, CurrentEnergy);
		}
	}

	/** \brief Accept the most recently selected sub-move.
		\param sys  System to update.
		\param add  Additional information (e.g. the interaction potential). */
	virtual void Accept(System & sys, AdditionalInfo add)
	{
		if(this->IsMoving1==true)
			this->m1->Accept(sys, add);
		else
			this->m2->Accept(sys, add);
	}

	/** \brief Return a heap-allocated copy of this move.
		\return A new clone of this CombinedMove. */
	virtual MCMove<System, AdditionalInfo> * clone()
	{
		return new CombinedMove(*this->m1, *this->m2, this->Move1Probability, this->gen.RandomInt());//full-width child seed (was RandomDouble()*10000, only 10^4 values)
	}
};


#endif
