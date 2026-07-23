/**	Author	: Ge Zhang
 *	Email	:
 *	Date	:
 *  Annotated by Jaeuk Kim (2019). */

/** \file CoolingSchedule.h
 *	\brief Cooling schedules for simulated annealing: the abstract
 *	CoolingSchedule plus the ExpCool / ZeroCool / ThermoCool implementations.
 *	Split out of MonteCarlo.h (2026-07-07) so that annealing protocols can be
 *	used (e.g. by MD_System) without pulling in the Monte Carlo drivers. */
#ifndef COOLINGSCHEDULE_INCLUDED
#define COOLINGSCHEDULE_INCLUDED

#include <cmath>
#include <cstddef>
#include <vector>

#include "AnalyzeData.h"

/**	\class	An abstract class to define and control the cooling schedule. */
class CoolingSchedule
{
public:
	bool Continue;	//< true = continue cooling the system.
	double Temperature;	//< Current temperature.
	size_t NumTry;		//<

	/** \brief Report the current energy so the schedule can advance to the next cooling stage.
		\param Energy  Current (average) energy reported by the Monte Carlo driver. */
	virtual void Report(double Energy) = 0;
	//when initialized, the variables are set to the 1st stage cooling.
	//after class SimulatedAnnealing performs these MC trials, it will call Report to report it's current energy
	//the class need to decide if another stage is needed, and update the variables
};


const size_t StageSize=1000;

/** \brief A class to define and control the exponential cooling schedule. Inherited from the CoolingSchedule class. */
class ExpCool : public CoolingSchedule
{
	double StageAlpha;	//<
	double Tcut;		//< The final temperature, below which the system is cooled to the zero temperature.
	signed long ZeroMove;
public:
	/** \brief Construct an exponential cooling schedule.
		\param Alpha     Per-step cooling factor (temperature multiplier raised to the stage size).
		\param Tinit     Initial temperature.
		\param Tcut      Final temperature below which the system is cooled to zero temperature.
		\param ZeroMove  Number of trial moves to perform at zero temperature. */
	ExpCool(double Alpha, double Tinit, double Tcut, long ZeroMove) : Tcut(Tcut), ZeroMove(ZeroMove)
	{
		if(this->ZeroMove<0)
			this->ZeroMove=0;//a negative count would wrap to a huge value in the (size_t) comparison below
		this->StageAlpha=std::pow(Alpha, static_cast<double>(StageSize));
		this->Continue=true;
		this->Temperature=Tinit;
		this->NumTry=StageSize;
	}

	/** \brief Advance the temperature exponentially after a cooling stage.
		\param Energy  Current energy (unused by this schedule). */
	virtual void Report(double /*Energy*/)
	{
		//////////////////////////
		if(this->Temperature==0)
		{
			this->ZeroMove-=this->NumTry;
			if(this->ZeroMove<=0)
				this->Continue=false;
		}
		else
		{
			this->Temperature*=this->StageAlpha;
			if(this->Temperature<this->Tcut)
			{
				this->Temperature=0;
				if(this->NumTry>(size_t)this->ZeroMove)
					this->NumTry=this->ZeroMove;
			}
		}
	}
};



const size_t ReportTime=100;
/** \brief A cooling schedule that holds the temperature at zero for a fixed number of reports. */
class ZeroCool : public CoolingSchedule
{
	size_t RemainingMoves;//zero-temperature trial moves still to perform

	/** \brief Number of moves for the next stage: the remaining total spread over ReportTime stages,
		rounded up so small totals still perform moves and no remainder is dropped.
		\return Moves to attempt in the next stage (never more than RemainingMoves). */
	size_t NextStage(void) const
	{
		if(this->RemainingMoves==0)
			return 0;
		size_t n=(this->RemainingMoves + ReportTime - 1)/ReportTime;
		return (n>this->RemainingMoves) ? this->RemainingMoves : n;
	}

public:
	/** \brief Construct a zero-temperature schedule.
		\param ZeroMoveTime  Total number of zero-temperature trial moves, split across the reports. */
	ZeroCool(long ZeroMoveTime)
	{
		this->Temperature=0.0;
		if(ZeroMoveTime<0)
			ZeroMoveTime=0;
		this->RemainingMoves=static_cast<size_t>(ZeroMoveTime);
		this->NumTry=this->NextStage();
		this->Continue=(this->RemainingMoves>0);
	}

	/** \brief Retire the moves just performed and schedule the next stage; stop once the total is exhausted.
		\param Energy  Current energy (unused by this schedule). */
	virtual void Report(double /*Energy*/)
	{
		this->RemainingMoves-=this->NumTry;
		this->NumTry=this->NextStage();
		if(this->RemainingMoves==0)
			this->Continue=false;
	}
};


const size_t DefaultSampleInterval = 100;
const double DefaultVsCoeff=500;
const size_t DefaultInitAnalyzeSize=200;
class ThermoCool : public CoolingSchedule
{
	size_t SampleInterval;
	double VsCoeff;
	size_t InitAnalyzeSize;

	double Vs;
	double Tcut;
	signed long ZeroMove;//deprecated parameter

	size_t AnalyzeSize;

	std::vector<double> data;
public:
	/** \brief Construct an adaptive (thermodynamic) cooling schedule driven by the measured heat capacity.
		\param Alpha            Reference per-step cooling factor used to set the entropy step Vs.
		\param Tinit            Initial temperature.
		\param Tcut             Final temperature below which the system is cooled to zero temperature.
		\param ZeroMove         Number of zero-temperature trial moves (deprecated parameter).
		\param SampleInterval   Number of trial moves between energy samples.
		\param VsCoeff          Coefficient setting the target entropy decrease per cooling step.
		\param InitAnalyzeSize  Initial number of samples collected before analyzing the energy series. */
	ThermoCool(double Alpha, double Tinit, double Tcut, long ZeroMove,
		size_t SampleInterval=DefaultSampleInterval, double VsCoeff=DefaultVsCoeff,
		size_t InitAnalyzeSize=DefaultInitAnalyzeSize) : SampleInterval(SampleInterval),
		VsCoeff(VsCoeff), InitAnalyzeSize(InitAnalyzeSize), Tcut(Tcut), ZeroMove(ZeroMove)
	{
		if(this->ZeroMove<0)
			this->ZeroMove=0;//a negative count would wrap huge when assigned to the size_t NumTry
		this->Vs=(1.0-Alpha)*SampleInterval*VsCoeff;
		this->Continue=true;
		this->Temperature=Tinit;
		this->NumTry=SampleInterval;

		this->AnalyzeSize=InitAnalyzeSize;
	}

	/** \brief Sample the energy and adaptively lower the temperature based on the estimated heat capacity.
		\param Energy  Current energy sample to add to the analysis buffer. */
	virtual void Report(double Energy)
	{
		//////////////////////////
		if(this->Temperature==0)
		{
			this->Continue=false;
		}
		else
		{
			data.push_back(Energy);

			if(data.size()>=AnalyzeSize)
			{
				double Average, error, Corrtime;
				AnalyzeData(data, Average, error, Corrtime);
				double H2c=error*error*data.size()/(1+2*Corrtime);//this is T^2*C_v
				double dlnT=this->Vs*this->Temperature/Corrtime/std::sqrt(H2c);

					//debug temp
					//std::fstream file("ErrorLog.txt", std::fstream::out);
					//for(auto iter=data.begin(); iter!=data.end(); iter++)
					//	file<<(*iter)<<'\n';
					//file.close();
					//exit(0);
				if( (!(dlnT>0)) )
				{
					if(this->AnalyzeSize < 5*InitAnalyzeSize)
					{
						this->Vs*=2;
						this->AnalyzeSize*=2;
					}
					this->Temperature/=std::exp(0.000001);//cool very cautiously when can't analyze data
				}
				else if(dlnT>0.1)
				{
					//std::cout<<"Warning in ThermoCool : Temperature dropping too quickly, decrease sample size.\n";
					if(this->AnalyzeSize>InitAnalyzeSize/2 && this->AnalyzeSize>1)
					{
						this->Vs/=2;
						this->AnalyzeSize/=2;
					}
					this->Temperature/=std::exp(0.1);
				}
				else
				{
					this->Temperature/=std::exp(dlnT);
				}

				data.clear();
			}

			if(this->Temperature<this->Tcut)
			{
				this->Temperature=0;
				if(this->ZeroMove==0)
					this->Continue=false;
				else
					this->NumTry=this->ZeroMove;
			}
		}
	}
};

#endif
