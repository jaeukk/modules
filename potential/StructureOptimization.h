/** author: Ge Zhnag
	date:
	email:		
	comments: Jaeuk Kim (June 2019).	*/

/**	\file StructureOptimization.h
	\brief Header file to perform energy minimization... */


#ifndef STRUCTUREOPTIMIZATION_INCLUDED
#define STRUCTUREOPTIMIZATION_INCLUDED

#include <cassert>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <nlopt.h>
/* in $(cores) */
#include <PeriodicCellList.h>

#include "Potential.h"

//default: 0
//if not zero, then partially optimized configurations are saved every these steps
//will save to two files: "StructureOptimization_Progress_a.configuration" 
//and "StructureOptimization_Progress_b.configuration" alternate
extern size_t StructureOptimization_SaveConfigurationInterval;
extern bool StructureOptimization_SaveNumberedConfiguration;

/** \brief Exception thrown when a NaN is encountered during structure optimization. (CLAUDE) */
class NotANumberFound : public std::exception
{
public:
	/** \brief Construct the exception. (CLAUDE) */
	NotANumberFound()
	{
	}
	/** \brief Return the exception's descriptive message. (CLAUDE)
		\return A C-string describing the error. */
	virtual const char * what() const throw()
	{
		return "Caught a NaN!\n";
	}
};

/** \brief Minimize the potential energy of a given initial configuration via MINOP algorithm.
 *	\param[in,out]	List	An initial configuration, and it becomes optimized.
 *	\param[in]		pot		A refernece of the potential energy.
 *	\param[in]		Pressure	A constant pressure in the isobaric ensembles.
 *	\param[in]		Switch	Determine which one is being optimized;
					0:	particle positions.
					1:	fundamental cell
					2:	both particle positions and fundamental cell
 *	\param[in]		MinDistance	Particle displacement is smaller than this value, the movement is rejected.
 *	\param[in]		MaxStep	The maximal number of optimization steps.
 *	\param[in]		Emin	Target energy below which optimization stops early. (CLAUDE) */
void RelaxStructure_MINOP_withoutPertubation(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=20000, double Emin=-1e10);

/** \brief do MINOP relax without pertubation, then with pertubation, return the best one.
	\param[in,out] List       Initial configuration; replaced by the optimized one. (CLAUDE)
	\param pot                Interaction potential to minimize. (CLAUDE)
	\param Pressure           Constant pressure for the isobaric ensemble. (CLAUDE)
	\param Switch             0: positions, 1: cell, 2: both. (CLAUDE)
	\param MinDistance        Reject moves that bring particles closer than this distance. (CLAUDE)
	\param MaxStep            Maximum number of optimization steps. (CLAUDE) */
void RelaxStructure_MINOP(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=20000);

//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
/** \brief Minimize the energy using an NLopt gradient-based optimizer. (CLAUDE)
	\param[in,out] List       Initial configuration; replaced by the optimized one.
	\param pot                Interaction potential to minimize.
	\param Pressure           Constant pressure for the isobaric ensemble.
	\param Switch             0: move atoms, 1: move basis vectors, 2: move both.
	\param MinDistance        Reject moves that bring particles closer than this distance.
	\param MaxStep            Maximum number of optimization steps. */
void RelaxStructure_NLOPT(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);
//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
/** \brief Minimize the energy using an NLopt optimizer, stopping when the energy falls below Emin. (CLAUDE)
	\param[in,out] List       Initial configuration; replaced by the optimized one.
	\param pot                Interaction potential to minimize.
	\param Pressure           Constant pressure for the isobaric ensemble.
	\param Switch             0: move atoms, 1: move basis vectors, 2: move both.
	\param MinDistance        Reject moves that bring particles closer than this distance.
	\param MaxStep            Maximum number of optimization steps.
	\param Emin               Target energy below which optimization stops early. */
void RelaxStructure_NLOPT_Emin(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double Emin);
//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
/** \brief Minimize the energy using steepest-descent. (CLAUDE)
	\param[in,out] List       Initial configuration; replaced by the optimized one.
	\param pot                Interaction potential to minimize.
	\param Pressure           Constant pressure for the isobaric ensemble.
	\param Switch             0: move atoms, 1: move basis vectors, 2: move both.
	\param MinDistance        Reject moves that bring particles closer than this distance.
	\param MaxStep            Maximum number of optimization steps. */
void RelaxStructure_SteepestDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);
/** \brief Minimize the energy using conjugate-gradient descent. (CLAUDE)
	\param[in,out] List       Initial configuration; replaced by the optimized one.
	\param pot                Interaction potential to minimize.
	\param Pressure           Constant pressure for the isobaric ensemble.
	\param Switch             0: move atoms, 1: move basis vectors, 2: move both.
	\param MinDistance        Reject moves that bring particles closer than this distance.
	\param MaxStep            Maximum number of optimization steps.
	\param minG               Gradient-norm threshold for convergence. */
void RelaxStructure_ConjugateGradient(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000, double minG=1e-13);
/** \brief Minimize the energy using a local gradient-descent scheme. (CLAUDE)
	\param[in,out] List       Initial configuration; replaced by the optimized one.
	\param pot                Interaction potential to minimize.
	\param Pressure           Constant pressure for the isobaric ensemble.
	\param Switch             0: move atoms, 1: move basis vectors, 2: move both.
	\param MinDistance        Reject moves that bring particles closer than this distance.
	\param MaxStep            Maximum number of optimization steps. */
void RelaxStructure_LocalGradientDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);

//slower but does not need derivative information
/** \brief Minimize the energy using a derivative-free NLopt optimizer (slower but needs no gradients). (CLAUDE)
	\param[in,out] List       Initial configuration; replaced by the optimized one.
	\param pot                Interaction potential to minimize.
	\param Pressure           Constant pressure for the isobaric ensemble.
	\param Switch             0: move atoms, 1: move basis vectors, 2: move both.
	\param MinDistance        Reject moves that bring particles closer than this distance.
	\param rescale            Scale factor applied to the optimization variables.
	\param MaxStep            Maximum number of optimization steps.
	\param fmin               Target energy below which optimization stops early. */
void RelaxStructure_NLOPT_NoDerivative(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, double rescale=1.0, size_t MaxStep=100000, double fmin=(-1.0)*MaxEnergy);

/** \brief Relax a structure (positions and cell), falling back to derivative-free relaxation if cell derivatives are unavailable. (CLAUDE)
	\param[in,out] List   Initial configuration; replaced by the relaxed one.
	\param pot            Interaction potential to minimize.
	\param Pressure       Constant pressure for the isobaric ensemble.
	\param MinDistance    Reject moves that bring particles closer than this distance. */
inline void RelaxStructure(Configuration & List, Potential & pot, double Pressure, double MinDistance)
{
	try
	{
		RelaxStructure_NLOPT(List, pot, Pressure, 2, MinDistance);
	}
	catch(EnergyDerivativeToBasisVectors_NotImplemented a)
	{
		for(int i=0; i<5; i++)
		{
			RelaxStructure_NLOPT(List, pot, Pressure, 0, MinDistance);
			RelaxStructure_NLOPT_NoDerivative(List, pot, Pressure, 1, MinDistance);
		}
		RelaxStructure_NLOPT_NoDerivative(List, pot, Pressure, 2, MinDistance);
	}
}

//calculate Elastic Constants
//relax the structure before calling this function
//\epsilon_{ij}=C_{ijkl}*e_{kl}
/** \brief calculate Elastic Constants; relax the structure before calling this function.
	\param pot          Interaction potential. (CLAUDE)
	\param structure    Relaxed configuration to probe. (CLAUDE)
	\param i            First index of the elastic-constant tensor component. (CLAUDE)
	\param j            Second index of the elastic-constant tensor component. (CLAUDE)
	\param k            Third index of the elastic-constant tensor component. (CLAUDE)
	\param l            Fourth index of the elastic-constant tensor component. (CLAUDE)
	\param Pressure     Constant pressure for the isobaric ensemble. (CLAUDE)
	\param InfiniteTime  If true, use the (infinite-time / fully relaxed) elastic response. (CLAUDE)
	\return The elastic constant C_{ijkl}. (CLAUDE) */
double ElasticConstant(Potential & pot, Configuration structure, DimensionType i, DimensionType j, DimensionType k, DimensionType l, double Pressure, bool InfiniteTime=false);

//if presult is not nullptr, then presult[0..dimension^4] will be filled with C_{ijkl}
/** \brief Print all elastic constants; if presult is not nullptr, then presult[0..dimension^4] will be filled with C_{ijkl}.
	\param out          Output stream to print the elastic constants to. (CLAUDE)
	\param stru         Relaxed configuration to probe. (CLAUDE)
	\param pot          Interaction potential. (CLAUDE)
	\param Pressure     Constant pressure for the isobaric ensemble. (CLAUDE)
	\param InfiniteTime  If true, use the (infinite-time / fully relaxed) elastic response. (CLAUDE)
	\param[out] presult  If non-null, filled with the full C_{ijkl} tensor (dimension^4 entries). (CLAUDE) */
void PrintAllElasticConstants(std::ostream & out, Configuration stru, Potential & pot, double Pressure, bool InfiniteTime = false, double * presult=nullptr);


//optimize the parameters for a elastic model and rotations so that it fits the target
/** \brief optimize the parameters for a elastic model and rotations so that it fits the target.
	\param getElasticFunc  Function mapping model parameters to the predicted elastic response. (CLAUDE)
	\param nparam          Number of model parameters. (CLAUDE)
	\param[in,out] param   Initial parameters; replaced by the optimized values. (CLAUDE)
	\param ubound          Upper bounds for the parameters. (CLAUDE)
	\param lbound          Lower bounds for the parameters. (CLAUDE)
	\param c               Reference configuration. (CLAUDE)
	\param pot             Interaction potential. (CLAUDE)
	\param pressure        Constant pressure for the isobaric ensemble. (CLAUDE) */
void ElasticOptimize(std::function<BigVector(const std::vector<double> & param)> getElasticFunc, int nparam, std::vector<double> & param, std::vector<double> ubound, std::vector<double> lbound, const Configuration & c, Potential & pot, double pressure);
#endif
