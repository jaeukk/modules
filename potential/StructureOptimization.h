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

/** \brief Exception thrown when a NaN is encountered during structure optimization. */
class NotANumberFound : public std::exception
{
public:
	/** \brief Construct the exception. */
	NotANumberFound()
	{
	}

	/** \brief Return the exception's descriptive message.
		\return A C-string describing the error. */
	virtual const char * what() const throw()
	{
		return "Caught a NaN!\n";
	}
};

/** \brief Minimize the configuration enthalpy with MINOP without adding an initial perturbation.
	\param[in,out] List Initial configuration, replaced by the optimized configuration.
	\param pot Interaction potential.
	\param Pressure Pressure used in the enthalpy.
	\param Switch Zero optimizes positions, one the cell, and two both.
	\param MinDistance Minimum permitted interparticle distance.
	\param MaxStep Maximum number of MINOP iterations.
	\param Emin Objective value below which optimization stops early. */
void RelaxStructure_MINOP_withoutPertubation(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=20000, double Emin=-1e10);

/** \brief Run MINOP from the original configuration and from a slightly perturbed copy, then retain the configuration with lower enthalpy.
	\param[in,out] List Initial configuration, replaced by the better optimized configuration.
	\param pot Interaction potential.
	\param Pressure Pressure used in the enthalpy.
	\param Switch Zero optimizes positions, one the cell, and two both.
	\param MinDistance Minimum permitted interparticle distance.
	\param MaxStep Maximum iterations for each MINOP run. */
void RelaxStructure_MINOP(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=20000);

//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
/** \brief Minimize the enthalpy using three NLopt L-BFGS passes with different relative-coordinate scales.
	\param[in,out] List Initial configuration, replaced by the optimized configuration.
	\param pot Interaction potential.
	\param Pressure Pressure used in the enthalpy.
	\param Switch Zero optimizes positions, one basis vectors, and two both.
	\param MinDistance Minimum permitted interparticle distance.
	\param MaxStep Maximum objective evaluations per pass. */
void RelaxStructure_NLOPT(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);
//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
/** \brief Run scaled NLopt L-BFGS passes using Emin as the objective stop value, skipping later passes when the potential energy is already below Emin.
	\param[in,out] List Initial configuration, replaced by the optimized configuration.
	\param pot Interaction potential.
	\param Pressure Pressure used in the objective enthalpy.
	\param Switch Zero optimizes positions, one basis vectors, and two both.
	\param MinDistance Minimum permitted interparticle distance.
	\param MaxStep Maximum objective evaluations per pass.
	\param Emin Stop value supplied to each optimization pass. */
void RelaxStructure_NLOPT_Emin(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double Emin);
//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
/** \brief Minimize the enthalpy using GSL steepest descent.
	\param[in,out] List Initial configuration, replaced by the optimized configuration.
	\param pot Interaction potential.
	\param Pressure Pressure used in the enthalpy.
	\param Switch Zero optimizes positions, one basis vectors, and two both.
	\param MinDistance Minimum permitted interparticle distance.
	\param MaxStep Maximum number of iterations. */
void RelaxStructure_SteepestDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);

/** \brief Minimize the enthalpy using GSL Polak-Ribiere conjugate-gradient descent.
	\param[in,out] List Initial configuration, replaced by the optimized configuration.
	\param pot Interaction potential.
	\param Pressure Pressure used in the enthalpy.
	\param Switch Zero optimizes positions, one basis vectors, and two both.
	\param MinDistance Minimum permitted interparticle distance.
	\param MaxStep Maximum number of iterations.
	\param minG Gradient-norm convergence threshold. */
void RelaxStructure_ConjugateGradient(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000, double minG=1e-13);

/** \brief Minimize the enthalpy using the local normalized-gradient descent implementation.
	\param[in,out] List Initial configuration, replaced by the optimized configuration.
	\param pot Interaction potential.
	\param Pressure Pressure used in the enthalpy.
	\param Switch Zero optimizes positions, one basis vectors, and two both.
	\param MinDistance Minimum permitted interparticle distance.
	\param MaxStep Objective-evaluation budget. */
void RelaxStructure_LocalGradientDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);

//slower but does not need derivative information
/** \brief Minimize the enthalpy using NLopt's derivative-free Subplex optimizer.
	\param[in,out] List Initial configuration, replaced by the optimized configuration.
	\param pot Interaction potential.
	\param Pressure Pressure used in the enthalpy.
	\param Switch Zero optimizes positions, one basis vectors, and two both.
	\param MinDistance Minimum permitted interparticle distance.
	\param rescale Scale applied to relative-coordinate variables.
	\param MaxStep Maximum number of objective evaluations.
	\param fmin Objective stop value. */
void RelaxStructure_NLOPT_NoDerivative(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, double rescale=1.0, size_t MaxStep=100000, double fmin=(-1.0)*MaxEnergy);

/** \brief Relax a structure (positions and cell), falling back to derivative-free relaxation if cell derivatives are unavailable.
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
	catch(EnergyDerivativeToBasisVectors_NotImplemented & a)
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
	\param pot          Interaction potential.
	\param structure    Relaxed configuration to probe.
	\param i            First index of the elastic-constant tensor component.
	\param j            Second index of the elastic-constant tensor component.
	\param k            Third index of the elastic-constant tensor component.
	\param l            Fourth index of the elastic-constant tensor component.
	\param Pressure     Constant pressure for the isobaric ensemble.
	\param InfiniteTime  If true, use the (infinite-time / fully relaxed) elastic response.
	\return The elastic constant C_{ijkl}. */
double ElasticConstant(Potential & pot, Configuration structure, DimensionType i, DimensionType j, DimensionType k, DimensionType l, double Pressure, bool InfiniteTime=false);

//if presult is not nullptr, then presult[0..dimension^4] will be filled with C_{ijkl}
/** \brief Print every elastic-constant component and optionally store the dimension^4 values in a flat array.
	\param out Output stream.
	\param stru Relaxed configuration to probe.
	\param pot Interaction potential.
	\param Pressure Pressure used in the enthalpy.
	\param InfiniteTime If true, relax internal particle coordinates after applying the strain.
	\param[out] presult If non-null, receives entries 0 through dimension^4-1 in i,j,k,l order. */
void PrintAllElasticConstants(std::ostream & out, Configuration stru, Potential & pot, double Pressure, bool InfiniteTime = false, double * presult=nullptr);


//optimize the parameters for a elastic model and rotations so that it fits the target
/** \brief optimize the parameters for a elastic model and rotations so that it fits the target.
	\param getElasticFunc  Function mapping model parameters to the predicted elastic response.
	\param nparam          Number of model parameters.
	\param[in,out] param   Initial parameters; replaced by the optimized values.
	\param ubound          Upper bounds for the parameters.
	\param lbound          Lower bounds for the parameters.
	\param c               Reference configuration.
	\param pot             Interaction potential.
	\param pressure        Constant pressure for the isobaric ensemble. */
void ElasticOptimize(std::function<BigVector(const std::vector<double> & param)> getElasticFunc, int nparam, std::vector<double> & param, std::vector<double> ubound, std::vector<double> lbound, const Configuration & c, Potential & pot, double pressure);
#endif
