/** Author: Ge Zhang
	email:	?
	Date:	?
	Comments:	Jaeuk Kim (June 2019).*/

/** \file LJPotential.h
 *	\brief Header file to define the Lennard-Jone potential. */

#ifndef LJPOTENTIAL_INCLUDED
#define LJPOTENTIAL_INCLUDED

#include "Potential.h"

/** \brief The Lennard-Jones potential. 
	It is derived from IsotropicPairPotential class. */
class LJPotential: public IsotropicPairPotential
{
public:
	double epsilon4;//< 4*epsilon
	double sigma;	//< Length scale.

	/** \brief Destructor. */
	virtual ~LJPotential();

	/** \brief Construct a Lennard-Jones potential.
		\param Dimension Spatial dimension of the system.
		\param epsilon Potential-well depth.
		\param sigma Lennard-Jones length scale.
		\param Rcut Pair-interaction cutoff distance. */
	LJPotential(DimensionType Dimension, double epsilon, double sigma, double Rcut);

	/** \brief LJ potential.
		\param distance	Separation.
		\return energy. */
	virtual double IsotropicEnergyFormula(double distance);

	/** \brief Force via the LJ potential.
		\param distance	Separation.
		\return the magnitude of force.	 */
	virtual double IsotropicForceFormula(double distance);
};

#endif