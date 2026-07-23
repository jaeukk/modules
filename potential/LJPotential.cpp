/** Author: Ge Zhang
	email:	?
	Date:	?
Comments:	Jaeuk Kim (June 2019).*/

/** \file LJPotential.cpp
 *	\brief Functional impelementations of the LJPotential class. */

#include "LJPotential.h"
#include <cmath>
#include <cstring>
	LJPotential::~LJPotential()
	{
	}

	LJPotential::LJPotential(DimensionType Dimension, double epsilon, double sigma, double Rcut) :IsotropicPairPotential(Dimension, Rcut)
	{
		this->epsilon4 = 4 * epsilon;
		this->sigma = sigma;
	}

	inline double LJPotential::IsotropicEnergyFormula(double distance)
	{
		double a=this->sigma*this->sigma/distance/distance;
		a=a*a*a;
		return this->epsilon4*(a*a-a);
	}

	inline double LJPotential::IsotropicForceFormula(double distance)
	{
		double term6 = this->sigma / distance ; 
		term6 = term6 * term6 * term6 * term6 * term6 * term6; //term6 = (sigma / distance)^6
		return (-1.0)*this->epsilon4* (6.0*term6 - 12.0*term6*term6) / distance;
	}
