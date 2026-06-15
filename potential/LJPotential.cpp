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
		double temp=this->sigma/distance;
		return (-1.0)*this->epsilon4*(6/std::pow(distance, (double)(7))-12/std::pow(distance, (double)(13)));
	}
