#include "PhononFrequency.h"

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_eigen.h>

#include <cmath>
#include <algorithm>
#include <iostream>

void PhononFrequency(const GeometryVector & k, const Configuration & crystal, PairPotential & pot, std::vector<double> & result2, double Mass)
{
	result2.clear();
	size_t dim=crystal.GetDimension();
	size_t nbr=crystal.NumParticle();
	if(nbr==0)   // gsl_matrix_complex_calloc(0,0) and the eigen-solvers cannot handle an empty crystal
	{
		std::cerr<<"Error in PhononFrequency: empty crystal (no particles); returning no frequencies.\n";
		return;
	}
	if(!(Mass>0.0) || !std::isfinite(Mass))   // Mass divides every eigenvalue below
	{
		std::cerr<<"Error in PhononFrequency: Mass must be a finite positive number; returning no frequencies.\n";
		return;
	}
	//Configuration::particle * start = crystal.GetParticle(0);

	gsl_matrix_complex * mat = gsl_matrix_complex_calloc(nbr*dim, nbr*dim);

	for(size_t i=0; i<nbr; i++)
	{
		//Configuration::particle * p1=crystal.GetParticle(i);
		crystal.IterateThroughNeighbors(crystal.GetRelativeCoordinates(i), pot.Rcut,
			[&k, &dim, &mat, &i, &pot] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long *, const size_t srcAtom) ->void
			{
				if(shift.Modulus2()< ::LengthPrecision* ::LengthPrecision)
					return;
				size_t j=srcAtom;
				double ex=std::cos(k.Dot(LatticeShift));
				double iex=std::sin(k.Dot(LatticeShift));
				for(size_t m=0; m<dim; m++)
				{
					for(size_t n=0; n<dim; n++)
					{
						(*gsl_matrix_complex_ptr(mat, i*dim+m, j*dim+n)).dat[0]+=ex* pot.SecondDerivative(shift, m, n);
						(*gsl_matrix_complex_ptr(mat, i*dim+m, j*dim+n)).dat[1]+=iex* pot.SecondDerivative(shift, m, n);
						(*gsl_matrix_complex_ptr(mat, i*dim+m, i*dim+n)).dat[0]-= pot.SecondDerivative(shift, m, n);
					}
				}
			}
		);
	}

	gsl_eigen_herm_workspace * wor = gsl_eigen_herm_alloc(dim*nbr);
	gsl_vector * eval = gsl_vector_alloc(dim*nbr);

	gsl_eigen_herm(mat, eval, wor);
	result2.clear();
	for(size_t i=0; i<eval->size; i++)
		result2.push_back((-1)*gsl_vector_get(eval, i)/Mass);
	std::sort(result2.begin(), result2.end());

	gsl_vector_free(eval);
	gsl_eigen_herm_free(wor);
	gsl_matrix_complex_free(mat);
}

