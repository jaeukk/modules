/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	Feb. 2019 */

/** \file FormFactors.cpp
 * \brief Function implementations of computing form factor (Fourier transform) of polygon shapes.
 * Definitions of GetFormFactors(const GeometryVector &)
 * for	1. PlanarFigure
 *		2. Polyhedron
 * 	In our purpose, the analytic expressions are accurate enough.
 *	We didn't implement the numerical calculations involving with series expansions.
 *	
 *	Algorithm is based on
 *	J. Wuttke, Form factor (Fourier shape transform) of polygon and polyhedron, arXiv:1703.00255, (2017). */

#include "ParticleShapes.h"
#include <GeometryVector.h>
#include <complex>

const double c_parallel	= 1e-10;	/// \var	When k.r_parallel < c, a series expansion formula (PlanarFigure) should be implemented.
const double c			= 1e-10;	/// \var	When k < c, a series expansion formula (Polyhedron) should be implmented.	
const double MachinePrecision = 1e-16;	/// \var	Define a machine precision.

/** Sinc function. 
 * \param x a real number
 * \return a value of sinc function */
inline double sinc(double x) {
	return ((x == 0.0) ? 1.0 : sin(x) / x);
}

/** \brief	Compute the form factor via an analytic expression. 
 * Fourier transform by exp(-i k x).
 * \param k	A wavevector in 2D or 3D.
 * \return A complex form factor of this PlanarFigure.	
 */
std::complex<double> PlanarFigure::FormFactor_analytic(const GeometryVector &k)
{
	GeometryVector k_perp = (this->Dimension == 2) ? GeometryVector((int)2) : k.Dot(this->normal)*this->normal;
	GeometryVector k_parallel = k - k_perp;
	GeometryVector Kx;
	if (this->Dimension == 2) {
		Kx = GeometryVector(-k_parallel.x[1], k_parallel.x[0]);
		Kx = (this->normal.x[2] > 0) ? Kx : (-1.0)*Kx;
	}
	else if (this->Dimension == 3) {
		Kx = this->normal.Cross(k_parallel);
	}

	
	double K_perp = sqrt(k_perp.Modulus2());
	double K_parallel = sqrt(k_parallel.Modulus2());

	GeometryVector Sum_Real((int)this->Dimension), Sum_Imag((int)this->Dimension);

	for (size_t j = 0; j < this->E.size(); j++) {
		Sum_Real = Sum_Real + sin(k.Dot(this->R[j]))*sinc(k.Dot(this->E[j]))*this->E[j];
		Sum_Imag = Sum_Imag + cos(k.Dot(this->R[j]))*sinc(k.Dot(this->E[j]))*this->E[j];
		//Sum_Imag = Sum_Imag + (cos(k.Dot(this->R[j]))*sinc(k.Dot(this->E[j]))-1.0)*this->E[j];
	}

	std::complex<double> result(Kx.Dot(Sum_Real), Kx.Dot(Sum_Imag));
	result *= 2.0 / (K_parallel*K_parallel);
	return result;
}

/** \brief	Compute the form factor via a series expasion formula.
* It is accurate and easy to compute when |k| is very small, but not implemented here.
* \param k	A wavevector in 2D or 3D.
* \return A complex form factor of this PlanarFigure.
*/
std::complex<double> PlanarFigure::FormFactor_series(const GeometryVector & k)
{
	std::complex<double> result;
	std::cerr << "PlanarFigure::FormFactor_series\tNot defined!\n";
	return result;
}

/** \brief	Compute the form factor.
* \param k	A wavevector in 2D or 3D.
* \return A complex form factor of this PlanarFigure.
*/
std::complex<double> PlanarFigure::GetFormFactor(const GeometryVector & k){
	std::complex<double> result;

	double K = k.Modulus2();
	double K_perp = (this->Dimension == 2)? 0:k.Dot(this->normal);	
	double K_parallel = sqrt(K - K_perp*K_perp);	K = sqrt(K);
	if (K_parallel < c_parallel) {
		result = this->GetVolume()* std::complex<double>(cos(K_perp*this->r_perp), -1.0*sin(K_perp*this->r_perp));
	}
	/*
	else if (K*this->GetMaxRadius() < c) {
		
	}
	else if (K_parallel*this->GetMaxRadius() < c_parallel) {
		
	}
	*/
	else {
		//The reference suggested to use series expansion for small k's
		//However, I found that the analytic expression is accurate enough down to ka < 1e-40;
		result = this->FormFactor_analytic(k);
	}
	return result;
}

/** \brief	Compute the form factor of this Polyhedron with respect to its centroid.
* \param k	A wavevector in 3D.
* \return A complex form factor of this Polyhedron.
*/
std::complex<double> Polyhedron::GetFormFactor(const GeometryVector &k) {
	std::complex<double> result;
	double K = sqrt(k.Modulus2());
	if (K == 0.0) {
		result = std::complex<double>(this->GetVolume(), 0.0);
	}
	/*
	else if (K*this->GetMaxRadius() < C) {
		result = this->FormFactor_series(k);
	}*/
	else {
		result = this->FormFactor_analytic(k);
	}
	return result;
}

/** \brief	Compute the form factor of this Polyhedron via an analytic expression.
* \param k	A wavevector in 3D.
* \return A complex form factor of this Polyhedron.
*/
std::complex<double> Polyhedron::FormFactor_analytic(const GeometryVector & k)
{
	std::complex<double> result(0.0, 0.0);
	GeometryVector real_part((int)3), imag_part((int)3);
	for (auto gamma = this->Gamma.begin(); gamma != this->Gamma.end(); gamma++) {
		std::complex<double> formfactor = gamma->GetFormFactor(k);
		real_part = real_part - std::imag(formfactor)*gamma->normal;
		imag_part = imag_part + std::real(formfactor)*gamma->normal;
	}

	return (1.0)/k.Modulus2() * std::complex<double>(k.Dot(real_part), k.Dot(imag_part));
}

/** \brief	Compute the form factor of this Polyhedron via a series expansion formula.
* It is accurate and easy to compute when |k| is small, but not implemented here. 
* \param k	A wavevector in 3D.
* \return A complex form factor of this Polyhedron.
*/
std::complex<double> Polyhedron::FormFactor_series(const GeometryVector & k)
{
	std::cerr << "Polyhedron::FormFactor_series is not implemented \n";
	return std::complex<double>();
}

