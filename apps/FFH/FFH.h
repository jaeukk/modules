#ifndef FFH_INCLUDED
#define FFH_INCLUDED

/** Author	:	Jaeuk Kim
*	Email	:	phy000.kim@gmail.com
*	Date	:	November	2019 */

/** \file FFH.h
*	\brief Header file to define classes and functions \
		for solving the homogenization problems with  */

#include <complex>	/* this should be compiled earlier than fftw3!! */
#include <cmath>
#include <fftw3.h>	/* library for computing the discerete fast Fourier transform. */
#include <vector>
#include <omp.h>
#include <fstream>
#include <boost/math/special_functions/hankel.hpp>

#include "Tensors.h"
#include "TwoPhaseDispersion.h"
#include <PeriodicCellList.h>
#include <RandomGenerator.h>

//const size_t size_fftw = sizeof(fftw_complex);
extern size_t FFH_verbosity;	/* determine the amount of information displayed on the screen. 
								 * 0 = display the minimal information.
								 * 1 = display progress.
								 * 2 = display details. */
const size_t num_limit = 200;	/* the number of maximum iterations. */


inline double Modulus2(const fftw_complex & x) {
	return (x[0] * x[0] + x[1] * x[1]);
}


inline std::complex<double> TruncatedGreenFn_Helmholtz(int d, double xi, double L, double k1);

struct FFH_base {
	//friend int main();
	friend int main(int, char**);

	bool optimized_reference = true;
	double Tol_Field = -1.0,	/* the average of the difference in the fields at the current and previous stages. */ 
		Tol_div = -1.0;			/* the average of the "divergence" of the computed fields. */
	size_t m;
	/* contructors */
	FFH_base();
	FFH_base(const TwoPhaseDispersion & src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);

	/* destructor. */
	~FFH_base(){
		for (int i = 0; i < FourierTransform_P.size(); i++) {
			if (FourierTransform_P[i]!=nullptr);
				fftw_destroy_plan(FourierTransform_P[i]);
			if (InverseFourierTransform_Response[i]!=nullptr)
				fftw_destroy_plan(InverseFourierTransform_Response[i]);
		}
		fftw_cleanup_threads();
		// Propagator_Fourier.~TensorField();
		// P_Fourier.~TensorField();	Response_Fourier.~TensorField();
		// P_direct.~TensorField();	Response_direct.~TensorField();
	};

	

	/* redefine problem, including composites ... */
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);

	//virtual void DefineProperties(const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);
	void SetMethod(bool strong_contrast);

	/* determine whehter we use the optimal reference phase or not. */
	void UseOptimalReference(bool use) { this->optimized_reference = use; };

	/* Initialize the status of the problem.
	 * @param[in] external field. */
	virtual void Prepare(const ConstantTensor & applied);
	virtual void Iterate(size_t steps = 1);
	/* repeate Iterate() until achieving the prescribed tolerances... */
	void Run(double target_tol_field = 1e-6, double target_tol_div = 1e-10, int num_iters = -1);

	/* check convergence... */
	/* compute the change in the response field at two consecutive stages.
	 @ return false = error > target_tolerenace. */
	bool CheckConvergence(double target_tolerence = 1e-7);
	/* compute the divergence of fields (supposedly to be divergence-free) in the Fourier space.
	 *@ return false erorr > target_tolerance. */
	virtual bool CheckDivergenceFree(double target_tolereance = 1e-10);

	bool IsCovergeningSolution() const { return this->is_converging; };

	void Print(const std::string & prefix) const;
	void Print_Fourier(const std::string & prefix) const;
	void Print_Grid(const std::string & prefix) const;
	void Print_GreenFn_2D(const std::string & prefix) const;

protected:
	/* member variables that must  be defined in each problem!!*/
	std::string name_div_free_field, name_field;
	size_t rank_of_fields;
	/* member variables: system paramters */
	DimensionType d;
	TwoPhaseDispersion Geometry;

	ConstantTensor Appl_field, 	/* Applied constant field.*/
					L1,			/* A constant tensor to define the constitutive relation from P_direct to Response_direct in the matrix phase (phase1). This tensor depends on the scheme that we use. */
					L2;			/* A constant tensor to define the constitutive relation from P_direct to Response_direct in the inclusion phase (phase2)*/
								//	L			reference_phase,	/* A constant tensor describing the properties of the reference phase. */
	//				polarized_phase;	/* A constant tensor describing the properties of the polarized phase. */
	/* member variables: tensor fields */
	TensorField P_direct, Response_direct, Response_direct_prev;
	TensorField P_Fourier, Response_Fourier;
	TensorField ks;	/* wavevectors associated with the DFT. */

	/* variables about the states */
	bool update_propagator, update_plan, use_strong_contrast = false;
	bool is_converging;
	TensorField Propagator_Fourier;
	/* fftw_plans */
	std::vector<fftw_plan> FourierTransform_P, InverseFourierTransform_Response;
	


	/* update the geomtry of the dispersion.
	 @return false = have to redefine the wavevectors assocate with DFT. */
	bool DefineDispersion(const TwoPhaseDispersion & src);
	


	virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const = 0;
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const = 0;
	virtual void ComputeField_Fourier(TensorField & out) const = 0;
	virtual void ComputeField_Direct(TensorField & out) const = 0;

	virtual void ComputePropagator_weak() = 0;
	virtual void Iterate_weak(size_t steps = 1) { this->Iterate_common(steps); };
	virtual void ComputeL_weak() = 0;

	virtual void ComputePropagator_strong() = 0;
	virtual void Iterate_strong(size_t steps = 1) { this->Iterate_common(steps); };
	virtual void ComputeL_strong() = 0;

	void Iterate_common(size_t steps = 1);
};

/* A class derived from FFH_base to solve the homogenization problem for the static dielectric problem. */
struct StaticDielectricProblem : public FFH_base {
	//friend int main();
	friend int main(int, char**);

	friend double curl(const StaticDielectricProblem & src);

	StaticDielectricProblem();
	StaticDielectricProblem(const TwoPhaseDispersion & src, 
		const std::complex<double> & epsilon1, const std::complex<double> & epsilon2);

	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);
	virtual void DefineProperties(const std::complex<double> & epsilon1, const std::complex<double> & epsilon2);
	//virtual void Prepare(const ConstantTensor & applied_electric_field);
	
	virtual std::complex<double> ComputeEffectiveEpsilon() const;
protected:
	std::complex<double> epsilon1, epsilon2, epsilon0;
	ConstantTensor D;

	/* Compute the displacement field D. */
	virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const;
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const;
	/* Compute the electric field E. */
	virtual void ComputeField_Fourier(TensorField & out) const;
	virtual void ComputeField_Direct(TensorField & out) const;

	/* Fourier transform of the Green's function.
	 * G_{ij} = -k_i k_j /k^2 /\epsilon_1. */
	virtual void ComputePropagator_weak();
	virtual void ComputeL_weak();

	/* Fourier transform of the H-tensor*/
	virtual void ComputePropagator_strong();
	virtual void ComputeL_strong();
};

struct DynamicDielectricProblem : public StaticDielectricProblem {
	RandomGenerator rng;


	DynamicDielectricProblem();
	DynamicDielectricProblem(const TwoPhaseDispersion & src,
		const std::complex<double> & epsilon1, const std::complex<double> & epsilon2);

	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector & k1);
	virtual void DefineProperties(const std::complex<double> & epsilon1, const std::complex<double> & epsilon2, const GeometryVector & k1);
	void SetWavevector(const GeometryVector & k1);
	
	//virtual void Prepare(const ConstantTensor & applied_electric_field);

	/* with cut-off... -> this is not rigorous... */
	void ComputeG_Fourier_test();
	/* compute the average values of ??*/
	void ComputeG_Fourier_test2();
	/* estimate the Green's function via the Monte Carlo method... */
	void ComputeG_Fourier_test3();

	void ComputeG_Fourier_test4(double regularized_part = 1e-10);
protected:
	GeometryVector k1;	/* wavevector of the electric field in the reference phase!! */

	/* Compute the displacement field D. */
	//virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const;
	/* It should multiply e^{ik1*x} */
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const {
		StaticDielectricProblem::ComputeDivergenceFreeField_Direct(out);
		ShiftFieldByk1(out);
	};
	/* Compute the electric field E. */
	//virtual void ComputeField_Fourier(TensorField & out) const;
	virtual void ComputeField_Direct(TensorField & out) const {
		StaticDielectricProblem::ComputeField_Direct(out);
		ShiftFieldByk1(out);
	};

	/* Fourier transform of the Green's function.
	 * G_{ij}(q+k_1) of the free-space Green's function . */
	virtual void ComputePropagator_weak();
	//virtual void  ComputeL_weak();

	/* Fourier transform of the H-tensor*/
	virtual void ComputePropagator_strong();
	//virtual void ComputeL_strong();


	ConstantTensor Hfunction_direct(const GeometryVector & position, const GeometryVector & k1, int d, double dx);
	
	inline ConstantTensor DyadicGreenFunction_Fourier(const GeometryVector & q, double qq, double kk);
	inline ConstantTensor DyadicGreenFunction_Fourier_regularized(const GeometryVector & q, double qq, double kk, double epsilon);
private:
	/* Shift the tensor field in the direct space by k1.*/
	void ShiftFieldByk1(TensorField & out) const;
};






/* \brief class to solve the problem for elastostatic problem via iterative solvers.
 *	Only isotropic phases are implemented. */
struct Homo_ElasticProblem_static{
	DimensionType d;
	TwoPhaseDispersion Geometry;
	double Tol_Field = -1, Tol_div = -1;
	std::complex<double> K1, K2, /* bulk moduli of the reference and polarized phases, respectively. */
		G1, G2, /* shear moduli of the reference and polarized phases, respectively. */
		kappa, /* bulk polarizability. */
		mu;	/* shear polarizability. */
	ConstantTensor L;	/* A constant rank-4 tensor defining the constitutive relation:
						P_direct(x) = I(x) L: Response_direct (x). 
							L I(x) = C(x) - C1 if USE_Strong_Contrast = false.
							L I(x) = L^(1) (x) if USE_Strong_Contrast = true. */
	ConstantTensor D, delta_C;
	bool update_propagator, update_plan;
	
	ConstantTensor applied_strain;	/* direction of the external strain. */
	TensorField Propagator_Fourier; /*< Static Green function in Fourier space*/
	TensorField P_Fourier, Response_Fourier;	/*< Polarization and electric field in Fourier space representations */
	TensorField P_direct, Response_direct;		/*< Polarization and electric field in direct space representations */
	std::vector<fftw_plan> FourierTransform_P, InverseFourierTransform_Response;
	size_t m;	/*< the number of stages */
	
	/* ------------------------- 
		Constructors.
	   ------------------------- */
	Homo_ElasticProblem_static() {
		d = 0; K1 = 0; G1 = 0; K2 = 0; G2 = 0; kappa = 0, mu=0;
		update_propagator = false; update_plan = false; m = 0;
		Use_StrongContrast = false; L = ConstantTensor(); D = ConstantTensor(), delta_C = ConstantTensor();
	};
	Homo_ElasticProblem_static(const TwoPhaseDispersion& src, const std::complex<double> & K1, const std::complex<double> & G1,
		const std::complex<double> & K2, const std::complex<double> & G2);
	~Homo_ElasticProblem_static() {
		for (int i = 0; i < FourierTransform_P.size(); i++) {
			fftw_destroy_plan(FourierTransform_P[i]);
			fftw_destroy_plan(InverseFourierTransform_Response[i]);
		}
		fftw_cleanup_threads();
		Propagator_Fourier.~TensorField();
		P_Fourier.~TensorField();	Response_Fourier.~TensorField();
		P_direct.~TensorField();	Response_direct.~TensorField();
	};


	/* -------------------------
		Member functions: define the problem.
	------------------------- */
	void RedefineProblem(const TwoPhaseDispersion& src, const std::complex<double> & K1, const std::complex<double> & G1,
		const std::complex<double> & K2, const std::complex<double> & G2);
	/* Initialize the problem.
	* @param[in] comp_epsilon0	The direction of external static field. */
	void SetStrongContrast(bool use_strong_contrast);

	/* -------------------------
	Member functions: solve the problem.
	------------------------- */
	void Prepare(const ConstantTensor & external_strain);
	/* @param[in] use_strong_contrast	true = use the strong-contrast expansion. */
	void Iterate(size_t steps = 1);

	void Run(double tol_field = 1e-10, double tol_div = 1e-14);

	bool CheckConvergence(double tol_field);
	bool CheckDivergenceFree(double tol_div);
	/* compute the effective bulk and shear moduli assuming the isotropy. 
	 * @return [0] bulk modulus, [1] shear modulus*/
	void ComputeEffectiveProperties(std::vector<std::complex<double>> & moduli) const;
	ConstantTensor ComputeEffectiveTensor() const;
	void Print(const std::string & prefix) const;
	void Print_Fourier(const std::string & prefix) const;
	void Print_Grid_2D(const std::string & prefix) const;
	void Print_GreenFn_2D(const std::string & prefix) const;
	
private:
	bool Use_StrongContrast;	/* true = use the strong-contrast expansion*/
	TensorField k, Response_direct_prev;	/* temporary storage */
	
	/* member function for weak-contrast expansions. */
	void Iterate_weak_contrast(size_t steps);
	/*\brief compute the Fourier transform of the Green's function G^{(q)}(x) tensor field.
	* G^{(q)}(x) is the rank-4 Green's function that solves the excess "strain" field
	from the polarization field p(x); see Eq. (20.100) in Ref. [3].
	* This function is identical to $-\Gamma_0$ in Ref. [1]. */
	void ComputeG_Fourier();

	/* member function for strong-contrast expansions. */
	void Iterate_strong_contrast(size_t steps);
	/*\brief compute the Fourier transform of the H^{(q)}(x) tensor field.
	 * H^{(q)}(x) is the (symmetric) double graident 
	 of the diadic Green's function $g_{ij}(x)$ that solves the excess displacement field
	 from the polarization field p(x); see Eq. (20.100) in Ref. [3].	 */
	void ComputeH_Fourier();
};

ConstantTensor HydrostaticStrain(DimensionType d);
ConstantTensor DeviatoricStrain(DimensionType d);


//TODO:: define the StaticElasticProblem struct and then inherit it to DynamicElasticProblem.

/* \brief class to solve the problem for elastodynamic problem via iterative solvers.
 *	Only isotropic phases are implemented. 
 *  mass densities of both phases are identical. */
struct DynamicElasticProblem : FFH_base {
	/* mass density = 1*/
	DynamicElasticProblem();
	DynamicElasticProblem(const TwoPhaseDispersion & src,
		const std::complex<double> & K1, const std::complex<double> & G1, 
		const std::complex<double> & K2, const std::complex<double> & G2);

	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);
	/* @param[in] kL	Wavevector of the longitudinal wave. */
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector kL);
	virtual void DefineProperties(const std::complex<double> & K1, const std::complex<double> & G1,
		const std::complex<double> & K2, const std::complex<double> & G2);
	virtual void SetWavevector(const GeometryVector & kL);
	/* check [\omega^2 \rho_e \FT[strain] - q_i q_j \FT[stress]_{jk} = 0. */
	virtual bool CheckDivergenceFree(double target_tolerance = 1e-10);
	/* strain is decomposed into the longitudial and transverse contributions. */
	virtual void Prepare(const ConstantTensor & applied);

	virtual void ComputeEffectiveModuli(std::vector<std::complex<double>> & moduli) const;

protected:
	std::complex<double> K1, K2, G1, G2, kappa1, mu1, kappa2, mu2;
	std::complex<double> K0, G0;
	std::complex<double> rho;
	ConstantTensor D;

	TensorField qq;		/* a tensor field to store the tensor products of spatial frequencies. */
	GeometryVector k;	/* the unit vector representing the direction of propagation. */
	GeometryVector k0;	/* temp: wavevector of the incident waves.*/
	ConstantTensor kk;	/* A projection operator onto the subspace spanned by k. */
	ConstantTensor strain_L, strain_T;
	std::complex<double> kL, kT;	/* the longitudinal / transverse wavenumbers */
	double cL, cT;		/* phase speeds of the longitudinal and transverse waves in the reference phase. */

	/* Compute the stress tensor. */
	virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const;
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const;
	/* Compute the strain tensor. */
	virtual void ComputeField_Fourier(TensorField & out) const;
	virtual void ComputeField_Direct(TensorField & out) const;

	/* Fourier transform of the 4th rank Green's function for the elastodynamic problem.
	 * G_{ijkl}. */
	virtual void ComputePropagator_weak();
	virtual void ComputeL_weak();

	/* Fourier transform of the H-tensor*/
	virtual void ComputePropagator_strong();
	virtual void ComputeL_strong();

	/* @param[in] spatial frequency (shifted by a wavevector of the incident waves). */
	inline ConstantTensor GreenFunction_Fourier(const GeometryVector & q);
private:
	/* Shift the tensor field in the direct space by k1.*/
	void ShiftFieldByk0(TensorField & out) const;
};

/* return the symmetrized tensor product. q*\delta*q */
inline ConstantTensor T2(const GeometryVector & q);

/* References:
[1] H. Moulinec and P. Suquet, A Numerical Method for Computing the Overall Response of Nonlinear Composites with Complex Microstructure,
Computer Methods in Applied Mechanics and Engineering, 157(1), 69-94 (1998),
doi = {10.1016/S0045-7825(97)00218-1}.

[2] D.J. Eyre and G.W. Milton, A Fast Numerical Scheme for Computing the Response of Composites Using Grid Refinement,
The European Physical Journal Applied Physics, 6(1), 41-47 (1999),
doi = {10.1051/epjap:1999150}.

[3] S. Torquat, Random Heterogeneous Materials: Microstructure and Macroscopic Properties, (2002).
*/



#endif