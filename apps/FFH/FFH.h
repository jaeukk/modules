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


/** \brief Squared modulus of a complex number stored as fftw_complex. (CLAUDE)
 *  \param x  The complex value (real, imaginary).
 *  \return   |x|^2 = re^2 + im^2. */
inline double Modulus2(const fftw_complex & x) {
	return (x[0] * x[0] + x[1] * x[1]);
}


/** \brief Truncated (finite-box) Green's function of the Helmholtz equation. (CLAUDE)
 *  \param d   Spatial dimension (2 or 3).
 *  \param xi  Wavenumber magnitude of the evaluation point.
 *  \param L   Linear size of the truncation box.
 *  \param k1  Background wavenumber.
 *  \return    The truncated Green's function value. */
inline std::complex<double> TruncatedGreenFn_Helmholtz(int d, double xi, double L, double k1);

struct FFH_base {
	//friend int main();
	friend int main(int, char**);

	bool optimized_reference = true;
	double Tol_Field = -1.0,	/* the average of the difference in the fields at the current and previous stages. */ 
		Tol_div = -1.0;			/* the average of the "divergence" of the computed fields. */
	size_t m;
	/* contructors */
	/** \brief Construct an empty homogenization problem. (CLAUDE) */
	FFH_base();

	/** \brief Construct a homogenization problem from a geometry and two phase tensors. (CLAUDE)
	 *  \param src              The two-phase dispersion (geometry).
	 *  \param matrix_phase     Constant tensor describing the matrix phase (phase 1).
	 *  \param inclusion_phase  Constant tensor describing the inclusion phase (phase 2). */
	FFH_base(const TwoPhaseDispersion & src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);

	/** \brief destructor. Destroys the fftw plans and cleans up the fftw threads. (CLAUDE) */
	~FFH_base(){
		for (size_t i = 0; i < FourierTransform_P.size(); i++) {
			if (FourierTransform_P[i]!=nullptr)
				fftw_destroy_plan(FourierTransform_P[i]);
			if (InverseFourierTransform_Response[i]!=nullptr)
				fftw_destroy_plan(InverseFourierTransform_Response[i]);
		}
		fftw_cleanup_threads();
		// Propagator_Fourier.~TensorField();
		// P_Fourier.~TensorField();	Response_Fourier.~TensorField();
		// P_direct.~TensorField();	Response_direct.~TensorField();
	};

	

	/** \brief redefine problem, including composites ...
	 *  \param src              The two-phase dispersion (geometry).
	 *  \param matrix_phase     Constant tensor describing the matrix phase (phase 1).
	 *  \param inclusion_phase  Constant tensor describing the inclusion phase (phase 2). */
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);

	//virtual void DefineProperties(const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);

	/** \brief Select the weak- or strong-contrast expansion. (CLAUDE)
	 *  \param strong_contrast  true = use the strong-contrast expansion. */
	void SetMethod(bool strong_contrast);

	/** \brief determine whehter we use the optimal reference phase or not.
	 *  \param use  true = use the optimal reference phase. */
	void UseOptimalReference(bool use) { this->optimized_reference = use; };

	/** \brief Initialize the status of the problem.
	 *  \param[in] applied  external field. */
	virtual void Prepare(const ConstantTensor & applied);

	/** \brief Advance the iterative solver by a number of stages. (CLAUDE)
	 *  \param steps  Number of iteration stages to perform. */
	virtual void Iterate(size_t steps = 1);

	/** \brief repeate Iterate() until achieving the prescribed tolerances...
	 *  \param target_tol_field  Target tolerance on the change in the field.
	 *  \param target_tol_div    Target tolerance on the divergence of the field.
	 *  \param num_iters         Maximum number of iterations (-1 for the built-in limit). */
	void Run(double target_tol_field = 1e-6, double target_tol_div = 1e-10, int num_iters = -1);

	/* check convergence... */
	/** \brief compute the change in the response field at two consecutive stages.
	 *  \param target_tolerence  Target tolerance on the field change.
	 *  \return false = error > target_tolerenace. */
	bool CheckConvergence(double target_tolerence = 1e-7);

	/** \brief compute the divergence of fields (supposedly to be divergence-free) in the Fourier space.
	 *  \param target_tolereance  Target tolerance on the divergence.
	 *  \return false erorr > target_tolerance. */
	virtual bool CheckDivergenceFree(double target_tolereance = 1e-10);

	/** \brief Report whether the current solution appears to be converging. (CLAUDE)
	 *  \return true if the solver still considers the solution to be converging. */
	bool IsCovergeningSolution() const { return this->is_converging; };

	/** \brief Print the divergence-free field and the conjugate field in direct space. (CLAUDE)
	 *  \param prefix  Output filename prefix. */
	void Print(const std::string & prefix) const;

	/** \brief Print the divergence-free field and the conjugate field in Fourier space. (CLAUDE)
	 *  \param prefix  Output filename prefix. */
	void Print_Fourier(const std::string & prefix) const;

	/** \brief Print the discretized geometry grid. (CLAUDE)
	 *  \param prefix  Output filename prefix. */
	void Print_Grid(const std::string & prefix) const;

	/** \brief Print the 2D Fourier-space Green's function components. (CLAUDE)
	 *  \param prefix  Output filename prefix. */
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
	


	/** \brief update the geomtry of the dispersion.
	 *  \param src  The new two-phase dispersion (geometry).
	 *  \return false = have to redefine the wavevectors assocate with DFT. */
	bool DefineDispersion(const TwoPhaseDispersion & src);



	/** \brief Compute the divergence-free field in Fourier space. (CLAUDE)
	 *  \param[out] out  The resulting field. */
	virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const = 0;

	/** \brief Compute the divergence-free field in direct space. (CLAUDE)
	 *  \param[out] out  The resulting field. */
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const = 0;

	/** \brief Compute the conjugate (gradient-type) field in Fourier space. (CLAUDE)
	 *  \param[out] out  The resulting field. */
	virtual void ComputeField_Fourier(TensorField & out) const = 0;

	/** \brief Compute the conjugate (gradient-type) field in direct space. (CLAUDE)
	 *  \param[out] out  The resulting field. */
	virtual void ComputeField_Direct(TensorField & out) const = 0;

	/** \brief Compute the Fourier-space propagator for the weak-contrast expansion. (CLAUDE) */
	virtual void ComputePropagator_weak() = 0;

	/** \brief Iterate the weak-contrast solver. (CLAUDE)
	 *  \param steps  Number of iteration stages. */
	virtual void Iterate_weak(size_t steps = 1) { this->Iterate_common(steps); };
	/** \brief Compute the constitutive L tensors for the weak-contrast expansion. (CLAUDE) */
	virtual void ComputeL_weak() = 0;

	/** \brief Compute the Fourier-space propagator for the strong-contrast expansion. (CLAUDE) */
	virtual void ComputePropagator_strong() = 0;

	/** \brief Iterate the strong-contrast solver. (CLAUDE)
	 *  \param steps  Number of iteration stages. */
	virtual void Iterate_strong(size_t steps = 1) { this->Iterate_common(steps); };
	/** \brief Compute the constitutive L tensors for the strong-contrast expansion. (CLAUDE) */
	virtual void ComputeL_strong() = 0;

	/** \brief Shared iteration kernel used by both expansions. (CLAUDE)
	 *  \param steps  Number of iteration stages. */
	void Iterate_common(size_t steps = 1);
};

/* A class derived from FFH_base to solve the homogenization problem for the static dielectric problem. */
struct StaticDielectricProblem : public FFH_base {
	//friend int main();
	friend int main(int, char**);

	friend double curl(const StaticDielectricProblem & src);

	/** \brief Construct an empty static dielectric problem. (CLAUDE) */
	StaticDielectricProblem();

	/** \brief Construct a static dielectric problem from a geometry and two permittivities. (CLAUDE)
	 *  \param src       The two-phase dispersion (geometry).
	 *  \param epsilon1  Permittivity of the matrix phase.
	 *  \param epsilon2  Permittivity of the inclusion phase. */
	StaticDielectricProblem(const TwoPhaseDispersion & src,
		const std::complex<double> & epsilon1, const std::complex<double> & epsilon2);

	/** \brief Define the static dielectric problem geometry. (CLAUDE)
	 *  \param src              The two-phase dispersion (geometry).
	 *  \param matrix_phase     Constant tensor for the matrix phase.
	 *  \param inclusion_phase  Constant tensor for the inclusion phase. */
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);

	/** \brief Set the two permittivities and the derived reference phase. (CLAUDE)
	 *  \param epsilon1  Permittivity of the matrix phase.
	 *  \param epsilon2  Permittivity of the inclusion phase. */
	virtual void DefineProperties(const std::complex<double> & epsilon1, const std::complex<double> & epsilon2);
	//virtual void Prepare(const ConstantTensor & applied_electric_field);

	/** \brief Compute the effective permittivity from the converged fields. (CLAUDE)
	 *  \return The effective (complex) permittivity. */
	virtual std::complex<double> ComputeEffectiveEpsilon() const;
protected:
	std::complex<double> epsilon1, epsilon2, epsilon0;
	ConstantTensor D;

	/** \brief Compute the displacement field D in Fourier space.
	 *  \param[out] out  The resulting displacement field. */
	virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const;

	/** \brief Compute the displacement field D in direct space.
	 *  \param[out] out  The resulting displacement field. */
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const;

	/** \brief Compute the electric field E in Fourier space.
	 *  \param[out] out  The resulting electric field. */
	virtual void ComputeField_Fourier(TensorField & out) const;

	/** \brief Compute the electric field E in direct space.
	 *  \param[out] out  The resulting electric field. */
	virtual void ComputeField_Direct(TensorField & out) const;

	/** \brief Fourier transform of the Green's function.
	 *  G_{ij} = -k_i k_j /k^2 /\epsilon_1. */
	virtual void ComputePropagator_weak();

	/** \brief Compute the weak-contrast L tensors from the permittivities. (CLAUDE) */
	virtual void ComputeL_weak();

	/** \brief Fourier transform of the H-tensor*/
	virtual void ComputePropagator_strong();

	/** \brief Compute the strong-contrast L tensors from the permittivities. (CLAUDE) */
	virtual void ComputeL_strong();
};

struct DynamicDielectricProblem : public StaticDielectricProblem {
	RandomGenerator rng;


	/** \brief Construct an empty dynamic dielectric problem. (CLAUDE) */
	DynamicDielectricProblem();

	/** \brief Construct a dynamic dielectric problem from a geometry and two permittivities. (CLAUDE)
	 *  \param src       The two-phase dispersion (geometry).
	 *  \param epsilon1  Permittivity of the matrix phase.
	 *  \param epsilon2  Permittivity of the inclusion phase. */
	DynamicDielectricProblem(const TwoPhaseDispersion & src,
		const std::complex<double> & epsilon1, const std::complex<double> & epsilon2);

	/** \brief Define the dynamic dielectric problem geometry and incident wavevector. (CLAUDE)
	 *  \param src              The two-phase dispersion (geometry).
	 *  \param matrix_phase     Constant tensor for the matrix phase.
	 *  \param inclusion_phase  Constant tensor for the inclusion phase.
	 *  \param k1               Wavevector of the field in the reference phase. */
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector & k1);

	/** \brief Make the two-argument base-class overload visible (fixes -Woverloaded-virtual). (CLAUDE) */
	using StaticDielectricProblem::DefineProperties;

	/** \brief Set the permittivities and the incident wavevector. (CLAUDE)
	 *  \param epsilon1  Permittivity of the matrix phase.
	 *  \param epsilon2  Permittivity of the inclusion phase.
	 *  \param k1        Wavevector of the field in the reference phase. */
	virtual void DefineProperties(const std::complex<double> & epsilon1, const std::complex<double> & epsilon2, const GeometryVector & k1);

	/** \brief Set the incident wavevector and shift the DFT wavevectors accordingly. (CLAUDE)
	 *  \param k1  Wavevector of the field in the reference phase. */
	void SetWavevector(const GeometryVector & k1);

	//virtual void Prepare(const ConstantTensor & applied_electric_field);

	/** \brief with cut-off... -> this is not rigorous... */
	void ComputeG_Fourier_test();

	/** \brief compute the average values of ??*/
	void ComputeG_Fourier_test2();

	/** \brief estimate the Green's function via the Monte Carlo method... */
	void ComputeG_Fourier_test3();

	/** \brief Compute the Fourier-space propagator with a regularized dyadic Green's function. (CLAUDE)
	 *  \param regularized_part  Imaginary regularization added to the denominator. */
	void ComputeG_Fourier_test4(double regularized_part = 1e-10);
protected:
	GeometryVector k1;	/* wavevector of the electric field in the reference phase!! */

	/* Compute the displacement field D. */
	//virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const;

	/** \brief Compute the displacement field D in direct space, including the e^{ik1*x} phase factor. (CLAUDE)
	 *  \param[out] out  The resulting displacement field. */
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const {
		StaticDielectricProblem::ComputeDivergenceFreeField_Direct(out);
		ShiftFieldByk1(out);
	};
	/** \brief Compute the electric field E in direct space, including the e^{ik1*x} phase factor. (CLAUDE)
	 *  \param[out] out  The resulting electric field. */
	virtual void ComputeField_Direct(TensorField & out) const {
		StaticDielectricProblem::ComputeField_Direct(out);
		ShiftFieldByk1(out);
	};

	/** \brief Fourier transform of the Green's function.
	 *  G_{ij}(q+k_1) of the free-space Green's function . */
	virtual void ComputePropagator_weak();
	//virtual void  ComputeL_weak();

	/** \brief Fourier transform of the H-tensor*/
	virtual void ComputePropagator_strong();
	//virtual void ComputeL_strong();


	/** \brief Evaluate the direct-space H-function (double-gradient Green's function). (CLAUDE)
	 *  \param position  Direct-space position.
	 *  \param k1        Wavevector of the field in the reference phase.
	 *  \param d         Spatial dimension.
	 *  \param dx        Grid spacing (used for the self-interaction at the origin).
	 *  \return          The rank-2 H-tensor at the position. */
	ConstantTensor Hfunction_direct(const GeometryVector & position, const GeometryVector & k1, int d, double dx);

	/** \brief Fourier-space dyadic Green's function at a wavevector. (CLAUDE)
	 *  \param q   The wavevector.
	 *  \param qq  Squared magnitude of q.
	 *  \param kk  Squared magnitude of the reference wavevector.
	 *  \return    The rank-2 dyadic Green's function. */
	inline ConstantTensor DyadicGreenFunction_Fourier(const GeometryVector & q, double qq, double kk);

	/** \brief Regularized Fourier-space dyadic Green's function at a wavevector. (CLAUDE)
	 *  \param q        The wavevector.
	 *  \param qq       Squared magnitude of q.
	 *  \param kk       Squared magnitude of the reference wavevector.
	 *  \param epsilon  Imaginary regularization in the denominator.
	 *  \return         The rank-2 regularized dyadic Green's function. */
	inline ConstantTensor DyadicGreenFunction_Fourier_regularized(const GeometryVector & q, double qq, double kk, double epsilon);
private:
	/** \brief Shift the tensor field in the direct space by k1.
	 *  \param[in,out] out  The tensor field multiplied by the e^{ik1*x} phase factor. */
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
	/** \brief Construct an empty static elastic problem. (CLAUDE) */
	Homo_ElasticProblem_static() {
		d = 0; K1 = 0; G1 = 0; K2 = 0; G2 = 0; kappa = 0, mu=0;
		update_propagator = false; update_plan = false; m = 0;
		Use_StrongContrast = false; L = ConstantTensor(); D = ConstantTensor(), delta_C = ConstantTensor();
	};
	/** \brief Construct a static elastic problem from a geometry and phase moduli. (CLAUDE)
	 *  \param src  The two-phase dispersion (geometry).
	 *  \param K1   Bulk modulus of the reference (matrix) phase.
	 *  \param G1   Shear modulus of the reference (matrix) phase.
	 *  \param K2   Bulk modulus of the polarized (inclusion) phase.
	 *  \param G2   Shear modulus of the polarized (inclusion) phase. */
	Homo_ElasticProblem_static(const TwoPhaseDispersion& src, const std::complex<double> & K1, const std::complex<double> & G1,
		const std::complex<double> & K2, const std::complex<double> & G2);

	/** \brief Destructor; destroys the fftw plans and tensor fields. (CLAUDE) */
	~Homo_ElasticProblem_static() {
		for (size_t i = 0; i < FourierTransform_P.size(); i++) {
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
	/** \brief Redefine the geometry and phase moduli of the problem. (CLAUDE)
	 *  \param src  The two-phase dispersion (geometry).
	 *  \param K1   Bulk modulus of the reference (matrix) phase.
	 *  \param G1   Shear modulus of the reference (matrix) phase.
	 *  \param K2   Bulk modulus of the polarized (inclusion) phase.
	 *  \param G2   Shear modulus of the polarized (inclusion) phase. */
	void RedefineProblem(const TwoPhaseDispersion& src, const std::complex<double> & K1, const std::complex<double> & G1,
		const std::complex<double> & K2, const std::complex<double> & G2);

	/** \brief Select the weak- or strong-contrast expansion and rebuild the L tensor. (CLAUDE)
	 *  \param use_strong_contrast  true = use the strong-contrast expansion. */
	void SetStrongContrast(bool use_strong_contrast);

	/* -------------------------
	Member functions: solve the problem.
	------------------------- */
	/** \brief Initialize the solver for a given applied strain and build the fftw plans. (CLAUDE)
	 *  \param external_strain  The applied (constant) strain tensor. */
	void Prepare(const ConstantTensor & external_strain);

	/** \brief Advance the iterative solver by a number of stages. (CLAUDE)
	 *  \param steps  Number of iteration stages. */
	void Iterate(size_t steps = 1);

	/** \brief Iterate until the field and divergence tolerances are met. (CLAUDE)
	 *  \param tol_field  Target tolerance on the change in the field.
	 *  \param tol_div    Target tolerance on the divergence of the stress. */
	void Run(double tol_field = 1e-10, double tol_div = 1e-14);

	/** \brief Check the change in the response field between two stages. (CLAUDE)
	 *  \param tol_field  Target tolerance.
	 *  \return true if the change is below tol_field. */
	bool CheckConvergence(double tol_field);

	/** \brief Check the divergence of the stress field in Fourier space. (CLAUDE)
	 *  \param tol_div  Target tolerance.
	 *  \return true if the divergence is below tol_div. */
	bool CheckDivergenceFree(double tol_div);

	/** \brief compute the effective bulk and shear moduli assuming the isotropy.
	 *  \param[out] moduli  [0] bulk modulus, [1] shear modulus. */
	void ComputeEffectiveProperties(std::vector<std::complex<double>> & moduli) const;

	/** \brief Compute the full effective stiffness tensor. (CLAUDE)
	 *  \return The effective rank-4 stiffness tensor. */
	ConstantTensor ComputeEffectiveTensor() const;

	/** \brief Print the response and polarization fields in direct space. (CLAUDE)
	 *  \param prefix  Output filename prefix. */
	void Print(const std::string & prefix) const;

	/** \brief Print the response and polarization fields in Fourier space. (CLAUDE)
	 *  \param prefix  Output filename prefix. */
	void Print_Fourier(const std::string & prefix) const;

	/** \brief Print the 2D discretized geometry grid. (CLAUDE)
	 *  \param prefix  Output filename prefix. */
	void Print_Grid_2D(const std::string & prefix) const;

	/** \brief Print the 2D Fourier-space Green's function components. (CLAUDE)
	 *  \param prefix  Output filename prefix. */
	void Print_GreenFn_2D(const std::string & prefix) const;

private:
	bool Use_StrongContrast;	/* true = use the strong-contrast expansion*/
	TensorField k, Response_direct_prev;	/* temporary storage */

	/** \brief member function for weak-contrast expansions.
	 *  \param steps  Number of iteration stages. */
	void Iterate_weak_contrast(size_t steps);
	/*\brief compute the Fourier transform of the Green's function G^{(q)}(x) tensor field.
	* G^{(q)}(x) is the rank-4 Green's function that solves the excess "strain" field
	from the polarization field p(x); see Eq. (20.100) in Ref. [3].
	* This function is identical to $-\Gamma_0$ in Ref. [1]. */
	void ComputeG_Fourier();

	/** \brief member function for strong-contrast expansions.
	 *  \param steps  Number of iteration stages. */
	void Iterate_strong_contrast(size_t steps);
	/*\brief compute the Fourier transform of the H^{(q)}(x) tensor field.
	 * H^{(q)}(x) is the (symmetric) double graident 
	 of the diadic Green's function $g_{ij}(x)$ that solves the excess displacement field
	 from the polarization field p(x); see Eq. (20.100) in Ref. [3].	 */
	void ComputeH_Fourier();
};

/** \brief Return the unit hydrostatic strain tensor in d dimensions. (CLAUDE)
 *  \param d  Spatial dimension.
 *  \return   The rank-2 hydrostatic strain tensor. */
ConstantTensor HydrostaticStrain(DimensionType d);

/** \brief Return a unit deviatoric (shear) strain tensor in d dimensions. (CLAUDE)
 *  \param d  Spatial dimension.
 *  \return   The rank-2 deviatoric strain tensor. */
ConstantTensor DeviatoricStrain(DimensionType d);


//TODO:: define the StaticElasticProblem struct and then inherit it to DynamicElasticProblem.

/* \brief class to solve the problem for elastodynamic problem via iterative solvers.
 *	Only isotropic phases are implemented. 
 *  mass densities of both phases are identical. */
struct DynamicElasticProblem : FFH_base {
	/** \brief Construct an empty dynamic elastic problem (mass density = 1). (CLAUDE) */
	DynamicElasticProblem();

	/** \brief Construct a dynamic elastic problem from a geometry and phase moduli. (CLAUDE)
	 *  \param src  The two-phase dispersion (geometry).
	 *  \param K1   Bulk modulus of the reference (matrix) phase.
	 *  \param G1   Shear modulus of the reference (matrix) phase.
	 *  \param K2   Bulk modulus of the polarized (inclusion) phase.
	 *  \param G2   Shear modulus of the polarized (inclusion) phase. */
	DynamicElasticProblem(const TwoPhaseDispersion & src,
		const std::complex<double> & K1, const std::complex<double> & G1,
		const std::complex<double> & K2, const std::complex<double> & G2);

	/** \brief Define the dynamic elastic problem geometry. (CLAUDE)
	 *  \param src              The two-phase dispersion (geometry).
	 *  \param matrix_phase     Stiffness tensor of the matrix phase.
	 *  \param inclusion_phase  Stiffness tensor of the inclusion phase. */
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase);

	/** \brief Define the dynamic elastic problem geometry and incident wavevector.
	 *  \param src              The two-phase dispersion (geometry).
	 *  \param matrix_phase     Stiffness tensor of the matrix phase.
	 *  \param inclusion_phase  Stiffness tensor of the inclusion phase.
	 *  \param kL	Wavevector of the longitudinal wave. */
	virtual void DefineProblem(const TwoPhaseDispersion& src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector kL);

	/** \brief Set the phase moduli and the derived reference/polarizability parameters. (CLAUDE)
	 *  \param K1  Bulk modulus of the reference (matrix) phase.
	 *  \param G1  Shear modulus of the reference (matrix) phase.
	 *  \param K2  Bulk modulus of the polarized (inclusion) phase.
	 *  \param G2  Shear modulus of the polarized (inclusion) phase. */
	virtual void DefineProperties(const std::complex<double> & K1, const std::complex<double> & G1,
		const std::complex<double> & K2, const std::complex<double> & G2);

	/** \brief Set the propagation direction and the longitudinal/transverse wavenumbers. (CLAUDE)
	 *  \param kL  Wavevector of the longitudinal wave. */
	virtual void SetWavevector(const GeometryVector & kL);

	/** \brief check [\omega^2 \rho_e \FT[strain] - q_i q_j \FT[stress]_{jk} = 0.
	 *  \param target_tolerance  Target tolerance.
	 *  \return true if the equation of motion residual is below the tolerance. */
	virtual bool CheckDivergenceFree(double target_tolerance = 1e-10);

	/** \brief Prepare the solver; strain is decomposed into the longitudial and transverse contributions.
	 *  \param applied  The applied (constant) strain tensor. */
	virtual void Prepare(const ConstantTensor & applied);

	/** \brief Compute the effective bulk and shear moduli from the converged fields. (CLAUDE)
	 *  \param[out] moduli  [0] effective bulk modulus, [1] effective shear modulus. */
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

	/** \brief Compute the stress tensor in Fourier space.
	 *  \param[out] out  The resulting stress field. */
	virtual void ComputeDivergenceFreeField_Fourier(TensorField & out) const;

	/** \brief Compute the stress tensor in direct space.
	 *  \param[out] out  The resulting stress field. */
	virtual void ComputeDivergenceFreeField_Direct(TensorField & out) const;

	/** \brief Compute the strain tensor in Fourier space.
	 *  \param[out] out  The resulting strain field. */
	virtual void ComputeField_Fourier(TensorField & out) const;

	/** \brief Compute the strain tensor in direct space.
	 *  \param[out] out  The resulting strain field. */
	virtual void ComputeField_Direct(TensorField & out) const;

	/** \brief Fourier transform of the 4th rank Green's function for the elastodynamic problem.
	 *  G_{ijkl}. */
	virtual void ComputePropagator_weak();

	/** \brief Compute the weak-contrast L tensors from the phase moduli. (CLAUDE) */
	virtual void ComputeL_weak();

	/** \brief Fourier transform of the H-tensor*/
	virtual void ComputePropagator_strong();

	/** \brief Compute the strong-contrast L tensors from the phase moduli. (CLAUDE) */
	virtual void ComputeL_strong();

	/** \brief Evaluate the Fourier-space rank-4 Green's function at a wavevector.
	 *  \param q  spatial frequency (shifted by a wavevector of the incident waves).
	 *  \return   The rank-4 Green's function tensor. */
	inline ConstantTensor GreenFunction_Fourier(const GeometryVector & q);
private:
	/** \brief Shift the tensor field in the direct space by k1.
	 *  \param[in,out] out  The tensor field multiplied by the e^{ik0*x} phase factor. */
	void ShiftFieldByk0(TensorField & out) const;
};

/** \brief return the symmetrized tensor product. q*\delta*q
 *  \param q  The wavevector.
 *  \return   The rank-4 symmetrized tensor product. */
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