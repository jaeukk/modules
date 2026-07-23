/** Author	:	Jaeuk Kim
*	Email	:	phy000.kim@gmail.com
*	Date	:	November	2019 */

/** \file FFH.cpp
*	\brief Function implementations of performing the fast Fourier transform homogenization */

#include "FFH.h"
#include <assert.h>
#include <algorithm>	/* library for using the sort() function. */
#include <cmath>





/* -------------------Member functions of FFH_base */
FFH_base::FFH_base()
{	d = 0;	update_propagator = false; update_plan = false; m = 0;}

FFH_base::FFH_base(const TwoPhaseDispersion & src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase) {
	this->DefineProblem(src, matrix_phase, inclusion_phase);
}

bool FFH_base::DefineDispersion(const TwoPhaseDispersion & src)
{
	bool SameSimulationBox = (d == src.packing.GetDimension());
	d = src.packing.GetDimension();
	/* compare simulation box */
	if (SameSimulationBox) {
		std::vector<size_t> prev_resolutions = Geometry.GetResolutions();
		std::vector<size_t> new_resolutions = src.GetResolutions();
		for (DimensionType i = 0; i < d; i++) {
			/* compare basis vectors */
			SameSimulationBox &= (Geometry.packing.GetBasisVector(i) == src.packing.GetBasisVector(i));
			/* compare resolutions in each side */
			SameSimulationBox &= (prev_resolutions[i] == new_resolutions[i]);
		}
	}
	this->Geometry = TwoPhaseDispersion(src);
	return SameSimulationBox;
}

/* default iterate method for weak-contrast expansions. */
void FFH_base::Iterate_common(size_t steps)
{
	for (size_t i = 0; i < steps; i++) {
		m++;
		if (FFH_verbosity > 1) {
			std::cout << "------------------" << m << " th stage ---------------\n";
			std::cout << "\t\t Compute F field \n";
		}

		/*==== 3. Compute
			 Response(k) <- Response_0(k) + Propagator(k)* P(k).*/
		InnerProduct(Response_Fourier, Propagator_Fourier, P_Fourier, rank_of_fields);
		/* Add Response[0] = V * applied */
		for (size_t i = 0; i < Appl_field.total_comp; i++) {
			Response_Fourier[i][0][0] = Appl_field[i].real() * Response_Fourier.num_total_pixels;
			Response_Fourier[i][0][1] = Appl_field[i].imag() * Response_Fourier.num_total_pixels;
		}

		/*==== 4. Compute the inverse FFT of Respose_Fourier. */
		Response_direct_prev = Response_direct;
		if (FFH_verbosity > 1)
			std::cout << "\t\t Inverse Fourier transform of F field\n";
		for (size_t j = 0; j < InverseFourierTransform_Response.size(); j++)
			fftw_execute(InverseFourierTransform_Response[j]);
		Response_direct *= (1.0 / P_Fourier.num_total_pixels);	//This step is required because fftw3 do not consider the 1/N factor.

		/*==== 5. Compute P(x) = L(x)* F(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Compute P field\n";
		/*#pragma omp parallel for
				for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
					size_t x = Geometry.prt_indicator[i];
					for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
						std::complex<double> P = InnerProduct(x, comp, this->L, Response_direct, rank_of_fields);
						P *= Geometry.WeightedPhaseIndicator(i);
						P_direct[comp][x][0] = P.real();	P_direct[comp][x][1] = P.imag();
					}
				}*/
		/* @ phase 1*/
		InnerProduct(P_direct, L1, Response_direct, rank_of_fields);
#pragma omp parallel for
		for (size_t i = 0; i < Geometry.prt_indicator.size(); i++) {
			/* @ phase 2*/
			size_t x = Geometry.prt_indicator[i];
			for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
				std::complex<double> P = InnerProduct(x, comp, this->L2, Response_direct, rank_of_fields);
				P *= Geometry.WeightedPhaseIndicator(i);
//#pragma omp atomic
				P_direct[comp][x][0] = P.real();
//#pragma omp atomic
				P_direct[comp][x][1] = P.imag();
			}
		}

		/*==== 1. Compute FFT of the polarization field P(x). ===*/
		if (FFH_verbosity > 0)
			std::cout << "\t\t Fourier transform of P (" << m << " th stage)\n";
		for (size_t j = 0; j < FourierTransform_P.size(); j++)
			fftw_execute(FourierTransform_P[j]);

	}
}

void FFH_base::DefineProblem(const TwoPhaseDispersion & src, const ConstantTensor & /*matrix_phase*/, const ConstantTensor & /*inclusion_phase*/)
{
	bool SameSimulationBox = this->DefineDispersion(src);
	
	if (!SameSimulationBox) {
		update_propagator = false;	update_plan = false;
		Propagator_Fourier = TensorField(d, 2*rank_of_fields, Geometry.GetResolutions());
		ks = src.GetKs();
	}
	/* define fields */
	P_Fourier = TensorField(d, rank_of_fields, Geometry.GetResolutions());
	Response_Fourier = TensorField(d, rank_of_fields, Geometry.GetResolutions());
	P_direct = TensorField(d, rank_of_fields, Geometry.GetResolutions());
	Response_direct = TensorField(d, rank_of_fields, Geometry.GetResolutions());
	m = 0;
}


void FFH_base::SetMethod(bool strong_contrast)
{
	if (!(use_strong_contrast = strong_contrast)) {
		/* redefine L tensor. */
		this->update_propagator = false;
	}
	this->use_strong_contrast = strong_contrast;
}

void FFH_base::Prepare(const ConstantTensor & applied)
{
	if(applied.d == d && applied.rank == rank_of_fields)
		this->Appl_field = applied;
	else {
		std::cerr << "FFH_base::Prepare() The rank of applied should be identical to the another tensor fields.\n";
	}

	/* Compute the plan for FFT */
	if (!update_plan) {
		std::cout << "fftw_plan is being prepared\n";
		FourierTransform_P.clear();
		InverseFourierTransform_Response.clear();

		FourierTransform_P.resize((size_t)pow(d, rank_of_fields));
		InverseFourierTransform_Response.resize((size_t)pow(d, rank_of_fields));
		
		if (d == 2) {
			for (size_t i = 0; i < FourierTransform_P.size(); i++) {
				/* Fast Fourier transform plan for the ith component. */
				fftw_plan_with_nthreads(omp_get_num_threads());
				FourierTransform_P[i] = fftw_plan_dft_2d(P_Fourier.num_points[0], P_Fourier.num_points[1],
					P_direct[i], P_Fourier[i], -1, FFTW_ESTIMATE);

				fftw_plan_with_nthreads(omp_get_num_threads());
				InverseFourierTransform_Response[i] = fftw_plan_dft_2d(Response_Fourier.num_points[0], Response_Fourier.num_points[1],
					Response_Fourier[i], Response_direct[i], 1, FFTW_ESTIMATE);
			}
		}
		else if (d == 3) {
			for (size_t i = 0; i < FourierTransform_P.size(); i++) {
				/* Fast Fourier transform plan for the ith component. */
				fftw_plan_with_nthreads(omp_get_num_threads());
				FourierTransform_P[i] = fftw_plan_dft_3d(P_Fourier.num_points[0], P_Fourier.num_points[1], P_Fourier.num_points[2],
					P_direct[i], P_Fourier[i], -1, FFTW_ESTIMATE);

				fftw_plan_with_nthreads(omp_get_num_threads());
				InverseFourierTransform_Response[i] = fftw_plan_dft_3d(Response_Fourier.num_points[0], Response_Fourier.num_points[1], Response_Fourier.num_points[2],
					Response_Fourier[i], Response_direct[i], 1, FFTW_ESTIMATE);
			}
		}
		else
			std::cerr << "unsupported dimension!!\n";

		update_plan = true;
	}

	/* Compute the Fourier transform of the Green's function */
	if (!update_propagator) {
		if (use_strong_contrast)
			this->ComputePropagator_strong();
		else
			this->ComputePropagator_weak();

		update_propagator = true;
	
		/*test*/
		//TensorField test;
		//InnerProduct(test, ks, Propagator_Fourier, 1);
		//test.Print_shifted_origin("k.G_untruncated");
	}
	/* Compute the constant tensor L for the constitutive relation. */
	if (use_strong_contrast)
		this->ComputeL_strong();
	else
		this->ComputeL_weak();

	m = 0;
	/* Initialize E field
	*	E(x) = e_direction. */
	Response_direct = TensorField(d, rank_of_fields, Propagator_Fourier.num_points);
	for (size_t x = 0; x < Response_direct.num_total_pixels; x++) {
		for (size_t i = 0; i < Appl_field.total_comp; i++) {
			Response_direct.field[i][x][0] = Appl_field[i].real();
			Response_direct.field[i][x][1] = Appl_field[i].imag();
		}
	}

	/* Initialize P field
	*	P(x) = (epsilon_2 - epsilon1) * I(x) * E(x). */
	P_direct = TensorField(d, rank_of_fields, Propagator_Fourier.num_points);
	if (!Geometry.update_indicator)
		Geometry.ComputeIndicator();
	/*
#pragma omp parallel for 
	for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
		size_t x_idx_domain = Geometry.prt_indicator[i], x_idx_P_direct;
		std::vector<size_t> comp_in_packing;// , comp_in_padded;
		Geometry.Index2Components(x_idx_domain, comp_in_packing);
		x_idx_P_direct = P_direct.at(comp_in_packing);
		for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
			std::complex<double> P = InnerProduct(x_idx_P_direct, comp, this->L, Response_direct, rank_of_fields);
			P *= Geometry.WeightedPhaseIndicator(i);
			P_direct[comp][x_idx_P_direct][0] = P.real();	P_direct[comp][x_idx_P_direct][1] = P.imag();
		}
	}*/

	/* @ phase 1*/
	InnerProduct(P_direct, L1, Response_direct, rank_of_fields);
#pragma omp parallel for
	for (size_t i = 0; i < Geometry.prt_indicator.size(); i++) {
		/* @ phase 2*/
		size_t x = Geometry.prt_indicator[i];
		for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
			std::complex<double> P = InnerProduct(x, comp, this->L2, Response_direct, rank_of_fields);
			P *= Geometry.WeightedPhaseIndicator(i);
			//#pragma omp atomic
			P_direct[comp][x][0] = P.real();
			//#pragma omp atomic
			P_direct[comp][x][1] = P.imag();
		}
	}
	/* compute the Fourier transform of the polarizability. P^1 */
	for (size_t i = 0; i < FourierTransform_P.size(); i++) {
		fftw_execute(FourierTransform_P[i]);
	}
}

void FFH_base::Iterate(size_t steps)
{
	if (use_strong_contrast)
		Iterate_strong(steps);
	else
		Iterate_weak(steps);
}

void FFH_base::Run(double target_tol_field, double target_tol_div, int num_iters)
{
	bool good_solution = false;
	this->is_converging = true;
	double tol_prev = 0.0L;
	/* num_iters < 0 is the sentinel for "iterate up to the built-in limit". */
	int cap = (num_iters < 0) ? static_cast<int>(num_limit) : num_iters;
	int count = 0, iter = 0;
	while (!good_solution && (iter < cap) && is_converging) {
		Iterate();
		tol_prev = Tol_Field;
		if (CheckConvergence(target_tol_field)) {
			good_solution = CheckDivergenceFree(target_tol_div);
		}

		/* check whether the solution is really converging. */
		if ((iter > 30)) {
			if((tol_prev < Tol_Field)&&(Tol_Field > target_tol_field))
				count++;

			if (Tol_Field < Tol_div) {
				std::cout << "The divergence of this solution will not decrease anymore...\n";
				is_converging = true;
				return;
			}
		}
		if (count > 3)
			is_converging = false;
		iter ++;
	}
	if (iter == cap)
		std::cout << "solutions failed to converge!!\n";
	else if (!is_converging)
		std::cout << "solutions don't look like converging: quit at the " << iter << "-th stage\n";
	if (FFH_verbosity > 0) {
		std::cout << "Error in field = " << Tol_Field << ",\tError in div = " << Tol_div << "\n";
		std::cout << "Routine ends at the " << iter << "-th stage.\n";
	}
}


bool FFH_base::CheckConvergence(double target_tolerence)
{
	Response_direct_prev -= Response_direct;
	
	double error = 0.0L, mean = 0.0L;
	/* compute the squared sum of the difference between response fields... */
#pragma omp parallel for
	for (size_t x = 0; x < Response_direct_prev.num_total_pixels; x++) {
		double sum_over_comp = 0.0L;
		double sum_over_comp2 = 0.0L;
		for (size_t i = 0; i < Response_direct_prev.total_comp; i++) {
			sum_over_comp += Modulus2(Response_direct_prev[i][x]);
			sum_over_comp2 += Modulus2(Response_direct[i][x]);
		}
#pragma omp atomic
		error += sum_over_comp;
#pragma omp atomic
		mean += sum_over_comp2;
	}
	
	error = sqrt(error/mean);
	Tol_Field = error;
	if (FFH_verbosity == 2) {
		std::cout << "\t\terrors in the local cavity field = " << error << "\n";
	}
	//F_direct_prev += F_direct;	/**/
	return error < target_tolerence;
}

bool FFH_base::CheckDivergenceFree(double target_tolereance)
{
	TensorField div_free_field, div;
	double error_in_div = 0.0L;
	double mean = 0.0L;

	this->ComputeDivergenceFreeField_Fourier(div_free_field);
	InnerProduct(div, this->ks, div_free_field, rank_of_fields);
	/* compute the squared sum of the divergence in the Fourier space. */
#pragma omp parallel for schedule(guided)
	for (size_t x = 0; x < div.num_total_pixels; x++) {
		double k2 = 0.0L;
		for (size_t i = 0; i < ks.total_comp; i++)
			k2 += Modulus2(ks[i][x]);

		for (size_t i = 0; i < div.total_comp; i++) {
			if (k2 > 1e-10) {
				double mod2 = Modulus2(div[i][x]) / sqrt(k2);
#pragma omp atomic
				error_in_div += mod2;
			}
		}

		for (size_t i = 0; i < div_free_field.total_comp; i++) {
			double mod = Modulus2(div_free_field[i][x]);
#pragma omp atomic
			mean += mod;
		}
	}
	mean = sqrt(mean);

	error_in_div = sqrt(error_in_div) / mean;
	Tol_div = error_in_div;
	if (FFH_verbosity == 2) {
		std::cout << "\t\t Error in divergence free field = " << error_in_div << "\n";
	}

	return (error_in_div < target_tolereance);
}


void FFH_base::Print(const std::string & prefix) const
{
	TensorField div_free, field;
	ComputeDivergenceFreeField_Direct(div_free);
	ComputeField_Direct(field);

	div_free.Print(prefix + name_div_free_field);
	field.Print(prefix + name_field);
}

void FFH_base::Print_Fourier(const std::string & prefix) const
{
	TensorField div_free, field;
	ComputeDivergenceFreeField_Fourier(div_free);
	ComputeField_Fourier(field);

	div_free.Print(prefix + name_div_free_field + "_Fourier");
	field.Print(prefix + name_field + "_Fourier");
}

void FFH_base::Print_Grid(const std::string & prefix) const
{
	this->Geometry.Print(prefix);
}

void FFH_base::Print_GreenFn_2D(const std::string & prefix) const
{
	for (size_t i = 0; i < Propagator_Fourier.total_comp; i++) {
		Propagator_Fourier.Print_shifted_origin(prefix + "_GreenFunction_", i, 0);
		Propagator_Fourier.Print_shifted_origin(prefix + "_GreenFunction_", i, 1);
	}
}


/* StaticDielectricProblem. */
StaticDielectricProblem::StaticDielectricProblem() : FFH_base()
{
	name_div_free_field = "D";	name_field = "E"; rank_of_fields = 1;
}

StaticDielectricProblem::StaticDielectricProblem(const TwoPhaseDispersion & src, const std::complex<double> &epsilon1, const std::complex<double> & epsilon2)
	: StaticDielectricProblem()
{
	this->epsilon1 = epsilon1;	this->epsilon2 = epsilon2;
	this->DefineProperties(epsilon1, epsilon2);
	this->DefineProblem(src, epsilon1*Identity(src.packing.GetDimension()), epsilon1*Identity(src.packing.GetDimension()));
}

void StaticDielectricProblem::DefineProblem(const TwoPhaseDispersion & src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase)
{
	FFH_base::DefineProblem(src, matrix_phase, inclusion_phase);
	//TODO: tensor expression for D-tensor
}

void StaticDielectricProblem::DefineProperties(const std::complex<double>& epsilon1, const std::complex<double>& epsilon2)
{
	if (this->epsilon1 == epsilon1) {
		this->update_propagator = false;
	}
	this->epsilon1 = epsilon1;
	this->epsilon2 = epsilon2;

	if (this->optimized_reference) {
		/* define the reference phase s.t. the median of L(x) is 0.*/
		if (use_strong_contrast) {
			double d = this->d;
			std::complex<double> temp = ((d - 2.0) * (epsilon1 + epsilon2) + sqrt(-16.0 * (1.0 - d) * epsilon1 * epsilon2 + pow((d - 2.0) * (epsilon1 + epsilon2), 2))) / (4.0 * (d - 1.0));
			if (std::abs(temp) < 1e-10 ){
				this->epsilon0 = 0.5*(epsilon2 + epsilon1);
			}
			else{
				this->epsilon0 = temp;
			}
		}
		else {
			this->epsilon0 = 0.5*(epsilon2 + epsilon1);
		}
	}
	else {
		/* reference phase = matrix phase. */
		this->epsilon0 = epsilon1;
	}


	this->D = 1.0 / ((double)d * epsilon0) * Identity(d);
}

std::complex<double> StaticDielectricProblem::ComputeEffectiveEpsilon() const
{
	
	if (use_strong_contrast) {
		ConstantTensor F_av = MeanField(Response_direct),
			P_av = MeanField(P_direct);
		if (FFH_verbosity > 1) {
			std::cout << "F_av = ";
			for (int i = 0; i < d; i++)
				std::cout << "[" << F_av[i] << "]\t";
			std::cout << "\n";

			std::cout << "P_av = ";
			for (int i = 0; i < d; i++)
				std::cout << "[" << P_av[i] << "]\t";
			std::cout << "\n";
		}

		std::complex<double> F_dot_P(0.0L, 0.0L), P2(0.0L, 0.0L);
		for (int i = 0; i < d; i++) {
			F_dot_P += (F_av[i] * conj(P_av[i]));
			P2 += abs(P_av[i])*abs(P_av[i]);
		}
		std::complex<double> Le = P2 / F_dot_P, epsilon_eff;
		epsilon_eff = epsilon0 * (1.0 - (double)d + (double)d / (1.0 - Le / ((double)d*epsilon0)));
		return epsilon_eff;
	}
	else {
		ConstantTensor E_av = MeanField(Response_direct),
			P_av = MeanField(P_direct);
		if (FFH_verbosity > 1) {
			std::cout << "E_av = ";
			for (int i = 0; i < d; i++)
				std::cout << "[" << E_av[i] << "]\t";
			std::cout << "\n";

			std::cout << "P_av = ";
			for (int i = 0; i < d; i++)
				std::cout << "[" << P_av[i] << "]\t";
			std::cout << "\n";
		}
		
		std::complex<double> P_dot_E(0.0L, 0.0L), E2(0.0L, 0.0L);
		for (int i = 0; i < d; i++) {
			P_dot_E += (P_av[i] * conj(E_av[i]));
			E2 += abs(E_av[i])*abs(E_av[i]);
		}

		std::complex<double> epsilon_eff = (P_dot_E / E2) + epsilon0;
		return epsilon_eff;
	}

}

void StaticDielectricProblem::ComputeDivergenceFreeField_Fourier(TensorField & out) const
{
	out = TensorField(this->P_Fourier);	//Polarization!!
	if (use_strong_contrast) {
		//Response = electric field + D:P
		//D_Fourier = epsilon_1*F_Fourier + (1.0-1.0/d)
		out *= (1.0 - 1.0 / (double)d) / epsilon0;
		out += Response_Fourier;
		out *= epsilon0;
	}
	else {
		//Response = electric field
		out += epsilon0*Response_Fourier;
	}
}

void StaticDielectricProblem::ComputeDivergenceFreeField_Direct(TensorField & out) const
{
	out = TensorField(this->P_direct);	//Polarization!!
	if (use_strong_contrast) {
		//Response = electric field + D:P
		//D_Fourier = epsilon_1*F_Fourier + (1.0-1.0/d)
		out *= (1.0 - 1.0 / (double)d) / epsilon0;
		out += Response_direct;
		out *= epsilon0;
	}
	else {
		//Response = electric field
		out += epsilon0 * Response_direct;
	}
}

void StaticDielectricProblem::ComputeField_Fourier(TensorField & out) const
{
	if (use_strong_contrast) {
		//Response = F = E + D*P
		// E = Response - D*P
		out = TensorField(this->P_Fourier);
		out *= (-1.0 / ((double)d*epsilon0));
		out += Response_Fourier;
	}
	else {
		//Response = electric field -> done!
		out = TensorField(this->Response_Fourier);
	}
}

void StaticDielectricProblem::ComputeField_Direct(TensorField & out) const
{
	if (use_strong_contrast) {
		//Response = F = E + D*P
		// E = Response - D*P
		out = TensorField(this->P_direct);
		out *= (-1.0 / ((double)d*epsilon0));
		out += Response_direct;
	}
	else {
		//Response = electric field -> done!
		out = TensorField(this->Response_direct);
	}
}

void StaticDielectricProblem::ComputePropagator_weak()
{
	this->Propagator_Fourier = TensorField(d, 2 * rank_of_fields, Geometry.GetResolutions());

#pragma omp parallel for
	for (size_t k_index = 0; k_index < Propagator_Fourier.num_total_pixels; k_index++) {
		double k2 = 0.0L;
		for (int i = 0; i < d; i++)
			k2 += Modulus2(ks[i][k_index]);
		
		if (k2 > 1e-10) {
			for (int i = 0; i < d; i++) {
				for (int j = 0; j < d; j++) {
					std::complex<double> val = (-1.0)* ks[i][k_index][0] * ks[j][k_index][0] / k2 / epsilon0;

					Propagator_Fourier[j + d * i][k_index][0] = val.real();
					Propagator_Fourier[j + d * i][k_index][1] = val.imag();
				}
			}
		}
	}

}

void StaticDielectricProblem::ComputeL_weak()
{
	this->L1 = (this->epsilon1 - this->epsilon0)*Identity(d);
	this->L2 = (this->epsilon2 - this->epsilon0)*Identity(d);
}

void StaticDielectricProblem::ComputePropagator_strong()
{
	this->ComputePropagator_weak();
	/* add the D constant tensor... */
	//std::complex<double> D = 1.0 / ((double)d * epsilon1);

#pragma omp parallel for
	for (size_t k_index = 0; k_index < Propagator_Fourier.num_total_pixels; k_index++) {
		for (int i = 0; i < d; i++) {
#pragma omp atomic
			this->Propagator_Fourier[i*(1 + d)][k_index][0] += D[i*(1 + d)].real();
#pragma omp atomic
			this->Propagator_Fourier[i*(1 + d)][k_index][1] += D[i*(1 + d)].imag();
		}
		
	}
}

void StaticDielectricProblem::ComputeL_strong()
{
	double d = this->d;
	//this->L = this->epsilon1* d*(this->epsilon2 - this->epsilon1) / (this->epsilon2 + (d - 1.0)*this->epsilon1) *Identity(this->d);
	this->L1 = this->epsilon0* d*(this->epsilon1 - this->epsilon0) / (this->epsilon1 + (d - 1.0)*this->epsilon0) *Identity(this->d);
	this->L2 = this->epsilon0* d*(this->epsilon2 - this->epsilon0) / (this->epsilon2 + (d - 1.0)*this->epsilon0) *Identity(this->d);
}

DynamicDielectricProblem::DynamicDielectricProblem(): StaticDielectricProblem(){}

DynamicDielectricProblem::DynamicDielectricProblem(const TwoPhaseDispersion & src, const std::complex<double>& epsilon1, const std::complex<double>& epsilon2):
	StaticDielectricProblem(src, epsilon1, epsilon2)
{}

void DynamicDielectricProblem::DefineProblem(const TwoPhaseDispersion & src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector & k1)
{
	StaticDielectricProblem::DefineProblem(src, matrix_phase, inclusion_phase);
	this->SetWavevector(k1);
}

void DynamicDielectricProblem::DefineProperties(const std::complex<double>& epsilon1, const std::complex<double>& epsilon2, const GeometryVector & k1)
{
	StaticDielectricProblem::DefineProperties(epsilon1, epsilon2);
	SetWavevector(k1);
}

void DynamicDielectricProblem::SetWavevector(const GeometryVector & k1)
{
	GeometryVector k0 = sqrt(abs(epsilon0 / epsilon1))*k1;
	if (k1.Dimension == d) {
		if (this->k1.Dimension != d || this->k1 != k1) {
			this->update_propagator = false;
		}
		this->k1 = k0;
		/* shift ks by k1. */
		ks = Geometry.GetKs();
#pragma omp parallel for
		for (size_t k_indx = 0; k_indx < ks.num_total_pixels; k_indx++) {
			for (int i = 0; i < d; i++) {
#pragma omp atomic
				ks[i][k_indx][0] += k0.x[i];
			}
		}
		this->update_propagator = false;
	}
	else {
		std::cerr << "k1 should have the same dimensions as the system.\n";
	}

}


ConstantTensor DynamicDielectricProblem::Hfunction_direct(const GeometryVector & position, const GeometryVector & k1, int d, double dx) {
	ConstantTensor H_ij(d, 2);
	double k = sqrt(k1.Modulus2());
	double kr = sqrt(k1.Modulus2()*position.Modulus2());
	double r = sqrt(position.Modulus2());

	std::complex<double> I(0.0L, 1.0L);
	if (r < 0.5*dx) {
		std::complex<double> G1;
		double dv = 1.0, dr;
		for (int i = 0; i < d; i++)
			dv *= dx;

		if (d == 2) {
			dr = dx/sqrt(pi);
			G1 = (-0.5 + I * 0.25*k*dr*pi *boost::math::cyl_hankel_1(1, k*dr))/this->epsilon1/dv;
		}
		else if (d == 3) {
			dr = dx / cbrt(4.0*pi / 3.0);
			G1 = 2.0 / 3.0*(-1.0 + exp(I*k*dr)*(1.0 - I * k*dr)) /(this->epsilon1)/dv;
		}
		H_ij = G1 * Identity(d);
	}
	else {
		std::complex<double> G1, G2;
		ConstantTensor rr(d, 2);
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				rr[j + d * i] = position.x[i] * position.x[j] / (r*r);
			}
		}

		if (d == 2) {
			G1 = 0.25*k*I*(-boost::math::cyl_hankel_1(1, kr) / r + k * boost::math::cyl_hankel_1(0, kr));
			G2 = 0.25*k*k*I*boost::math::cyl_hankel_1(2, kr);
		}
		else if (d == 3) {
			G1 = (-1.0 + I * kr + kr * kr)*exp(I*kr) / (4.0*pi*r*r*r);
			G2 = (3.0 - 3.0*I*kr - kr * kr)*exp(I*kr) / (4.0*pi*r*r*r);
		}

		H_ij = 1.0/this->epsilon1*(G1 * Identity(d) + G2 * rr);
	}
	return exp(-I*k1.Dot(position))*H_ij;
}

inline ConstantTensor DynamicDielectricProblem::DyadicGreenFunction_Fourier(const GeometryVector & q, double qq, double kk)
{
	//double qq = q.Modulus2(), kk = k1.Modulus2();
	ConstantTensor out = qq*Identity(d);
	for (int i = 0; i < d; i++)
		for (int j = 0; j < d; j++)
			out[j+d*i] -= (q.x[i] * q.x[j]);
	
	
	out = 1.0 / (this->epsilon0 * (qq - kk)) * out;
	out = out -1.0/(this->epsilon0)*Identity(d);
	return out;
}

inline ConstantTensor DynamicDielectricProblem::DyadicGreenFunction_Fourier_regularized(const GeometryVector & q, double qq, double kk, double epsilon)
{
	//double qq = q.Modulus2(), kk = k1.Modulus2();
	ConstantTensor out = std::complex<double>(qq, 0.0L) * Identity(d);
	for (int i = 0; i < d; i++)
		for (int j = 0; j < d; j++)
			out[j + d * i] -= (q.x[i] * q.x[j]);

	out = 1.0 / (this->epsilon1 * (qq - kk - std::complex<double>(0.0L, epsilon))) * out;
	out = out - 1.0 / (this->epsilon1)*Identity(d);

	//if (abs(kk - qq) < 4 * pi*pi / pow(Geometry.packing.PeriodicVolume(), 2.0 / (double)this->d)) {
	//	out = 1.0 / this->epsilon1 * ((double)this->P_Fourier.num_total_pixels * pi *0.5/ sqrt(kk)) * out;
	//}
	//else {
	//	out = 1.0 / (this->epsilon1 * (qq - kk))*out;
	//}
	//out = out - 1.0 / (this->epsilon1)*Identity(d);
	return out;
}


void DynamicDielectricProblem::ComputeG_Fourier_test()
{
	/* Compute the Green's function in the direct space. */
	std::vector<size_t> pixels;
	for (int i = 0; i < d; i++)
		pixels.emplace_back(2*Propagator_Fourier.num_points[i]);
	TensorField H(d, 2*rank_of_fields, pixels);

	std::vector<double> dx;
	{
		std::vector<size_t> resol = Geometry.GetResolutions();
		for (int i = 0; i < d; i++)
			dx.push_back(sqrt(Geometry.packing.GetBasisVector(i).Modulus2()) / (double)resol[i]);

	}
	double dv = 1.0;
	for (int i = 0; i < d; i++)
		dv *= dx[i];
	double min_dist = *std::min_element(dx.begin(), dx.end());
#pragma omp parallel for
	for (size_t x = 0; x < H.num_total_pixels; x++) {
		std::vector<size_t> comp;			H.GetPositionComponent(x, comp);

		GeometryVector pos(static_cast<int>(d));

		//Geometry.Index2Components(x, comp);
		for (int i = 0; i < d; i++) {
			if (H.num_points[i] < 2 * comp[i])
				pos.x[i] = dx[i] * ((double)comp[i] - (double)H.num_points[i]);
			else
				pos.x[i] = dx[i] * (double)(comp[i]);
		}

		ConstantTensor H_ij = Hfunction_direct(pos, k1, d, min_dist);
		if(x==0)
			H_ij = H_ij - 1.0/dv * D ;
	
		for (size_t i = 0; i < H_ij.total_comp; i++) {
			H[i][x][0] = H_ij[i].real();
			H[i][x][1] = H_ij[i].imag();
		}
	}


	H.Print_shifted_origin("H_function_real", -1, 0);
	H.Print_shifted_origin("H_function_real", -1, 1);



	TensorField Propagator_double(H);
	std::vector<fftw_plan> fft_G;
	fft_G.resize(H.total_comp);
	if (d == 2) {
		for (size_t i = 0; i < H.total_comp; i++) {
			fft_G[i] = fftw_plan_dft_2d(H.num_points[0], H.num_points[1], 
				H[i], Propagator_double[i], -1, FFTW_ESTIMATE);
		}
	}
	else if (d == 3) {
		for (size_t i = 0; i < H.total_comp; i++) {
			fft_G[i] = fftw_plan_dft_3d(H.num_points[0], H.num_points[1], H.num_points[2], 
				H[i], Propagator_double[i], -1, FFTW_ESTIMATE);
		}
	}
	
	/* Compute its Fourier transform. */
	for (size_t i = 0; i < fft_G.size(); i++)
		fftw_execute(fft_G[i]);
	for (size_t i = 0; i < fft_G.size(); i++)
		fftw_destroy_plan(fft_G[i]);

#pragma omp parallel for
	for (size_t q_idx = 0; q_idx < Propagator_Fourier.num_total_pixels; q_idx++) {
		std::vector<size_t> q_comp;
		Geometry.Index2Components(q_idx, q_comp);
		for (int i = 0; i < d; i++) {
			q_comp[i] *= 2;
			//if (2 * q_comp[i] > Propagator_Fourier.num_points[i])
			//	q_comp[i] += Propagator_Fourier.num_points[i];
		}
		size_t q_idx_doubled = Propagator_double.at(q_comp);
		for (size_t j = 0; j < Propagator_double.total_comp; j++) {
			Propagator_Fourier[j][q_idx][0] = Propagator_double[j][q_idx_doubled][0];
			Propagator_Fourier[j][q_idx][1] = Propagator_double[j][q_idx_doubled][1];
		}
	}


	Propagator_Fourier *= dv;
	
	TensorField div;
	InnerProduct(div, ks, Propagator_Fourier, 1);
	div.Print_shifted_origin("test_div");
}

void DynamicDielectricProblem::ComputeG_Fourier_test2()
{

	Propagator_Fourier = TensorField(d, 2 * rank_of_fields, Geometry.GetResolutions());
	//TODO:: generalize it
	double L = 0, k_wave = sqrt(k1.Modulus2());
	double factor = Geometry.packing.PeriodicVolume() / Propagator_Fourier.num_total_pixels;
	for (int i = 0; i < d; i++) {
		L += Geometry.packing.GetBasisVector(i).Modulus2();
	}
	L = sqrt(L);

	for (size_t k_idx = 0; k_idx < ks.num_total_pixels; k_idx++) {
		double q = 0;
		for (int i = 0; i < d; i++)
			q += Modulus2(ks[i][k_idx]);
		q = sqrt(q);
		
		ConstantTensor val_at_k = (q * q)*Identity(d);
		for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++)
				val_at_k[j+d*i] -= (ks[i][k_idx][0] * ks[j][k_idx][0]);

		//val_at_k = - D*Identity(d) + TruncatedGreenFn_Helmholtz(d, q, L, k_wave)*val_at_k;
		std::complex<double> val = factor*TruncatedGreenFn_Helmholtz(d, q, L, k_wave);
		val_at_k = -1.0/this->epsilon0 * Identity(d) + val*val_at_k;
		for (size_t i = 0; i < Propagator_Fourier.total_comp; i++) {
			Propagator_Fourier[i][k_idx][0] = val_at_k[i].real();
			Propagator_Fourier[i][k_idx][1] = val_at_k[i].imag();
		}
	}

	//TensorField h(d, 0, Geometry.GetResolutions());
	//for (int k_idx = 0; k_idx < ks.num_total_pixels; k_idx++) {
	//	double q = 0;
	//	double k = sqrt(k1.Modulus2());
	//	for (int i = 0; i < d; i++)
	//		q += (ks[i][k_idx][0] - k)*(ks[i][k_idx][0] - k);//Modulus2(ks[i][k_idx]);
	//	q = sqrt(q);
	//	std::complex<double> val = TruncatedGreenFn_Helmholtz(d, q, L, k_wave);
	//	val *= (q*q);
	//	h[0][k_idx][0] = val.real();
	//	h[0][k_idx][1] = val.imag();
	//}
	//h.Print_shifted_origin("truncated_gH_real",0,0);
	//h.Print_shifted_origin("truncated_gH_imag",0,1);
	//ks.Print_shifted_origin("test");
}

void DynamicDielectricProblem::ComputeG_Fourier_test3()
{
	std::vector<double> dk(d, 0.0L);
	for (int i = 0; i < d; i++)
		dk[i] = 1.0*pi / sqrt(Geometry.packing.GetBasisVector(i).Modulus2());

	double kk = k1.Modulus2(), dk2 = dk[0]*dk[0];
	Propagator_Fourier = TensorField(d, 2 * rank_of_fields, Geometry.GetResolutions());
	/* Note that ks is already shifted by k1. */
#pragma omp parallel for
	for (size_t q_indx = 0; q_indx < ks.num_total_pixels; q_indx++) {
		GeometryVector q(static_cast<int> (d));
		for (int i = 0; i < d; i++)
			q.x[i] = ks[i][q_indx][0];

		double qq = q.Modulus2();
		ConstantTensor val;
		/*if (abs(qq - kk) < 1e-15) {
			val = -1.0/this->epsilon0 * Identity(d);
		}*/
		if(abs(qq - kk) < dk2) {
			val = qq * Identity(d);
			for (int i = 0; i < d; i++)
				for (int j = 0; j < d; j++)
					val[j + d * i] -= (kk/qq*q.x[i] * q.x[j]);

			val = ((double)ks.num_total_pixels * std::complex<double>(0.0, pi / sqrt(kk) / 2.0)) * val;
			val = val -1.0 / (this->epsilon0)*Identity(d);
			//double count = 0, qq_sample;
			//GeometryVector q_sample(static_cast<int>(d));

			///* Compute the average via the Monte Carlo method */
			//for (int i = 0; i < num_sampling_pts; i++) {
			//	/* */
			//	for (int j = 0; j < d; j++)
			//		q_sample.x[j] = dk[j] * (rng.RandomDouble() - 0.5);

			//	q_sample = q_sample + q;
			//	qq_sample = q_sample.Modulus2();
			//	if (abs(q_sample.Modulus2() - kk) > 1e-10) {
			//		val = val + DyadicGreenFunction_Fourier(q_sample, qq_sample, kk);
			//		count++;
			//	}
			//}
			//
			//val = 1.0 / (double)count * val;

		}
		else {
			val = DyadicGreenFunction_Fourier(q, qq, kk);
		}

		for (size_t i = 0; i < val.total_comp; i++) {
			Propagator_Fourier[i][q_indx][0] = val[i].real() + this->D[i].real();
			Propagator_Fourier[i][q_indx][1] = val[i].imag() + this->D[i].imag();
		}	
	}

	//TensorField div;
	//InnerProduct(div, ks, Propagator_Fourier, 1);
	//div.Print("test_div");
}

void DynamicDielectricProblem::ComputeG_Fourier_test4(double epsilon)
{
	std::vector<double> dk(d, 0.0L);
	for (int i = 0; i < d; i++)
		dk[i] = 1.0*pi / sqrt(Geometry.packing.GetBasisVector(i).Modulus2());

	double kk = k1.Modulus2();
	Propagator_Fourier = TensorField(d, 2 * rank_of_fields, Geometry.GetResolutions());
	/* Note that ks is already shifted by k1. */
#pragma omp parallel for
	for (size_t q_indx = 0; q_indx < ks.num_total_pixels; q_indx++) {
		GeometryVector q(static_cast<int> (d));
		for (int i = 0; i < d; i++)
			q.x[i] = ks[i][q_indx][0];

		double qq = q.Modulus2();
		ConstantTensor val;
		
		val = DyadicGreenFunction_Fourier_regularized(q, qq, kk, epsilon);
		
		for (size_t i = 0; i < val.total_comp; i++) {
			Propagator_Fourier[i][q_indx][0] = val[i].real();
			Propagator_Fourier[i][q_indx][1] = val[i].imag();
		}
	}

	TensorField div;
	InnerProduct(div, ks, Propagator_Fourier, 1);
	div.Print("test_div");
}

/* initial solution. */
void DynamicDielectricProblem::ComputePropagator_weak()
{
	double kk = k1.Modulus2();
	Propagator_Fourier = TensorField(d, 2 * rank_of_fields, Geometry.GetResolutions());
	/* Note that ks is already shifted by k1. */
#pragma omp parallel for
	for (size_t q_indx = 0; q_indx < ks.num_total_pixels; q_indx++) {
		GeometryVector q(static_cast<int> (d));
		for (int i = 0; i < d; i++)
			q.x[i] = ks[i][q_indx][0];

		double qq = q.Modulus2();
		ConstantTensor val;
		if (abs(qq - kk) < 1e-15) {
			val = -1.0 / this->epsilon0 * Identity(d);
		}
		else {
			val = DyadicGreenFunction_Fourier(q, qq, kk);
		}

		for (size_t i = 0; i < val.total_comp; i++) {
			Propagator_Fourier[i][q_indx][0] = val[i].real();
			Propagator_Fourier[i][q_indx][1] = val[i].imag();
		}
	}
	//Propagator_Fourier *= (1.0/epsilon0);
	//TensorField test;
	//InnerProduct(test, ks, Propagator_Fourier, 1);
	//test.Print("test_k.G");
	//ks.Print("test_ks");
}

void DynamicDielectricProblem::ComputePropagator_strong()
{
	this->ComputePropagator_weak();
#pragma omp parallel for
	for (size_t x = 0; x < Propagator_Fourier.num_total_pixels; x++) {
		for (size_t i = 0; i < Propagator_Fourier.total_comp; i++) {
#pragma omp atomic
			Propagator_Fourier[i][x][0] += this->D[i].real();
#pragma omp atomic
			Propagator_Fourier[i][x][1] += this->D[i].imag();
		}
	}
}

void DynamicDielectricProblem::ShiftFieldByk1(TensorField & out) const
{
	std::vector<GeometryVector> normalized_basis;
	std::vector<size_t> pixels = Geometry.GetResolutions();
	for (int i = 0; i < d; i++)
		normalized_basis.emplace_back((1.0 / (double)pixels[i]) *Geometry.packing.GetBasisVector(i));

	/* multiply the phase factor. */
#pragma omp parallel for
	for (size_t x = 0; x < out.num_total_pixels; x++) {
		std::vector<size_t> n_x;
		Geometry.Index2Components(x, n_x);

		GeometryVector pos(static_cast<int>(d));
		for (int i = 0; i < d; i++)
			pos = pos + (double)n_x[i] * normalized_basis[i];

		std::complex<double> phase_factor = std::exp(std::complex<double>(0.0L, k1.Dot(pos)));

		for (int i = 0; i < d; i++) {
			std::complex<double> value = std::complex<double>(out[i][x][0], out[i][x][1]);
			value *= phase_factor;
			out[i][x][0] = value.real();
			out[i][x][1] = value.imag();
		}

	}
}










inline std::complex<double> TruncatedGreenFn_Helmholtz(int d, double xi, double L, double k1)
{
	double sL = L * xi, kL = L * k1;
	std::complex<double> out(0.0L, 0.0L);
	if (d == 2) {
		if (abs(xi - k1) < 1e-3) {
			out = 0.25 * L*L *pi *std::complex<double>(0.0L, 1.0L) 
				* (boost::math::cyl_bessel_j(0, kL)*boost::math::cyl_hankel_1(0, kL)
					+ boost::math::cyl_bessel_j(1, kL)*boost::math::cyl_hankel_1(1,kL));
		}
		else {
			out = 1.0 + std::complex<double>(0.0, 1.0L)*0.5*pi*
				(sL* boost::math::cyl_bessel_j(1,sL)* boost::math::cyl_hankel_1(0, kL) 
					- kL*boost::math::cyl_bessel_j(0, sL) * boost::math::cyl_hankel_1(1, kL));
			out /= (xi*xi- k1*k1);
		}
	}
	else if (d == 3) {
		if (abs(xi - k1) < 1e-3) {
			if (k1 < 1e-3) {
				out = 0.5*L*L;
			}
			else {
				out = (1.0 - exp(std::complex<double>(0.0L, 2.0*kL)) + std::complex<double>(0.0L, 2.0*kL));
				out /= (4.0*k1*k1);
			}
		}
		else {
			out = (-1.0 + exp(std::complex<double>(0, kL))*(std::complex<double>(cos(sL), -k1/xi *sin(sL))));
			out /= (xi*xi - k1 * k1);
		}
	}

	return out;
}

ConstantTensor HydrostaticStrain(DimensionType d)
{
	ConstantTensor out(d, 2);
	for (int i = 0; i < d; i++)
		out[i*(d + 1)] = 1.0 / (double)d;
	return out;
}

ConstantTensor DeviatoricStrain(DimensionType d)
{
	ConstantTensor out(d, 2);
	for(int i=0; i<d; i++)
		for (int j = 0; j < d; j++) {
			if (i != j) {
				out[i + d * j] = 1.0L;
			}
		}
	return out;
}






Homo_ElasticProblem_static::Homo_ElasticProblem_static(const TwoPhaseDispersion & src, const std::complex<double>& K1, const std::complex<double>& G1, const std::complex<double>& K2, const std::complex<double>& G2)
{
	this->RedefineProblem(src, K1, G1, K2, G2);
}

void Homo_ElasticProblem_static::RedefineProblem(const TwoPhaseDispersion & src, const std::complex<double>& K1, const std::complex<double>& G1, const std::complex<double>& K2, const std::complex<double>& G2)
{
	bool SameSimulationBox = (Geometry.packing.GetDimension() == src.packing.GetDimension());
	d = src.packing.GetDimension();
	/* compare simulation box */
	if (SameSimulationBox) {
		std::vector<size_t> prev_resolutions = Geometry.GetResolutions();
		std::vector<size_t> new_resolutions = src.GetResolutions();
		for (DimensionType i = 0; i < Geometry.packing.GetDimension(); i++) {
			/* compare basis vectors */
			SameSimulationBox &= (Geometry.packing.GetBasisVector(i) == src.packing.GetBasisVector(i));
			/* compare resolutions in each side */
			SameSimulationBox &= (prev_resolutions[i] == new_resolutions[i]);
		}
	}
	if (!SameSimulationBox) {
		this->Geometry = TwoPhaseDispersion(src);
	}
	else {
		this->Geometry.packing = SpherePacking(src.packing);
	}
	double dimension = d;
	if ((this->K1 != K1) || (this->G1 != G1)) {
		update_propagator = false;
	}
	this->K1 = K1;	this->G1 = G1;	
	this->K2 = K2;	this->G2 = G2;
	this->kappa = (K2 - K1) / (K2 + 2.0*(dimension - 1.0) / dimension * G1);
	this->mu = (G2 - G1) / (G2 + G1 * (0.5*dimension*K1 + (dimension +1.0)*(dimension -2.0)/ dimension * G1 ) / (K1 + 2.0*G1));
	
	if (!SameSimulationBox) {
		update_propagator = false;	update_plan = false;
		Propagator_Fourier = TensorField(d, 4, Geometry.GetResolutions());
	}
	/* defin constant tensors */
	delta_C = dimension * (K2 - K1)*Lambda_h(d) + 2.0*(G2 - G1)*Lambda_s(d);

	D = Lambda_h(d) + dimension * (K1 + 2.0*G1) / ((dimension + 2.0)*G1)*Lambda_s(d);
	D = 1.0 / (dimension*K1 + 2.0*(dimension - 1.0)*G1)*D;

	L = kappa * Lambda_h(d) + (dimension + 2.0)*G1 / (dimension*(K1 + 2.0*G1)) *mu*Lambda_s(d);
	L = (dimension*K1 + 2.0*(dimension - 1.0)*G1)*L;
	
	/* define fields */
	P_Fourier = TensorField(d, 2, Geometry.GetResolutions());
	Response_Fourier = TensorField(d, 2, Geometry.GetResolutions());
	P_direct = TensorField(d, 2, Geometry.GetResolutions());
	Response_direct = TensorField(d, 2, Geometry.GetResolutions());
	m = 0;
	
}

void Homo_ElasticProblem_static::SetStrongContrast(bool use_strong_contrast)
{
	if (!(Use_StrongContrast = use_strong_contrast && L.total_comp == (size_t)pow(d, 4)) ){
		/* redefine L tensor. */
		L = ConstantTensor(d, 4);

		double d = this->d;
		if (use_strong_contrast) {
			L = (d * K1 + 2.0*(d - 1.0)*G1)*(kappa *Lambda_h(this->d) + (d + 2.0)*G1 / (d*(K1 + 2.0*G1))*mu*Lambda_s(this->d));
			if (FFH_verbosity == 2)
				std::cout << "L-tensor is defined for the strong-contrast expansions\n";
		}
		else {
			L = d * (K2 - K1)*Lambda_h(this->d) + 2.0*(G2 - G1)*Lambda_s(this->d);
			if (FFH_verbosity == 2)
				std::cout << "L-tensor is defined for the weak-contrast expansions\n";
		}
		this->update_propagator = false;
	}

	this->Use_StrongContrast = use_strong_contrast;
}

void Homo_ElasticProblem_static::Prepare(const ConstantTensor & external_strain)
{
	if (external_strain.d != d) {
		std::cerr << "external_strain has a different dimension!!\n";
		return;
	}

	applied_strain = ConstantTensor(external_strain);

	/* Compute the plan for FFT */
	if (!update_plan) {
		if(FFH_verbosity == 2)
			std::cout << "fftw_plan is being prepared\n";
		
		if (d == 2) {
			// P and Response are rank-2 tensors.
			for (size_t i = 0; i < d*d; i++) {
				/* Fast Fourier transform plan for the ith component. */
				fftw_plan_with_nthreads(omp_get_num_threads());
				FourierTransform_P.emplace_back(fftw_plan_dft_2d(P_Fourier.num_points[0], P_Fourier.num_points[1],
					P_direct.field[i], P_Fourier.field[i], -1, FFTW_ESTIMATE));
				fftw_plan_with_nthreads(omp_get_num_threads());
				InverseFourierTransform_Response.emplace_back(fftw_plan_dft_2d(Response_Fourier.num_points[0], Response_Fourier.num_points[1],
					Response_Fourier.field[i], Response_direct.field[i], 1, FFTW_ESTIMATE));
			}
			update_plan = true;
		}
		else if (d == 3) {
			for (size_t i = 0; i < d*d; i++) {
				/* Fast Fourier transform plan for the ith component. */
				fftw_plan_with_nthreads(omp_get_num_threads());
				FourierTransform_P.emplace_back(fftw_plan_dft_3d(P_Fourier.num_points[0], P_Fourier.num_points[1], P_Fourier.num_points[2],
					P_direct.field[i], P_Fourier.field[i], -1, FFTW_ESTIMATE));
				fftw_plan_with_nthreads(omp_get_num_threads());
				InverseFourierTransform_Response.emplace_back(fftw_plan_dft_3d(Response_Fourier.num_points[0], Response_Fourier.num_points[1], Response_Fourier.num_points[2],
					Response_Fourier.field[i], Response_direct.field[i], 1, FFTW_ESTIMATE));
			}
			update_plan = true;
		}
		else
			std::cerr << "unsupported dimension!!\n";
	}

	/* Compute the Fourier transform of the Green's function */
	if (!update_propagator) {
		k = Geometry.GetKs();

		if (Use_StrongContrast)
			ComputeH_Fourier();
		else
			ComputeG_Fourier();
	}
	m = 0;
	/* Initialize strain field
	*	E(x) = e_direction. */
	Response_direct = TensorField(d, 2, Geometry.GetResolutions());
#pragma omp parallel for
	for (size_t x = 0; x < Response_direct.num_total_pixels; x++) {
		for (size_t i = 0; i < applied_strain.total_comp; i++) {
			Response_direct[i][x][0] = applied_strain[i].real();
			Response_direct[i][x][1] = applied_strain[i].imag();
		}
	}

	/* Initialize P field
	*	P(x) = (epsilon_2 - epsilon1) * I(x) * E(x). */
	P_direct = TensorField(d, 2, Geometry.GetResolutions());
	if (!Geometry.update_indicator)
		Geometry.ComputeIndicator();

	//#pragma omp parallel for 
	for (size_t i = 0; i < Geometry.prt_indicator.size(); i++) {
		size_t x = Geometry.prt_indicator[i];
		for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
			std::complex<double> P = InnerProduct(x, comp, this->L, Response_direct, 2);
			P *= Geometry.WeightedPhaseIndicator(i);
			P_direct[comp][x][0] = P.real();	P_direct[comp][x][1] = P.imag();
		}
	}
}

void Homo_ElasticProblem_static::Iterate(size_t steps) {
	if (Use_StrongContrast)
		Iterate_strong_contrast(steps);
	else
		Iterate_weak_contrast(steps);
};

void Homo_ElasticProblem_static::Run(double tol_field, double tol_div)
{
	bool good_solution = false;
	while (!good_solution) {
		Iterate();
		if (CheckConvergence(tol_field)) {
			good_solution = CheckDivergenceFree(tol_div);
		}
	}
	std::cout << "Routine ends at the " << m << "-th stage.\n";
}

bool Homo_ElasticProblem_static::CheckConvergence(double tol_field)
{
	Response_direct_prev -= Response_direct;
	double error = 0.0L;
	for (size_t x = 0; x < Response_direct_prev.num_total_pixels; x++) {
		double sum_over_comp = 0.0L;
		for (size_t i = 0; i < Response_direct_prev.total_comp; i++) {
			sum_over_comp += Modulus2(Response_direct_prev.field[i][x]);
		}
		error += sqrt(sum_over_comp);
	}
	error /= (double)Response_direct_prev.num_total_pixels;
	if (FFH_verbosity == 2) {
		std::cout << "errors in the local cavity field = " << error << "\n";
	}
	Tol_Field = error;
	//F_direct_prev += F_direct;	/**/
	return (error < tol_field);
}

bool Homo_ElasticProblem_static::CheckDivergenceFree(double tol_div)
{
	double d = this->d;
	/* Compute the Fourier transform of the stress tensor. */
	TensorField sigma_Fourier(P_Fourier);
	ConstantTensor C1 = d * K1*Lambda_h(this->d) + 2.0*G1*Lambda_s(this->d);
	
	if (Use_StrongContrast) {
		/* sigma(x) = C1:f(x) + (I - C1: D ):P(x) */
		ConstantTensor I_C1D = Lambda_h(this->d) + Lambda_s(this->d) - InnerProduct(C1, D, 2);

		for (size_t x = 0; x < sigma_Fourier.num_total_pixels; x++) {
			for (size_t i = 0; i < sigma_Fourier.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, C1, Response_Fourier, 2) + InnerProduct(x, i, I_C1D, P_Fourier, 2);
				sigma_Fourier[i][x][0] = sigma.real();
				sigma_Fourier[i][x][1] = sigma.imag();
			}
		}
	}
	else {
		/* sigma(x) = C1:epsilon(x) + P(x) */
		//TODO:: simplify it...
		for (size_t x = 0; x < sigma_Fourier.num_total_pixels; x++) {
			for (size_t i = 0; i < sigma_Fourier.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, C1, Response_Fourier, 2);
				sigma_Fourier[i][x][0] += sigma.real();
				sigma_Fourier[i][x][1] += sigma.imag();
			}
		}
	}

	TensorField Div_sigma;
	InnerProduct(Div_sigma, k, sigma_Fourier, 1);
	double error = 0.0L;
	for (size_t x = 0; x < Div_sigma.num_total_pixels; x++) {
		double abs_div = 0.0L;
		for (int i = 0; i < d; i++) {
			abs_div += Modulus2(Div_sigma[i][x]);
		}
		error += abs_div;
	}
	double sigma0 = 0.0L;
	for (size_t i = 0; i < sigma_Fourier.total_comp; i++)
		sigma0 += (sigma_Fourier[i][0][0] * sigma_Fourier[i][0][0] + sigma_Fourier[i][0][1] * sigma_Fourier[i][0][1]);

	error = sqrt(error/(double)sigma_Fourier.num_total_pixels) / sqrt(sigma0);
	Tol_div = error;
	if (FFH_verbosity == 2)
		std::cout << "Error in divergence-free field (stress) = " << Tol_div << "\n";
	return error < tol_div;
}

void Homo_ElasticProblem_static::ComputeEffectiveProperties(std::vector<std::complex<double>> & result) const
{
	double d = this->d;
	result.resize(0);
	if (Use_StrongContrast) {
		ConstantTensor f_av = MeanField(Response_direct),
			P_av = MeanField(P_direct);

		std::complex<double> Ke;
		{
			std::complex<double> A = Trace(P_av) / Trace(f_av);
			std::complex<double> kappa_e = A / (d*K1 + 2.0*(d - 1.0)*G1);
			Ke = (d*K1 + 2.0*G1*(d - 1.0)*kappa_e) / (d*(1.0-kappa_e));
		}
		/* deviatoric tensors */
		f_av = f_av - Trace(f_av)*HydrostaticStrain(this->d);
		P_av = P_av - Trace(P_av)*HydrostaticStrain(this->d);
		std::complex<double> Ge;
		{
			std::complex<double> B(0.0L, 0.0L);
			double count = 0;
			for (size_t i = 0; i < f_av.total_comp; i++) {
				if (abs(f_av[i]) > Tol_Field && abs(P_av[i]) > Tol_Field) {
					B += P_av[i] / f_av[i];
					count++;
				}
			}
			if(count > 0)
				B /= count;
	
			std::complex<double> mu_e = B / (d*K1 + 2.0*(d - 1.0)*G1);
			mu_e /= ((d + 2.0)*G1 / d / (K1 + 2.0*G1));


			Ge = -G1*((2.0* d*(2.0* G1 + K1) + 2.0* (-2.0 + d) *(1.0 + d)* G1* mu_e + d*d* K1* mu_e)) / (2.0* d*(2.0* G1 + K1)* (-1.0 + mu_e));
			//Ge = G1 * (2.0*d*(2.0*G1 + K1) + (2.0*(d-2.0)*(d+1.0)*G1 + d*d*K1)*mu_e) / (2.0*d*(2.0*G1+K1)*(1.0-mu_e));
		}
		//TODO: implement for anisotropic case.

		result.push_back(Ke);
		result.push_back(Ge);
	}
	else {
		/* Response = strain! */
		ConstantTensor stress_av(this->d, 2), strain_av(this->d, 2);
		TensorField stress(P_direct);
		
		ConstantTensor C1 = d * K1*Lambda_h(this->d) + 2.0*G1*Lambda_s(this->d);
		//TODO:: simplify it...
		for (size_t x = 0; x < stress.num_total_pixels; x++) {
			for (size_t i = 0; i < stress.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, C1, Response_direct, 2);
				stress[i][x][0] += sigma.real();
				stress[i][x][1] += sigma.imag();
			}
		}

		/* compute the means of stress and strain. */
#pragma omp parallel for 
		for (size_t idx = 0; idx < stress.total_comp; idx++) {
			for (size_t x_idx = 0; x_idx < stress.num_total_pixels; x_idx++) {
#pragma omp critical
				{//TODO:: check whether this is thread-safe.
					stress_av[idx] += std::complex<double>(stress[idx][x_idx][0], stress[idx][x_idx][1]);
					strain_av[idx] += std::complex<double>(Response_direct[idx][x_idx][0], Response_direct[idx][x_idx][1]);
				}
			}
		}
		
		stress_av = 1.0 / (double)stress.num_total_pixels * stress_av;
		strain_av = 1.0 / (double)stress.num_total_pixels * strain_av;
	
		std::complex<double> Ke(0.0L, 0.0L), Ge(0.0L, 0.0L),
			stress_h(0.0L, 0.0L), strain_h(0.0L, 0.0L);

		for (int i = 0; i < d; i++) {
			stress_h += stress_av[i];
			strain_h += strain_av[i];
		}
		Ke = 1.0 / d * stress_h / strain_h;
		
		stress_av = stress_av - stress_h / d * HydrostaticStrain(this->d);
		strain_av = strain_av - strain_h / d * HydrostaticStrain(this->d);
		
		double count = 0;
		for (size_t i = 0; i < strain_av.total_comp; i++) {
			if (abs(stress_av[i]) > 1e-10) {
				Ge += (stress_av[i] / strain_av[i]);
				count++;
			}
		}
		Ge *= (0.5/count);
		//TODO::implement anisotropic effective tensor...
		result.push_back(Ke);
		result.push_back(Ge);
	}
}

void Homo_ElasticProblem_static::Print(const std::string & prefix) const
{
		for (size_t i = 0; i < Response_direct.total_comp; i++) {
			if (Use_StrongContrast) 
				Response_direct.Print((prefix + "_f").c_str(), i);
			else
				Response_direct.Print((prefix + "_strain").c_str(), i);
			P_direct.Print((prefix + "_P").c_str(), i);
		}
}

void Homo_ElasticProblem_static::Print_Fourier(const std::string & prefix) const
{
	for (size_t i = 0; i < Response_Fourier.total_comp; i++) {
		if (Use_StrongContrast)
			Response_Fourier.Print_shifted_origin((prefix + "_f").c_str(), i);
		else
			Response_Fourier.Print_shifted_origin((prefix + "_strain").c_str(), i);
		P_Fourier.Print_shifted_origin((prefix + "_P").c_str(), i);
	}
}

void Homo_ElasticProblem_static::Print_GreenFn_2D(const std::string & prefix) const
{
	for (size_t i = 0; i < Propagator_Fourier.total_comp; i++) {
		Propagator_Fourier.Print_shifted_origin(prefix + "_GreenFunction", i);
	}
}

void Homo_ElasticProblem_static::Iterate_weak_contrast(size_t steps)
{
	const TensorField & G_Fourier = Propagator_Fourier;
	TensorField & epsilon_Fourier = Response_Fourier;
	TensorField & epsilon_direct = Response_direct;

	for (size_t i = 0; i < steps; i++) {
		/*1. Compute FFT of P(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Fourier transform of P (" << m << " th stage)\n";
		for (size_t j = 0; j<P_direct.total_comp; j++)
			fftw_execute(FourierTransform_P[j]);

		m++;
		if (FFH_verbosity > 0) {
			std::cout << "------------------" << m << " th stage ---------------\n";
			std::cout << "\t\t Compute F field \n";
		}
		/* epsilon(k) <- epsilon_0(k) + G(k).P(k)*/
		/*3. Compute \tilde{E}(k) = \tilde{G}(k). \tilde{P}(k) if k \neq 0
		\tilde{E}(k) = \tilde{E_0} */
		InnerProduct(epsilon_Fourier, G_Fourier, P_Fourier, 2);
		/* Add \tilde{epsilon_0} */
		for (size_t i = 0; i < applied_strain.total_comp; i++) {
			Response_Fourier[i][0][0] = applied_strain[i].real() * epsilon_Fourier.num_total_pixels;
			Response_Fourier[i][0][1] = applied_strain[i].imag() * epsilon_Fourier.num_total_pixels;
		}

		/*4. Compute the inverse FFT of \tilde{E}(k). */
		Response_direct_prev = Response_direct;
		if (FFH_verbosity > 0)
			std::cout << "\t\t Inverse Fourier transform of F field\n";
		for (size_t j = 0; j<epsilon_direct.total_comp; j++)
			fftw_execute(InverseFourierTransform_Response[j]);
		epsilon_direct *= (1.0 / P_Fourier.num_total_pixels);	//This step is required because fftw3 do not consider the 1/N factor.

		/*5. Compute P(x) = L(x)* F(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Compute P field\n";
#pragma omp parallel for 
		for (size_t i = 0; i < Geometry.prt_indicator.size(); i++) {
			size_t x = Geometry.prt_indicator[i];
			for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
				std::complex<double> P = InnerProduct(x, comp, this->delta_C, Response_direct, 2);
				P *= Geometry.WeightedPhaseIndicator(i);
				P_direct[comp][x][0] = P.real();	P_direct[comp][x][1] = P.imag();
			}
		}
	}
}

void Homo_ElasticProblem_static::ComputeG_Fourier()
{
	std::complex<double> lambda = K1 - 2.0*G1 / (double)d;	/*Lame parameter */
	std::complex<double> coeff1 = 1.0 / (4.0*G1),
		coeff2 = -(lambda + G1) / (G1*(lambda + 2.0*G1));
	Propagator_Fourier = TensorField(d, 4, k.num_points);

#pragma omp parallel for
	for (size_t k_idx = 0; k_idx < k.num_total_pixels; k_idx++) {	//Wavevectork
		std::vector<size_t> H_comp;
		double k2 = 0.0L;
		for (int i = 0; i < d; i++)
			k2 += Modulus2(k[i][k_idx]);

		for (size_t H_idx = 0; H_idx < Propagator_Fourier.total_comp; H_idx++) {	//component of the H-tensor
			std::complex<double> value = 0;
			Propagator_Fourier.GetComponent(H_idx, H_comp);
			size_t i = H_comp[0], j = H_comp[1], m = H_comp[2], n = H_comp[3];
			if (k2 > 1e-6) {
				value = (i == m) ? k[j][k_idx][0] * k[n][k_idx][0] : 0.0L;
				value += (i == n) ? k[j][k_idx][0] * k[m][k_idx][0] : 0.0L;
				value += (j == m) ? k[i][k_idx][0] * k[n][k_idx][0] : 0.0L;
				value += (j == n) ? k[i][k_idx][0] * k[m][k_idx][0] : 0.0L;
				value *= (coeff1 / k2);

				value += coeff2 * (k[i][k_idx][0] * k[j][k_idx][0] * k[m][k_idx][0] * k[n][k_idx][0]) / (k2*k2);
			}
			Propagator_Fourier[H_idx][k_idx][0] = value.real();
			Propagator_Fourier[H_idx][k_idx][1] = value.imag();
		}
	}

	Propagator_Fourier *= -1.0;

	update_propagator = true;
	if(FFH_verbosity == 2)
		std::cout << "The Fourier transform of G tensor (rank 4) is computed\n";

}

void Homo_ElasticProblem_static::Iterate_strong_contrast(size_t steps)
{
	const TensorField & H_Fourier = Propagator_Fourier;
	TensorField & f_Fourier = Response_Fourier;
	TensorField & f_direct = Response_direct;

	for (size_t i = 0; i < steps; i++) {
		/*1. Compute FFT of P(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Fourier transform of P (" << m << " th stage)\n";
		for (size_t j = 0; j<P_direct.total_comp; j++)
			fftw_execute(FourierTransform_P[j]);

		m++;
		if (FFH_verbosity > 0) {
			std::cout << "------------------" << m << " th stage ---------------\n";
			std::cout << "\t\t Compute F field \n";
		}
		/* f(k) <- E_0(k) + H(k).P(k)*/
		/*3. Compute \tilde{E}(k) = \tilde{G}(k). \tilde{P}(k) if k \neq 0
		\tilde{E}(k) = \tilde{E_0} */
		InnerProduct(f_Fourier, H_Fourier, P_Fourier, 2);
		/* Add \tilde{E_0} */
		for (size_t i = 0; i < applied_strain.total_comp; i++) {
			Response_Fourier[i][0][0] = applied_strain[i].real() * Response_Fourier.num_total_pixels;
			Response_Fourier[i][0][1] = applied_strain[i].imag() * Response_Fourier.num_total_pixels;
		}

		/*4. Compute the inverse FFT of \tilde{E}(k). */
		Response_direct_prev = Response_direct;
		if (FFH_verbosity > 0)
			std::cout << "\t\t Inverse Fourier transform of F field\n";
		for (size_t j = 0; j<f_direct.total_comp; j++)
			fftw_execute(InverseFourierTransform_Response[j]);
		f_direct *= (1.0 / P_Fourier.num_total_pixels);	//This step is required because fftw3 do not consider the 1/N factor.

																/*5. Compute P(x) = L(x)* F(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Compute P field\n";
#pragma omp parallel for 
		for (size_t i = 0; i < Geometry.prt_indicator.size(); i++) {
			size_t x = Geometry.prt_indicator[i];
			for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
				std::complex<double> P = InnerProduct(x, comp, this->L, Response_direct, 2);
				P *= Geometry.WeightedPhaseIndicator(i);
				P_direct[comp][x][0] = P.real();	P_direct[comp][x][1] = P.imag();
			}
		}
	}
}

void Homo_ElasticProblem_static::ComputeH_Fourier()
{
	ComputeG_Fourier();
	
	/* Exclude the contribution from polarized phase at the position. */
	for (size_t x = 0; x < Propagator_Fourier.num_total_pixels; x++) {
		for (size_t i = 0; i < Propagator_Fourier.total_comp; i++) {
			Propagator_Fourier[i][x][0] += D[i].real();
			Propagator_Fourier[i][x][1] += D[i].imag();
		}
	}

	update_propagator = true;
	if (FFH_verbosity == 2)
		std::cout << "The Fourier transform of H tensor (rank 4) is computed\n";
}

DynamicElasticProblem::DynamicElasticProblem()
{
	name_div_free_field = "stress"; name_field = "strain"; rank_of_fields = 2;
	rho = 1.0L;
}

DynamicElasticProblem::DynamicElasticProblem(const TwoPhaseDispersion & src, const std::complex<double>& K1, const std::complex<double>& G1, const std::complex<double>& K2, const std::complex<double>& G2)
	:DynamicElasticProblem()
{
	this->K1 = K1;	this->G1 = G1;
	this->K2 = K2;	this->G2 = G2;
	int d = src.packing.GetDimension();
	ConstantTensor C1 = (double)d*K1*Lambda_h(d) + 2.0*G1*Lambda_s(d) ;//stiffness tensor for phase 1
	ConstantTensor C2 = (double)d*K2*Lambda_h(d) + 2.0*G2*Lambda_s(d) ;//stiffness tensor for phase 2
	this->DefineProblem(src, C1, C2);
	this->DefineProperties(K1, G1, K2, G2);
}

void DynamicElasticProblem::DefineProblem(const TwoPhaseDispersion & src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase)
{
	FFH_base::DefineProblem(src, matrix_phase, inclusion_phase);
	DefineProperties(K1, G1, K2, G2);
	/* define L tensor */
	if (use_strong_contrast)
		this->ComputeL_strong();
	else
		this->ComputeL_weak();
}

void DynamicElasticProblem::DefineProblem(const TwoPhaseDispersion & src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector kL)
{
	this->DefineProblem(src, matrix_phase, inclusion_phase);
	this->SetWavevector(kL);
}

void DynamicElasticProblem::DefineProperties(const std::complex<double>& K1, const std::complex<double>& G1, const std::complex<double>& K2, const std::complex<double>& G2)
{
	double d = this->d;

	if (this->optimized_reference) {
		if (use_strong_contrast) {
			if(FFH_verbosity==2)
				std::cout << "the optimized reference is not defined in the strong-contrast expansion yet. We use those of the weak-contrast expansions\n";
			K0 = 0.5*(K1 + K2);	G0 = 0.5*(G1 + G2);
		}
		else {
			K0 = 0.5*(K1 + K2);	G0 = 0.5*(G1 + G2);
		}
	}
	else {
		this->K0 = K1;	this->G0 = G1;
	}

	cL = sqrt((K0 + 2.0*(1.0-1.0/d)*G0)/rho ).real();
	cT = sqrt((G0)/rho ).real();

	kappa1 = (K1 - K0) / (K1 + 2.0*(1.0 - 1.0 / d)*G0);
	mu1 = (G1 - G0) / (G1 + (d*K0 / 2.0 + (d + 1.0)*(d - 2.0)*G0 / d )*G0 / (K0 + 2.0*G0));

	kappa2 = (K2 - K0) / (K2 + 2.0*(1.0 - 1.0 / d)*G0);
	mu2 = (G2 - G0) / (G2 + (d*K0 / 2.0 + (d + 1.0)*(d - 2.0)*G0 / d)*G0 / (K0 + 2.0*G0));

	this->D = 1.0 / (d*cL*cL)*Lambda_h(this->d) + (2.0 / (d*cL*cL) + 1.0 / (cT*cT)) / (d + 2)*Lambda_s(d);
}

void DynamicElasticProblem::SetWavevector(const GeometryVector & kL)
{
	if (kL.Dimension == d) {
		double cL_phase1 = sqrt(abs(K1 + 2.0*(1.0-1.0/(double)d)*G1));

		this->kL = cL_phase1/cL * sqrt(kL.Modulus2());
		this->kT = (cL / cT) * this->kL;
		GeometryVector new_k = 1.0 / sqrt(kL.Modulus2())*kL;
		if (this->k.Dimension != d || new_k != this->k) {
			this->update_propagator = false;
		}
		
		this->k = new_k;
		this->kk = ConstantTensor(d, 2);
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				kk[j + d * i] = new_k.x[i] * new_k.x[j];
			}
		}

	}
	else {
		std::cerr << "kL should have the same dimensions as the system.\n";
	}
}

bool DynamicElasticProblem::CheckDivergenceFree(double target_tolerance)
{
	TensorField stress, projection, strain;
	double omega2 = pow(cL*abs(kL),2),/* the square of angular frequencies. */
		error_in_div = 0.0L;


	/* compute the Fourier transform of the stress field*/
	this->ComputeDivergenceFreeField_Fourier(stress);

	this->ComputeField_Fourier(strain);

	double mean = 0.0L;
	for (size_t i = 0; i < strain.total_comp; i++)
		mean += Modulus2(strain[i][0]);
	mean = sqrt(mean);

	strain *= (omega2 * rho);
	
	
	InnerProduct(projection, this->qq, stress, 1);
	/* Symmetrize projection. */
#pragma omp parallel for
	for (size_t q_indx = 0; q_indx < projection.num_total_pixels; q_indx++) {
		for (int i = 0; i < d; i++) {
			for (int j = i; j < d; j++) {
				int idx1 = i + d * j, idx2 = j + d * i;
				std::complex<double> val(projection[idx1][q_indx][0]+projection[idx2][q_indx][0],
					projection[idx1][q_indx][1]+projection[idx2][q_indx][1]);
				val *= 0.5;
//#pragma omp atomic
				projection[idx1][q_indx][0] = val.real();
//#pragma omp atomic
				projection[idx1][q_indx][1] = val.imag();
//#pragma omp atomic
				projection[idx2][q_indx][0] = val.real();
//#pragma omp atomic
				projection[idx2][q_indx][1] = val.imag();
			}
		}
	}
//	strain -= projection;

	/* compute the squared sum of the divergence in the Fourier space. */
#pragma omp parallel for schedule(guided)
	for (size_t x = 0; x < strain.num_total_pixels; x++) {
		double diff = 0.0L, value = 0.0L;

		for (size_t i = 0; i < strain.total_comp; i++) {
			double mod2 = 0.0L;

			for (int part = 0; part < 2; part++) {
				value = strain[i][x][part];
				diff = value - projection[i][x][part];
				/*check the significant digits. */
				if ((log10(abs(value)) - log10(abs(diff))) > 15) {
					/* 0 in numerical precision */
					diff = 0.0L;
				}
			
				mod2 += (diff*diff);
			}

#pragma omp atomic
			error_in_div += mod2;
		}
	}
	

	error_in_div = sqrt(error_in_div / (double)strain.num_total_pixels) / mean;
	Tol_div = error_in_div;
	if (FFH_verbosity == 2) {
		std::cout << "\t\t Error in the equation of motion = " << error_in_div << "\n";
	}

	{
		// compute the divergence of the stress
		TensorField div;
		InnerProduct(div, this->ks, stress, 1);
		double error = 0.0L;
		for (size_t x = 0; x < div.num_total_pixels; x++) {
			double abs_div = 0.0L;
			for (int i = 0; i < d; i++) {
				abs_div += Modulus2(div[i][x]);
			}
			error += abs_div;
		}
		double sigma0 = 0.0L;
		for (size_t i = 0; i < stress.total_comp; i++)
			sigma0 += Modulus2(stress[i][0]);

		error = sqrt(error / (double)stress.num_total_pixels) / sqrt(sigma0);
		//Tol_div = error;
		if (FFH_verbosity == 2)
			std::cout << "Error in the divergence of stress tensor = " << error << "\n";

	}

	return (error_in_div < target_tolerance);

}

void DynamicElasticProblem::Prepare(const ConstantTensor & applied)
{
	double d = this->d;
	{/* decompose the applied strain into the longitudinal and transverse contributions. */
		std::complex<double> epsilon_L = InnerProduct(this->kk, applied, 2)[0],
			epsilon_T = InnerProduct(Identity(this->d) - this->kk, applied, 2)[0];
		strain_L = epsilon_L * this->kk;
		//strain_T = epsilon_T * (1.0 / d * Identity(this->d) - this->kk);


		this->ks = Geometry.GetKs();
		
		if (abs(epsilon_L) < 1e-10) {
			std::cout << "Incident waves are transverse\n";
			this->k0 = kT.real()*this->k;
		}
		else if (abs(epsilon_T) < 1e-10) {
			std::cout << "Incident waves are longitudinal\n";
			this->k0 = kL.real()*this->k;
		}
		else {
			std::cerr << "The applied field has both longitudinal and transverse components.\n"
				<< "Such cases are not implemented yet!!\n";
			return;
		}
		/* shift ks by k0. */
#pragma omp parallel for
		for (size_t k_indx = 0; k_indx < ks.num_total_pixels; k_indx++) {
			for (int i = 0; i < d; i++) {
#pragma omp atomic
				ks[i][k_indx][0] += k0.x[i];
			}
		}
		if (FFH_verbosity == 2)
			ks.Print_shifted_origin("test_wavevector");

		this->qq = TensorField(d, 2, Geometry.GetResolutions());
#pragma omp parallel for 
		for (size_t k_indx = 0; k_indx < ks.num_total_pixels; k_indx++) {
			for (int i = 0; i < d; i++) {
				for (int j = 0; j < d; j++) {
					qq[i+d*j][k_indx][0] = ks[i][k_indx][0]*ks[j][k_indx][0];
				}
			}
		}
		if (FFH_verbosity == 2)
			this->qq.Print_shifted_origin("test_qq");
		this->update_propagator = false;
	}
	//TODO:: generalize the code to the cases where the longitudinal and transverse components are mixed.

	FFH_base::Prepare(applied);

}

void DynamicElasticProblem::ComputeEffectiveModuli(std::vector<std::complex<double>>& result) const
{
	double d = this->d;
	result.resize(0);
	if (use_strong_contrast) {
		ConstantTensor f_av = MeanField(Response_direct),
			P_av = MeanField(P_direct);

		std::complex<double> Ke;
		{
			std::complex<double> A = Trace(P_av) / Trace(f_av);
			std::complex<double> kappa_e = A / (d*K1 + 2.0*(d - 1.0)*G1);
			Ke = (d*K1 + 2.0*G1*(d - 1.0)*kappa_e) / (d*(1.0 - kappa_e));
		}
		/* deviatoric tensors */
		f_av = f_av - Trace(f_av)*HydrostaticStrain(this->d);
		P_av = P_av - Trace(P_av)*HydrostaticStrain(this->d);
		std::complex<double> Ge;
		{
			std::complex<double> B(0.0L, 0.0L);
			double count = 0;
			for (size_t i = 0; i < f_av.total_comp; i++) {
				if (abs(f_av[i]) > Tol_Field && abs(P_av[i]) > Tol_Field) {
					B += P_av[i] / f_av[i];
					count++;
				}
			}
			if (count > 0)
				B /= count;

			std::complex<double> mu_e = B / (d*K1 + 2.0*(d - 1.0)*G1);
			mu_e /= ((d + 2.0)*G1 / d / (K1 + 2.0*G1));


			Ge = -G1 * ((2.0* d*(2.0* G1 + K1) + 2.0* (-2.0 + d) *(1.0 + d)* G1* mu_e + d * d* K1* mu_e)) / (2.0* d*(2.0* G1 + K1)* (-1.0 + mu_e));
			//Ge = G1 * (2.0*d*(2.0*G1 + K1) + (2.0*(d-2.0)*(d+1.0)*G1 + d*d*K1)*mu_e) / (2.0*d*(2.0*G1+K1)*(1.0-mu_e));
		}
		//TODO: implement for anisotropic case.

		result.push_back(Ke);
		result.push_back(Ge);
	}
	else {
		/* Response = strain! */
		ConstantTensor stress_av(this->d, 2), strain_av(this->d, 2);
		TensorField stress(P_direct);

		ConstantTensor C1 = d * K1*Lambda_h(this->d) + 2.0*G1*Lambda_s(this->d);
		//TODO:: simplify it...
		for (size_t x = 0; x < stress.num_total_pixels; x++) {
			for (size_t i = 0; i < stress.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, C1, Response_direct, 2);
				stress[i][x][0] += sigma.real();
				stress[i][x][1] += sigma.imag();
			}
		}

		/* compute the means of stress and strain. */
#pragma omp parallel for 
		for (size_t idx = 0; idx < stress.total_comp; idx++) {
			for (size_t x_idx = 0; x_idx < stress.num_total_pixels; x_idx++) {
#pragma omp critical
				{//TODO:: check whether this is thread-safe.
					stress_av[idx] += std::complex<double>(stress[idx][x_idx][0], stress[idx][x_idx][1]);
					strain_av[idx] += std::complex<double>(Response_direct[idx][x_idx][0], Response_direct[idx][x_idx][1]);
				}
			}
		}

		stress_av = 1.0 / (double)stress.num_total_pixels * stress_av;
		strain_av = 1.0 / (double)stress.num_total_pixels * strain_av;

		std::complex<double> Ke(0.0L, 0.0L), Ge(0.0L, 0.0L),
			stress_h(0.0L, 0.0L), strain_h(0.0L, 0.0L);

		for (int i = 0; i < d; i++) {
			/*stress_h += stress_av[i];
			strain_h += strain_av[i];*/
			stress_h += stress_av[i*(this->d + 1)];
			strain_h += strain_av[i*(this->d + 1)];
		}
		Ke = 1.0 / d * stress_h / strain_h;

		stress_av = stress_av - stress_h / d * HydrostaticStrain(this->d);
		strain_av = strain_av - strain_h / d * HydrostaticStrain(this->d);

		double count = 0;
		for (size_t i = 0; i < strain_av.total_comp; i++) {
			if (abs(stress_av[i]) > 1e-10) {
				Ge += (stress_av[i] / strain_av[i]);
				count++;
			}
		}
		Ge *= (0.5 / count);
		//TODO::implement anisotropic effective tensor...
		result.push_back(Ke);
		result.push_back(Ge);
	}

}

void DynamicElasticProblem::ComputeDivergenceFreeField_Fourier(TensorField & out) const
{
	double d = this->d;
	out = TensorField(P_Fourier);
	ConstantTensor C0 = d * K0*Lambda_h(this->d) + 2.0*G0*Lambda_s(this->d);

	if (use_strong_contrast) {
		/* stress = C1:f(x) + (I - C1:D): P(x)*/
		ConstantTensor I_C1D = Lambda_h(this->d) + Lambda_s(this->d) - InnerProduct(C0, D, rank_of_fields);

#pragma omp parallel for
		for (size_t x = 0; x < out.num_total_pixels; x++) {
			for (size_t i = 0; i < out.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, C0, Response_Fourier, rank_of_fields) + InnerProduct(x, i, I_C1D, P_Fourier, rank_of_fields);
				out[i][x][0] = sigma.real();
				out[i][x][1] = sigma.imag();
			}
		}
	}
	else {
		/* sigma(x) = C1:epsilon(x) + P(x) */
#pragma omp parallel for
		for (size_t x = 0; x < out.num_total_pixels; x++) {
			for (size_t i = 0; i < out.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, C0, Response_Fourier, rank_of_fields);
#pragma omp atomic
				out[i][x][0] += sigma.real();
#pragma omp atomic
				out[i][x][1] += sigma.imag();
			}
		}
	}
}

void DynamicElasticProblem::ComputeDivergenceFreeField_Direct(TensorField & out) const
{
	double d = this->d;
	out = TensorField(P_direct);
	ConstantTensor C0 = d * K0*Lambda_h(this->d) + 2.0*G0*Lambda_s(this->d);

	if (use_strong_contrast) {
		/* stress = C1:f(x) + (I - C1:D): P(x)*/
		ConstantTensor I_C1D = Lambda_h(this->d) + Lambda_s(this->d) - InnerProduct(C0, D, rank_of_fields);

#pragma omp parallel for
		for (size_t x = 0; x < out.num_total_pixels; x++) {
			for (size_t i = 0; i < out.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, C0, Response_direct, 2) + InnerProduct(x, i, I_C1D, P_direct, 2);
				out[i][x][0] = sigma.real();
				out[i][x][1] = sigma.imag();
			}
		}
	}
	else {
		/* sigma(x) = C1:epsilon(x) + P(x) */
#pragma omp parallel for
		for (size_t x = 0; x < out.num_total_pixels; x++) {
			for (size_t i = 0; i < out.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, C0, Response_direct, 2);
#pragma omp atomic
				out[i][x][0] += sigma.real();
#pragma omp atomic
				out[i][x][1] += sigma.imag();
			}
		}
	}

	ShiftFieldByk0(out);
}

void DynamicElasticProblem::ComputeField_Fourier(TensorField & out) const
{
	out = TensorField(Response_Fourier);

	if (use_strong_contrast) {
		/*	Response = f, where f = strain + D : Polarization and out = strain. */
#pragma omp parallel for
		for (size_t x = 0; x < out.num_total_pixels; x++) {
			std::complex<double> sigma;
			for (size_t i = 0; i < out.total_comp; i++) {
				sigma = InnerProduct(x, i, D, P_Fourier, rank_of_fields);
#pragma omp atomic
				out[i][x][0] -= sigma.real();
#pragma omp atomic
				out[i][x][1] -= sigma.imag();
			}
		}
	}
	else {
		/* Response = strain. */
		//out = TensorField(Response_Fourier);
	}
}

void DynamicElasticProblem::ComputeField_Direct(TensorField & out) const
{
	out = TensorField(Response_direct);

	if (use_strong_contrast) {
		/* f = strain + D:P, where out = strain. */
#pragma omp parallel for
		for (size_t x = 0; x < out.num_total_pixels; x++) {
			for (size_t i = 0; i < out.total_comp; i++) {
				std::complex<double> sigma;
				sigma = InnerProduct(x, i, D, P_direct, 2);
#pragma omp atomic
				out[i][x][0] -= sigma.real();
#pragma omp atomic
				out[i][x][1] -= sigma.imag();
			}
		}
	}
	else {
		/* Response = strain. */
		//out = TensorField(Response_Fourier);
	}
	ShiftFieldByk0(out);
}

void DynamicElasticProblem::ComputePropagator_weak()
{
	Propagator_Fourier = TensorField(d, 2 * rank_of_fields, Geometry.GetResolutions());
	/* Note that ks is already shifted by k0. */
#pragma omp parallel for
	for (size_t q_indx = 0; q_indx < ks.num_total_pixels; q_indx++) {
		GeometryVector q(static_cast<int>(d));
		for (int i = 0; i < d; i++)
			q.x[i] = ks[i][q_indx][0];

		ConstantTensor val = GreenFunction_Fourier(q);
		

		for (size_t i = 0; i < val.total_comp; i++) {
			Propagator_Fourier[i][q_indx][0] = val[i].real();
			Propagator_Fourier[i][q_indx][1] = val[i].imag();
		}
	}

}

void DynamicElasticProblem::ComputeL_weak()
{
	//TODO::implement!!
	//L = (double)d * (K2 - K1) *Lambda_h(d) + 2.0*(G2 - G1)*Lambda_s(d);
	L1 = (double)d * (K1 - K0) *Lambda_h(d) + 2.0*(G1 - G0)*Lambda_s(d);
	L2 = (double)d * (K2 - K0) *Lambda_h(d) + 2.0*(G2 - G0)*Lambda_s(d);
}

void DynamicElasticProblem::ComputePropagator_strong()
{
	this->ComputePropagator_weak();
#pragma omp parallel for
	for (size_t x = 0; x < Propagator_Fourier.num_total_pixels; x++) {
		for (size_t i = 0; i < Propagator_Fourier.total_comp; i++) {
#pragma omp atomic
			Propagator_Fourier[i][x][0] += this->D[i].real();
#pragma omp atomic
			Propagator_Fourier[i][x][1] += this->D[i].imag();
		}
	}
}

void DynamicElasticProblem::ComputeL_strong()
{
	double d = this->d;

	L1 = kappa1 * Lambda_h(d) + (d + 2.0)*cT*cT / (d*cL*cL + 2.0*cT*cT) *mu1*Lambda_s(d);
	L1 = d * cL*cL*L1;

	L2 = kappa2 * Lambda_h(d) + (d + 2.0)*cT*cT / (d*cL*cL + 2.0*cT*cT) *mu2*Lambda_s(d);
	L2 = d * cL*cL*L2;
}


ConstantTensor DynamicElasticProblem::GreenFunction_Fourier(const GeometryVector & q)
{
	int indx = 0;
	double tol = 1e-15;
	double qq = q.Modulus2();
	ConstantTensor T3(d, 2 * rank_of_fields);
	ConstantTensor G(T3);
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			for (int k = 0; k < d; k++) {
				for (int l = 0; l < d; l++) {
					indx = i + d * (j + d * (k + d * l));
				// indx = l + d * (k + d * (j + d * i));
					T3[indx] = q.x[i] * q.x[j] * q.x[k] * q.x[l];
				}
			}
		}
	}



	std::complex<double> I(0.0L, 1.0L);
	
	std::complex<double> contrib_long(d, 2 * rank_of_fields);
	std::complex<double> contrib_trans(d, 2 * rank_of_fields);
	if (abs(qq - kL * kL) < tol) {
		//princial value
		contrib_long = I * pi*(double)Response_Fourier.num_total_pixels;
	}
	else {//
		contrib_long = 1.0 / (qq - kL * kL);
	}

	if (abs(qq - kT * kT) < tol) {
		contrib_trans = I * pi*(double)Response_Fourier.num_total_pixels;
	}
	else {
		contrib_trans = 1.0 / (qq - kT * kT);
	}

	G = -1.0 / (cT*cT)*contrib_trans * T2(q) + (-1.0/(cL*cL) + 1.0/(cT*cT))*contrib_long*contrib_trans*T3;
	
	return G;
}

void DynamicElasticProblem::ShiftFieldByk0(TensorField & out) const
{
	std::vector<GeometryVector> normalized_basis;
	std::vector<size_t> pixels = Geometry.GetResolutions();
	for (int i = 0; i < d; i++)
		normalized_basis.emplace_back((1.0 / (double)pixels[i]) *Geometry.packing.GetBasisVector(i));

	/* multiply the phase factor. */
#pragma omp parallel for
	for (size_t x = 0; x < out.num_total_pixels; x++) {
		std::vector<size_t> n_x;
		Geometry.Index2Components(x, n_x);

		GeometryVector pos(static_cast<int>(d));
		for (int i = 0; i < d; i++)
			pos = pos + (double)n_x[i] * normalized_basis[i];

		std::complex<double> phase_factor = std::exp(std::complex<double>(0.0L, k0.Dot(pos)));

		for (int i = 0; i < d; i++) {
			std::complex<double> value = std::complex<double>(out[i][x][0], out[i][x][1]);
			value *= phase_factor;
			out[i][x][0] = value.real();
			out[i][x][1] = value.imag();
		}

	}
}

ConstantTensor T2(const GeometryVector & q) {
	int d = q.Dimension, indx;
	ConstantTensor out(d, 4);
	
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			for (int k = 0; k < d; k++) {
				for (int l = 0; l < d; l++) {
					indx = i + d * (j + d * (k + d * l));
					if (i == k)
						out[indx] += (q.x[j] * q.x[l]);
					if (j == k)
						out[indx] += (q.x[i] * q.x[l]);
					if (i == l)
						out[indx] += (q.x[j] * q.x[k]);
					if (j == l)
						out[indx] += (q.x[i] * q.x[k]);
				}
			}
		}
	}

	return	0.25*out;
}