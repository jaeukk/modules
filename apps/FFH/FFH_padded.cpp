#include "FFH_padded.h"

DynamicDielectricProblem_padded::DynamicDielectricProblem_padded() : StaticDielectricProblem() {}

DynamicDielectricProblem_padded::DynamicDielectricProblem_padded(const TwoPhaseDispersion & src, const std::complex<double>& epsilon1, const std::complex<double>& epsilon2) :
	StaticDielectricProblem()
{
	this->epsilon1 = epsilon1;	this->epsilon2 = epsilon2;
	this->DefineProblem(src, epsilon1*Identity(d), epsilon2*Identity(d), GeometryVector(static_cast<int>(src.packing.GetDimension())));
	this->D = 1.0 / ((double)src.packing.GetDimension() * epsilon1) * Identity(src.packing.GetDimension());
}

void DynamicDielectricProblem_padded::DefineProblem(const TwoPhaseDispersion & src, const ConstantTensor & matrix_phase, const ConstantTensor & inclusion_phase, const GeometryVector & k1)
{
	bool SameSimulationBox = this->DefineDispersion(src);

	if (!SameSimulationBox) {
		update_propagator = false;	update_plan = false;
		//Propagator_Fourier = TensorField(d, 2 * rank_of_fields, Geometry.GetResolutions());
		//ks = src.GetKs();
	}

	this->SetWavevector(k1);

	std::vector<size_t> padded_size(Geometry.GetResolutions());
	for (int i = 0; i < d; i++)
		padded_size[i] *= 2;

	Response_direct = TensorField(d, rank_of_fields, padded_size);
	Response_Fourier = TensorField(d, rank_of_fields, padded_size);
	P_direct = TensorField(d, rank_of_fields, padded_size);
	P_Fourier = TensorField(d, rank_of_fields, padded_size);
	Propagator_Fourier = TensorField(d, 2*rank_of_fields, padded_size);
	Propagator_direct = TensorField(d, 2*rank_of_fields, padded_size);

}

void DynamicDielectricProblem_padded::DefineProperties(const std::complex<double>& epsilon1, const std::complex<double>& epsilon2, const GeometryVector & k1)
{
	StaticDielectricProblem::DefineProperties(epsilon1, epsilon2);
	SetWavevector(k1);
}

void DynamicDielectricProblem_padded::SetWavevector(const GeometryVector & k1)
{
	if (k1.Dimension == d) {
		if (this->k1.Dimension != d || this->k1 != k1) {
			this->update_propagator = false;
		}
		GeometryVector k_prev(this->k1);
		this->k1 = k1;
		/* shift ks by k1. */
		ks = Geometry.GetKs();
#pragma omp parallel for
		for (int k_indx = 0; k_indx < ks.num_total_pixels; k_indx++) {
			for (int i = 0; i < d; i++) {
#pragma omp atomic
				ks[i][k_indx][0] += k1.x[i];
			}
		}
		/* We already have the data. */
		this->update_propagator = false;
	}
	else {
		std::cerr << "k1 should have the same dimensions as the system.\n";
	}
}

void DynamicDielectricProblem_padded::Prepare(const ConstantTensor & applied_electric_field)
{
	std::complex<double> error = 0;
	for (int i = 0; i < d; i++)
		error += applied_electric_field[i] * k1.x[i];
	
	if (abs(error) < 1e-20) {
		FFH_base::Prepare(applied_electric_field);
	}
	else {
		std::cerr << "DynamicDielectricProblem_padded::applied field must be perpendicular to k1.\n";
	}
	
	
	{
		P_direct.Print("test_P_direct");
	}
}

void DynamicDielectricProblem_padded::ComputeDivergenceFreeField_Fourier(TensorField & out) const
{
	TensorField aux;
	StaticDielectricProblem::ComputeDivergenceFreeField_Fourier(aux);
	CropTensorField_Fourier(out, aux);
}

void DynamicDielectricProblem_padded::ComputeField_Fourier(TensorField & out) const
{
	TensorField aux;
	StaticDielectricProblem::ComputeField_Fourier(aux);
	CropTensorField_Fourier(out, aux);
}

void DynamicDielectricProblem_padded::ComputePropagator_weak()
{
	this->ComputePropagator_strong();
#pragma omp parallel for
	for (int x = 0; x < Propagator_Fourier.num_total_pixels; x++) {
		for (int i = 0; i < Propagator_Fourier.total_comp; i++) {
#pragma omp atomic
			Propagator_Fourier[i][x][0] -= this->D[i].real();
#pragma omp atomic
			Propagator_Fourier[i][x][1] -= this->D[i].imag();
		}
	}

}

void DynamicDielectricProblem_padded::ComputePropagator_strong()
{
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
	for (int x_idx = 0; x_idx < Propagator_direct.num_total_pixels; x_idx++) {
		std::vector<size_t> x_comp;		Propagator_direct.GetPositionComponent(x_idx, x_comp);
		GeometryVector pos(static_cast<int>(d));

		for (int i = 0; i < d; i++) {
			if (Propagator_direct.num_points[i] < 2 * x_comp[i])
				pos.x[i] = dx[i] * ((double)x_comp[i] - (double)Propagator_direct.num_points[i]);
			else
				pos.x[i] = dx[i] * (double)(x_comp[i]);
		}
		ConstantTensor H_ij = Hfunction_direct(pos, k1, d, min_dist);

		for (int i = 0; i < H_ij.total_comp; i++) {
			Propagator_direct[i][x_idx][0] = H_ij[i].real();
			Propagator_direct[i][x_idx][1] = H_ij[i].imag();
		}
	}


	{
		Propagator_direct.Print_shifted_origin("test_H_direct");
	}
	
	/* prepare for DFT */
	std::vector<fftw_plan> fft_H;
	fft_H.resize(Propagator_direct.total_comp);
	if (d == 2) {
		for (int i = 0; i < Propagator_direct.total_comp; i++) {
			fft_H[i] = fftw_plan_dft_2d(Propagator_direct.num_points[0], Propagator_direct.num_points[1],
				Propagator_direct[i], Propagator_Fourier[i], -1, FFTW_ESTIMATE);
		}
	}
	else if (d == 3) {
		for (int i = 0; i < Propagator_direct.total_comp; i++) {
			fft_H[i] = fftw_plan_dft_3d(Propagator_direct.num_points[0], Propagator_direct.num_points[1], Propagator_direct.num_points[2],
				Propagator_direct[i], Propagator_Fourier[i], -1, FFTW_ESTIMATE);
		}
	}

	/* Compute its Fourier transform. */
	for (int i = 0; i < fft_H.size(); i++)
		fftw_execute(fft_H[i]);

	Propagator_Fourier *= dv;

	{
		TensorField H_crop;
		CropTensorField_Fourier(H_crop, Propagator_Fourier);
		TensorField div;
		InnerProduct(div, ks, H_crop, 1);
		div.Print_shifted_origin("test_div");
	}
}

inline ConstantTensor DynamicDielectricProblem_padded::Hfunction_direct(const GeometryVector & position, const GeometryVector & k1, int d, double dx)
{
	ConstantTensor H_ij(d, 2);
	double k = sqrt(k1.Modulus2());
	double r = sqrt(position.Modulus2());
	double kr = k * r;

	std::complex<double> I(0.0L, 1.0L);
	if (r < 0.5*dx) {
		std::complex<double> G1;
		double dv = 1.0, dr;
		for (int i = 0; i < d; i++)
			dv *= dx;

		if (d == 2) {
			dr = dx / sqrt(pi);
			G1 = (-0.5 + I * 0.25*k*dr*pi *boost::math::cyl_hankel_1(1, k*dr)) / this->epsilon1 / dv;
		}
		else if (d == 3) {
			dr = dx / cbrt(4.0*pi / 3.0);
			G1 = 2.0 / 3.0*(-1.0 + exp(I*k*dr)*(1.0 - I * k*dr)) / (this->epsilon1) / dv;
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
			G1 = 0.25*k*k*I*(boost::math::cyl_hankel_1(0, kr) - boost::math::cyl_hankel_1(1, kr) / (kr));
			G2 = 0.25*k*k*I*boost::math::cyl_hankel_1(2, kr);
		}
		else if (d == 3) {
			G1 = (-1.0 + I * kr + kr * kr)*exp(I*kr) / (4.0*pi*r*r*r);
			G2 = (3.0 - 3.0*I*kr - kr * kr)*exp(I*kr) / (4.0*pi*r*r*r);
		}

		H_ij = 1.0 / this->epsilon1*(G1 * Identity(d) + G2 * rr);
	}
	return exp(-I * k1.Dot(position))*H_ij;
}

void DynamicDielectricProblem_padded::ShiftFieldByk1(TensorField & out) const
{
	std::vector<GeometryVector> normalized_basis;
	std::vector<size_t> pixels = Geometry.GetResolutions();
	for (int i = 0; i < d; i++)
		normalized_basis.emplace_back((1.0 / (double)pixels[i]) *Geometry.packing.GetBasisVector(i));

	/* multiply the phase factor. */
#pragma omp parallel for
	for (int x = 0; x < out.num_total_pixels; x++) {
		std::vector<size_t> n_x;
		Geometry.Index2Components(x, n_x);

		GeometryVector pos(static_cast<int>(d));
		for (int i = 0; i < d; i++)
			pos = pos + (double)n_x[i] * normalized_basis[i];

		std::complex<double> phase_factor = std::exp(std::complex<double>(0.0L, k1.Dot(pos)));

		for (int i = 0; i < d; i++) {
			std::complex<double> value = std::complex<double>(out[i][x][0], out[i][x][0]);
			value *= phase_factor;
			out[i][x][0] = value.real();
			out[i][x][1] = value.imag();
		}

	}
}

void DynamicDielectricProblem_padded::CropTensorField_direct(TensorField & out, const TensorField & in) const
{
	out = TensorField(d, in.rank, Geometry.GetResolutions());
	for (int x_idx = 0; x_idx < out.num_total_pixels; x_idx++) {
		std::vector<size_t> out_position_comps; out.GetPositionComponent(x_idx, out_position_comps);
		size_t in_position_idx = in.at(out_position_comps);
		for (int field_comp = 0; field_comp < in.total_comp; field_comp++) {
			out[field_comp][x_idx][0] = in[field_comp][in_position_idx][0];
			out[field_comp][x_idx][1] = in[field_comp][in_position_idx][1];
		}
	}
}

void DynamicDielectricProblem_padded::CropTensorField_Fourier(TensorField & out, const TensorField & in) const
{
	out = TensorField(d, in.rank, Geometry.GetResolutions());
	for (int x_idx = 0; x_idx < out.num_total_pixels; x_idx++) {
		std::vector<size_t> out_position_comps;  out.GetPositionComponent(x_idx, out_position_comps);
		for (int i = 0; i < d; i++)
			out_position_comps[i] *= 2;
		size_t in_position_idx = in.at(out_position_comps);
		for (int field_comp = 0; field_comp < in.total_comp; field_comp++) {
			out[field_comp][x_idx][0] = in[field_comp][in_position_idx][0];
			out[field_comp][x_idx][1] = in[field_comp][in_position_idx][1];
		}
	}
}
