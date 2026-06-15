#pragma once

/* temporary ... */

/* Class for homogenizing the static dielectric constant of two-phase sphere packings. */
struct Homo_DielectricConstant_static {
	DimensionType d;
	TwoPhaseDispersion Geometry;
	double Tol = 1e-10;
	double epsilon1, epsilon2, delta_epsilon;
	bool update_G, update_plan;
	size_t direction = 0;
	TensorField G_Fourier; /*< Static Green function in Fourier space*/
	TensorField P_Fourier, E_Fourier;	/*< Polarization and electric field in Fourier space representations */
	TensorField P_direct, E_direct;		/*< Polarization and electric field in direct space representations */
	std::vector<fftw_plan> FourierTransform_P, InverseFourierTransform_E;
	size_t m;	/*< the number of stages */

	Homo_DielectricConstant_static() {
		d = 0; epsilon1 = 0; epsilon2 = 0; delta_epsilon = 0;
		update_G = false; update_plan = false; m = 0;
	};
	Homo_DielectricConstant_static(const TwoPhaseDispersion& src, double epsilon1, double epsilon2);
	~Homo_DielectricConstant_static() {
		for (int i = 0; i < FourierTransform_P.size(); i++) {
			fftw_destroy_plan(FourierTransform_P[i]);
			fftw_destroy_plan(InverseFourierTransform_E[i]);
		}
		fftw_cleanup_threads();
		G_Fourier.~TensorField();
		P_Fourier.~TensorField();	E_Fourier.~TensorField();
		P_direct.~TensorField();	E_direct.~TensorField();
	};

	void RedefineProblem(const TwoPhaseDispersion& src, double epsilon1, double epsilon2);
	/* Initialize the problem.
	 * @param[in] direction	The direction of external static field. */
	void Prepare(size_t direction = 0);
	void Iterate(size_t steps = 1);
	void Iterate_test(size_t steps = 1);
	void Iterate_test2(size_t steps = 1);

	void Run();

	bool ConvergenceTest(double tolerance = 1e-14);
	std::complex<double> ComputeEffectiveDielectric();
	void Print(const std::string & prefix) const;
	void Print_Fourier(const std::string & prefix) const;
	void Print_Grid_2D(const std::string & prefix) const;
	void Print_GreenFn_2D(const std::string & prefix) const;

private:
	TensorField k, temp1, temp2;	/* temporary storage */
};

/* Class for homogenizing the static dielectric constant of two-phase sphere packings (Strong-contrast approach). */
struct Homo_DielectricConstant_static_SCA {
	DimensionType d;
	TwoPhaseDispersion Geometry;
	double Tol_F = -1, Tol_div = -1;
	double epsilon1, epsilon2, delta_epsilon;
	bool update_H, update_plan;
	double L;
	size_t direction = 0;
	TensorField H_Fourier; /*< Static Green function in Fourier space*/
	TensorField P_Fourier, F_Fourier;	/*< Polarization and electric field in Fourier space representations */
	TensorField P_direct, F_direct;		/*< Polarization and electric field in direct space representations */
	std::vector<fftw_plan> FourierTransform_P, InverseFourierTransform_F;
	size_t m;	/*< the number of stages */

	Homo_DielectricConstant_static_SCA() {
		d = 0; epsilon1 = 0; epsilon2 = 0; delta_epsilon = 0, L = 0;
		update_H = false; update_plan = false; m = 0;
	};
	Homo_DielectricConstant_static_SCA(const TwoPhaseDispersion& src, double epsilon1, double epsilon2);
	~Homo_DielectricConstant_static_SCA() {
		for (int i = 0; i < FourierTransform_P.size(); i++) {
			fftw_destroy_plan(FourierTransform_P[i]);
			fftw_destroy_plan(InverseFourierTransform_F[i]);
		}
		fftw_cleanup_threads();
		H_Fourier.~TensorField();
		P_Fourier.~TensorField();	F_Fourier.~TensorField();
		P_direct.~TensorField();	F_direct.~TensorField();
	};

	virtual void RedefineProblem(const TwoPhaseDispersion& src, double epsilon1, double epsilon2);
	/* Initialize the problem.
	* @param[in] direction	The direction of external static field. */
	virtual void Prepare(size_t direction = 0);
	virtual void Iterate(size_t steps = 1);

	virtual void Run(double tol_field = 1e-10, double tol_div = 1e-14);

	bool CheckConvergence(double tol_field);
	bool CheckDivergenceFree(double tol_div);
	std::complex<double> ComputeEffectiveDielectric();
	void Print(const std::string & prefix) const;
	void Print_Fourier(const std::string & prefix) const;
	void Print_Grid_2D(const std::string & prefix) const;
	void Print_GreenFn_2D(const std::string & prefix) const;

protected:
	TensorField k, F_direct_prev;	/* temporary storage */

	virtual void ComputeH_Fourier();
};


/* Class for homogenizing the static dielectric constant of two-phase sphere packings (Strong-contrast approach). */
struct Homo_DielectricConstant_dynamic_SCA : public Homo_DielectricConstant_static_SCA {

	Homo_DielectricConstant_dynamic_SCA() : Homo_DielectricConstant_static_SCA() {};
	Homo_DielectricConstant_dynamic_SCA(const TwoPhaseDispersion& src, double epsilon1, double epsilon2, const ConstantTensor & k1);
	~Homo_DielectricConstant_dynamic_SCA() {};

	void RedefineDynamicProblem(const TwoPhaseDispersion& src, double epsilon1, double epsilon2, const ConstantTensor & k1);
	void RedefineWavevector(const ConstantTensor & k1);
	/* Initialize the problem.
	* @param[in] direction	The direction of external static field. */
	virtual void Prepare(size_t direction = 0);
	virtual void Iterate(size_t steps = 1);
	virtual void Use_Strong_Contrast(bool use_sc);
	//	virtual void Run(double tol_field = 1e-10, double tol_div = 1e-14);

	//	bool CheckConvergence(double tol_field);
	//	bool CheckDivergenceFree(double tol_div);
	//	std::complex<double> ComputeEffectiveDielectric();
	//	void Print(const std::string & prefix) const;
	//	void Print_Fourier(const std::string & prefix) const;
	//	void Print_Grid_2D(const std::string & prefix) const;
	//	void Print_GreenFn_2D(const std::string & prefix) const;

protected:
	virtual void ComputeG_Fourier();
	virtual void ComputeH_Fourier();
	ConstantTensor k1;
	bool use_strong_contrast = false;
};


void DFT(TensorField & out, const TensorField &in);
void IDFT(TensorField & out, const TensorField &in);





Homo_DielectricConstant_static::Homo_DielectricConstant_static(const TwoPhaseDispersion & src, double epsilon1, double epsilon2)
{
	this->RedefineProblem(src, epsilon1, epsilon2);
}

void Homo_DielectricConstant_static::RedefineProblem(const TwoPhaseDispersion & src, double epsilon1, double epsilon2)
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
	this->epsilon1 = epsilon1;	this->epsilon2 = epsilon2;	this->delta_epsilon = epsilon2 - epsilon1;

	if (!SameSimulationBox) {
		update_G = false;	update_plan = false;
		G_Fourier = TensorField(d, 2, Geometry.GetResolutions());
	}

	/* define fields */
	P_Fourier = TensorField(d, 1, Geometry.GetResolutions());
	E_Fourier = TensorField(d, 1, Geometry.GetResolutions());
	P_direct = TensorField(d, 1, Geometry.GetResolutions());
	E_direct = TensorField(d, 1, Geometry.GetResolutions());
	m = 0;
}

void Homo_DielectricConstant_static::Prepare(size_t direction)
{
	this->direction = direction;
	/* Compute the plan for FFT */
	if (!update_plan) {
		std::cout << "fftw_plan is being prepared\n";

		if (d == 2) {
			for (int i = 0; i < d; i++) {
				/* Fast Fourier transform plan for the ith component. */
				fftw_plan_with_nthreads(omp_get_num_threads());
				FourierTransform_P.emplace_back(fftw_plan_dft_2d(P_Fourier.num_points[0], P_Fourier.num_points[1],
					P_direct.field[i], P_Fourier.field[i], -1, FFTW_ESTIMATE));
				fftw_plan_with_nthreads(omp_get_num_threads());
				InverseFourierTransform_E.emplace_back(fftw_plan_dft_2d(E_Fourier.num_points[0], E_Fourier.num_points[1],
					E_Fourier.field[i], E_direct.field[i], 1, FFTW_ESTIMATE));
			}
		}
		else if (d == 3) {
			for (size_t i = 0; i < d; i++) {
				/* Fast Fourier transform plan for the ith component. */
				fftw_plan_with_nthreads(omp_get_num_threads());
				FourierTransform_P.emplace_back(fftw_plan_dft_3d(P_Fourier.num_points[0], P_Fourier.num_points[1], P_Fourier.num_points[2],
					P_direct.field[i], P_Fourier.field[i], -1, FFTW_ESTIMATE));
				fftw_plan_with_nthreads(omp_get_num_threads());
				InverseFourierTransform_E.emplace_back(fftw_plan_dft_3d(E_Fourier.num_points[0], E_Fourier.num_points[1], E_Fourier.num_points[2],
					E_Fourier.field[i], E_direct.field[i], 1, FFTW_ESTIMATE));
			}
		}
		else
			std::cerr << "unsupported dimension!!\n";
	}

	/* Compute the Fourier transform of the Green's function */
	if (!update_G) {
		/* The Fourier transform of the Green's function of electric field is
		 * G_{ij} = -k_i k_j /k^2 */
		std::vector<size_t> SystemSize = Geometry.GetResolutions();
		G_Fourier = TensorField(d, 2, SystemSize);
		k = TensorField(d, 1, SystemSize);
		std::cout << "Green's function is being initialized!!\n";
		std::vector<double> a;
		for (int i = 0; i < d; i++)
			a.push_back(sqrt(Geometry.packing.GetBasisVector(i).Modulus2()) / SystemSize[i]);

		if (d == 2) {
#pragma omp parallel for num_threads(omp_get_num_threads())
			for (int k_index = 0; k_index < G_Fourier.num_total_pixels; k_index++) {
				std::vector<size_t> comp;
				Geometry.Index2Components(k_index, comp);
				double kx = (double)comp[0] / SystemSize[0],
					ky = (double)comp[1] / SystemSize[1];
				/* we didn't multiply this wavevector with 2*\pi because this factor will be canceled out. */
				kx = (kx - round(kx)) / a[0];
				ky = (ky - round(ky)) / a[1];
				double k2 = kx * kx + ky * ky;
				k.field[0][k_index][0] = kx;
				k.field[1][k_index][0] = ky;
				if (k2 > 0) {
					/*G_{0,0}*/
					G_Fourier.field[0][k_index][0] = -kx * kx / k2 / epsilon1;
					/*G_{1,0}*/
					G_Fourier.field[1][k_index][0] = -ky * kx / k2 / epsilon1;
					/*G_{0,1}*/
					G_Fourier.field[2][k_index][0] = -kx * ky / k2 / epsilon1;
					/*G_{1,1}*/
					G_Fourier.field[3][k_index][0] = -ky * ky / k2 / epsilon1;
				}
			}
		}
		else if (d == 3) {
#pragma omp parallel for num_threads(omp_get_num_threads())
			for (int k_index = 0; k_index < G_Fourier.num_total_pixels; k_index++) {
				std::vector<size_t> comp;
				Geometry.Index2Components(k_index, comp);
				double kx = (double)comp[0] / SystemSize[0],
					ky = (double)comp[1] / SystemSize[1],
					kz = (double)comp[2] / SystemSize[2];
				/* we didn't multiply this wavevector with 2*\pi because this factor will be canceled out. */
				kx = (kx - round(kx)) / a[0];
				ky = (ky - round(ky)) / a[1];
				kz = (kz - round(kz)) / a[2];

				k.field[0][k_index][0] = kx;
				k.field[1][k_index][0] = ky;
				k.field[2][k_index][0] = kz;

				double k2 = kx * kx + ky * ky + kz * kz;
				if (k2 > 0) {
					/*G_{0,0}*/
					G_Fourier.field[0][k_index][0] = -kx * kx / k2;
					/*G_{1,0}*/
					G_Fourier.field[1][k_index][0] = -ky * kx / k2;
					/*G_{2,0}*/
					G_Fourier.field[2][k_index][0] = -kz * kx / k2;

					/*G_{0,1}*/
					G_Fourier.field[3][k_index][0] = -kx * ky / k2;
					/*G_{1,1}*/
					G_Fourier.field[4][k_index][0] = -ky * ky / k2;
					/*G_{2,1}*/
					G_Fourier.field[5][k_index][0] = -kz * ky / k2;

					/*G_{0,2}*/
					G_Fourier.field[6][k_index][0] = -kx * kz / k2;
					/*G_{1,2}*/
					G_Fourier.field[7][k_index][0] = -ky * kz / k2;
					/*G_{2,2}*/
					G_Fourier.field[8][k_index][0] = -kz * kz / k2;
				}
			}
		}
		else
			std::cerr << "Green function is not defined in " << d << " dimensions\n";
	}

	m = 0;
	/* Initialize E field
	 *	E(x) = e_direction. */
	E_direct = TensorField(d, 1, Geometry.GetResolutions());
	for (size_t x = 0; x < E_direct.num_total_pixels; x++) {
		E_direct.field[direction][x][0] = 1.0L;
		//E_direct.field[direction][x][1] = 0.0L;
	}

	/* Initialize P field
	 *	P(x) = (epsilon_2 - epsilon1) * I(x) * E(x). */
	P_direct = TensorField(d, 1, Geometry.GetResolutions());
	if (!Geometry.update_indicator)
		Geometry.ComputeIndicator();

	//#pragma omp parallel for 
	for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
		size_t x = Geometry.prt_indicator[i];
		for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
			P_direct.field[comp][x][0] = delta_epsilon * Geometry.WeightedPhaseIndicator(i) * E_direct.field[comp][x][0];
			P_direct.field[comp][x][1] = delta_epsilon * Geometry.WeightedPhaseIndicator(i) * E_direct.field[comp][x][1];
		}
	}

}

void Homo_DielectricConstant_static::Iterate_test(size_t steps) {
	for (int i = 0; i < steps; i++) {
		std::cout << "\t\t Fourier transform of P (" << m << " th stage)\n";
		/* discrete Fourier transform */
		DFT(P_Fourier, P_direct);

		//fftw_execute(FourierTransform_P[j]);

		double interal_factor = 1.0 / P_Fourier.num_total_pixels;
		//for (int i = 0; i < P_Fourier.total_comp; i++) {
		//	for (int x = 0; x < P_Fourier.num_total_pixels; x++) {
		//		P_Fourier.field[i][x][0] *= interal_factor;
		//		P_Fourier.field[i][x][1] *= interal_factor;
		//	}
		//}

		if (ConvergenceTest())
			std::cout << "good!\n";
		else
			std::cout << "convergence test fails\n";

		m++;
		std::cout << "------------------" << m << " th stage ---------------\n";
		std::cout << "\t\t Compute E field \n";
		/* E(k) <- E(k) + G(k).P(k)*/
		//InnerProduct(temp1, G_Fourier, P_Fourier, 1);
		//E_Fourier += temp1;
		InnerProduct(E_Fourier, G_Fourier, P_Fourier, 1);
		for (int i = 0; i < E_Fourier.total_comp; i++) {
			for (int x = 0; x < E_Fourier.num_total_pixels; x++) {
				//E_Fourier.field[i][x][0] *= (interal_factor*interal_factor);
				//E_Fourier.field[i][x][1] *= (interal_factor*interal_factor);
			}
		}
		E_Fourier.field[direction][0][0] = 1.0L;// *E_Fourier.num_total_pixels;
		E_Fourier.field[direction][0][1] = 0.0L;

		std::cout << "\t\t Inverse Fourier transform of E field\n";
		IDFT(E_direct, E_Fourier);

		//fftw_execute(InverseFourierTransform_E[j]);
	/*for (int i = 0; i < E_direct.total_comp; i++) {
	for (int x = 0; x < E_direct.num_total_pixels; x++) {
	E_direct.field[i][x][0] *= interal_factor;
	E_direct.field[i][x][1] *= interal_factor;
	}
	}*/
		std::cout << "\t\t Compute P field\n";
#pragma omp parallel for 
		for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
			size_t x = Geometry.prt_indicator[i];
			for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
				P_direct.field[comp][x][0] = delta_epsilon * (Geometry.WeightedPhaseIndicator(i) + 1.0) * E_direct.field[comp][x][0];
				P_direct.field[comp][x][1] = delta_epsilon * (Geometry.WeightedPhaseIndicator(i) + 1.0) * E_direct.field[comp][x][1];
			}
		}
	}
}

void Homo_DielectricConstant_static::Iterate_test2(size_t steps) {
	for (int i = 0; i < steps; i++) {
		std::cout << "\t\t Fourier transform of P (" << m << " th stage)\n";
		for (int j = 0; j < P_direct.total_comp; j++)
			fftw_execute(FourierTransform_P[j]);

		//P_Fourier.Print("Test_tildeP_",0);
		double integral_factor = 1.0 / P_Fourier.num_total_pixels;
		//P_Fourier *= integral_factor;

		if (ConvergenceTest())
			std::cout << "good!\n";
		else
			std::cout << "convergence test fails\n";

		TensorField D_Fourier(P_Fourier);
		D_Fourier += (E_Fourier);
		m++;
		std::cout << "------------------" << m << " th stage ---------------\n";
		std::cout << "\t\t Compute E field \n";
		/* E(k) <- E(k) + G(k).P(k)*/
		TensorField E_prev(E_Fourier);
		InnerProduct(E_Fourier, G_Fourier, P_Fourier, 1);


		//E_Fourier *= sqrt(integral_factor);
		//E_Fourier += E_prev;
		E_Fourier.field[direction][0][0] = 1.0L *E_Fourier.num_total_pixels;
		E_Fourier.field[direction][0][1] = 0.0L;

		//E_Fourier.Print("Test_tildeE_", 0);
		std::cout << "\t\t Inverse Fourier transform of E field\n";
		for (int j = 0; j < E_direct.total_comp; j++)
			fftw_execute(InverseFourierTransform_E[j]);
		E_direct *= integral_factor;


		std::cout << "\t\t Compute P field\n";
#pragma omp parallel for 
		for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
			size_t x = Geometry.prt_indicator[i];
			for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
				P_direct.field[comp][x][0] = delta_epsilon * Geometry.WeightedPhaseIndicator(i) * E_direct.field[comp][x][0];
				P_direct.field[comp][x][1] = delta_epsilon * Geometry.WeightedPhaseIndicator(i) * E_direct.field[comp][x][1];
			}
		}
	}
}

void Homo_DielectricConstant_static::Iterate(size_t steps) {
	for (int i = 0; i < steps; i++) {
		/*1. Compute FFT of P(x) */
		std::cout << "\t\t Fourier transform of P (" << m << " th stage)\n";
		for (int j = 0; j < P_direct.total_comp; j++)
			fftw_execute(FourierTransform_P[j]);

		//P_Fourier.Print("Test_tildeP_",0);
		double integral_factor = 1.0 / P_Fourier.num_total_pixels;
		//P_Fourier *= integral_factor;

		/*2. Convergence test */
		if (ConvergenceTest(1e-4))
			std::cout << "good!\n";
		else
			std::cout << "convergence test fails\n";

		m++;
		std::cout << "------------------" << m << " th stage ---------------\n";
		std::cout << "\t\t Compute E field \n";
		/* E(k) <- E_0(k) + G(k).P(k)*/
		TensorField E_prev(E_Fourier);
		/*3. Compute \tilde{E}(k) = \tilde{G}(k). \tilde{P}(k) if k \neq 0
					 \tilde{E}(k) = \tilde{E_0} */
		InnerProduct(E_Fourier, G_Fourier, P_Fourier, 1);
		//E_Fourier *= integral_factor;	//!!!!!To  Ge: how do you think this step?
		//E_Fourier += E_prev;
		/* Add \tilde{E_0} */
		E_Fourier.field[direction][0][0] = 1.0L *E_Fourier.num_total_pixels;	//!!!!!To  Ge: how do you think this step?
		E_Fourier.field[direction][0][1] = 0.0L;

		//E_Fourier.Print("Test_tildeE_", 0);

		/*4. Compute the inverse FFT of \tilde{E}(k). */
		std::cout << "\t\t Inverse Fourier transform of E field\n";
		for (int j = 0; j < E_direct.total_comp; j++)
			fftw_execute(InverseFourierTransform_E[j]);
		E_direct *= integral_factor;	//This step is required because fftw3 do not consider the 1/N factor.

		/*5. Compute P(x) = (epsilon(x) - epsilon 1) E(x) */
		std::cout << "\t\t Compute P field\n";
#pragma omp parallel for 
		for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
			size_t x = Geometry.prt_indicator[i];
			for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
				P_direct.field[comp][x][0] = delta_epsilon * Geometry.WeightedPhaseIndicator(i) * E_direct.field[comp][x][0];
				P_direct.field[comp][x][1] = delta_epsilon * Geometry.WeightedPhaseIndicator(i) * E_direct.field[comp][x][1];
			}
		}
	}
}

bool Homo_DielectricConstant_static::ConvergenceTest(double tolerance)
{
	/* compute the average of |k. tilde{D}(k)|*/
	TensorField D_Fourier(P_Fourier), Div;
	D_Fourier += E_Fourier;

	double error = 0.0L;
	InnerProduct(Div, k, D_Fourier, 1);	/* Div is a scalar field...*/
	for (int x = 0; x < Div.num_total_pixels; x++) {
		error += sqrt(Div.field[0][x][0] * Div.field[0][x][0] + Div.field[0][x][1] * Div.field[0][x][1]);
	}
	//Div.Print(std::to_string(m));
	error /= (double)Div.num_total_pixels;

	std::cout << "\t\t||k. tilde{D}(k) || = " << error << "\n";

	return (error < tolerance);
}

std::complex<double> Homo_DielectricConstant_static::ComputeEffectiveDielectric()
{
	std::vector<std::complex<double>> E_av(d, std::complex<double>(0.0L, 0.0L));
	std::vector<std::complex<double>> P_av(d, std::complex<double>(0.0L, 0.0L));

	for (int i = 0; i < d; i++) {
		for (int x = 0; x < E_direct.num_total_pixels; x++) {
			E_av[i] += std::complex<double>(E_direct.field[i][x][0], E_direct.field[i][x][1]);
			P_av[i] += std::complex<double>(P_direct.field[i][x][0], P_direct.field[i][x][1]);
		}
		E_av[i] /= (double)E_direct.num_total_pixels;
		P_av[i] /= (double)E_direct.num_total_pixels;
	}

	std::cout << "E_av = ";
	for (int i = 0; i < d; i++)
		std::cout << "[" << E_av[i] << "] \t";
	std::cout << "\n";
	std::cout << "P_av = ";
	for (int i = 0; i < d; i++)
		std::cout << "[" << P_av[i] << "] \t";
	std::cout << "\n------------------------\n";

	double E2 = 0.0L;
	for (int i = 0; i < d; i++)
		E2 += abs(E_av[i])*abs(E_av[i]);

	std::complex<double> P_dot_E(0.0L, 0.0L);
	for (int i = 0; i < d; i++)
		P_dot_E += (P_av[i] * conj(E_av[i]));

	std::complex<double> result = P_dot_E / E2 + std::complex<double>(epsilon1, 0.0L);
	return result;
}

void Homo_DielectricConstant_static::Print(const std::string & prefix) const
{
	std::fstream outputEx((prefix + "_E_x.txt").c_str(), std::fstream::out);
	std::fstream outputEy((prefix + "_E_y.txt").c_str(), std::fstream::out);
	std::fstream outputPx((prefix + "_P_x.txt").c_str(), std::fstream::out);
	std::fstream outputPy((prefix + "_P_y.txt").c_str(), std::fstream::out);

	std::vector<size_t> px = Geometry.GetResolutions();
	std::vector<size_t> comp({ 0, 0 });
	for (size_t j = 0; j < px[1]; j++) {
		for (size_t i = 0; i < px[0]; i++) {
			comp[0] = i; comp[1] = j;
			outputEx << sqrt(E_direct.field[0][Geometry.Components2Index(comp)][0] * E_direct.field[0][Geometry.Components2Index(comp)][0]
				+ E_direct.field[0][Geometry.Components2Index(comp)][1] * E_direct.field[0][Geometry.Components2Index(comp)][1]) << "\t";
			outputEy << sqrt(E_direct.field[1][Geometry.Components2Index(comp)][0] * E_direct.field[1][Geometry.Components2Index(comp)][0]
				+ E_direct.field[1][Geometry.Components2Index(comp)][1] * E_direct.field[1][Geometry.Components2Index(comp)][1]) << "\t";
			outputPx << sqrt(P_direct.field[0][Geometry.Components2Index(comp)][0] * P_direct.field[0][Geometry.Components2Index(comp)][0]
				+ P_direct.field[0][Geometry.Components2Index(comp)][1] * P_direct.field[0][Geometry.Components2Index(comp)][1]) << "\t";
			outputPy << sqrt(P_direct.field[1][Geometry.Components2Index(comp)][0] * P_direct.field[1][Geometry.Components2Index(comp)][0]
				+ P_direct.field[1][Geometry.Components2Index(comp)][1] * P_direct.field[1][Geometry.Components2Index(comp)][1]) << "\t";
		}
		outputEx << "\n";
		outputEy << "\n";
		outputPx << "\n";
		outputPy << "\n";
	}
	outputEx.close();
	outputEy.close();
	outputPx.close();
	outputPy.close();
}

void Homo_DielectricConstant_static::Print_Fourier(const std::string & prefix) const
{
	std::fstream outputEx((prefix + "_tildeE_x.txt").c_str(), std::fstream::out);
	std::fstream outputEy((prefix + "_tildeE_y.txt").c_str(), std::fstream::out);
	std::fstream outputPx((prefix + "_tildeP_x.txt").c_str(), std::fstream::out);
	std::fstream outputPy((prefix + "_tildeP_y.txt").c_str(), std::fstream::out);

	std::vector<size_t> px = Geometry.GetResolutions();
	std::vector<size_t> comp({ 0, 0 });
	for (size_t j = 0; j < px[1]; j++) {
		for (size_t i = 0; i < px[0]; i++) {
			comp[0] = i; comp[1] = j;
			outputEx << sqrt(E_Fourier.field[0][Geometry.Components2Index(comp)][0] * E_Fourier.field[0][Geometry.Components2Index(comp)][0]
				+ E_Fourier.field[0][Geometry.Components2Index(comp)][1] * E_Fourier.field[0][Geometry.Components2Index(comp)][1]) << "\t";
			outputEy << sqrt(E_Fourier.field[1][Geometry.Components2Index(comp)][0] * E_Fourier.field[1][Geometry.Components2Index(comp)][0]
				+ E_Fourier.field[1][Geometry.Components2Index(comp)][1] * E_Fourier.field[1][Geometry.Components2Index(comp)][1]) << "\t";
			outputPx << sqrt(P_Fourier.field[0][Geometry.Components2Index(comp)][0] * P_Fourier.field[0][Geometry.Components2Index(comp)][0]
				+ P_Fourier.field[0][Geometry.Components2Index(comp)][1] * P_Fourier.field[0][Geometry.Components2Index(comp)][1]) << "\t";
			outputPy << sqrt(P_Fourier.field[1][Geometry.Components2Index(comp)][0] * P_Fourier.field[1][Geometry.Components2Index(comp)][0]
				+ P_Fourier.field[1][Geometry.Components2Index(comp)][1] * P_Fourier.field[1][Geometry.Components2Index(comp)][1]) << "\t";
		}
		outputEx << "\n";
		outputEy << "\n";
		outputPx << "\n";
		outputPy << "\n";
	}
	outputEx.close();
	outputEy.close();
	outputPx.close();
	outputPy.close();
}


void Homo_DielectricConstant_static::Print_Grid_2D(const std::string & prefix) const
{
	std::fstream output((prefix + "_grid.txt").c_str(), std::fstream::out);
	//TODO:: implement it.

}

void Homo_DielectricConstant_static::Print_GreenFn_2D(const std::string & prefix) const
{
	std::fstream output0((prefix + "_G_xx.txt").c_str(), std::fstream::out);
	std::fstream output1((prefix + "_G_xy.txt").c_str(), std::fstream::out);
	std::fstream output2((prefix + "_G_yx.txt").c_str(), std::fstream::out);
	std::fstream output3((prefix + "_G_yy.txt").c_str(), std::fstream::out);

	std::vector<size_t> px = Geometry.GetResolutions();
	std::vector<size_t> comp({ 0, 0 });
	for (size_t j = 0; j < px[1]; j++) {
		for (size_t i = 0; i < px[0]; i++) {
			comp[0] = i; comp[1] = j;
			output0 << sqrt(G_Fourier.field[0][Geometry.Components2Index(comp)][0] * G_Fourier.field[0][Geometry.Components2Index(comp)][0]
				+ G_Fourier.field[0][Geometry.Components2Index(comp)][1] * G_Fourier.field[0][Geometry.Components2Index(comp)][1]) << "\t";
			output1 << sqrt(G_Fourier.field[1][Geometry.Components2Index(comp)][0] * G_Fourier.field[1][Geometry.Components2Index(comp)][0]
				+ G_Fourier.field[1][Geometry.Components2Index(comp)][1] * G_Fourier.field[1][Geometry.Components2Index(comp)][1]) << "\t";
			output2 << sqrt(G_Fourier.field[2][Geometry.Components2Index(comp)][0] * G_Fourier.field[2][Geometry.Components2Index(comp)][0]
				+ G_Fourier.field[2][Geometry.Components2Index(comp)][1] * G_Fourier.field[2][Geometry.Components2Index(comp)][1]) << "\t";
			output3 << sqrt(G_Fourier.field[3][Geometry.Components2Index(comp)][0] * G_Fourier.field[3][Geometry.Components2Index(comp)][0]
				+ G_Fourier.field[3][Geometry.Components2Index(comp)][1] * G_Fourier.field[3][Geometry.Components2Index(comp)][1]) << "\t";
		}
		output0 << "\n";
		output1 << "\n";
		output2 << "\n";
		output3 << "\n";
	}
	output0.close();
	output1.close();
	output2.close();
	output3.close();
}


void DFT(TensorField & out, const TensorField &in) {
	out = TensorField(in.d, in.rank, in.num_points);
	if (in.d == 1) {
		for (int j = 0; j < out.total_comp; j++) {
			for (size_t k1 = 0; k1 < out.num_points[0]; k1++) {
				out.field[j][k1][0] = 0.0L;
				out.field[j][k1][1] = 0.0L;

				for (int x1 = 0; x1 < in.num_points[0]; x1++) {
					double real = in.field[j][x1][0],
						imag = in.field[j][x1][1];
					double phase1 = 2.0*pi*(double)(x1 * k1) / (double)in.num_points[0];

					std::complex<double> factor = std::exp(std::complex<double>(0.0, -phase1));
					out.field[j][k1][0] +=
						(real * factor.real() + imag * factor.imag());
					out.field[j][k1][1] +=
						(imag * factor.real() - real * factor.imag());
				}
			}
		}
	}


	else if (in.d == 2) {
		for (int j = 0; j < out.total_comp; j++) {
			for (size_t k1 = 0; k1 < out.num_points[0]; k1++) {
				for (size_t k2 = 0; k2 < out.num_points[1]; k2++) {
					std::vector<size_t> k({ k1,k2 });
					std::vector<size_t> x({ 0,0 });
					out.field[j][k2 + out.num_points[1] * k1][0] = 0.0L;
					out.field[j][k2 + out.num_points[1] * k1][1] = 0.0L;

					for (int x1 = 0; x1 < in.num_points[0]; x1++) {
						x[0] = x1;
						for (int x2 = 0; x2 < in.num_points[1]; x2++) {
							x[1] = x2;
							double real = in.field[j][x2 + in.num_points[1] * x1][0],
								imag = in.field[j][x2 + in.num_points[1] * x1][1];
							double phase1 = 2.0*pi*(double)(x1 * k1) / (double)in.num_points[0];
							double phase2 = 2.0*pi*(double)(x2 * k2) / (double)in.num_points[1];

							std::complex<double> factor = std::exp(std::complex<double>(0.0, -phase1)) * std::exp(std::complex<double>(0.0, -phase2));
							out.field[j][k2 + out.num_points[1] * k1][0] +=
								(real * factor.real() - imag * factor.imag());
							out.field[j][k2 + out.num_points[1] * k1][1] +=
								(imag * factor.real() + real * factor.imag());
						}
					}
				}
			}
		}
	}

}

void IDFT(TensorField & out, const TensorField &in) {
	out = TensorField(in.d, in.rank, in.num_points);
	if (in.d == 1) {
		for (int j = 0; j < out.total_comp; j++) {
			for (size_t k1 = 0; k1 < out.num_points[0]; k1++) {
				out.field[j][k1][0] = 0.0L;
				out.field[j][k1][1] = 0.0L;

				for (int x1 = 0; x1 < in.num_points[0]; x1++) {
					double real = in.field[j][x1][0],
						imag = in.field[j][x1][1];
					double phase1 = 2.0*pi*(double)(x1 * k1) / (double)in.num_points[0];

					std::complex<double> factor = std::exp(std::complex<double>(0.0, -phase1));
					out.field[j][k1][0] +=
						(real * factor.real() - imag * factor.imag());
					out.field[j][k1][1] +=
						(imag * factor.real() + real * factor.imag());
				}
			}
		}
	}


	else if (in.d == 2) {
		for (int j = 0; j < out.total_comp; j++) {
			for (size_t k1 = 0; k1 < out.num_points[0]; k1++) {
				for (size_t k2 = 0; k2 < out.num_points[1]; k2++) {
					std::vector<size_t> k({ k1,k2 });
					std::vector<size_t> x({ 0,0 });
					out.field[j][k2 + out.num_points[1] * k1][0] = 0.0L;
					out.field[j][k2 + out.num_points[1] * k1][1] = 0.0L;

					for (int x1 = 0; x1 < in.num_points[0]; x1++) {
						x[0] = x1;
						for (int x2 = 0; x2 < in.num_points[1]; x2++) {
							x[1] = x2;
							double real = in.field[j][x2 + in.num_points[1] * x1][0],
								imag = in.field[j][x2 + in.num_points[1] * x1][1];
							double phase1 = 2.0*pi*(double)(x1 * k1) / (double)in.num_points[0];
							double phase2 = 2.0*pi*(double)(x2 * k2) / (double)in.num_points[1];

							std::complex<double> factor = std::exp(std::complex<double>(0.0, phase1)) * std::exp(std::complex<double>(0.0, phase2));
							out.field[j][k2 + out.num_points[1] * k1][0] +=
								(real * factor.real() - imag * factor.imag());
							out.field[j][k2 + out.num_points[1] * k1][1] +=
								(imag * factor.real() + real * factor.imag());
						}
					}
				}
			}
		}
	}


	for (int i = 0; i < out.total_comp; i++) {
		for (int x = 0; x < out.num_total_pixels; x++) {
			out.field[i][x][0] /= (double)out.num_total_pixels;
			out.field[i][x][1] /= (double)out.num_total_pixels;
		}
	}
}



Homo_DielectricConstant_static_SCA::Homo_DielectricConstant_static_SCA(const TwoPhaseDispersion & src, double epsilon1, double epsilon2)
{
	this->RedefineProblem(src, epsilon1, epsilon2);
}

void Homo_DielectricConstant_static_SCA::RedefineProblem(const TwoPhaseDispersion & src, double epsilon1, double epsilon2)
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

	this->epsilon1 = epsilon1;	this->epsilon2 = epsilon2;	this->delta_epsilon = epsilon2 - epsilon1;
	this->L = (double)d*epsilon1*(epsilon2 - epsilon1) / (epsilon2 + (d - 1.0) *epsilon1);

	if (!SameSimulationBox) {
		update_H = false;	update_plan = false;
		H_Fourier = TensorField(d, 2, Geometry.GetResolutions());
	}

	/* define fields */
	P_Fourier = TensorField(d, 1, Geometry.GetResolutions());
	F_Fourier = TensorField(d, 1, Geometry.GetResolutions());
	P_direct = TensorField(d, 1, Geometry.GetResolutions());
	F_direct = TensorField(d, 1, Geometry.GetResolutions());
	m = 0;
}

void Homo_DielectricConstant_static_SCA::Prepare(size_t direction)
{
	this->direction = direction;
	/* Compute the plan for FFT */
	if (!update_plan) {
		std::cout << "fftw_plan is being prepared\n";
		fftw_plan_with_nthreads(omp_get_num_threads());
		if (d == 2) {
			for (size_t i = 0; i < d; i++) {
				/* Fast Fourier transform plan for the ith component. */
				FourierTransform_P.emplace_back(fftw_plan_dft_2d(P_Fourier.num_points[0], P_Fourier.num_points[1],
					P_direct.field[i], P_Fourier.field[i], -1, FFTW_ESTIMATE));
				InverseFourierTransform_F.emplace_back(fftw_plan_dft_2d(F_Fourier.num_points[0], F_Fourier.num_points[1],
					F_Fourier.field[i], F_direct.field[i], 1, FFTW_ESTIMATE));
			}
		}
		else if (d == 3) {
			for (size_t i = 0; i < d; i++) {
				/* Fast Fourier transform plan for the ith component. */
				FourierTransform_P.emplace_back(fftw_plan_dft_3d(P_Fourier.num_points[0], P_Fourier.num_points[1], P_Fourier.num_points[2],
					P_direct.field[i], P_Fourier.field[i], -1, FFTW_ESTIMATE));
				InverseFourierTransform_F.emplace_back(fftw_plan_dft_3d(F_Fourier.num_points[0], F_Fourier.num_points[1], F_Fourier.num_points[2],
					F_Fourier.field[i], F_direct.field[i], 1, FFTW_ESTIMATE));
			}
		}
		else
			std::cerr << "unsupported dimension!!\n";
	}

	/* Compute the Fourier transform of the Green's function */
	if (!update_H) {
		this->ComputeH_Fourier();
	}

	m = 0;
	/* Initialize E field
	*	E(x) = e_direction. */
	F_direct = TensorField(d, 1, Geometry.GetResolutions());
	for (size_t x = 0; x < F_direct.num_total_pixels; x++) {
		F_direct.field[direction][x][0] = 1.0L;
		//E_direct.field[direction][x][1] = 0.0L;
	}

	/* Initialize P field
	*	P(x) = (epsilon_2 - epsilon1) * I(x) * E(x). */
	P_direct = TensorField(d, 1, Geometry.GetResolutions());
	if (!Geometry.update_indicator)
		Geometry.ComputeIndicator();

	//#pragma omp parallel for 
	for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
		size_t x = Geometry.prt_indicator[i];
		for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
			P_direct.field[comp][x][0] = L * Geometry.WeightedPhaseIndicator(i) * F_direct.field[comp][x][0];
			P_direct.field[comp][x][1] = L * Geometry.WeightedPhaseIndicator(i) * F_direct.field[comp][x][1];
		}
	}

}

void Homo_DielectricConstant_static_SCA::Iterate(size_t steps) {
	for (int i = 0; i < steps; i++) {
		/*1. Compute FFT of P(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Fourier transform of P (" << m << " th stage)\n";
		for (int j = 0; j < P_direct.total_comp; j++)
			fftw_execute(FourierTransform_P[j]);

		m++;
		if (FFH_verbosity > 0) {
			std::cout << "------------------" << m << " th stage ---------------\n";
			std::cout << "\t\t Compute F field \n";
		}
		/* F(k) <- E_0(k) + H(k).P(k)*/
		TensorField F_prev(F_Fourier);
		/*3. Compute \tilde{E}(k) = \tilde{G}(k). \tilde{P}(k) if k \neq 0
		\tilde{E}(k) = \tilde{E_0} */
		InnerProduct(F_Fourier, H_Fourier, P_Fourier, 1);
		/* Add \tilde{E_0} */
		F_Fourier.field[direction][0][0] = 1.0L *F_Fourier.num_total_pixels;	//!!!!!To  Ge: how do you think this step?
		F_Fourier.field[direction][0][1] = 0.0L;

		/*4. Compute the inverse FFT of \tilde{E}(k). */
		F_direct_prev = F_direct;
		if (FFH_verbosity > 0)
			std::cout << "\t\t Inverse Fourier transform of F field\n";
		for (int j = 0; j < F_direct.total_comp; j++)
			fftw_execute(InverseFourierTransform_F[j]);
		F_direct *= (1.0 / P_Fourier.num_total_pixels);	//This step is required because fftw3 do not consider the 1/N factor.

		/*5. Compute P(x) = L(x)* F(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Compute P field\n";
#pragma omp parallel for 
		for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
			size_t x = Geometry.prt_indicator[i];
			for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
				P_direct.field[comp][x][0] = L * Geometry.WeightedPhaseIndicator(i) * F_direct.field[comp][x][0];
				P_direct.field[comp][x][1] = L * Geometry.WeightedPhaseIndicator(i) * F_direct.field[comp][x][1];
			}
		}
	}
}

/* Iteratively solve the equation until achieving the mean change in field is smaller than [tol_field]
 and the mean divergence of the D field is smaller than [tol_div].
 @param[in] tol_field	The target tolerance in the mean change in the local cavity field
 @param[in] tol_div		The target tolerance in the mean divergence of the D field. */
void Homo_DielectricConstant_static_SCA::Run(double tol_field, double tol_div)
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

/* Compute the divergence of the local cavity field in the Fourier space */
bool Homo_DielectricConstant_static_SCA::CheckDivergenceFree(double tolerance)
{
	/* compute the Fourier transform of the D field.
	 * D = epsilon_1*E + P = epsilon1* F + (1-1/d)*P */
	TensorField D_Fourier(P_Fourier), Div, epsilon1F_Fourier(F_Fourier);
	D_Fourier *= (1.0 - 1.0 / d);
	epsilon1F_Fourier *= epsilon1;
	D_Fourier += epsilon1F_Fourier;

	/* compute the average of sum |k. tilde{D}(k)|/N*/
	double error = 0.0L;
	InnerProduct(Div, k, D_Fourier, 1);	/* Div is a scalar field...*/
	for (int x = 0; x < Div.num_total_pixels; x++) {
		error += sqrt(Div.field[0][x][0] * Div.field[0][x][0] + Div.field[0][x][1] * Div.field[0][x][1]);
	}
	//Div.Print(std::to_string(m));
	error /= (double)Div.num_total_pixels;

	double mean = 0.0L;
	for (int i = 0; i < d; i++)
		mean += (D_Fourier.field[i][0][0] * D_Fourier.field[i][0][0] + D_Fourier.field[i][0][1] * D_Fourier.field[i][0][1]);
	mean /= (double)D_Fourier.num_total_pixels;
	error = sqrt(error) / mean;

	std::cout << "\t\t||k. tilde{D}(k) || = " << error << "\n";
	Tol_div = error;
	return (error < tolerance);
}

/* Compute the difference in F_direct at the current stage and the previous stage.
 * @param[in] tol_field	The tolerance in */
bool Homo_DielectricConstant_static_SCA::CheckConvergence(double tol_field)
{
	F_direct_prev -= F_direct;
	double error = 0.0L;
	for (int x = 0; x < F_direct_prev.num_total_pixels; x++) {
		double sum_over_comp = 0.0L;
		for (int i = 0; i < F_direct_prev.total_comp; i++) {
			sum_over_comp += (F_direct_prev.field[i][x][0] * F_direct_prev.field[i][x][0] + F_direct_prev.field[i][x][1] * F_direct_prev.field[i][x][1]);
		}
		error += sqrt(sum_over_comp);
	}
	if (FFH_verbosity == 2) {
		std::cout << "errors in the local cavity field = " << error << "\n";
	}
	Tol_F = error;
	//F_direct_prev += F_direct;	/**/
	return error < tol_field;
}

std::complex<double> Homo_DielectricConstant_static_SCA::ComputeEffectiveDielectric()
{
	std::vector<std::complex<double>> F_av(d, std::complex<double>(0.0L, 0.0L));
	std::vector<std::complex<double>> P_av(d, std::complex<double>(0.0L, 0.0L));

	for (int i = 0; i < d; i++) {
		for (int x = 0; x < F_direct.num_total_pixels; x++) {
			F_av[i] += std::complex<double>(F_direct.field[i][x][0], F_direct.field[i][x][1]);
			P_av[i] += std::complex<double>(P_direct.field[i][x][0], P_direct.field[i][x][1]);
		}
		F_av[i] /= (double)F_direct.num_total_pixels;
		P_av[i] /= (double)P_direct.num_total_pixels;
	}

	std::cout << "F_av = ";
	for (int i = 0; i < d; i++)
		std::cout << "[" << F_av[i] << "] \t";
	std::cout << "\n";
	std::cout << "P_av = ";
	for (int i = 0; i < d; i++)
		std::cout << "[" << P_av[i] << "] \t";
	std::cout << "\n------------------------\n";

	double P2 = 0.0L;
	for (int i = 0; i < d; i++)
		P2 += abs(P_av[i])*abs(P_av[i]);

	std::complex<double> F_dot_P(0.0L, 0.0L);
	for (int i = 0; i < d; i++)
		F_dot_P += (F_av[i] * conj(P_av[i]));

	std::complex<double> Le = P2 / F_dot_P, result;
	result = epsilon1 * (1.0 - (double)d + (double)d / (1.0 - Le / ((double)d*epsilon1)));
	return result;
}

void Homo_DielectricConstant_static_SCA::Print(const std::string & prefix) const
{
	F_direct.Print((prefix + "_F").c_str(), 0);
	F_direct.Print((prefix + "_F").c_str(), 1);
	P_direct.Print((prefix + "_P").c_str(), 0);
	P_direct.Print((prefix + "_P").c_str(), 1);
}

void Homo_DielectricConstant_static_SCA::Print_Fourier(const std::string & prefix) const
{
	F_Fourier.Print((prefix + "_tildeF").c_str(), 0);
	F_Fourier.Print((prefix + "_tildeF").c_str(), 1);
	P_Fourier.Print((prefix + "_tildeP").c_str(), 0);
	P_Fourier.Print((prefix + "_tildeP").c_str(), 1);
}


void Homo_DielectricConstant_static_SCA::Print_Grid_2D(const std::string & prefix) const
{
	std::fstream output((prefix + "_grid.txt").c_str(), std::fstream::out);
	//TODO:: implement it.

}

void Homo_DielectricConstant_static_SCA::Print_GreenFn_2D(const std::string & prefix) const
{
	H_Fourier.Print((prefix + "_tildeH_xx").c_str(), 0);
	H_Fourier.Print((prefix + "_tildeH_yx").c_str(), 1);
	H_Fourier.Print((prefix + "_tildeH_xy").c_str(), 2);
	H_Fourier.Print((prefix + "_tildeH_yy").c_str(), 3);
}

void Homo_DielectricConstant_static_SCA::ComputeH_Fourier()
{
	//TODO: treat the highest frequencies as discussed.
/* The Fourier transform of the Green's function of electric field is
* H_{ij} = ... -k_i k_j /k^2 */
	std::vector<size_t> SystemSize = Geometry.GetResolutions();
	H_Fourier = TensorField(d, 2, SystemSize);
	k = TensorField(d, 1, SystemSize);
	std::cout << "Green's function is being initialized!!\n";
	std::vector<double> a;
	for (int i = 0; i < d; i++)
		a.push_back(sqrt(Geometry.packing.GetBasisVector(i).Modulus2()) / SystemSize[i]);

	if (d == 2) {
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (int k_index = 0; k_index < H_Fourier.num_total_pixels; k_index++) {
			std::vector<size_t> comp;
			Geometry.Index2Components(k_index, comp);
			double kx = (double)comp[0] / SystemSize[0],
				ky = (double)comp[1] / SystemSize[1];
			/* we didn't multiply this wavevector with 2*\pi because this factor will be canceled out. */
			kx = (kx - round(kx)) / a[0];
			ky = (ky - round(ky)) / a[1];
			double k2 = kx * kx + ky * ky;
			k.field[0][k_index][0] = kx;
			k.field[1][k_index][0] = ky;
			if (k2 > 0) {
				/*H_{0,0}*/
				H_Fourier.field[0][k_index][0] = (0.5 - kx * kx / k2) / epsilon1;
				/*H_{1,0}*/
				H_Fourier.field[1][k_index][0] = -ky * kx / k2 / epsilon1;
				/*H_{0,1}*/
				H_Fourier.field[2][k_index][0] = -kx * ky / k2 / epsilon1;
				/*H_{1,1}*/
				H_Fourier.field[3][k_index][0] = (0.5 - ky * ky / k2) / epsilon1;
			}
		}
	}
	else if (d == 3) {
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (int k_index = 0; k_index < H_Fourier.num_total_pixels; k_index++) {
			std::vector<size_t> comp;
			Geometry.Index2Components(k_index, comp);
			double kx = (double)comp[0] / SystemSize[0],
				ky = (double)comp[1] / SystemSize[1],
				kz = (double)comp[2] / SystemSize[2];
			/* we didn't multiply this wavevector with 2*\pi because this factor will be canceled out. */
			kx = (kx - round(kx)) / a[0];
			ky = (ky - round(ky)) / a[1];
			kz = (kz - round(kz)) / a[2];

			k.field[0][k_index][0] = kx;
			k.field[1][k_index][0] = ky;
			k.field[2][k_index][0] = kz;

			double k2 = kx * kx + ky * ky + kz * kz;
			if (k2 > 0) {
				/*H_{0,0}*/
				H_Fourier.field[0][k_index][0] = (1.0 / 3.0 - kx * kx / k2) / epsilon1;
				/*H_{1,0}*/
				H_Fourier.field[1][k_index][0] = -ky * kx / k2;
				/*H_{2,0}*/
				H_Fourier.field[2][k_index][0] = -kz * kx / k2;

				/*H_{0,1}*/
				H_Fourier.field[3][k_index][0] = -kx * ky / k2;
				/*H_{1,1}*/
				H_Fourier.field[4][k_index][0] = (1.0 / 3.0 - ky * ky / k2) / epsilon1;
				/*H_{2,1}*/
				H_Fourier.field[5][k_index][0] = -kz * ky / k2;

				/*H_{0,2}*/
				H_Fourier.field[6][k_index][0] = -kx * kz / k2;
				/*H_{1,2}*/
				H_Fourier.field[7][k_index][0] = -ky * kz / k2;
				/*H_{2,2}*/
				H_Fourier.field[8][k_index][0] = (1.0 / 3.0 - kz * kz / k2) / epsilon1;
			}
		}
	}
	else
		std::cerr << "Green function is not defined in " << d << " dimensions\n";

}


Homo_DielectricConstant_dynamic_SCA::Homo_DielectricConstant_dynamic_SCA(const TwoPhaseDispersion & src, double epsilon1, double epsilon2, const ConstantTensor & k1)
{
	this->RedefineDynamicProblem(src, epsilon1, epsilon2, k1);
}

void Homo_DielectricConstant_dynamic_SCA::RedefineDynamicProblem(const TwoPhaseDispersion & src, double epsilon1, double epsilon2, const ConstantTensor & k1)
{
	this->Homo_DielectricConstant_dynamic_SCA::RedefineProblem(src, epsilon1, epsilon2);
	this->RedefineWavevector(k1);
}

void Homo_DielectricConstant_dynamic_SCA::RedefineWavevector(const ConstantTensor & k1)
{
	if (k1.d == d && k1.rank == 1) {
		this->k1 = k1;
		this->update_H = false;
	}
	else {
		std::cerr << "The wavevector k1 should be rank-1 tensor in d dimensions\n";
	}
}

void Homo_DielectricConstant_dynamic_SCA::Prepare(size_t direction)
{
	this->Homo_DielectricConstant_static_SCA::Prepare(direction);
}

void Homo_DielectricConstant_dynamic_SCA::Iterate(size_t steps)
{
	for (int i = 0; i < steps; i++) {
		/*1. Compute FFT of P(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Fourier transform of P (" << m << " th stage)\n";
		for (int j = 0; j < P_direct.total_comp; j++)
			fftw_execute(FourierTransform_P[j]);

		m++;
		if (FFH_verbosity > 0) {
			std::cout << "------------------" << m << " th stage ---------------\n";
			std::cout << "\t\t Compute F field \n";
		}
		/* F(k) <- E_0(k) + H(k).P(k)*/
		TensorField F_prev(F_Fourier);
		/*3. Compute \tilde{E}(k) = \tilde{G}(k). \tilde{P}(k) if k \neq 0
		\tilde{E}(k) = \tilde{E_0} */
		InnerProduct(F_Fourier, H_Fourier, P_Fourier, 1);
		/* Add \tilde{E_0} */
		F_Fourier.field[direction][0][0] = 1.0L *F_Fourier.num_total_pixels;	//!!!!!To  Ge: how do you think this step?
		F_Fourier.field[direction][0][1] = 0.0L;

		/*4. Compute the inverse FFT of \tilde{E}(k). */
		F_direct_prev = F_direct;
		if (FFH_verbosity > 0)
			std::cout << "\t\t Inverse Fourier transform of F field\n";
		for (int j = 0; j < F_direct.total_comp; j++)
			fftw_execute(InverseFourierTransform_F[j]);
		F_direct *= (1.0 / P_Fourier.num_total_pixels);	//This step is required because fftw3 do not consider the 1/N factor.

		/*5. Compute P(x) = L(x)* F(x) */
		if (FFH_verbosity > 0)
			std::cout << "\t\t Compute P field\n";
#pragma omp parallel for 
		for (int i = 0; i < Geometry.prt_indicator.size(); i++) {
			size_t x = Geometry.prt_indicator[i];
			for (size_t comp = 0; comp < P_direct.total_comp; comp++) {
				P_direct.field[comp][x][0] = L * Geometry.WeightedPhaseIndicator(i) * F_direct.field[comp][x][0];
				P_direct.field[comp][x][1] = L * Geometry.WeightedPhaseIndicator(i) * F_direct.field[comp][x][1];
			}
		}
	}
}

void Homo_DielectricConstant_dynamic_SCA::Use_Strong_Contrast(bool use_sc)
{
}

void Homo_DielectricConstant_dynamic_SCA::ComputeG_Fourier()
{
}

/* Compute the Fourier transform of the dyadic Green's function... */
void Homo_DielectricConstant_dynamic_SCA::ComputeH_Fourier()
{

}

{


	{	/* Square array */
		double epsilon1 = 1.0L, phi = 0.6;
		std::vector<GeometryVector> epsilon_eff, errors;



		for (double epsilon2 = 10.0L; epsilon2 <= 10.0L; epsilon2 += 0.1) {
			//ifile >> epsilon2;
			ofile << "epsilon2 = " << epsilon2 << "\n";
			SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
			c.Insert(pow(phi / HyperSphere_Volume(2, 1.0), 0.5), GeometryVector(0.5, 0.5));

			TwoPhaseDispersion packing(c, 201);	//It gives a discretized version of a sphere packing.
			packing.ComputeIndicator();

			// define the problem.
			Homo_DielectricConstant_static_SCA test(packing, epsilon1, epsilon2);

			test.Prepare(0);

			// Iteratively compute the resulting electric field
			test.Run(1e-6, 1e-10);

			{
				TensorField D_Fourier(test.P_Fourier), Div, epsilon1F_Fourier(test.F_Fourier);
				D_Fourier *= (1.0 - 1.0 / (double)test.d);
				epsilon1F_Fourier *= epsilon1;
				D_Fourier += epsilon1F_Fourier;

				TensorField psi(2, 2, D_Fourier.num_points);
				for (int x = 0; x < D_Fourier.num_total_pixels; x++) {
					for (int i = 0; i < 2; i++) {
						for (int j = 0; j < 2; j++) {
							std::complex<double> tildeD_i(D_Fourier[i][x][0], D_Fourier[i][x][1]),
								tildeD_j(D_Fourier[j][x][0], -D_Fourier[j][x][1]),
								result = tildeD_i * tildeD_j;

							psi[i + 2 * j][x][0] = result.real();
							psi[i + 2 * j][x][1] = result.imag();

						}
					}
				}
				psi.Print("SpectralDensity-xx-real_2Dsquare_array_0.6_10.0", 0, 0);
				psi.Print("SpectralDensity-xx-imag_2Dsquare_array_0.6_10.0", 0, 1);
				psi.Print("SpectralDensity-yx-real_2Dsquare_array_0.6_10.0", 1, 0);
				psi.Print("SpectralDensity-yx-imag_2Dsquare_array_0.6_10.0", 1, 1);
				psi.Print("SpectralDensity-xy-real_2Dsquare_array_0.6_10.0", 2, 0);
				psi.Print("SpectralDensity-xy-imag_2Dsquare_array_0.6_10.0", 2, 1);
				psi.Print("SpectralDensity-yy-real_2Dsquare_array_0.6_10.0", 3, 0);
				psi.Print("SpectralDensity-yy-imag_2Dsquare_array_0.6_10.0", 3, 1);

				TensorField kSpectral, k = test.Geometry.GetKs();
				InnerProduct(kSpectral, k, psi, 1);
				kSpectral.Print("k.Spectral_x-real_2Dsquare_array_0.6_10.0", 0, 0);
				kSpectral.Print("k.Spectral_x-imag_2Dsquare_array_0.6_10.0", 0, 1);
				kSpectral.Print("k.Spectral_y-real_2Dsquare_array_0.6_10.0", 1, 0);
				kSpectral.Print("k.Spectral_y-imag_2Dsquare_array_0.6_10.0", 1, 1);
			}

			std::complex<double> effective_property = test.ComputeEffectiveDielectric();
			epsilon_eff.emplace_back(epsilon2, effective_property.real());
			errors.emplace_back(epsilon2, test.m, test.Tol_F, test.Tol_div);
			ofile << "epsilon_eff = " << effective_property << "\n";
			sprintf(name, "Strong-contrast_%0.3f_%0.3f", phi, epsilon2);
			test.Print(name);
		}
		sprintf(name, "epsilon_eff_square_array_%0.3f", phi);
		WriteFunction(epsilon_eff, name);
		sprintf(name, "errors_square_array_%0.3f", phi);
		WriteFunction(errors, name);
	}
	return 0;





	{

		double G1 = 1.0L, G2 = 1.3L;
		double K1 = 1.0L, K2 = 1.3L;
		double phi = 0.05;	//packing fraction
		//ifile >> K2;
		ofile << "K2 = " << K2 << "\n";
		//ifile >> G2;
		ofile << "G2 = " << G2 << "\n";
		SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
		c.Insert(pow(phi / HyperSphere_Volume(2, 1.0), 0.5), GeometryVector(0.5, 0.5));

		TwoPhaseDispersion packing(c, 101);	//It gives a discretized version of a sphere packing.
		packing.ComputeIndicator();

		// define the problem.
		Homo_ElasticProblem_static test(packing, K1, G1, K2, G2);
		test.SetStrongContrast(true);
		test.Prepare(HydrostaticStrain(2));
		test.Print_GreenFn_2D("test_G");	// Print each component of the Fourier transform of the Green's function.
		//test.Print("FFTW_" + std::to_string(phi) + "_0__" + std::to_string(epsilon2 / epsilon1));
		sprintf(name, "Test_strong-%0.2f_%d__%.2f_%0.2f", phi, 0, K2 / K1, G2 / G1);
		test.Print(name);
		//Iteratively compute the resulting electric field
		for (int i = 1; i < 20; i++) {
			test.Iterate();
			test.CheckConvergence(1e-6);
			ofile << "changes in field = " << test.Tol_Field << "\n";
			test.CheckDivergenceFree(1e-10);
			ofile << "divergence in field = " << test.Tol_div << "\n";
			sprintf(name, "Test_strong-%0.2f_%d__%.2f_%0.2f", phi, i, K2 / K1, G2 / G1);
			test.Print(name);
		}
		ofile << "prediction of the strain inside the inclusion = " << (1.0 - test.kappa) << "\n";
		ofile << "prediction of the polarization inside the inclusion = " << (double)test.d *(K1 + 2.0*(1.0 - 1.0 / test.d)*G1) * test.kappa << "\n";

		std::vector<std::complex<double>> moduli;
		test.ComputeEffectiveProperties(moduli);
		ofile << "computed Ke = " << moduli[0] << "\n";
		ofile << "computed Ge = " << moduli[1] << "\n";
		//		std::vector<std::complex<double>> result = test.ComputeEffective
//		ofile << "epsilon_eff = " << test.ComputeEffectiveDielectric() << "\n";


	}
	return 0;

	{
		TensorField A, B, C;
		{
			A = TensorField(2, 2, std::vector<size_t>({ 10,10 }));
		}
		ofile << A[2][1][1] << "\n";

		{
			B = TensorField(2, 3, std::vector<size_t>({ 10,10 }));
		}
		ofile << B[2][1][1] << "\n";
		{
			C = TensorField(2, 4, std::vector<size_t>({ 10,10 }));
		}
		ofile << C[2][1][1] << "\n";

		for (int d = 1; d < 10; d++) {
			ofile << "Dimension = " << d << "\n";
			ofile << "quadruple inner product of hydrostatic: "
				<< InnerProduct(Lambda_h(d), Lambda_h(d), 4)[0] << ",\t"
				<< 1 << "\n";
			ofile << "quadruple inner product of shear: "
				<< InnerProduct(Lambda_s(d), Lambda_s(d), 4)[0] << ",\t"
				<< (d - 1)*(d + 2) / 2 << "\n";
			ofile << "\n\n\n";
		}
	}

	{
		ConstantTensor A2 = Lambda_h(2), B2 = Lambda_s(2);
		ConstantTensor t1, t2, t3;

		t1 = InnerProduct(A2, A2, 2);

		t1 = InnerProduct(A2, B2, 2);

		t1 = InnerProduct(B2, B2, 2);
		std::complex<double> trace1 = 0, trace2 = 0;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++) {
				trace1 += A2[i + 2 * (i + 2 * (j + 2 * j))];
				trace2 += B2[i + 2 * (i + 2 * (j + 2 * j))];
			}

		ConstantTensor A3 = Lambda_h(3), B3 = Lambda_s(3);
		trace1 = 0; trace2 = 0;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++) {
				trace1 += A3[i + 3 * (i + 3 * (j + 3 * j))];
				trace2 += B3[i + 3 * (i + 3 * (j + 3 * j))];
			}
		ofile << "end?\n";
	}

	{
		TensorField test(2, 2, std::vector<size_t>({ 2,2 }));
		for (int i = 0; i < test.total_comp; i++)
			for (int x = 0; x < test.num_total_pixels; x++) {
				test.field[i][x][0] = rng.RandomDouble();
				test.field[i][x][1] = rng.RandomDouble();
			}
		{
			TensorField test1 = 2.0*test;
			ofile << "compare\n";
			for (int i = 0; i < test.total_comp; i++)
				for (int x = 0; x < test.num_total_pixels; x++) {
					ofile << "test: " << test.field[i][x][0] << " + i" << test.field[i][x][1] << "\n";
					ofile << "test1: " << test1.field[i][x][0] << " + i" << test1.field[i][x][1] << "\n";
				}
			TensorField t = test1 - test;
			for (int i = 0; i < test.total_comp; i++)
				for (int x = 0; x < test.num_total_pixels; x++) {
					ofile << "test: " << test.field[i][x][0] << " + i" << test.field[i][x][1] << "\n";
					ofile << "x: " << t.field[i][x][0] << " + i" << t.field[i][x][1] << "\n";
				}

			ofile << "\n\n\n\nchange\n";
			for (int i = 0; i < test.total_comp; i++)
				for (int x = 0; x < test.num_total_pixels; x++) {
					test1.field[i][x][0] = rng.RandomDouble();
					ofile << "test: " << test.field[i][x][0] << " + i" << test.field[i][x][1] << "\n";
					ofile << "test1: " << test1.field[i][x][0] << " + i" << test1.field[i][x][1] << "\n";
				}
		}
		ofile << "\n\n\n\nDelete\n";

		for (int i = 0; i < test.total_comp; i++)
			for (int x = 0; x < test.num_total_pixels; x++) {
				ofile << "test: " << test.field[i][x][0] << " + i" << test.field[i][x][1] << "\n";
			}

	}
	return 0;

	{
		TensorField test(2, 2, std::vector<size_t>({ 3,3 }));
		for (int i = 0; i < test.total_comp; i++)
			for (int x = 0; x < test.num_total_pixels; x++) {
				test.field[i][x][0] = rng.RandomDouble();
				test.field[i][x][1] = rng.RandomDouble();
			}
		{
			TensorField & test1 = test;
			ofile << "compare\n";
			for (int i = 0; i < test.total_comp; i++)
				for (int x = 0; x < test.num_total_pixels; x++) {
					ofile << "test: " << test.field[i][x][0] << " + i" << test.field[i][x][1] << "\n";
					ofile << "test1: " << test1.field[i][x][0] << " + i" << test1.field[i][x][1] << "\n";
				}
			ofile << "\n\n\n\nchange\n";
			for (int i = 0; i < test.total_comp; i++)
				for (int x = 0; x < test.num_total_pixels; x++) {
					test1.field[i][x][0] = rng.RandomDouble();
					ofile << "test: " << test.field[i][x][0] << " + i" << test.field[i][x][1] << "\n";
					ofile << "test1: " << test1.field[i][x][0] << " + i" << test1.field[i][x][1] << "\n";
				}

		}
		ofile << "\n\n\n\nDelete\n";

		for (int i = 0; i < test.total_comp; i++)
			for (int x = 0; x < test.num_total_pixels; x++) {
				ofile << "test: " << test.field[i][x][0] << " + i" << test.field[i][x][1] << "\n";
			}

	}
	return 0;

	{
		std::complex<double> x(1.0, 1.0);
		ofile << x << "\n";

		double y = 0;
		x = y;

		ofile << x << "\n";


	}
	return 0;


	{/* periodic checker board. */
		int num_domains = 4, num_pixels = 200, num_configs = 1;


		std::vector<GeometryVector> epsilon_eff, errors;
		std::vector < double> list;
		for (double t = 0.1; t < 13; t += 0.1)
			list.push_back(t);

		int index = 1; double factor = 0.5;
		for (auto i = 0; i < list.size(); i++) {
			double epsilon2 = list[i];

			GeometryVector e_eff(0, 0, 0, 0);
			GeometryVector error(0, 0, 0, 0);
			TwoPhaseDispersion Checkerboard;
			PeriodicCheckerBoard(Checkerboard, num_domains, num_pixels);

			Homo_DielectricConstant_static_SCA test(Checkerboard, 1.0, epsilon2);
			test.Prepare(0);
			//test.P_direct.Print("test", 0);

			test.Run(1e-6, 1e-10);
			std::complex<double> effective_property = test.ComputeEffectiveDielectric();

			if (abs(epsilon2 - factor * index) < 0.001) {
				sprintf(name, "4x4_checkerboard__epsilon2_%0.2f", epsilon2);
				test.Print(name);
				index++;
			}

			ofile << "epsilon_eff = " << effective_property << "\n";
			e_eff.x[0] += effective_property.real();
			e_eff.x[1] += effective_property.real()*effective_property.real();

			error.x[0] += test.m;
			error.x[1] += test.Tol_F;
			error.x[2] += test.Tol_div;

			double mean = e_eff.x[0] / (double)num_configs;
			double std_err = e_eff.x[1] / (double)num_configs - mean * mean;
			std_err = (num_configs > 1) ? sqrt(abs(std_err) / (double)(num_configs - 1.0)) : 0.0L;
			epsilon_eff.emplace_back(epsilon2, mean, std_err);
			errors.emplace_back(epsilon2, error.x[0] / (double)num_configs, error.x[1] / (double)num_configs, error.x[2] / (double)num_configs);
			//sprintf(name, "Strong-contrast_%dx%d_checkerboard_%0.3f", num_domains, num_domains, epsilon2);
			//test.Print(name);
		}
		sprintf(name, "epsilon_e_%dx%d_periodic_checkerboard", num_domains, num_domains);
		WriteFunction(epsilon_eff, name);
		sprintf(name, "errors_%dx%d_periodic_checkerboard", num_domains, num_domains);
		WriteFunction(errors, name);
	}
	return 0;


	{/* random checker board. */
		int num_domains = 4, num_pixels = 200, num_configs = 100;

		std::vector<GeometryVector> epsilon_eff, errors;
		std::vector < double> list({ 0.1, 0.5, 1.1, 2, 4, 6, 12 });
		for (auto i = 0; i < list.size(); i++) {
			double epsilon2 = list[i];

			GeometryVector e_eff(0, 0, 0, 0);
			GeometryVector error(0, 0, 0, 0);
			for (int n = 0; n < num_configs; n++) {
				TwoPhaseDispersion Checkerboard;
				RandomCheckerBoard(Checkerboard, num_domains, num_pixels, rng);
				//generate a random checkboard;

				Homo_DielectricConstant_static_SCA test(Checkerboard, 1.0, epsilon2);
				test.Prepare(0);
				//test.P_direct.Print("test", 0);

				test.Run(1e-6, 1e-10);
				std::complex<double> effective_property = test.ComputeEffectiveDielectric();

				ofile << "epsilon_eff = " << effective_property << "\n";
				e_eff.x[0] += effective_property.real();
				e_eff.x[1] += effective_property.real()*effective_property.real();

				error.x[0] += test.m;
				error.x[1] += test.Tol_F;
				error.x[2] += test.Tol_div;
			}
			double mean = e_eff.x[0] / (double)num_configs;
			double std_err = e_eff.x[1] / (double)num_configs - mean * mean;
			std_err = sqrt(abs(std_err) / (double)(num_configs - 1.0));
			epsilon_eff.emplace_back(epsilon2, mean, std_err);
			errors.emplace_back(epsilon2, error.x[0] / (double)num_configs, error.x[1] / (double)num_configs, error.x[2] / (double)num_configs);
			//sprintf(name, "Strong-contrast_%dx%d_checkerboard_%0.3f", num_domains, num_domains, epsilon2);
			//test.Print(name);
		}
		sprintf(name, "epsilon_e_%dx%d_checkerboard", num_domains, num_domains);
		WriteFunction(epsilon_eff, name);
		sprintf(name, "errors_%dx%d_checkerboard", num_domains, num_domains);
		WriteFunction(errors, name);
	}
	return 0;




	for (;;) {
		double epsilon1 = 1.0L, epsilon2 = 10.0L;	//epsilon1 and epsilon2 are for the matrix and inclusion phases, respectively.
		double phi = 0.05;	//packing fraction
		ifile >> epsilon2;
		ofile << "epsilon2 = " << epsilon2 << "\n";
		SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
		c.Insert(pow(phi / HyperSphere_Volume(2, 1.0), 0.5), GeometryVector(0.5, 0.5));

		TwoPhaseDispersion packing(c, 101);	//It gives a discretized version of a sphere packing.
		packing.ComputeIndicator();
		/*	for (int i = 0; i < 10; i++) {
		std::vector<size_t> comp, n;
		comp.push_back((size_t)round(101 * rng.RandomDouble()));
		comp.push_back((size_t)round(101 * rng.RandomDouble()));

		ofile << i << ": ";
		ofile << comp[0] << ", " << comp[1] << ">>\t" << packing.Components2Index(comp) << "\n";
		packing.Index2Components(packing.Components2Index(comp), n);
		ofile << n[0] <<", "<< n[1] << "\n\n\n";
		}*/
		std::vector<size_t> v = packing.GetResolutions();


		// define the problem.
		Homo_DielectricConstant_static_SCA test(packing, epsilon1, epsilon2);

		test.Prepare(0);
		test.Print_GreenFn_2D("test");	// Print each component of the Fourier transform of the Green's function.
		ofile << "epsilon_eff = " << test.ComputeEffectiveDielectric() << "\n";
		//test.Print("FFTW_" + std::to_string(phi) + "_0__" + std::to_string(epsilon2 / epsilon1));

		//Iteratively compute the resulting electric field
		for (int i = 0; i < 100; i++) {
			test.Iterate();
			//	test.P_Fourier.Print("FFTW_" + std::to_string(phi) + "_" + std::to_string(i) + "__" + std::to_string(epsilon2 / epsilon1) + "_tildeP",0);
			//	test.E_direct.Print("FFTW_" + std::to_string(phi) + "_" + std::to_string(i + 1) + "__" + std::to_string(epsilon2 / epsilon1)+ "_E");
			test.Print("FFTW-strong_" + std::to_string(phi) + "_" + std::to_string(i + 1) + "__" + std::to_string(epsilon2 / epsilon1));
			//	test.Print_Fourier("FFTW_" + std::to_string(phi) + "_" + std::to_string(i + 1) + "__" + std::to_string(epsilon2 / epsilon1));
			//ofile << "epsilon_eff = " << test.ComputeEffectiveDielectric() << "\n";
		}
		ofile << "epsilon_eff = " << test.ComputeEffectiveDielectric() << "\n";
	}

	return 0;

	{
		double epsilon1 = 1.0L, epsilon2 = 2.0L;	//epsilon1 and epsilon2 are for the matrix and inclusion phases, respectively.
		double phi = 0.05;	//packing fraction
		ifile >> epsilon2;
		ofile << "epsilon2 = " << epsilon2 << "\n";
		SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
		c.Insert(pow(phi / HyperSphere_Volume(2, 1.0), 0.5), GeometryVector(0.5, 0.5));

		TwoPhaseDispersion packing(c, 101);	//It gives a discretized version of a sphere packing.
		packing.ComputeIndicator();
		/*	for (int i = 0; i < 10; i++) {
		std::vector<size_t> comp, n;
		comp.push_back((size_t)round(101 * rng.RandomDouble()));
		comp.push_back((size_t)round(101 * rng.RandomDouble()));

		ofile << i << ": ";
		ofile << comp[0] << ", " << comp[1] << ">>\t" << packing.Components2Index(comp) << "\n";
		packing.Index2Components(packing.Components2Index(comp), n);
		ofile << n[0] <<", "<< n[1] << "\n\n\n";
		}*/
		std::vector<size_t> v = packing.GetResolutions();


		// define the problem.
		Homo_DielectricConstant_static test(packing, epsilon1, epsilon2);

		test.Prepare(0);
		test.Print_GreenFn_2D("test");	// Print each component of the Fourier transform of the Green's function.
		ofile << "epsilon_eff = " << test.ComputeEffectiveDielectric() << "\n";
		//test.Print("FFTW_" + std::to_string(phi) + "_0__" + std::to_string(epsilon2 / epsilon1));

		//Iteratively compute the resulting electric field
		for (int i = 0; i < 100; i++) {
			test.Iterate();
			//	test.P_Fourier.Print("FFTW_" + std::to_string(phi) + "_" + std::to_string(i) + "__" + std::to_string(epsilon2 / epsilon1) + "_tildeP",0);
			//	test.E_direct.Print("FFTW_" + std::to_string(phi) + "_" + std::to_string(i + 1) + "__" + std::to_string(epsilon2 / epsilon1)+ "_E");
			//	test.Print("FFTW_" + std::to_string(phi) + "_" + std::to_string(i + 1) + "__" + std::to_string(epsilon2 / epsilon1));
			//	test.Print_Fourier("FFTW_" + std::to_string(phi) + "_" + std::to_string(i + 1) + "__" + std::to_string(epsilon2 / epsilon1));
			test.Print("FFTW-" + std::to_string(phi) + "_" + std::to_string(i + 1) + "__" + std::to_string(epsilon2 / epsilon1));

		}
		ofile << "epsilon_eff = " << test.ComputeEffectiveDielectric() << "\n";


	}
	return 0;






	{
		TensorField G(2, 2, std::vector<size_t>({ 2, 2 }));
		TensorField E(2, 1, G.num_points), P(2, 1, G.num_points);

		for (int i = 0; i < G.total_comp; i++) {
			for (int x = 0; x < G.num_total_pixels; x++) {
				G.field[i][x][0] = rng.RandomDouble();
			}
		}
		for (int i = 0; i < P.total_comp; i++) {
			for (int x = 0; x < G.num_total_pixels; x++) {
				P.field[i][x][0] = rng.RandomDouble();
				P.field[i][x][1] = rng.RandomDouble();
			}
		}

		//InnerProduct(E, G, P, 1);

		//for (int x = 0; x < E.num_total_pixels; x++) {
		//	ofile << E.field[0][x][0] << " +i " << E.field[0][x][1] << "=" << G.field[0][x][0] << "," << G.field[2][x][0] << "* " << P.field[0][x][0] << " + i " << P.field[0][x][1] << "\n";
		//	ofile << E.field[1][x][0] << " +i " << E.field[1][x][1] << "=" << G.field[1][x][0] << "," << G.field[3][x][0] << "* " << P.field[1][x][0] << " + i " << P.field[1][x][1] << "\n";
		//	ofile << "\n\n";
		//}
		E += P;
		for (int x = 0; x < E.num_total_pixels; x++) {
			for (int i = 0; i < E.total_comp; i++) {
				ofile << E.field[i][x][0] << "+i " << E.field[i][x][1] << "\t";
			}
			ofile << "\n";
		}
		ofile << "\n\n\n";
		E += P;
		for (int x = 0; x < E.num_total_pixels; x++) {
			for (int i = 0; i < E.total_comp; i++) {
				ofile << E.field[i][x][0] << "+i " << E.field[i][x][1] << "\t";
			}
			ofile << "\n";
		}
	}
	return 0;





	{
		int num_px = 100;
		double phi = 0.6;
		SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
		c.Insert(pow(phi / HyperSphere_Volume(2, 1.0), 0.5), GeometryVector(0.5, 0.5));

		TwoPhaseDispersion packing(c, num_px);
		packing.ComputeIndicator();

		TensorField x(2, 0, std::vector<size_t>({ (size_t)num_px, (size_t)num_px })), X;
		for (auto iter = packing.prt_indicator.begin(); iter != packing.prt_indicator.end(); iter++)
			x.field[0][*iter][0] = 1.0;

		x.Print("manual_x");
		//	DFT(X, x);
		X.Print("manual_XX");

		TensorField y;
		//	IDFT(y, X);
		y.Print("manual_y");


		x.Print("fftw_x");
		//DFT(X, x);
		fftw_plan a = fftw_plan_dft_2d(num_px, num_px, x.field[0], X.field[0], -1, FFTW_ESTIMATE);
		fftw_plan b = fftw_plan_dft_2d(num_px, num_px, X.field[0], y.field[0], 1, FFTW_ESTIMATE);
		fftw_execute(a);
		X.Print("fftw_XX");

		fftw_execute(b);
		y *= (1.0 / (num_px*num_px));
		//IDFT(y, X);
		y.Print("fftw_y");
	}
	return 0;





	{
		TensorField x(2, 0, std::vector<size_t>({ 10, 10 })), X(x), y(x);
		for (int i = 0; i < x.num_total_pixels; i++)
			x.field[0][i][0] = rng.RandomDouble();

		ofile << x.field[0] << "\n";
		ofile << X.field[0] << "\n";
		ofile << y.field[0] << "\n";
		x.Print("test_x");
		//DFT(X, x);
		fftw_plan a = fftw_plan_dft_2d(10, 10, x.field[0], X.field[0], -1, FFTW_ESTIMATE);
		fftw_plan b = fftw_plan_dft_2d(10, 10, X.field[0], y.field[0], 1, FFTW_ESTIMATE);
		fftw_execute(a);
		ofile << "after...\n";
		ofile << x.field[0] << "\n";
		ofile << X.field[0] << "\n";
		ofile << y.field[0] << "\n";
		X.Print("test_XX");

		fftw_execute(b);
		ofile << "after...\n";
		ofile << x.field[0] << "\n";
		ofile << X.field[0] << "\n";
		ofile << y.field[0] << "\n";
		y *= (1.0 / 100);
		//IDFT(y, X);
		y.Print("test_y");

	}
	return 0;





	{
		TensorField x(1, 0, std::vector<size_t>({ 10 })), X;
		for (int i = 0; i < 10; i++)
			x.field[0][i][0] = rng.RandomDouble();

		ofile << "x = \n";
		for (int i = 0; i < 10; i++) {
			ofile << x.field[0][i][0] << " +i " << x.field[0][i][1] << "\n";
		}

		DFT(X, x);
		ofile << "X = \n";
		for (int i = 0; i < 10; i++) {
			ofile << X.field[0][i][0] << " +i " << X.field[0][i][1] << "\n";
		}

		TensorField y;
		IDFT(y, X);
		ofile << "y = \n";
		for (int i = 0; i < 10; i++) {
			ofile << y.field[0][i][0] << " +i " << y.field[0][i][1] << "\n";
		}

	}
	return 0;

}