#include "FFH.h"
#include <GeometryVector.h>
#include <vector>
#include <StructureFactor.h>
#include <vector>
#include <functional>

size_t FFH_verbosity = 2;

/** \brief Sample driver: homogenize a simple-cubic sphere packing (static dielectric). (CLAUDE)
 *  \return 0 on success. */
int main() {
	char tempstring[1000] = {};
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[200] = {};
	char dir[] = "./results/n1/";
	omp_set_num_threads(4);
	RandomGenerator rng(0);
	fftw_init_threads();
	std::string dump;

	{       /* Simple cubic array of sphere packing (static dielectric) */
		int d = 3;
		FFH_verbosity = 2;
		double epsilon1 = 1.0L, phi = 0.50;
		double epsilon2 = 10.0L;
		int num_pixels = 101;// std::atoi(argv[1]);

		std::vector<GeometryVector> epsilon_eff, errors;

		sprintf(tempstring, "Z3-packing-%0.2f-%.2f", phi, epsilon2);
		std::string system(tempstring), method("strong");
		double a = pow(phi / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);
		SpherePacking c(GetUnitCubicBox(d, 1), 0.0);
		{
			GeometryVector pos(static_cast<int>(d));
			for (int i = 0; i < d; i++) {
				pos.x[i] = 0.5;
			}
			c.Insert(a, pos);
		}

		TwoPhaseDispersion packing(c, num_pixels);      //It gives a discretized version of a sphere packing.
		packing.ComputeIndicator();


		StaticDielectricProblem test(packing, epsilon1, epsilon2);
		test.SetMethod(true);
		ConstantTensor E0(d, 1);
		E0[0] = 1.0L;
		//	test.SetMethod(false);
		test.DefineProperties(epsilon1, epsilon2);
		test.Prepare(E0);

		sprintf(tempstring, "Dimension = %d, Packing fraction = %0.3f, Particle number = %d, \n# resolutions = %d^%d", d, phi, packing.packing.NumParticle(), num_pixels, d);
		std::string sys_descrip(tempstring);
		sprintf(tempstring, "epsilon1 = %0.3f, epsilon2 = %0.3f", epsilon1, epsilon2);
		std::string phase_properties(tempstring);
		sprintf(tempstring, "method = %s", method.c_str());
		std::string solver(tempstring);

		std::vector<std::string> label1 = { sys_descrip, phase_properties, tempstring, "The effective dielectric constant according to the tolerance in the divergence-free condition.\n",
				"tolerance", "Re[epsilon_e]", "Im[epsilon_e]" };

		std::vector<std::string> label2 = { sys_descrip, phase_properties, tempstring, "The effective dielectric constant according to the tolerance in the divergence-free condition.\n",
"tolerance", "m", "tol_field", "tol_curl" };
		int m_prev = 0;

		test.Run(1e-5, 1e-9);
		std::complex<double> epsilone = test.ComputeEffectiveEpsilon();
		ofile << "epsilon2/epsilon1 = " << epsilon2 << "\n";
		ofile << "Effective dielectric constant = " << epsilone << "\n";
		ofile << "done\n";

	}
	return 0;
}