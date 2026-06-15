#include "FFH.h"
#include <GeometryVector.h>
#include <vector>
#include <StructureFactor.h>
#include <vector>
#include <functional>

size_t FFH_verbosity = 2;

void PeriodicCheckerBoard(TwoPhaseDispersion & out, int num_domains, int num_pixels)
{
	SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
	out = TwoPhaseDispersion(c, num_pixels);
	if ((num_domains % 2) == 0)
	{
		std::vector<std::vector<size_t>> domains(num_domains*num_domains);
		int index = 0;
		for (int i = 0; i < num_domains; i++) {
			for (int j = 0; j < num_domains; j++) {
				domains[index].push_back(i);
				domains[index].push_back(j);
				index++;
			}
		}

		for (int i = 0; i < domains.size(); i++) {
			if ((domains[i][0] + domains[i][1]) % 2 == 0) {
				domains.erase(domains.begin() + i);
				i--;
			}
		}


		std::vector<GeometryVector> basis;
		basis.emplace_back(out.packing.GetBasisVector(0));
		basis.emplace_back(out.packing.GetBasisVector(1));
		basis[0] = 1.0 / (double)num_pixels * basis[0];
		basis[1] = 1.0 / (double)num_pixels * basis[1];

		Configuration cell(2, &basis[0], 1);
		cell.Insert("a", GeometryVector(0.5, 0.5));

		Configuration grid = MultiplicateStructure(cell, std::vector<size_t>({ (size_t)num_pixels, (size_t)num_pixels }));
		out.prt_indicator.resize(grid.NumParticle(), 0);

		index = 0;
		for (int i = 0; i < grid.NumParticle(); i++) {
			GeometryVector pos = grid.GetRelativeCoordinates(i);
			int x = (int)floor(num_domains*pos.x[0]);
			int y = (int)floor(num_domains*pos.x[1]);
			bool inside = false;
			for (int j = 0; j < domains.size(); j++) {
				inside = inside || (domains[j][0] == x && domains[j][1] == y);
			}
			if (inside) {
				std::vector<size_t> comp({ (size_t)floor(pos.x[0] * num_pixels), (size_t)floor(pos.x[1] * num_pixels) });
				out.prt_indicator[index] = out.Components2Index(comp);
				index++;
			}
		}
		out.prt_indicator.resize(index);
		out.update_indicator = true;
	}

	else {
		std::cout << "num_domains should be an even number\n";
	}
}

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

	{       /* Periodic Checkerboard array */
		int d = 2;
		FFH_verbosity = 2;
		double epsilon1 = 1.0L, phi = 0.5;
		double epsilon2 = 10.0L;
		int num_pixels = 1001;
		std::vector<GeometryVector> epsilon_eff, errors;


		SpherePacking c(GetUnitCubicBox(d, 1), 0.0);
		TwoPhaseDispersion packing(c, num_pixels);      //It gives a discretized version of a sphere packing.
		PeriodicCheckerBoard(packing, 4, num_pixels);


		StaticDielectricProblem test(packing, epsilon1, epsilon2);
		test.SetMethod(true);
		for (epsilon2 = 10.0L; epsilon2 <= 100L; epsilon2 += 10) {
			//ifile >> epsilon2;
			ofile << "epsilon2 = " << epsilon2 << "\n";


			// define the problem.
			ConstantTensor E0(d, 1);
			E0[0] = 1.0L;
			test.DefineProperties(epsilon1, epsilon2);
			test.Prepare(E0);

			// Iteratively compute the resulting electric field
			test.Run(1e-4, 1e-10);
			std::complex<double> epsilone = test.ComputeEffectiveEpsilon();
			ofile << "epsilon2/epsilon1 = " << epsilon2 << "\n";
			ofile << "Effective dielectric constant = " << epsilone << "\n";
			ofile << "next...\n";
			epsilon_eff.emplace_back(epsilon2, epsilone.real(), epsilone.imag());
			errors.emplace_back(epsilon2, test.Tol_Field, test.Tol_div, test.m);

			sprintf(name, "strong_%.2f", epsilon2);
			test.Print(name);
			test.Print_Fourier(name);
		}
		sprintf(name, "PeriodicCheckerBoard_epsilone_strong");
		WriteFunction(epsilon_eff, name);
		sprintf(name, "PeriodicCheckerBoard_errors_strong");
		WriteFunction(errors, name);
	}
	return 0;
}