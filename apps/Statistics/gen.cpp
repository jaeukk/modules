/**
 * @file gen.cpp
 * @author Jaeuk Kim (phy000.kim@gmail.com)
 * @brief Generate ConfigPacks of 1D/2D/3D perturbed lattices ("PL") or Poisson
 *        patterns ("Poisson"); parameters are read from stdin.
 * @version 0.1
 * @date 2023-01-09
 *
 * PL reads the GenPerturb block (dimension, distribution, parameters), then the
 * side length L (cells per side), the number of configurations, and the output
 * prefix.  The source lattice is integer/square/cubic with unit spacing, and the
 * random seed is fixed to 0.  With "uniform <R>" each point moves by a vector
 * uniform in a ball of radius R; in 1D this is the uniformly randomized lattice
 * (URL) with a = 2R, which is cloaked (no Bragg peaks) for integer a.
 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

#include <etc.h>
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>

#include <GenerateConfigs.h>


/** @brief CLI entry point.
 *  @param argc  Argument count (must be 2).
 *  @param argv  argv[1] = "PL" or "Poisson".
 *  @return 0. */
int main(int argc, char ** argv){
	
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	RandomGenerator rng(0);
	
	std::string savename;
	
	
	if (argc == 2){
		if (strcmp(argv[1], "PL") == 0){
			/* Perturbed lattice with uniformly distributed displacements. */

			GenPerturb sys(ifile, ofile);
			sys.SetSeed(0);
			size_t L;
			size_t Nc;
			ofile << "side length = ";
			ifile >> L;
			ofile << "num configs = ";
			ifile >> Nc;
			std::string savename;
			ofile << "saveas = ";
			ifile >> savename;

			std::string lattice;
			if (sys.GetDimension() == 1){
				lattice = std::string("integer");
			}
			else if (sys.GetDimension() == 2){
				lattice = std::string("square");
			}
			else if (sys.GetDimension() == 3){
				lattice = std::string("cubic");
			}
				
			Configuration Z0;
			{
				LatticeGen Zn(lattice);
				Zn.SetNumCells(L);
				Zn.GenerateC(Z0);
			}

			ofile << "Generating Perturbed lattices\n";
			ConfigurationPack CP(savename);
			progress_display pd(Nc, ofile);
			for(size_t i=0; i<Nc; i++){
				Configuration c(Z0);
				sys.GenerateC(c);
				CP.AddConfig(c);
				pd ++ ;
			}
			ofile << "Done\n";

		}
		else if (strcmp(argv[1], "Poisson") == 0){
			size_t vol, Nc;
			DimensionType d;
			std::string savename;
			PoissonPtGen sys(ifile, ofile);
			ofile << "dimension = ";
			ifile >> d;
			ofile << "volume = ";
			ifile >> vol;
			ofile << "num configs = ";
			ifile >> Nc;
			ofile << "saveas = ";
			ifile >> savename;

			Configuration C0 = GetUnitCubicBox(d);
			C0.Resize(vol);
			ofile << "Generating spatially uncorrelated point patterns\n";
			ConfigurationPack CP(savename);
			progress_display pd(Nc, ofile);
			for(size_t i=0; i<Nc; i++){
				Configuration c(C0);
				sys.GenerateC(c);
				CP.AddConfig(c);
				pd ++ ;
			}
			ofile << "Done\n";

			ofile << "\n\nnum. particles\n";
			for (long long i = 0; i<CP.NumConfig(); i++){
				Configuration cc = CP.GetConfig(i);
				ofile << i << ": " << cc.NumParticle() <<"\n";  
			}
		}

		
	}
	else{
		ofile << "It needs a command line argument.\n";
	}

	return 0;
}
