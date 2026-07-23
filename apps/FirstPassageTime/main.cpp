#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>



#include <etc.h>
size_t Verbosity = 2;  // each executable defines it (new_c++ convention)
#include <GeometryVector.h>
#include <RandomGenerator.h>

#include <ParticleShapes.h>

#include <PeriodicCellList.h>
#include <Dispersions.h>

#include "FirstPassageTime.h"
#include <chrono>

/** \brief Write a list of GeometryVectors as columns to an output stream, with header fields. (CLAUDE)
 *  \param result  Data rows to write.
 *  \param ofile   Output stream to write to.
 *  \param fields  Header/comment lines emitted before the data. */
void WriteFunction(const std::vector<GeometryVector> & result, std::ostream & ofile, const std::vector<std::string> & fields);

/** \brief Write a list of GeometryVectors to a file given by prefix, with header fields. (CLAUDE)
 *  \param result         Data rows to write.
 *  \param OutFilePrefix  Output file path prefix.
 *  \param fields         Header/comment lines emitted before the data. */
void WriteFunction(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::vector<std::string> & fields);

/** \brief Read a list of GeometryVectors from an input stream, taking a fixed number of columns. (CLAUDE)
 *  \param result      Output vector populated with the parsed rows.
 *  \param ifile       Input stream to read from.
 *  \param NumColumns  Number of columns to read per row. */
void ReadFunction(std::vector<GeometryVector> & result, std::istream & ifile, size_t NumColumns);

/** \brief Read a list of GeometryVectors from a file given by prefix, taking a fixed number of columns. (CLAUDE)
 *  \param result       Output vector populated with the parsed rows.
 *  \param InFilePrefix  Input file path prefix.
 *  \param NumColumns   Number of columns to read per row. */
void ReadFunction(std::vector<GeometryVector> & result, const std::string & InFilePrefix, size_t NumColumns);

/** \brief Run the first-passage-time time-evolution simulation on a rescaled honeycomb network. (CLAUDE)
 *  \param time_max       Total simulation time.
 *  \param time_interval  Time increment between recorded slices.
 *  \param phi2           Target area fraction of phase 2 to rescale to.
 *  \param sigma1         Conductivity of the host phase.
 *  \param sigma2         Conductivity of the particle phase. */
void HoneycombNetwork(double, double, double, double, double);

/** \brief Run the first-passage-time time-evolution simulation on a checkerboard configuration. (CLAUDE)
 *  \param time_max       Total simulation time.
 *  \param time_interval  Time increment between recorded slices.
 *  \param sigma1         Conductivity of the host phase.
 *  \param sigma2         Conductivity of the particle phase. */
void Checkerboard(double, double, double, double);


/** \brief Program entry point; dispatches to a simulation mode selected by the first command-line argument. (CLAUDE)
 *  \param argc  Number of command-line arguments.
 *  \param argv  Command-line argument values; argv[1] selects the mode ("0"-"4"), remaining arguments are mode-specific.
 *  \return      Exit status (0 on success). */
int main(int argc, char ** argv) {
	std::string tempStr;
	char tempstring[1000];
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	RandomGenerator rng(0);
	int num_threads = 6;
	omp_set_num_threads(num_threads);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

//	Trajectory test;
//	ofile << test;
//	double x = std::atof(argv[1]);
//	double y = std::atof(argv[2]);
//	Test_simulation2();
	if (argc < 2) {
		ofile << "Usage: FirstPassageTime.out <mode 0-4> [args...]\n";
		return 1;
	}
	if(strcmp(argv[1],"0")==0){
		if (argc < 6) { ofile << "Usage: mode 0 requires 4 arguments\n"; return 1; }
		Checkerboard(std::stof(argv[2]),std::stof(argv[3]),std::stof(argv[4]),std::stof(argv[5]));
	}
	else if(strcmp(argv[1],"1")==0){
		if (argc < 7) { ofile << "Usage: mode 1 requires 5 arguments\n"; return 1; }
		HoneycombNetwork(std::stof(argv[2]),std::stof(argv[3]),std::stof(argv[4]),std::stof(argv[5]), std::stof(argv[6]));
	}
	else if(strcmp(argv[1],"2")==0){
		if (argc < 6) { ofile << "Usage: mode 2 requires <load> <start> <num> <save> [log]\n"; return 1; }
		for (int i=0; i<argc; i++)
			ofile << argv[i] <<std::endl;

		std::string loadname (argv[2]);
		size_t start_idx = std::atoi(argv[3]);
		size_t num_configs = std::atoi(argv[4]);
		std::string savename (argv[5]);
		std::string logname("");
		if (argc == 7){
			logname = std::string(argv[6]);
		}
		
		std::vector<std::string> header;
		header.emplace_back("First-passage-time calculation for 2D packings \n");

		snprintf(tempstring, sizeof(tempstring), "Load %s\n", loadname.c_str());
		header.emplace_back(tempstring);
		ofile << tempstring;

		snprintf(tempstring, sizeof(tempstring), "\tfrom %zu to %zu\n", start_idx, start_idx + num_configs);
		header.emplace_back(tempstring);
		ofile << tempstring;

		snprintf(tempstring, sizeof(tempstring), "save to %s\n", savename.c_str());
		header.emplace_back(tempstring);
		ofile << tempstring;

		if (argc == 7){
			snprintf(tempstring, sizeof(tempstring), "save log files as %s\n", logname.c_str());
			header.emplace_back(tempstring);
			ofile << tempstring;
		}

		double sigma1, sigma2;
		double time_max;
		size_t num_test_prts_per_config, num_threads;


		ofile << "sigma in host phase = ";
		ifile >> sigma1;
		ofile << sigma1 <<"\n";
		ofile << "sigma in particle phase = ";
		ifile >> sigma2;
		ofile << sigma2 << "\n";
		ofile << "target simulation time = ";
		ifile >> time_max;
		ofile << time_max <<"\n";
		ofile << "num. test particles / config = ";
		ifile >> num_test_prts_per_config;
		ofile << num_test_prts_per_config <<"\n";
		ofile << "num. threads = ";
		ifile >> num_threads;
		ofile << num_threads <<"\n";


		auto getc = [&loadname, &start_idx](size_t i)->Dispersion{
			Dispersion c;
			char name[400] = {};
			// loadname is a user printf pattern; it must use an int conversion (e.g. %d/%04d).
			snprintf(name, sizeof(name), loadname.c_str(), static_cast<int>(start_idx + i));
			ReadDispersion(c, name);
			return c;
		};

		FirstPassageTimeParameters p(sigma1, sigma2, getc(0).GetDimension());
		p.num_test_prts = num_test_prts_per_config;
		p.num_threads = num_threads;

		std::vector<GeometryVector> result;
		p.PerformSimulation2D_parallel(getc, num_configs, time_max, result, logname);

		
		header.emplace_back("simga1 = " + std::to_string(sigma1) + "\n");
		header.emplace_back("simga2 = " + std::to_string(sigma2) + "\n");
		header.emplace_back("num. random walkers / configs= " + std::to_string(p.num_test_prts) + "\n");
		header.emplace_back("<time>\t<X2>\terr. <time>\terr. <X2>");
		WriteFunction(result, savename, header);
	}
	else if(strcmp(argv[1],"3")==0){
		if (argc < 4) { ofile << "Usage: mode 3 requires <x> <y>\n"; return 1; }
		int num = 500;
		double x = atof(argv[2]);
		double y = atof(argv[3]);
		SpherePacking c;
		ReadConfiguration(c, "./config");
		ofile <<"read\n";
		double count_total_trials = (double) num, count2host = 0.0;
		double sigma1 = 1.0, sigma2 = 0.0;
		ofile << std::setprecision(8);
		ofile << std::scientific;
		// for (double X = 0.3-0.0001; X < 0.6 ; X += 0.05){
		// 	GeometryVector v (X, 0.0), x0_Relative, n2host;
		// 	std::vector<GeometryVector> list;
		// 	c.FindNearestNeighbors(v, list, x0_Relative, n2host);
		
		// 	ofile << "v = \n";
		// 	ofile << v;
		// 	ofile << "list = \n";
		// 	for (auto u = list.begin(); u!=list.end(); u++){
		// 		ofile << * u;
		// 	}
		// 	ofile << "x0 = \n";
		// 	ofile << x0_Relative << "\n";
		// 	ofile << "n2 host = \n";
		// 	ofile << n2host;

		// 	ofile << "\n\n\n";
		// }

		
		FirstPassageTimeParameters param(sigma1, sigma2, 2);
		Trajectory one;
		RandomGenerator rng(0);
		one.sigma_host = sigma1;	one.sigma_particle = sigma2;
		one.dimension = 2;
		one.delta_pos = GeometryVector(0.0,0.0);
		one.pos_curr = GeometryVector(x,y);
		std::cout << "initial position = " << one.pos_curr;
		std::vector<GeometryVector> debug;
		for(size_t i=0; i<1; i++){
			param.PerformSimulation2D_test(c, num, one, rng);
			//if(c.IsInsideParticle(one.pos_curr) == -1){
			//	count2host ++;
			//}

			//one.pos_curr = GeometryVector(x,y);
			//std::cout<<one.pos_curr;
			//param.JumpOneStep_2D(c, one, rng, debug, true);

		}
		ofile << one;
		
		std::cout << one.GetTimeTotal() << "\t" << one.GetX2()<<"\n";
		std::cout<<"p_1 = "<< count2host / count_total_trials <<"\n";
	}
	else if(strcmp(argv[1],"4")==0){
		if (argc < 8) { ofile << "Usage: mode 4 requires <geom> <phi2> <sigma_host> <sigma_particle> <time_max> <time_interval> [log]\n"; return 1; }
		for (int i=0; i<argc; i++)
			ofile << argv[i] <<std::endl;

		SpherePacking c;
		if(strcmp(argv[2],"square")==0){
			ReadConfiguration(c, "./square");
		}
		else if(strcmp(argv[2],"triangular")==0){
			ReadConfiguration(c, "./triangular");
		}

		std::string name(argv[2]);
		double phi2 = atof(argv[3]);
		double sigma1 = atof(argv[4]);	// sigma_host
		double sigma2 = atof(argv[5]);	// sigma_particle
		double time_max = atof(argv[6]);
		double time_interval = atof(argv[7]);
		std::string logname("");
		if (argc == 9){
			logname=std::string(argv[8]);
		}

		double phi_curr = c.PackingFraction();
		double factor = sqrt(phi2/phi_curr);
		c.Rescale_prtOnly(factor);
		std::cout << "particles are rescaled: packing fraction = "<< c.PackingFraction()<<"\n";
		
		FirstPassageTimeParameters param(sigma1, sigma2, 2);
		param.num_test_prts = 10; //00;
		param.num_threads = 4;

		std::vector<GeometryVector> result;
		param.PerformSimulation2D_TimeEvolution(c, time_max, time_interval, result, logname);

		char filename[400]={};
		std::vector<std::string> header;

		if(strcmp(argv[2],"square")==0){
			snprintf(filename, sizeof(filename), "/tigress/jaeukk/FPT/square/square-%0.2f-%0.1f-%0.1f", phi2, sigma1, sigma2);
			header.emplace_back("First-passage-time calculation for 2D periodic square lattice packing of identical circles\n");
		}
		else if(strcmp(argv[2],"triangular")==0){
			snprintf(filename, sizeof(filename), "/tigress/jaeukk/FPT/triangular/triangular-%0.2f-%0.1f-%0.1f", phi2, sigma1, sigma2);
			header.emplace_back("First-passage-time calculation for 2D periodic triangular lattice packing of identical circles\n");
		}
				
		header.emplace_back("area fraction of phase 2 = " + std::to_string(phi2) + "\n");
		header.emplace_back("simga1 = " + std::to_string(sigma1) + "\n");
		header.emplace_back("simga2 = " + std::to_string(sigma2) + "\n");
		header.emplace_back("num. random walkers = " + std::to_string(param.num_test_prts) + "\n");
		header.emplace_back("mean time\tmean X2\terr. time\terr. X2");

		WriteFunction(result, filename, header);
	}

//	Test_simulation(std::atof(argv[1]),std::atof(argv[2]), std::atoi(argv[3]));

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	ofile << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

//	RandomGenerator RNG;
//	Test_simulation(RNG.RandomDouble(),RNG.RandomDouble(), num_steps);

//	Test_interface(std::atof(argv[1]),std::atof(argv[2]), num_steps);
//	Test_simulation3(std::atof(argv[1]),std::atof(argv[2]), num_steps);

	//Test_p1();
	// std::vector<GeometryVector> test;
	// for (int i=0; i<1000; i++){
	// 	test.push_back(RandomUnitVector(2, rng));
	// }
	// WriteFunction(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/test");

	return 0;
}

void Checkerboard(double time_max, double time_interval, double sigma1, double sigma2){
	
	Dispersion c;
	ReadDispersion(c, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/checkerboard/regular2");
	
	FirstPassageTimeParameters param(sigma1, sigma2, 2);
	param.num_test_prts = 100; //00;
	param.num_threads = 1;

	char name[400]={};
	snprintf(name, sizeof(name), "./checkerboard-%0.2f-%0.2f", sigma1, sigma2);
	std::vector<GeometryVector> result;
	param.PerformSimulation2D_TimeEvolution(c, time_max, time_interval, result, name);

	snprintf(name, sizeof(name), "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/checkerboard-%0.1f-%0.1f", sigma1, sigma2);
	std::vector<std::string> header;
	header.emplace_back("First-passage-time calculation for 2D periodic checkerboard\n");
	header.emplace_back("simga1 = " + std::to_string(sigma1) + "\n");
	header.emplace_back("simga2 = " + std::to_string(sigma2) + "\n");
	header.emplace_back("num. random walkers = " + std::to_string(param.num_test_prts) + "\n");
	header.emplace_back("mean time\tmean X2\terr. time\terr. X2");

	WriteFunction(result, name, header);
}

void HoneycombNetwork(double time_max, double time_interval, double phi2, double sigma1, double sigma2){
	
	Dispersion c;
	ReadDispersion(c, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/honeycomb-0.50");
	double phi_curr = c.PackingFraction();
	double factor = sqrt(phi2/phi_curr);
	c.Rescale_prtOnly(factor);
	std::cout << "particles are rescaled: packing fraction = "<< c.PackingFraction()<<"\n";
	
	FirstPassageTimeParameters param(sigma1, sigma2, 2);
	param.num_test_prts = 1000000; //00;
	param.num_threads = 20;

	std::vector<GeometryVector> result;
	param.PerformSimulation2D_TimeEvolution(c, time_max, time_interval, result);

	char name[400]={};
	snprintf(name, sizeof(name), "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/honeycomb-%0.2f-%0.1f-%0.1f", phi2, sigma1, sigma2);
	std::vector<std::string> header;
	header.emplace_back("First-passage-time calculation for 2D periodic honeycomb network\n");
	header.emplace_back("area fraction of phase 2 = " + std::to_string(phi2) + "\n");
	header.emplace_back("simga1 = " + std::to_string(sigma1) + "\n");
	header.emplace_back("simga2 = " + std::to_string(sigma2) + "\n");
	header.emplace_back("num. random walkers = " + std::to_string(param.num_test_prts) + "\n");
	header.emplace_back("mean time\tmean X2\terr. time\terr. X2");

	WriteFunction(result, name, header);
}


