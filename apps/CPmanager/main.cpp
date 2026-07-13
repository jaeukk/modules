/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Created	: September 2020 
 *	Modified: November 2025 */

/** \file main.cpp
 *	\brief Implementations for converting ConfigPack files to txt formats 
 *	and vice versa.	 */


#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

//-------Core part--------
//#include <etc.h>
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>





//#include <Windows.h>
#include <algorithm>

/** \brief Entry point dispatching ConfigurationPack <-> txt conversion utilities by option code. (CLAUDE)
 *  \param argc  Argument count; argv[1] selects the mode and the remaining entries are mode-specific parameters.
 *  \param argv  Argument vector: argv[1] is the option (0=CheckCP statistics, 1=CP2txt, 2=txt2CP), followed by the parameters required for that option.
 *  \return      Process exit status (0 on normal completion). */
int main(int argc, char ** argv) {

	//char tempstring[1000];
	//std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[200] = {};
	omp_set_num_threads(4);
	


	if (argc < 2) {
		ofile << "Usage: CPmanager <option: 0=CheckCP, 1=CP2txt, 2=txt2CP> [parameters]\n";
		return 1;
	}
	int option = std::atoi(argv[1]);

	switch (option){
	case 0:
	{
		/* CheckCP */
		int num_param = 1;
		if(argc == num_param+2){
			/* see statistics in the ConfigurationPack */
			std::string cp_name(argv[2]);

			ofile 	<< "===================================\n"
					<< "\t\tStatistics in " << cp_name <<"\n"
					<< "\n\n";

			long long num_config = 0;
			int dim = 0;
			double num_prts = 0, num_density = 0;
			
			ConfigurationPack cp(cp_name);
			num_config = cp.NumConfig();

			ofile << "# configs = "<< num_config <<"\n";
			if (num_config > 0){
				dim = cp.GetConfig(0).GetDimension();
				ofile << "dimension = "<< dim <<"\n";

	#pragma omp parallel for
				for(int i=0; i<num_config; i++){
					Configuration c_temp = cp.GetConfig(i);
					double n = (double)c_temp.NumParticle(), rho = n / c_temp.PeriodicVolume();

	#pragma omp atomic
					num_prts+=n;
	#pragma omp atomic
					num_density+=rho;
				}

				num_prts /= (double)num_config;
				num_density /= (double)num_config;

				ofile << "av. # particles = "<< num_prts<<"\n";
				ofile << "av. # density = "<< num_density<<"\n";
				ofile << "==================================="<<std::endl;
			}

		}
	
		else{
			ofile << num_param << " parameters are required\n";
			ofile << "e.g., [Name of a ConfigPack file]\n";
		}
	}
		break;

	case 1:
	{
		/* CP2txt */
		int num_param = 3;
		ofile	<< "===================================\n"
				<< "\tConvert ConfigPack to txt formats\n\n";
	
		if(argc == num_param+2){
			std::string cp_name(argv[2]);
			ConfigurationPack cp(cp_name);
			long long num_tot_configs = cp.NumConfig();
			int start_idx = std::atoi(argv[3]), num_configs = std::atoi(argv[4]),
				i = 0;

			ofile << "prefix of a target ConfigPack file = "<< cp_name<<"\n";
			ofile << "start index = "<< start_idx<<"\n";
			ofile << "configuration number = "<< num_configs<<"\n";

			for(i=0; (i<num_configs) && (i+start_idx < num_tot_configs); i++){
				Configuration c = cp.GetConfig(i+start_idx);
				snprintf(name, sizeof(name), "%s__%04d", cp_name.c_str(), i);
				WriteConfiguration(c, name);
			}

			ofile << i <<" (out of " << num_tot_configs <<  ") configurations are sucessfully converted\n";
			

		}
		else{
			ofile << num_param << " parameters are required\n";
			ofile << "e.g., [Name of a ConfigPack file] [start index] [number of configs]"<<std::endl;
		}
		
	
	}
	break;

	case 2:
	{
		/* txt2CP */
		int num_param = 4;
		ofile	<< "===================================\n"
				<< "\tConvert txt formats to ConfigPack\n\n";

		if(argc == num_param+2){

			int start_index = std::atoi(argv[3]), num_configs = std::atoi(argv[4]);
			std::string prefix(argv[2]);
			snprintf(name, sizeof(name), "%s", argv[5]);

			ofile << "prefix = "<<  prefix <<"\n";
			ofile << "start index = "<< start_index<<"\n";
			ofile << "configuration number = "<< num_configs<<"\n";
			ofile << "save prefix = "	<< name<<"\n";

			ConfigurationPack save(name);
			int i=0;
			long long num_conf_prev = save.NumConfig();
			for (i=0; i<num_configs; i++){
					int idx = start_index + i;
					Configuration c;
					snprintf(name, sizeof(name), prefix.c_str(), idx);
					ReadConfiguration (c, name);
					if(c.GetDimension()>0){
						save.AddConfig(c);
					}
					else
						break;
			}
			ofile << save.NumConfig() - num_conf_prev << " configurations are converted successfully\n";			
		}
	
		else{
			ofile << num_param << " parameters are required\n";
			ofile << "e.g., [Name of txt files] [start index] [number of configs] [name of the saved ConfigPack file]"<<std::endl;
		}

	}
	break;
	
	default:
	break;
	}

}
