/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: September 2022 */

/** \file main.cpp 
 * \brief A simple CLI for computing various pair statistics from point configurations. 
 * This is a simplified version of the PairStatisticsCLI code written by Ge.*/

#include <PeriodicCellList.h>
#include <etc.h>
size_t Verbosity = 2;  // each executable defines it (new_c++ convention)
// Verbosity is owned by the core module (R5).

#include "ComputationPacking.h"



class ReadTextFiles
{
public:
	std::string filename;
	int start_idx, num_configs;
	double target_phi2 = -1.0;
	char tempstring[1000];
	/** \brief Construct by reading the text-file input parameters. (CLAUDE)
	 *  \param ifile  Input stream supplying filename, start_idx and num_configs.
	 *  \param ofile  Output stream prompting for each parameter. */
	ReadTextFiles(std::istream & ifile, std::ostream & ofile){
		ofile << "Input prefix = ";
		ifile >> filename;
		ofile << "Start index = ";
		ifile >> start_idx;
		ofile << "Num configurations = ";
		ifile >> num_configs;
	}


	/** \brief Load the i-th configuration from a text file, optionally rescaling to target_phi2. (CLAUDE)
	 *  \param i  Offset added to start_idx to select the configuration file.
	 *  \return   The loaded (and possibly rescaled) SpherePacking. */
	SpherePacking operator() (size_t i)
	{
		SpherePacking result;
		snprintf(tempstring, sizeof(tempstring), this->filename.c_str(), this->start_idx + i);
		ReadConfiguration(result, tempstring);
		
		if (this->target_phi2 > 0.0){
			double curr_phi2 = result.PackingFraction();
			int d = result.GetDimension();
			double scale = std::pow(target_phi2/curr_phi2,1.0/d);
			result.Rescale_prtOnly(scale);
		}

		return result;
	}
};


class ReadConfigPack
{
public:
	double target_phi2 = -1.0;
	double Rescale;
	size_t StartIdx;
	ConfigurationPack p;
	RandomGenerator rng;
	bool sequential_sampling = true;
	/** \brief Construct by opening a ConfigurationPack from the given prefix. (CLAUDE)
	 *  \param ifile  Input stream supplying the pack prefix.
	 *  \param ofile  Output stream prompting for the prefix. */
	ReadConfigPack(std::istream & ifile, std::ostream & ofile)
	{
//		this->Rescale=Rescale;
//		this->StartIdx = StartIdx;
//		this-> target_phi2 = targetphi2;
		ofile<<"Input Prefix:";
		std::string prefix;
		ifile>>prefix;
		p.Open(prefix);
	}

	/** \brief Load the i-th configuration from the pack as a SpherePacking at target_phi2. (CLAUDE)
	 *  \param i  Offset added to StartIdx to select the configuration in the pack.
	 *  \return   The SpherePacking whose radius is set to reach target_phi2 (empty if target_phi2<=0). */
	SpherePacking operator() (size_t i)
	{
		if (this->target_phi2 <= 0.0){
			std::cout << "Using ReadConfigPack, target packing fraction must be positive!!\n";
			return SpherePacking();
		}

		Configuration c = p.GetConfig(this->StartIdx + i);
		size_t d = c.GetDimension();
		double radius = std::pow( c.PeriodicVolume() * this->target_phi2 / c.NumParticle() / ::HyperSphere_Volume(d, 1.0), 1./(double)d);
		SpherePacking result(c,radius);
		return result;
	}
};



/** \brief Run the interactive command loop driving the pair-statistics computations. (CLAUDE)
 *  \return  Status code: 0 on "Exit", 2 on end-of-input. */
int PairStatisticsCLI();

/** \brief Program entry point; dispatches to PairStatisticsCLI(). (CLAUDE)
 *  \param argc  Number of command-line arguments (unused).
 *  \param argv  Command-line argument values (unused).
 *  \return      Always 0. */
int main(int, char **){
	PairStatisticsCLI();

	return 0;
}

int PairStatisticsCLI(){
	char tempstring[1000];
	std::istream & ifile=std::cin;
	std::ostream & ofile=std::cout;

	double Rescale=1.0;
	double target_phi2 = -1.0;
	size_t NumConfig =0;
	size_t NumThreads = 1;
	size_t StartIdx = 0;
	std::string OutputPrefix;
	std::function<const SpherePacking(size_t i)> GetConfigsFunction = nullptr;	

	std::vector<std::string> available_computations{
		"chivComputation","LocalVarianceComputation",
		"HoleProbability",
		"LinealPath", "ChordLength"
	};
	std::vector<std::string> todo_computations{
		"S2Computation"
	};

	std::vector<ComputationPacking *> vpComputations;
	for(;;)
	{
		ifile>>tempstring;
		/* Basic parameters */
		if(strcmp(tempstring, "Exit")==0)
		{
			for(auto iter=vpComputations.begin(); iter!=vpComputations.end(); iter++)
				delete *iter;
			return 0;
		}
		else if(strcmp(tempstring, "Available")==0){
			for(auto iter=available_computations.begin(); iter!=available_computations.end(); iter++)
				ofile << *iter <<";\n";
		}
		else if(strcmp(tempstring, "ToDo")==0){
			for(auto iter=todo_computations.begin(); iter!=todo_computations.end(); iter++)
				ofile << *iter <<";\n";
		}
		else if(strcmp(tempstring, "OutputPrefix")==0)
		{
			ifile>>OutputPrefix;
		}
		else if(strcmp(tempstring, "Rescale")==0)
		{
			std::cerr << "Warning : Rescale is deprecated\n";
			ifile>>Rescale;
		}
		else if(strcmp(tempstring, "phi2")==0)
		{
			ifile >> target_phi2;
			ofile << "This should be called before \"GetConfigsFunction\"\n";
		}
		else if(strcmp(tempstring, "NumConfig")==0)
		{
			ifile>>NumConfig;
		}
		else if(strcmp(tempstring, "NumThreads")==0)
		{
			ifile>>NumThreads;
		}

		/* define computations */
		else if(strcmp(tempstring, "chivComputation")==0)
		{
			vpComputations.push_back( new chivComputation(ifile, ofile) );
		}
		else if(strcmp(tempstring, "LocalVarianceComputation")==0)
		{
			//TODO vpComputations.push_back( new LocalVolumeFractionVariance(ifile, ofile) );
		}
		else if(strcmp(tempstring, "HoleProbability")==0)
		{
			//TODO vpComputations.push_back( new LocalVolumeFractionVariance(ifile, ofile) );
		}
		else if(strcmp(tempstring, "LinealPath")==0)
		{
			vpComputations.push_back( new 
			LinealPath(ifile, ofile) );
		}
		else if(strcmp(tempstring, "ChordLength")==0)
		{
			vpComputations.push_back( new 
			ChordLength(ifile, ofile) );
		}


		else if (strcmp(tempstring, "GetConfigsFunction") == 0)
		{
			ofile<<"Input Type:";
			ifile>>tempstring;
			if(strcmp(tempstring, "TextFiles")==0)
			{
				ReadTextFiles c(ifile, ofile);
				c.target_phi2 = target_phi2;
				GetConfigsFunction = c;
				NumConfig = c.num_configs;
			}
			else if (strcmp(tempstring, "ReadConfigPack")==0){
				ReadConfigPack c(ifile, ofile);
				c.target_phi2 = target_phi2;
				c.StartIdx = StartIdx;
				std::cout<<"The config pack contains "<<c.p.NumConfig()<<" configurations. ";
				if (c.p.NumConfig()-StartIdx < NumConfig){
					NumConfig=c.p.NumConfig()-StartIdx;
					std::cout << "Because this value exceeds NumConfig, re-set NumConfig to this value\n";
				}
				std::cout << "\n";
				GetConfigsFunction = c;
			}
			else
			{
				std::cout<<"Unrecognized Type!\n";
				std::cin.clear();
			}
		}
		else if (strcmp(tempstring, "StartIndex") == 0)
		{
			ifile>> StartIdx;
			ofile<< "Start Index = "<< StartIdx <<std::endl;
		}

		else if(strcmp(tempstring, "Calculation")==0)
		{
			if(GetConfigsFunction == nullptr)
				std::cerr<<"Specify GetConfigsFunction before compute!\n";
			else
			{
				for(auto iter=vpComputations.begin(); iter!=vpComputations.end(); iter++)
				{
					(*iter)->SetNumThreads(NumThreads);
					(*iter)->Compute(GetConfigsFunction, NumConfig);
					(*iter)->Write(OutputPrefix);
				}
			}
		}

/* 		else if (strcmp(tempstring, "NearestNeighborComputation") == 0)
		{//for historical reasons
			vpComputations.push_back( new HpComputation(ifile, ofile) );
		}
		else if(strcmp(tempstring, "HpComputation")==0)
		{
			vpComputations.push_back( new HpComputation(ifile, ofile) );
		}
		else if(strcmp(tempstring, "HvComputation")==0)
		{
			vpComputations.push_back( new HvComputation(ifile, ofile) );
		}
		else if (strcmp(tempstring, "VoronoiVolumeComputation") == 0)
		{
			vpComputations.push_back(new VoronoiVolumeComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "VoronoiNumSidesComputation") == 0)
		{
			vpComputations.push_back(new VoronoiNumSidesComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "CoveringRadiusDistributionComputation") == 0)
		{
			vpComputations.push_back(new CoveringRadiusDistributionComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "PackingRadiusDistributionComputation") == 0)
		{
			vpComputations.push_back(new PackingRadiusDistributionComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "AverageClusterSizeComputation") == 0)
		{
			vpComputations.push_back(new AverageClusterSizeComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "VoronoiVolumeCorrelationComputation") == 0)
		{
			vpComputations.push_back(new VoronoiVolumeCorrelationComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "Psi6CorrelationComputation") == 0)
		{
			vpComputations.push_back(new Psi6CorrelationComputation(ifile, ofile) );
		}
		else if (strcmp(tempstring, "LocalQ6DistributionComputation") == 0)
		{
			vpComputations.push_back(new LocalQ6DistributionComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "WindowNumberVarianceComputation") == 0)
		{
			vpComputations.push_back(new WindowNumberVarianceComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "AnalyticWindowNumberVarianceComputation") == 0)
		{
			vpComputations.push_back(new AnalyticWindowNumberVarianceComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "WindowNumberDistributionComputation") == 0)
		{
			vpComputations.push_back( new WindowNumberDistributionComputation(ifile, ofile) );
		}
		else if(strcmp(tempstring, "NeighborLinkTortuosityComputation")==0)
		{
			vpComputations.push_back( new NeighborLinkTortuosityComputation(ifile, ofile) );
		}
		else if (strcmp(tempstring, "VarianceOverRd_1Computation") == 0)
		{
			vpComputations.push_back(new VarianceOverRd_1Computation(ifile, ofile));
		}
		else if (strcmp(tempstring, "PercolationVolumeFractionDistribution") == 0)
		{
			vpComputations.push_back(new PercolationVolumeFractionDistribution(ifile, ofile));
		}
		else if (strcmp(tempstring, "PercolationP1") == 0)
		{
			vpComputations.push_back(new PercolationP1Calculation_v2(ifile, ofile));
		}
		else if (strcmp(tempstring, "PercolationP1Computation") == 0)
		{
			vpComputations.push_back(new PercolationP1Calculation(ifile, ofile));
		}
		else if (strcmp(tempstring, "PercolationP1Calculation") == 0)
		{
			vpComputations.push_back(new PercolationP1Calculation(ifile, ofile));
		}
		else if (strcmp(tempstring, "M2Computation") == 0)
		{
			vpComputations.push_back(new M2Computation(ifile, ofile));
		}
		else if (strcmp(tempstring, "DiscretizationVolumeFractionComputation") == 0)
		{
			vpComputations.push_back(new DiscretizationVolumeFractionComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "FirstPassageTimeComputation") == 0)
		{
			vpComputations.push_back(new FirstPassageTimeComputation(ifile, ofile));
		} */

		else
		{
			if (vpComputations.empty())
				std::cout << "Unrecognized command!\n";
			else
				vpComputations.back()->ProcessAdditionalOption(tempstring, std::cin, std::cout);

			std::cin.clear();
		}
		ifile.ignore(1000,'\n');
		if(ifile.eof()) return 2;
	}
	return 1;
}

