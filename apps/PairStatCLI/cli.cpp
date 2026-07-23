/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: March 2022 */

/**
 * \file main.cpp
 * \brief A simple CLI for computing various pair statistics from point configurations.
 * This is a simplified version of the PairStatisticsCLI code written by Ge.
 */

#include <memory>

#include "Computation.h"
//#include "ContactStatistics.h"

/* header fields in $(cores) */
#include <PeriodicCellList.h>
#include <etc.h>
size_t Verbosity = 2;  // each executable defines it (new_c++ convention)

// Verbosity is defined by the core module (modules/core/etc.cpp).

//helper type 1 : callable types to get configurations in various file formats, storing configurations in individual txt files.
class PrefixNumberSuffix
{
public:
	std::string path;
	std::string PostName;
	double Rescale;
	Configuration (*GetFunction)(std::istream & ifile);
	
	/**
	 * \brief Construct the loader, reading the file name prefix, suffix, and file type from the input stream.
	 * \param ifile    Input stream supplying the prefix, suffix, and file type.
	 * \param ofile    Output stream for prompts.
	 * \param Rescale  Linear scaling factor applied to each configuration when read.
	 */
	PrefixNumberSuffix(std::istream & ifile, std::ostream & ofile, double Rescale)
	{
		char tempstring[1000];
		this->Rescale=Rescale;
		ofile<<"Input Prefix:";
		ifile>>this->path;
		ofile<<"Input Suffix:";
		ifile>>this->PostName;

inputtype:
		ofile<<"Input File type:";
		ifile>>tempstring;
		if(strcmp(tempstring, "Pos")==0)
			this->GetFunction=ReadPos;
		else if(strcmp(tempstring, "Xml")==0)
			this->GetFunction=ReadHoomdXml;
		else if(strcmp(tempstring, "Stealth")==0)
			this->GetFunction=ReadStealthOutput;
		else if(strcmp(tempstring, "Coordinate")==0)
			this->GetFunction=ReadCoordinate;
		else
		{
			std::cout<<"Unrecognized file type!\n";
			goto inputtype;
		}
	}

	/**
	 * \brief Read the configuration whose file name is prefix + i + suffix, then rescale it.
	 * \param i  Index inserted between prefix and suffix to form the file name.
	 * \return   The loaded and rescaled configuration.
	 */
	Configuration operator() (size_t i)const
	{
		std::stringstream name;
		name<<path;
		name<<i;
		name<<PostName;
		std::string filename;
		name>>filename;
		std::fstream ifile(filename.c_str(), std::fstream::in);
		Configuration result=GetFunction(ifile);
		result.Rescale(Rescale);
		return result;
	}
};

//helper type 2 : callable types to get configurations
class PrefixNumberSuffix_withInitialNumber
{
public:
	std::string path;
	std::string PostName;
	double Rescale;
	signed long InitialNumber;
	Configuration (*GetFunction)(std::istream & ifile);

	/**
	 * \brief Construct the loader, reading the prefix, suffix, file type, and initial index offset from the input stream.
	 * \param ifile    Input stream supplying the prefix, suffix, file type, and initial number.
	 * \param ofile    Output stream for prompts.
	 * \param Rescale  Linear scaling factor applied to each configuration when read.
	 */
	PrefixNumberSuffix_withInitialNumber(std::istream & ifile, std::ostream & ofile, double Rescale)
	{
		char tempstring[1000];
		this->Rescale=Rescale;
		ofile<<"Input Prefix:";
		ifile>>this->path;
		ofile<<"Input Suffix:";
		ifile>>this->PostName;

inputtype:
		ofile<<"Input File type:";
		ifile>>tempstring;
		if(strcmp(tempstring, "Pos")==0)
			this->GetFunction=ReadPos;
		else if(strcmp(tempstring, "Xml")==0)
			this->GetFunction=ReadHoomdXml;
		else if(strcmp(tempstring, "Stealth")==0)
			this->GetFunction=ReadStealthOutput;
		else if(strcmp(tempstring, "Coordinate")==0)
			this->GetFunction=ReadCoordinate;
		else
		{
			std::cout<<"Unrecognized file type!\n";
			goto inputtype;
		}
		ofile<<"Input Initial Number:";
		ifile>>this->InitialNumber;
	}

	/**
	 * \brief Read the configuration whose file name is prefix + (i + InitialNumber) + suffix, then rescale it.
	 * \param i  Index offset, added to InitialNumber to form the file name.
	 * \return   The loaded and rescaled configuration.
	 */
	Configuration operator() (size_t i)const
	{
		std::stringstream name;
		name<<path;
		name<<i+InitialNumber;
		name<<PostName;
		std::string filename;
		name>>filename;
		std::fstream ifile(filename.c_str(), std::fstream::in);
		Configuration result=GetFunction(ifile);
		if(result.NumParticle()<1)
		{
			std::cerr<<"Error reading file:"<<filename<<", abort!\n";
			exit(0);
		}
		result.Rescale(Rescale);
		return result;
	}
};

//helper type 3 : callable types to get configurations from a ConfigPack file.
class ReadConfigPack
{
public:
	double Rescale;
	size_t StartIdx;
	ConfigurationPack p;
	RandomGenerator rng;
	bool sequential_sampling = true;

	/**
	 * \brief Construct the loader, reading the configuration-pack prefix and opening the pack.
	 * \param ifile     Input stream supplying the pack prefix.
	 * \param ofile     Output stream for prompts.
	 * \param Rescale   Linear scaling factor applied to each configuration when read.
	 * \param StartIdx  Index of the first configuration to read from the pack.
	 */
	ReadConfigPack(std::istream & ifile, std::ostream & ofile, double Rescale, size_t StartIdx = 0)
	{
		this->Rescale=Rescale;
		this->StartIdx = StartIdx;
		ofile<<"Input Prefix:";
		std::string prefix;
		ifile>>prefix;
		p.Open(prefix);
	}

	/**
	 * \brief Return the i-th configuration from the pack (sequentially from StartIdx, or randomly), then rescale it.
	 * \param i  Sequential index added to StartIdx; ignored when random sampling is enabled.
	 * \return   The loaded and rescaled configuration.
	 */
	Configuration operator() (size_t i)
	{
		Configuration result;
		if (sequential_sampling){
			result=p.GetConfig(i+this->StartIdx);
		}
		else{
			result=p.GetConfig((int)(std::floor(p.NumConfig() * rng.RandomDouble())));
		}
		result.Rescale(Rescale);
		return result;
	}
};

//take every configuration from GetConfigsFunction, then rescale it to unit number density, then output
class RescaleToUnitDensity
{
public:
	std::function<const Configuration(size_t i)> GetConfigsFunction;

	/**
	 * \brief Construct the wrapper around an existing configuration source.
	 * \param GetConfigsFunction  The underlying source returning the configuration of index i.
	 */
	RescaleToUnitDensity(std::function<const Configuration(size_t i)> GetConfigsFunction) : GetConfigsFunction(GetConfigsFunction)
	{}

	/**
	 * \brief Return the i-th configuration rescaled to unit number density.
	 * \param i  Index passed to the underlying source.
	 * \return   The configuration rescaled to unit number density.
	 */
	Configuration operator() (size_t i)
	{
		Configuration result = GetConfigsFunction(i);
		result.Resize(result.NumParticle());
		return result;
	}
};

class ReadConfigPack_withInitialNumber
{
public:
	double Rescale;
	ConfigurationPack p;
	size_t InitialNumber;
	
	/**
	 * \brief Construct the loader, reading the configuration-pack prefix and the initial index offset.
	 * \param ifile    Input stream supplying the pack prefix and initial number.
	 * \param ofile    Output stream for prompts.
	 * \param Rescale  Linear scaling factor applied to each configuration when read.
	 */
	ReadConfigPack_withInitialNumber(std::istream & ifile, std::ostream & ofile, double Rescale)
	{
		this->Rescale=Rescale;
		ofile<<"Input Prefix:";
		std::string prefix;
		ifile>>prefix;
		p.Open(prefix);
		ofile<<"Input Initial Number:";
		ifile>>InitialNumber;
	}

	/**
	 * \brief Return configuration (i + InitialNumber) from the pack, then rescale it.
	 * \param i  Index offset, added to InitialNumber.
	 * \return   The loaded and rescaled configuration.
	 */
	Configuration operator() (size_t i)
	{
		Configuration result=p.GetConfig(i+InitialNumber);
		result.Rescale(Rescale);
		return result;
	}
};

class SingleConfiguration
{
public:
	Configuration c;
	/**
	 * \brief Construct a source that always returns the single given configuration.
	 * \param C  The configuration to store and return.
	 */
	SingleConfiguration( Configuration C ) : c(C)
	{
	}

	/**
	 * \brief Return the stored configuration.
	 * \param i  Configuration index; must be 0.
	 * \return   The stored configuration.
	 */
	Configuration operator() (size_t i)
	{
		assert(i==0);
		return this->c;
	}
};

/**
 * \brief Run the interactive pair-statistics command-line interpreter.
 * \return  Exit status code of the interpreter.
 */
int PairStatisticsCLI();

/**
 * \brief Program entry point; launches the pair-statistics CLI. 
 * \param argc  Number of command-line arguments.
 * \param argv  Array of command-line argument strings.
 * \return      Process exit status (always 0).
 */
int main(int, char **){
	PairStatisticsCLI();

	return 0;
}

int PairStatisticsCLI(){
	char tempstring[1000];
	/* An extern variable defined in ${core}/etc.h */
	Verbosity = 3;

	std::istream & ifile=std::cin;
	std::ostream & ofile=std::cout;

	double Rescale=1.0;
	size_t NumConfig =0;
	size_t StartIdx = 0;
	size_t NumThreads = 1;
	std::string OutputPrefix;
	std::string GraceTitle;
	std::function<const Configuration(size_t i)> GetConfigsFunction = nullptr;	
	bool AlsoWriteGrace = false;
	bool RandomSampling = false;

	std::vector<std::string> available_computations{
		"g2Computation","Directionalg2Computation", "SkComputation", "NearestNeighborStatistics", "LocalNumberVariance", "CoordinationNumber"};
	std::vector<std::string> todo_computations{
		"HpComputation",
		"HvComputation",
		"VoronoiVolumeComputation",
		"VoronoiNumSidesComputation",
		"CoveringRadiusDistributionComputation",
		"PackingRadiusDistributionComputation",
		"PercolationVolumeFractionDistribution",
		"PercolationP1Calculation",
		"PercolationP1Calculation_v2",
		"VoronoiVolumeCorrelationComputation",
		"Psi6CorrelationComputation",
		"LocalQ6DistributionComputation",
		"AnalyticWindowNumberVarianceComputation",
		"WindowNumberDistributionComputation",
		"AverageClusterSizeComputation",
		"NeighborLinkTortuosityComputation",
		"M2Computation",
		"DiscretizationVolumeFractionComputation",
		"FirstPassageTimeComputation"
	};

	std::vector<std::unique_ptr<Computation>> vpComputations;
	for(;;)
	{
		ifile>>tempstring;
		/* Basic parameters */
		if(strcmp(tempstring, "Exit")==0)
		{
			// unique_ptr elements free their computations automatically.
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
		else if(strcmp(tempstring, "GraceTitle")==0)
		{
			ifile>>GraceTitle;
		}
		else if(strcmp(tempstring, "OutputPrefix")==0)
		{
			ifile>>OutputPrefix;
		}
		else if(strcmp(tempstring, "AlsoWriteEPS")==0)
		{
			ifile>>AlsoWriteEPS;
		}
		else if(strcmp(tempstring, "AlsoWriteGrace")==0)
		{
			ifile>>AlsoWriteGrace;
		}
		else if(strcmp(tempstring, "Rescale")==0)
		{
			std::cerr << "Warning : Rescale is deprecated\n";
			ifile>>Rescale;
		}
		else if(strcmp(tempstring, "NumConfig")==0)
		{
			ifile>>NumConfig;
		}
		else if (strcmp(tempstring, "StartIndex") == 0)
		{
			ifile>> StartIdx;
			ofile<< "Start Index = "<< StartIdx <<std::endl;
		}
		else if (strcmp(tempstring, "RandomSampling") == 0)
		{
			ifile >> RandomSampling;
			if (RandomSampling){
				ofile << "Turn on random sampling in GetConfigsFunction" << std::endl;
				ofile << "This option is effective only for ReadConfigPack" << std::endl;
				//TODO: It is applied only to the ReadConfigPack case. Please apply it ot other cases

			}
			else{
				ofile << "Turn off random sampling in GetConfigsFunction" << std::endl;				
			}
		}
		/* define computations */
		else if(strcmp(tempstring, "g2Computation")==0)
		{
			vpComputations.push_back( std::make_unique<g2Computation>(ifile, ofile) );
		}
		else if(strcmp(tempstring, "Directionalg2Computation")==0)
		{
			vpComputations.push_back( std::make_unique<Directionalg2Computation>(ifile, ofile) );
		}
		else if (strcmp(tempstring, "SkComputation") == 0)
		{
			vpComputations.push_back(std::make_unique<SkComputation>(ifile, ofile));
		}
		else if (strcmp(tempstring, "NearestNeighborStatistics") == 0)
		{
			vpComputations.push_back(std::make_unique<NearestNeighbors>(ifile, ofile));
		}
		else if (strcmp(tempstring, "CoordinationNumber") == 0)
		{
			vpComputations.push_back(std::make_unique<CoordinationNumberComputation>(ifile, ofile));
		}



		else if (strcmp(tempstring, "GetConfigsFunction") == 0)
		{
			ofile<<"Input Type:";
			ifile>>tempstring;
			if(strcmp(tempstring, "PrefixNumberSuffix")==0)
			{
				PrefixNumberSuffix c(ifile, ofile, Rescale);
				GetConfigsFunction = c;
			}
			else if (strcmp(tempstring, "PrefixNumberSuffix_withInitialNumber") == 0)
			{
				PrefixNumberSuffix_withInitialNumber c(ifile, ofile, Rescale);
				GetConfigsFunction = c;
			}
			else if (strcmp(tempstring, "ReadConfigPack") == 0)
			{
				ReadConfigPack c(ifile, ofile, Rescale, StartIdx);
				if (RandomSampling){
					c.sequential_sampling = false;
					NumConfig=c.p.NumConfig();
					std::cout<<"The config pack contains "<<NumConfig<<" configurations, set NumConfig to this value.\n";
				}
				//automatic sets NumConfig since Configuration Pack contains this information
				else{
					c.sequential_sampling = true;
					if (StartIdx >= static_cast<size_t>(c.p.NumConfig())){
						std::cout<<"Start index exceeds the number of configurations in the pack; set NumConfig to 0.\n";
						NumConfig = 0;
					}
					else{
						NumConfig=c.p.NumConfig()-StartIdx;
						std::cout<<"The config pack contains "<<NumConfig<<" configurations, set NumConfig to this value.\n";
					}
				}
				GetConfigsFunction = c;

			}
			else if (strcmp(tempstring, "ReadConfigPack_withInitialNumber") == 0)
			{
				ReadConfigPack_withInitialNumber c(ifile, ofile, Rescale);
				GetConfigsFunction = c;
				//automatic sets NumConfig since Configuration Pack contains this information
				if (c.InitialNumber >= static_cast<size_t>(c.p.NumConfig())){
					std::cout << "Initial number exceeds the number of configurations in the pack; set NumConfig to 0.\n";
					NumConfig = 0;
				}
				else{
					NumConfig = c.p.NumConfig() - c.InitialNumber;
					std::cout << "The config pack contains " << c.p.NumConfig() << " configurations, set NumConfig to " << NumConfig << "\n";
				}
			}
			else if (strcmp(tempstring, "SingleConfiguration") == 0)
			{
				std::string prefix;
				ofile << "Input prefix:";
				ifile >> prefix;
				SingleConfiguration c(ReadPos(prefix));
				GetConfigsFunction = c;
				//automatic sets NumConfig 
				NumConfig = 1;
				std::cout << "set NumConfig to " << 1 << "\n";
			}
			else
			{
				std::cout<<"Unrecognized Type!\n";
				std::cin.clear();
			}
		}
		else if (strcmp(tempstring, "PrintConfigs") == 0) { //print ConfigurationPack => readable txt files
			if (GetConfigsFunction == nullptr)
				std::cerr << "Specify GetConfigsFunction before printing!\n";
			else 
			{
				
				for (size_t index = 0; index < NumConfig; index++) {
					std::string file = OutputPrefix + std::string("__") + std::to_string(index);
					//sprintf(file, "%s__%d", OutputPrefix, index);
					WriteConfiguration(GetConfigsFunction(index), file);
				}
			}
		}

		else if (strcmp(tempstring, "RescaleToUnitDensity") == 0)
		{
			if (GetConfigsFunction == nullptr)
				std::cerr << "Specify GetConfigsFunction before specifying rescaling!\n";
			else
			{
				RescaleToUnitDensity c(GetConfigsFunction);
				GetConfigsFunction = c;
			}
		}
		else if(strcmp(tempstring, "NumThreads")==0){
			ifile >> NumThreads;
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
					if(AlsoWriteGrace)
						(*iter)->Plot(OutputPrefix, GraceTitle);
				}
			}
		}
		else if (strcmp(tempstring, "LocalNumberVariance") == 0)
		{
			vpComputations.push_back( std::make_unique<LocalNumberVariance>(ifile, ofile) );
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

