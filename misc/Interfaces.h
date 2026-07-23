#ifndef INTERFACES_INCLUDED
#define INTERFACES_INCLUDED
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <PeriodicCellList.h>
#include <Plots.h>
#include <DisplaySpheres.h>

/** \brief Run the inverse statistical-mechanics command-line interface.
 *  \return  Exit status code. */
int InvStatMechCLI();

/** \brief Run the pair-statistics command-line interface.
 *  \return  Exit status code. */
int PairStatisticsCLI();

/** \brief Run the collective-coordinate command-line interface.
 *  \return  Exit status code. */
int CollectiveCoordinateCLI();

/** \brief Run the tiling-problem command-line interface.
 *  \return  Exit status code. */
int TilingProblemCLI();

/** \brief Run the first collective-coordinate inversion problem.
 *  \param outConfig Whether to write all numerical solutions. */
void CCInversionProblem1(bool outConfig);

//defined in debug.cpp
//int Debug();

/** \brief Run the simple debug routine.
 *  \return  Exit status code. */
int DebugSimple();
//defined in debug.cpp
//int CalculateTau(void);
//int CalculateEvAndQuantizerError();

//How can I make a generic class for inputting parameters?
template<typename T> class InputParameters {
	T parameter;
	std::string question, response, keyword;

	/** \brief Construct an input-parameter record from an initial value, keyword, prompt, and response text.
	 * \param PARAMETER Initial stored parameter value.
	 * \param KEYWORD Keyword that selects this parameter.
	 * \param QUESTION Prompt text shown when requesting the value.
	 * \param RESPONSE Response or echo text associated with the parameter. */
	InputParameters(T & PARAMETER, const std::string & KEYWORD, const std::string & QUESTION, const std::string & RESPONSE) {
		keyword = std::string(KEYWORD);
		response = std::string(RESPONSE);
		question = std::string(QUESTION);
		parameter = PARAMETER;
	}

	/** \brief Test whether the given keyword matches this parameter's keyword.
	 *  \param KEYWORD  Keyword to compare against.
	 *  \return         True if the keywords are equal. */
	bool CompareKeyword(const std::string & KEYWORD) { 
		//return strcmp(this->keyword, KEYWORD) == 0; 
		return KEYWORD.compare(this->keyword) == 0;
		}

	/** \brief Input variables
	 *  \param ifile  Input stream to read the parameter value from.
	 *  \param ofile  Output stream for echoing prompts/responses. */
	void Input(std::istream & ifile, std::ofstream & ofile) { 

	}
};


/** \brief Read commands from standard input and dispatch them to program subsystems.
 * \return Status returned by PairStatisticsCLI or CollectiveCoordinateCLI, or 0 on Exit or end of input. */
inline int cli()
{
	char tempstring[1000] = {};
	std::istream & ifile=std::cin;
	std::ostream & ofile=std::cout;
	double BondLength = 0.0, Radius = 0.0;
	for(;;)
	{
		if(!(ifile >> std::setw(1000) >> tempstring))//setw bounds the read to the buffer; a failed read (EOF) stops the loop instead of re-dispatching the previous command forever
			return 0;
		if(strcmp(tempstring, "Exit")==0)
		{
			return 0;
		}
		else if (strcmp(tempstring, "Verbosity") == 0)
		{
			ifile >> Verbosity;
		}
		else if (strcmp(tempstring, "AlsoWriteEPS") == 0)
		{
			ifile >> AlsoWriteEPS;
		}
		else if (strcmp(tempstring, "TimeLimit") == 0)
		{
			size_t limit;
			ifile>>limit;
			::TimeLimit=::ProgramStart+limit;
		}
		/*else if(strcmp(tempstring, "InvStatMech")==0)
		{
			return InvStatMechCLI();
		}*/
		else if(strcmp(tempstring, "PairStatistics")==0)
		{
			return PairStatisticsCLI();
		}
		else if(strcmp(tempstring, "CollectiveCoordinate")==0)
		{
			return CollectiveCoordinateCLI();
		}
		else if(strcmp(tempstring, "Debug1")==0)
		{
			DebugSimple();
		}
		else if (strcmp(tempstring, "Tiling") == 0)
		{
			TilingProblemCLI();
		}
		else if (strcmp(tempstring, "CCInversion") == 0)
		{
			bool outConfig;
			ofile <<"Do we print all numerical solutions?\n";
			ifile >> outConfig;
			CCInversionProblem1(outConfig);
		}
		else if (strcmp(tempstring, "BondLength") == 0)
		{
			ifile >> BondLength;
		}
		else if (strcmp(tempstring, "Radius") == 0)
		{
			ifile >> Radius;
		}
		else if (strcmp(tempstring, "Draw") == 0)
		{
			std::string prefix;
			std::cin>>prefix;
			Configuration c=ReadPos(prefix);
			c.StandardizeOrientation();
			c.Resize(1.0);
			::AlsoWriteEPS=true;
			::Plot(prefix, c);
		}
		else if(strcmp(tempstring, "Display")==0)
		{
			std::string prefix;
			std::cin>>prefix;
			Configuration c=ReadPos(prefix);
			DisplayConfiguration(c, BondLength, Radius);
		}
		else if(strcmp(tempstring, "DisplayMultiple")==0)
		{
			std::string prefix;
			std::cin>>prefix;
			Configuration c=ReadPos(prefix);
			DisplayConfiguration(MultiplicateStructure(c, 5), BondLength, Radius);
		}
		else if(strcmp(tempstring, "DisplaySk")==0)
		{
			std::string prefix;
			std::cin>>prefix;
			Configuration c=ReadPos(prefix);
			double MaxK;
			std::cin>>MaxK;
			Display3DStructureFactor(c, MaxK);
		}
		else if(strcmp(tempstring, "Movie")==0)
		{
			std::string prefix;
			std::cin>>prefix;
			ConfigurationPack pk(prefix);
			Display3DConfigurationMovie(pk, BondLength, Radius);
		}
		else if (strcmp(tempstring, "Separate") == 0)
		{
			std::string prefix1, prefix2;
			std::cout << "Input ConfigPack=";
			std::cin >> prefix1;
			std::cout << "Output ConfigPack=";
			std::cin >> prefix2;
			size_t m, n;
			std::cout << "number of sub-config per configuration=";
			std::cin >> m;
			std::cout << "the sub-config to extract=";
			std::cin >> n;
			ConfigurationPack pk1(prefix1);
			ConfigurationPack pk2(prefix2);
			pk2.Clear();
			for (int i = 0; i < pk1.NumConfig(); i++)
			{
				Configuration c = pk1.GetConfig(i);
				c = SeparateList(c, m, n);
				pk2.AddConfig(c);
			}
		}
		/*else if (strcmp(tempstring, "CalculateTau") == 0)
		{
			CalculateTau();
		}
		else if (strcmp(tempstring, "CalculateEvAndQuantizerError") == 0)
		{
			CalculateEvAndQuantizerError();
		}*/
		else
			std::cerr << "Unrecognized command!\n";
	}
}
#endif