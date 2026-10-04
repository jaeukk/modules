/** \file main.cpp
 *  \brief Entry point for the CC-inversion application.
 *
 *  Reconstructs point configurations from a subset of their collective
 *  coordinates (structure-factor constraints) by minimizing the
 *  FrankInvertCCPotential family of objectives.
 *
 *  This driver is extracted from the legacy `cli()` dispatcher in
 *  Interfaces.h: of all the commands that dispatcher handled, only the
 *  "CCInversion" branch belongs to this application, so we keep just that one
 *  here (plus the shared "Verbosity"/"Exit" controls).  The heavy command set
 *  (PairStatistics, CollectiveCoordinate, Tiling, Draw/Display, ...) lives with
 *  its own modules and is intentionally not reachable from this binary.
 *
 *  Usage:
 *      ./CCInversion.out [masterDir]      # masterDir defaults to "."
 *  then on stdin:
 *      CCInversion <0|1>                  # 0/1 = print all numerical solutions?
 *      Exit
 *  Output is written under <masterDir>/results/Reconst/ (see CCInversionProblem1).
 */

#include <iostream>
#include <string>
#include <cstring>

#include "CCInversionProblem.h"

// Library globals owned by the executable (modular Verbosity convention: the
// core/etc.cpp object does NOT define these; each main provides its own).
size_t Verbosity = 2;
// Base directory for the results/Reconst/ output tree used by CCInversionProblem1.
std::string masterDir = ".";

int main(int argc, char ** argv)
{
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;

	// Optional first argument overrides the output base directory.
	if (argc > 1)
		masterDir = argv[1];

	char tempstring[1000] = {};
	for (;;)
	{
		ifile >> tempstring;
		if (!ifile)
			return 0;                       // EOF -> done
		if (strcmp(tempstring, "Exit") == 0)
		{
			return 0;
		}
		else if (strcmp(tempstring, "Verbosity") == 0)
		{
			ifile >> Verbosity;
		}
		else if (strcmp(tempstring, "CCInversion") == 0)
		{
			bool outConfig;
			ofile << "Do we print all numerical solutions?\n";
			ifile >> outConfig;
			CCInversionProblem1(outConfig);
		}
		else
			std::cerr << "Unrecognized command!\n";
	}
	return 0;
}
