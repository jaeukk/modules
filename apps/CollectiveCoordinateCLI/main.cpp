/** \file main.cpp
 *  \brief Thin entry point for the CollectiveCoordinateCLI application.
 *
 *  The whole command language lives in CollectiveCoordinateCLI() (adapted from
 *  Ge Zhang's legacy c++/CollectiveCoordinateCLI.cpp).  This wrapper provides
 *  the library globals the executable owns and the startup calls the legacy
 *  c++/main.cpp performed, then hands control to the dispatcher:
 *    - Verbosity  : declared extern in core/etc.h; the core object leaves it
 *                   undefined, so each main defines its own (PROGRESS.md R5 /
 *                   new_c++ convention).
 *    - ProgramStart / TimeLimit : core/etc.cpp zero-initializes these, so the
 *                   legacy startup set them before running the CLI.  We
 *                   replicate that: RelaxStructure_NLOPT and the MD loops compare
 *                   std::time(nullptr) > ::TimeLimit and would exit(0)
 *                   immediately if TimeLimit stayed 0.
 *    - nlopt_srand / omp_set_nested / stream precision : legacy c++/main.cpp
 *                   parity (deterministic NLopt RNG, nested OpenMP, 10-digit
 *                   numeric output).
 *
 *  Usage:
 *      ./CollectiveCoordinateCLI.out < commands.txt
 *  The command tokens are read from stdin; see README.md for the full table.
 */

#include <ctime>
#include <iostream>

#include <nlopt.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <etc.h>

/// Library verbosity global owned by this executable (declared extern in
/// core/etc.h; the core object does not define it).
size_t Verbosity = 2;

/** \brief Run the collective-coordinate command interpreter.
 *  \return Dispatcher exit status (0 on the "Exit" command, 2 on stdin EOF). */
int CollectiveCoordinateCLI();

/** \brief Program entry point; sets up the legacy startup state then delegates.
 *  \return The dispatcher's return value. */
int main()
{
	nlopt_srand(99999);                     // legacy c++/main.cpp parity: deterministic NLopt RNG
	ProgramStart = std::time(nullptr);
	TimeLimit = ProgramStart + 86400000;    // default time limit of ~1000 days
#ifdef _OPENMP
	omp_set_nested(true);                   // legacy parity: allow nested OpenMP regions
#endif
	std::cout.precision(10);                // legacy parity: 10-digit numeric output
	std::cerr.precision(10);
	return CollectiveCoordinateCLI();
}
