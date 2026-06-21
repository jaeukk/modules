#ifndef CCINVERSION_INCLUDED
#define CCINVERSION_INCLUDED

#include "Potential.h"
#include "StructureFactor.h"
#include "StructureOptimization.h"
#include "CCInversion.h"   // CC-inversion potentials + LatticeUtil (replaces etc1.h)
#include "Discrepancy.h"   // Merge_Sort / MergeSort (was reached via etc1.h)
#include "PeriodicCellList.h"

extern std::string masterDir;

//results [0]=N	[1]=Nk	[2]=max{dx}	[3]=max{ground state energy}	[4]=the number of tried initial configurations
//the default optimization method is NLOPT+MINOP
//possible methods = NLOP+MINOP
//					SteepestDescent
//					ConjugateGradient
//					LocalGradientDescent
/** \brief Test reconstruction (non-)uniqueness of a target from its lowest-|k| collective coordinates. (CLAUDE)
 *  \param target         Configuration to reconstruct (also the source of the target collective coordinates).
 *  \param Nc             Number of random initial configurations to relax.
 *  \param N_kPoints      Number of collective-coordinate constraints (interpreted via RealImagParts).
 *  \param[out] results   Summary vector: [0]=N, [1]=Nk, [2]=max{dx}, [3]=max{ground-state energy}, [4]=number of tried initial configurations.
 *  \param outFile        Output stream for diagnostics.
 *  \param Tolerance      Maximum allowed displacement of a reconstructed particle from the target before it is deemed non-unique.
 *  \param NKPointsSkip   Number of lowest-|k| points to skip before taking the constraint k-points.
 *  \param RealImagParts  If true, real and imaginary parts of each collective coordinate are constrained together (Nk=N_kPoints); if false, separately (Nk=N_kPoints/2).
 *  \param method         Optimization method ("NLOP+MINOP", "SteepestDescent", "ConjugateGradient", "LocalGradientDescent"; null defaults to NLOPT+MINOP).
 *  \return               True if every reconstruction stays within Tolerance of the target (unique), false otherwise. */
bool FrankInvCCUniqueTest(const Configuration & target, size_t Nc, size_t N_kPoints, std::vector<double> &results, std::fstream & outFile, double Tolerance = 1e-4, size_t NKPointsSkip = 0, bool RealImagParts = true, char *method = nullptr);
/** \brief Convenience overload of FrankInvCCUniqueTest that discards the results vector. (CLAUDE)
 *  \param target        Configuration to reconstruct.
 *  \param Nc            Number of random initial configurations to relax.
 *  \param N_kPoints     Number of collective-coordinate constraints.
 *  \param outFile       Output stream for diagnostics.
 *  \param Tolerance     Maximum allowed displacement before the target is deemed non-unique.
 *  \param NKPointsSkip  Number of lowest-|k| points to skip.
 *  \return              True if the target is uniquely reconstructed, false otherwise. */
bool FrankInvCCUniqueTest(const Configuration & target, size_t Nc, size_t N_kPoints, std::fstream & outFile, double Tolerance = 1e-4, size_t NKPointsSkip = 0);

//results [0]=N	[1]=Nk	[2]=max{dx}	[3]=max{ground state energy}	[4]=the number of tried initial configurations [5]=the number of correctly reconstructed configs.
//the default optimization method is NLOPT+MINOP
//possible methods = NLOP+MINOP
//					SteepestDescent
//					ConjugateGradient
//					LocalGradientDescent
/** \brief Like FrankInvCCUniqueTest but also counts how many random initial configurations reconstruct correctly. (CLAUDE)
 *  \param target         Configuration to reconstruct.
 *  \param Nc             Number of random initial configurations to relax.
 *  \param N_kPoints      Number of collective-coordinate constraints.
 *  \param[out] results   Summary vector: [0]=N, [1]=Nk, [2]=max{dx}, [3]=max{ground-state energy}, [4]=number of tried initial configurations, [5]=number of correctly reconstructed configurations.
 *  \param outFile        Output stream for diagnostics.
 *  \param Tolerance      Maximum allowed displacement before a reconstruction is deemed incorrect.
 *  \param NKPointsSkip   Number of lowest-|k| points to skip.
 *  \param RealImagParts  If true, real and imaginary parts are constrained together; if false, separately.
 *  \param method         Optimization method (see FrankInvCCUniqueTest; null defaults to NLOPT+MINOP).
 *  \return               True if the target is uniquely reconstructed, false otherwise. */
bool FrankInvCCUniqueTest2(const Configuration & target, size_t Nc, size_t N_kPoints, std::vector<double> &results, std::fstream & outFile, double Tolerance = 1e-4, size_t NKPointsSkip = 0, bool RealImagParts = true, char *method = nullptr);

//constraint types = both: (default) real and imaginary parts of collective coordinates
//					real: real parts first
//					imag: imaginary parts first
//if the "proximity" vector  is empty => don't calculate proximity
/** \brief Reconstruct a target from its collective coordinates, collecting all distinct solutions, their energy spreads, and optional proximity metrics. (CLAUDE)
 *  \param target           Configuration to reconstruct.
 *  \param numGroundStates  Number of random initial configurations to relax (trial count).
 *  \param numConstraints   Number of collective-coordinate constraints (real/imag counting depends on constType).
 *  \param[out] results     Summary vector: [0]=N, [1]=numConstraints, [2]=max{dx}, [3]=max energy, [4]=number of ground states found, [5]=number correctly reconstructed, followed by (minE,maxE,maxDx) per distinct solution.
 *  \param[out] sols        Distinct solution configurations found (sols[0] is the sorted target).
 *  \param outFile          Output stream for diagnostics.
 *  \param PositionTolerance  Maximum displacement at which two configurations are treated as identical.
 *  \param EnergyTolerance    Energy below which a relaxed configuration is accepted as a ground state.
 *  \param NKPointsSkip     Number of lowest-|k| points to skip.
 *  \param constType        Constraint type: "both" (default), "real" (real parts first), or "image" (imaginary parts first).
 *  \param method           Optimization method selector (0=NLOPT+MINOP, 1=SteepestDescent, 2=ConjugateGradient, 3=LocalGradientDescent).
 *  \param[out] energy      Per-solution lists of ground-state energies (energy[0]=target, energy[1]=secondary, energy[2]=others).
 *  \param initConfig       Optional initial configurations to relax from; if empty, random configurations are used.
 *  \param[in,out] proximity  If non-empty on entry, filled with proximity metrics: [0]=all solutions vs target, [n+1]=metrics for the n-th solution; left untouched if empty.
 *  \return                 True if the target is uniquely reconstructed within PositionTolerance, false otherwise. */
bool FrankInvCCUniqueTest3(const Configuration & target, size_t numGroundStates, size_t numConstraints, std::vector<double> &results, std::vector<Configuration> &sols, std::fstream & outFile, double PositionTolerance, double EnergyTolerance, size_t NKPointsSkip, char * constType, int method, std::vector<std::vector<double>> &energy, const std::vector<Configuration> & initConfig, std::vector<std::vector<double>> & proximity);

//perform inversion procedure for d=2,3 systems
//two constraints types are implemented; spherical and linear (k_i, k_i+k_j)
/** \brief Reconstruction test for d=2,3 systems using spherical or linear (k_i, k_i+k_j) constraint sets. (CLAUDE)
 *  \param target           Configuration to reconstruct.
 *  \param numGroundStates  Number of random or loaded initial configurations to relax.
 *  \param numConstraints   Constraint count: number per line for "linear", or total number for "spherical".
 *  \param[out] results     Summary vector: [0]=N, [1]=constraint density, [2]=max{dx}, [3]=max energy, [4]=number of ground states, [5]=number correctly reconstructed, then (minE,maxE,maxDx) per distinct solution.
 *  \param[out] sols        Distinct solution configurations found (sols[0] is the sorted target).
 *  \param outFile          Output stream for diagnostics.
 *  \param PositionTolerance  Maximum displacement at which two configurations are treated as identical.
 *  \param EnergyTolerance    Energy below which a relaxed configuration is accepted as a ground state.
 *  \param NKPointsSkip     Number of lowest-|k| points (or k-axis indices) to skip.
 *  \param constType        Constraint geometry: "linear" (k_i and k_i+k_j axes) or anything else for "spherical".
 *  \param method           Optimization method ("NLOP+MINOP", "SteepestDescent", "ConjugateGradient", "LocalGradientDescent"; null defaults to NLOPT+MINOP).
 *  \param initConfigs      Optional path to a ConfigurationPack of initial configurations; null or too few falls back to random configurations.
 *  \return                 True if the target is uniquely reconstructed within PositionTolerance, false otherwise. */
bool FrankInvCCUniqueTest4(const Configuration & target, size_t numGroundStates, size_t numConstraints, std::vector<double> &results, std::vector<Configuration> &sols, std::fstream & outFile, double PositionTolerance = 1e-4, double EnergyTolerance = 1e-20, size_t NKPointsSkip = 0, char * constType = nullptr, char *method = nullptr, char * initConfigs = nullptr);

//performed inversion problems from given collective coordinates using FrankInvCCUniqueTest...
// outConfig = true: output configurations
/** \brief Interactive driver: sweep system sizes/constraint counts and run the CC-inversion uniqueness tests from stdin parameters. (CLAUDE)
 *  \param outConfig  If true, write out the reconstructed solution configurations in addition to the summary statistics. */
void CCInversionProblem(bool outConfig);

//use the same target configurations!!
/** \brief Variant of CCInversionProblem that reuses the same target configurations and writes output under masterDir/results/Reconst/. (CLAUDE)
 *  \param outConfig  If true, write out the reconstructed solution configurations in addition to the summary statistics. */
void CCInversionProblem1(bool outConfig);

#endif