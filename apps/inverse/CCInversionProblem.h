#ifndef CCINVERSION_INCLUDED
#define CCINVERSION_INCLUDED

#include "Potential.h"
#include "StructureFactor.h"
#include "StructureOptimization.h"
#include "CCInversion.h"   // CC-inversion potentials (FrankInvertCCPotential*)
#include <LatticeUtil.h>   // GetLattice_config / Find_LeastCommonDivider (misc module)
#include "Discrepancy.h"   // Merge_Sort / MergeSort (was reached via etc1.h)
#include "PeriodicCellList.h"

extern std::string masterDir;

//results [0]=N	[1]=Nk	[2]=max{dx}	[3]=max{ground state energy}	[4]=the number of tried initial configurations
//the default optimization method is NLOPT+MINOP
//possible methods = NLOP+MINOP
//					SteepestDescent
//					ConjugateGradient
//					LocalGradientDescent
/** \brief Test reconstruction of a target from its lowest-magnitude collective coordinates.
 * \param target Target configuration.
 * \param Nc Number of random initial configurations.
 * \param N_kPoints Requested constraint count.
 * \param[out] results Summary: particle count, requested constraint count, maximum displacement, maximum energy, and relaxation-attempt count.
 * \param outFile Diagnostic stream.
 * \param Tolerance Displacement threshold for uniqueness.
 * \param NKPointsSkip Number of low-magnitude k-points to skip.
 * \param RealImagParts If true, use N_kPoints k-points with both components; if false, use round(N_kPoints/2 + 0.1) k-points with scalar components counted separately.
 * \param method Optimization method; null selects NLOPT followed by MINOP.
 * \return True if the maximum reconstructed displacement does not exceed Tolerance. */
bool FrankInvCCUniqueTest(const Configuration & target, size_t Nc, size_t N_kPoints, std::vector<double> &results, std::fstream & outFile, double Tolerance = 1e-4, size_t NKPointsSkip = 0, bool RealImagParts = true, char *method = nullptr);

/** \brief Convenience overload of FrankInvCCUniqueTest that discards the results vector.
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
/** \brief Like FrankInvCCUniqueTest but also counts how many random initial configurations reconstruct correctly.
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
/** \brief Reconstruct a target, collecting distinct accepted ground-state solutions, energy samples, and optional proximity metrics.
 * \param target Target configuration.
 * \param numGroundStates Number of initial configurations to relax.
 * \param numConstraints Requested constraint count.
 * \param[out] results Summary containing N, the constraint count, the maximum displacement recorded for target-matching states, accepted-state counts, and per-solution energy/displacement data.
 * \param[out] sols Distinct accepted solutions; the sorted target is sols[0].
 * \param outFile Diagnostic stream.
 * \param PositionTolerance Threshold used to classify two configurations as the same solution.
 * \param EnergyTolerance Energy threshold for accepting a relaxed state.
 * \param NKPointsSkip Number of low-magnitude k-points to skip.
 * \param constType "both", "real", or "image".
 * \param method Optimization-method selector from 0 through 3.
 * \param[out] energy Energy samples grouped by solution category.
 * \param initConfig Optional initial configurations; random configurations are used when empty.
 * \param[in,out] proximity Nonempty input enables collection of proximity metrics.
 * \return False only if the displacement recorded for target-matching states exceeds PositionTolerance; distinct non-target solutions do not affect the return value. */
bool FrankInvCCUniqueTest3(const Configuration & target, size_t numGroundStates, size_t numConstraints, std::vector<double> &results, std::vector<Configuration> &sols, std::fstream & outFile, double PositionTolerance, double EnergyTolerance, size_t NKPointsSkip, const char * constType, int method, std::vector<std::vector<double>> &energy, const std::vector<Configuration> & initConfig, std::vector<std::vector<double>> & proximity);

//perform inversion procedure for d=2,3 systems
//two constraints types are implemented; spherical and linear (k_i, k_i+k_j)
/** \brief Reconstruction test for two- and three-dimensional systems using linear or spherical k-point sets.
 * \param target Target configuration.
 * \param numGroundStates Number of initial configurations to relax.
 * \param numConstraints Number per reciprocal-space line for "linear", or total count for spherical mode.
 * \param[out] results Summary containing N, M/(d(N-1)), the maximum displacement recorded for target-matching states, accepted-state counts, and per-solution data.
 * \param[out] sols Distinct accepted solutions; the sorted target is sols[0].
 * \param outFile Diagnostic stream.
 * \param PositionTolerance Solution-matching threshold.
 * \param EnergyTolerance Ground-state acceptance threshold.
 * \param NKPointsSkip Number of reciprocal indices or low-magnitude k-points to skip.
 * \param constType Non-null constraint geometry; "linear" selects line constraints and other strings select spherical constraints.
 * \param method Non-null optimization-method string.
 * \param initConfigs Optional ConfigurationPack path; it is used only when the pack contains more than numGroundStates configurations.
 * \return False only if the displacement recorded for target-matching states exceeds PositionTolerance; distinct non-target solutions do not affect the return value. */
bool FrankInvCCUniqueTest4(const Configuration & target, size_t numGroundStates, size_t numConstraints, std::vector<double> &results, std::vector<Configuration> &sols, std::fstream & outFile, double PositionTolerance = 1e-4, double EnergyTolerance = 1e-20, size_t NKPointsSkip = 0, char * constType = nullptr, char *method = nullptr, char * initConfigs = nullptr);

//performed inversion problems from given collective coordinates using FrankInvCCUniqueTest...
// outConfig = true: output configurations
/** \brief Interactive driver: sweep system sizes/constraint counts and run the CC-inversion uniqueness tests from stdin parameters.
 *  \param outConfig  If true, write out the reconstructed solution configurations in addition to the summary statistics. */
void CCInversionProblem(bool outConfig);

//use the same target configurations!!
/** \brief Variant of CCInversionProblem that reuses the same target configurations and writes output under masterDir/results/Reconst/.
 *  \param outConfig  If true, write out the reconstructed solution configurations in addition to the summary statistics. */
void CCInversionProblem1(bool outConfig);

#endif