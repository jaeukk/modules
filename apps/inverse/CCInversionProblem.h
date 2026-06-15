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
bool FrankInvCCUniqueTest(const Configuration & target, size_t Nc, size_t N_kPoints, std::vector<double> &results, std::fstream & outFile, double Tolerance = 1e-4, size_t NKPointsSkip = 0, bool RealImagParts = true, char *method = nullptr);
bool FrankInvCCUniqueTest(const Configuration & target, size_t Nc, size_t N_kPoints, std::fstream & outFile, double Tolerance = 1e-4, size_t NKPointsSkip = 0);

//results [0]=N	[1]=Nk	[2]=max{dx}	[3]=max{ground state energy}	[4]=the number of tried initial configurations [5]=the number of correctly reconstructed configs.
//the default optimization method is NLOPT+MINOP
//possible methods = NLOP+MINOP
//					SteepestDescent
//					ConjugateGradient
//					LocalGradientDescent
bool FrankInvCCUniqueTest2(const Configuration & target, size_t Nc, size_t N_kPoints, std::vector<double> &results, std::fstream & outFile, double Tolerance = 1e-4, size_t NKPointsSkip = 0, bool RealImagParts = true, char *method = nullptr);

//constraint types = both: (default) real and imaginary parts of collective coordinates
//					real: real parts first
//					imag: imaginary parts first
//if the "proximity" vector  is empty => don't calculate proximity
bool FrankInvCCUniqueTest3(const Configuration & target, size_t numGroundStates, size_t numConstraints, std::vector<double> &results, std::vector<Configuration> &sols, std::fstream & outFile, double PositionTolerance, double EnergyTolerance, size_t NKPointsSkip, char * constType, int method, std::vector<std::vector<double>> &energy, const std::vector<Configuration> & initConfig, std::vector<std::vector<double>> & proximity);

//perform inversion procedure for d=2,3 systems
//two constraints types are implemented; spherical and linear (k_i, k_i+k_j)
bool FrankInvCCUniqueTest4(const Configuration & target, size_t numGroundStates, size_t numConstraints, std::vector<double> &results, std::vector<Configuration> &sols, std::fstream & outFile, double PositionTolerance = 1e-4, double EnergyTolerance = 1e-20, size_t NKPointsSkip = 0, char * constType = nullptr, char *method = nullptr, char * initConfigs = nullptr);

//performed inversion problems from given collective coordinates using FrankInvCCUniqueTest...
// outConfig = true: output configurations
void CCInversionProblem(bool outConfig);

//use the same target configurations!!
void CCInversionProblem1(bool outConfig);

#endif