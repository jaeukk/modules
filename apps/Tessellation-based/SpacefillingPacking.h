#ifndef TESS_SPACEFILLINGPACKING_H
#define TESS_SPACEFILLINGPACKING_H

/** \file SpacefillingPacking.h
 *  \brief Space-filling multiscale sphere packing.
 *
 *  Split out of the legacy monolithic TilingProblem.{h,cpp}: this unit owns the
 *  `SpacefillingPacking` class (an RSA that keeps adding ever-smaller spheres,
 *  one size per "stage", to fill space) and the free-function drivers
 *  `getSpacefillingPacking` / `getSpacefillingPacking2`.
 */

#include <boost/math/special_functions/zeta.hpp>

#include <etc.h>
#include <PeriodicCellList.h>
#include <RandomSequentialAddition.h>

// Library globals owned by the executable's main (see main.cpp).
extern std::string masterDir;
extern size_t Verbosity;

const int MAX_NUM_VOXEL = (int)5e8;

struct SpacefillingPacking : public RSA
{
protected:
	//Packing Fraction: it should be updated at each stage.
	double phi;
	double voxelSize;
public:
	double maxR;
	//scaling function for particle-volumes
	std::function <double(int n)> f;
	int N;

	/** \brief Construct from an initial packing; sets phi=0 and a default voxel size. (CLAUDE)
	 *  \param config  Initial sphere packing used as the substrate. */
	SpacefillingPacking(const SpherePacking& config);
	/** \brief Construct and derive N and maxR from the scaling-function parameters. (CLAUDE)
	 *  \param config  Initial sphere packing used as the substrate.
	 *  \param param   param[0] = maximum particle volume, param[1] = scaling exponent/ratio.
	 *  \param fnName  Scaling-function name ("Zeta" or "Power"). */
	SpacefillingPacking(const SpherePacking& config, const std::vector<double> & param, const char * fnName);

	/** \brief Replace the current packing and reset cell size and packing fraction. (CLAUDE)
	 *  \param newPacking  New sphere packing to copy in as the current configuration. */
	inline void ReplaceConfig(SpherePacking & newPacking) {
		this->pSpheres = new SpherePacking(newPacking);
		this->pSpheres->SetCellSize(2.0*this->maxR);
		this->phi = this->pSpheres->PackingFraction();
	}
	//Manual inputs
	/** \brief Set the number of spheres added per stage. (CLAUDE)
	 *  \param N  Number of spheres to insert at each stage. */
	void setN(int N);
	/** \brief Set the maximum (stage-1) particle radius. (CLAUDE)
	 *  \param Rmax  Maximum particle radius. */
	void setRmax(double Rmax);
	/** \brief set the scaling function of radius
	 *  \param param   param[0] = max volume, param[1] = parameter.
	 *  \param fnName  "Zeta": f(n) = 1/n^param[1]; "Power": f(n) = 1/param[1]^(n-1). */
	void setScalingFunction(const std::vector<double> & param, const char * fnName);
	/** \brief Set the per-stage volume scaling function directly. (CLAUDE)
	 *  \param function  Callable mapping stage index n to the relative particle volume. */
	void setScalingFunction(const std::function <double(int n)> & function);

	/** \brief Add (numAdd) spheres of radius (r); updates this->phi
	 *  \param r           Radius of the spheres to add.
	 *  \param numAdd      Number of spheres to insert.
	 *  \param numThreads  Number of threads used for the insertion. */
	void AddSpheres(double r, int numAdd, int numThreads = 1);

	/** \brief for (n = n0; n< n0+numSteps; n++) AddSpheres(r*f(n), N);
	 *  \param numSteps  Number of successive stages to run.
	 *  \param n0        Stage index of the first step. */
	void RepeatAddition(int numSteps, int n0);
	/** \brief Compute the packing fraction, when one can insert "N" particles at each stage from 1~stage.
	 *  \param stage  Final stage up to which the theoretical packing fraction is summed.
	 *  \return       Theoretical packing fraction accumulated through the given stage. */
	double GetTheoreticalPackingFraction(int stage);
	/** \brief Return the current packing fraction.
	 *  \return  Current packing fraction of the configuration. */
	double getPackingFraction();
};

/** \brief For a given scaling function, starting from a certain particle size, fill
 *  particles; if there is no space for particles of a certain size, go to the next stage.
 *  \param c            Initial packing; the result accumulates here.
 *  \param numStages    Number of stages (particle sizes) to fill.
 *  \param maxPrtVol    Maximum (stage-1) particle volume.
 *  \param parameter    Scaling-function parameter (exponent/ratio).
 *  \param fnName       Scaling-function name ("Zeta" or "Power").
 *  \param prefix       Output filename prefix for intermediate configurations.
 *  \param outputIndex  Stage indices at which to write a configuration.
 *  \return             Difference between the achieved and the expected covering fraction. */
double getSpacefillingPacking(SpherePacking & c, int numStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, const std::vector<int> & outputIndex = std::vector<int> ());

/** \brief for (n = 1; n< 1+ numStages; n++){ r = r(n); PHI_RESIDUE = PHI_THEORY - PHI_CURR; N = N0 + PHI_RESIDUE/vol(r); c.AddSpheres(r(n), N); }
 *  \param c            Initial packing; the result accumulates here.
 *  \param numStages    Number of stages (particle sizes) to fill.
 *  \param maxPrtVol    Maximum (stage-1) particle volume.
 *  \param parameter    Scaling-function parameter (exponent/ratio).
 *  \param fnName       Scaling-function name ("Zeta" or "Power").
 *  \param prefix       Output filename prefix for intermediate configurations.
 *  \param randomSeed   Seed for the random generator.
 *  \param numThreads   Number of threads used for the insertion.
 *  \param outputIndex  Stage indices at which to write a configuration.
 *  \return             Difference between the achieved and the theoretical packing fraction. */
double getSpacefillingPacking2(SpherePacking & c, int numStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, int randomSeed, int numThreads=1, const std::vector<int> & outputIndex = std::vector<int>());
/** \brief Stage-ranged variant of getSpacefillingPacking2 running stages initStage..finStages. (CLAUDE)
 *  \param c            Initial packing; the result accumulates here.
 *  \param initStage    First stage to fill.
 *  \param finStages    Final stage to fill.
 *  \param maxPrtVol    Maximum (stage-1) particle volume.
 *  \param parameter    Scaling-function parameter (exponent/ratio).
 *  \param fnName       Scaling-function name ("Zeta" or "Power").
 *  \param prefix       Output filename prefix for intermediate configurations.
 *  \param randomSeed   Seed for the random generator.
 *  \param numThreads   Number of threads used for the insertion.
 *  \param outputIndex  Stage indices at which to write a configuration.
 *  \return             Difference between the achieved and the theoretical packing fraction. */
double getSpacefillingPacking2(SpherePacking & c, int initStage, int finStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, int randomSeed, int numThreads = 1, const std::vector<int> & outputIndex = std::vector<int>());

#endif
