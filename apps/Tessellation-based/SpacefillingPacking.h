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

	SpacefillingPacking(const SpherePacking& config);
	SpacefillingPacking(const SpherePacking& config, const std::vector<double> & param, const char * fnName);

	inline void ReplaceConfig(SpherePacking & newPacking) {
		this->pSpheres = new SpherePacking(newPacking);
		this->pSpheres->SetCellSize(2.0*this->maxR);
		this->phi = this->pSpheres->PackingFraction();
	}
	//Manual inputs
	void setN(int N);
	void setRmax(double Rmax);
	/* set the scaling function of radius
	param[in]	param[0] = max volume
				param[1] = parameter
	param[in]	fnName	"Zeta": f(n) = 1/n^param[1]
						"Power":f(n) = 1/param[1]^(n-1) */
	void setScalingFunction(const std::vector<double> & param, const char * fnName);
	void setScalingFunction(const std::function <double(int n)> & function);

	/*Add (numAdd) spheres of radius (r); updates this->phi */
	void AddSpheres(double r, int numAdd, int numThreads = 1);

	/*for (n = n0; n< n0+numSteps; n++) AddSpheres(r*f(n), N); */
	void RepeatAddition(int numSteps, int n0);
	/* Compute the packing fraction, when one can insert "N" particles at each
	   stage from 1~stage. */
	double GetTheoreticalPackingFraction(int stage);
	//Return the current packing fraction.
	double getPackingFraction();
};

//For a given scaling function, starting from a certain particle size, fill
//particles; if there is no space for particles of a certain size, go to the
//next stage.
double getSpacefillingPacking(SpherePacking & c, int numStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, const std::vector<int> & outputIndex = std::vector<int> ());

/* for (n = 1; n< 1+ numStages; n++){ r = r(n); PHI_RESIDUE = PHI_THEORY -
   PHI_CURR; N = N0 + PHI_RESIDUE/vol(r); c.AddSpheres(r(n), N); } */
double getSpacefillingPacking2(SpherePacking & c, int numStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, int randomSeed, int numThreads=1, const std::vector<int> & outputIndex = std::vector<int>());
double getSpacefillingPacking2(SpherePacking & c, int initStage, int finStages, double & maxPrtVol, double parameter, const char * fnName, const char * prefix, int randomSeed, int numThreads = 1, const std::vector<int> & outputIndex = std::vector<int>());

#endif
