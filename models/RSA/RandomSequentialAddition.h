/**
 *	Author	: Ge Zhang
 *	Email	: - 
 *	Created Date	: 2014 
 */

#ifndef RandomSequentialAddition_Included
#define RandomSequentialAddition_Included 

#define MAX_MEM_BOX 40000000000 //Maximal memory space for voxel list 40GB

#include <vector>
#include <PeriodicCellList.h>
#include <cstdint>


class VoxelList
{
private:
	typedef unsigned short voxeltype;
	int dimension, level;
	double halfsize, originhalfsize;
	std::vector<voxeltype> * originindex;
	std::vector<unsigned short> * subindex;
	unsigned short nbrsplitvoxel;
	signed char * splitvoxellist;

	/** \brief convert double * to GeometryVector.
	 *  \param coord  Pointer to an array of `dimension` coordinates.
	 *  \return       The coordinates wrapped in a GeometryVector. */
	GeometryVector Coord2Vect(const double * coord);

	/** \brief Build the list of relative offsets used to split a voxel into sub-voxels. */
	void getsplitvoxellist(void);

public:
	/** \brief Default-initialize a voxel list with no allocated grid or index storage. */
	VoxelList();

	/** \brief Destructor: release the internal index and split-offset arrays. */
	~VoxelList();

	//Owns raw arrays; copying would alias them and double-free. Copying is disabled.
	VoxelList(const VoxelList &) = delete;
	VoxelList & operator=(const VoxelList &) = delete;

	/** \brief Initialize storage for a regular voxel grid over the unit box.
	 *  The voxel indices are populated later by GetOriginalList().
	 *  \param Dimension Spatial dimension.
	 *  \param VoxelRank Number of voxels per side; must be positive.
	 *  \param ExpectedSize Expected number of retained voxels, used to reserve storage. */
	VoxelList(int Dimension, voxeltype VoxelRank, size_t ExpectedSize=10000);

	/** \brief Populate the initial list with grid voxels not completely covered by an exclusion sphere.
	 *  Partially covered voxels are retained because they may contain an insertable point.
	 *  \param Config Existing sphere packing the new sphere must avoid.
	 *  \param Radius Radius of the spheres to be inserted. */
	void GetOriginalList(const SpherePacking & Config, double Radius);

	/** \brief Number of voxels currently in the list.
	 *  \return the voxel count. */
	unsigned long NumVoxel(void) const;

	/** \brief get the center coordinate of the nbr-th voxel.
	 *  \param[out] result  Array (size dimension) receiving the voxel-center coordinates.
	 *  \param nbr          Index of the voxel. */
	void GetCenter(double * result, size_t nbr) const;//get the center coordinate of the nbr-th voxel
	/** \brief generate a random coordinate which is inside a random voxel, write it into result.
	 *  \param[out] result  Array (size dimension) receiving the random coordinate.
	 *  \param gen          Random number generator. */
	void GetRandomCoordinate(double * result, RandomGenerator & gen) const;//generate a random coordinate which is inside a random voxel, write it into result
	/** \brief Recursively test whether a voxel region is fully covered (overlapped) by the spheres.
	 *  \param Config      Current sphere packing.
	 *  \param Radius      Radius of the spheres to be inserted.
	 *  \param Center      Center coordinate of the voxel region being tested.
	 *  \param QuaterSize  Quarter of the current voxel edge length.
	 *  \param checklevel  Remaining recursion depth for the deep overlap check.
	 *  \return            true if the whole region is covered (no insertable point), false otherwise. */
	bool CheckOverlapDeep(const SpherePacking & Config, const double & Radius, const double * Center, double QuaterSize, unsigned int checklevel);

	/** \brief Split each voxel and keep insertable sub-voxels, serial version.
	 *  \param Config      Current sphere packing.
	 *  \param Radius      Radius of the spheres to be inserted.
	 *  \param checklevel  Recursion depth used in the deep overlap check. */
	void SplitVoxelSerial(const SpherePacking & Config, double Radius, int checklevel=0);

	/** \brief Remove voxels that are fully covered by the spheres without splitting.
	 *  \param Config      Current sphere packing.
	 *  \param Radius      Radius of the spheres to be inserted.
	 *  \param checklevel  Recursion depth used in the deep overlap check. */
	void DistillVoxelList(const SpherePacking & Config, double Radius, int checklevel=0);

	/** \brief Split each voxel and keep insertable sub-voxels, parallel version.
	 *  \param Config      Current sphere packing.
	 *  \param Radius      Radius of the spheres to be inserted.
	 *  \param checklevel  Recursion depth used in the deep overlap check. */
	void SplitVoxel(const SpherePacking & Config, double Radius, int checklevel=0);

	/** \brief Current voxel edge length.
	 *  \return the voxel edge length. */
	double GetVoxelSize() const;
};

class RSA
{
public:
	SpherePacking * pSpheres;
	VoxelList * pVoxels;
	int dimension;
	RandomGenerator gen;
	//double diameter;
	double radius;
	double spherevolume;//the volume of a sphere

	/** \brief convert double * to GeometryVector.
	 *  \param coord  Pointer to an array of `dimension` coordinates.
	 *  \return       The coordinates wrapped in a GeometryVector. */
	GeometryVector Coord2Vect(const double * coord);

	/** \brief initialize square box with side length 1, sphere radii are set so that around NumSphere spheres exist in saturation limit.
	 *  \param Dimension  Spatial dimension.
	 *  \param NumSphere  Target number of spheres at the saturation limit. */
	RSA(int Dimension, unsigned long NumSphere);

	/** \brief Initialize from an existing packing; upcoming spheres have one radii.
	 *  pa should support spheres of different sizes, but this is not tested.
	 *  \param pa      initial configuration.
	 *  \param radius  Radius of the spheres to be inserted. */
	RSA(const SpherePacking & pa, double radius);

	/** \brief Destructor: free the owned sphere packing and voxel list. */
	~RSA();

	//Owns pSpheres/pVoxels; copying would alias them and double-free. Copying is disabled.
	RSA(const RSA &) = delete;
	RSA & operator=(const RSA &) = delete;

	/** \brief try to insert sphere, if less than NumInsertedLimit spheres inserted in TrialTime trials, stop.
	 *  \param NumInsertedLimit  Minimum number of insertions per round below which the loop stops.
	 *  \param TrialTime         Number of trial insertions per round. */
	void RSA_I(size_t NumInsertedLimit, size_t TrialTime);//try to insert sphere, if less than NumInsertedLimit spheres inserted in TrialTime trials, stop
	/** \brief try to insert sphere until there are NumCut spheres.
	 *  \param NumCut  Target number of spheres at which to stop. */
	void RSA_I(size_t NumCut);//try to insert sphere until there are NumCut spheres
	/** \brief Perform serial voxel-list insertion rounds until a round inserts no more than the limit.
	 *  \param NumInsertedLimit Maximum accepted insertion count for the terminating round.
	 *  \param TrialTime Number of trial insertions per round. */
	void RSA_II_Serial( size_t NumInsertedLimit, size_t TrialTime);

	/** \brief Generate voxel-list candidates in parallel and commit them serially until a round inserts no more than the limit.
	 *  \param NumInsertedLimit Maximum accepted insertion count for the terminating round.
	 *  \param TrialTime Number of candidate trials per round. */
	void RSA_II_Parallel( size_t NumInsertedLimit, size_t TrialTime);

	/** \brief Run voxel-list insertion rounds using the serial or parallel-candidate implementation.
	 *  Processing stops when a round inserts no more than the specified limit.
	 *  \param NumInsertedLimit Maximum accepted insertion count for the terminating round.
	 *  \param TrialTime Number of trials per round.
	 *  \param parallel If true, generate candidates in parallel before serial commitment. */
	void RSA_II(size_t NumInsertedLimit, size_t TrialTime, bool parallel);

	/** \brief Return the particle count times the volume of a sphere of the configured insertion radius.
	 *  This equals the packing fraction only for a unit-volume cell containing spheres of that radius.
	 *  \return The particle count multiplied by the configured single-sphere volume. */
	double Density(void) const;

	/** \brief get original voxel list.
	 *  \param VoxelDiameterRatio  Ratio of the sphere diameter to the initial voxel size. */
	void GetVoxelList(double VoxelDiameterRatio);//get original voxel list
	/** \brief Split the voxel list, dispatching to the serial or parallel routine.
	 *  \param parallel    If true, use the parallel splitting routine.
	 *  \param checklevel  Recursion depth used in the deep overlap check. */
	void SplitVoxel(bool parallel, int checklevel=0);

	/** \brief Write the sphere radius and configuration to an output stream.
	 *  \param out  Output stream to write to. */
	void OutputSpheres(std::ostream & out);
};

/** \brief Generate an RSA packing, returning exactly NumCut particles if that count is reached or a saturated packing otherwise.
 *  \param dimension Spatial dimension; a negative value suppresses writing SphereCoordinate.txt.
 *  \param nbrsphere Target saturated particle count used to choose the sphere radius.
 *  \param nbrinsertedlimit Terminating-round insertion limit; a negative value selects parallel candidate generation.
 *  \param trial1 Number of trial insertions per round in the initial voxel-free stage.
 *  \param trial2 Number of candidate trials per round in the voxel-list stage.
 *  \param voxelratio Approximate ratio of initial voxel edge length to sphere diameter; a negative value also re-tests voxel coverage after saturation.
 *  \param seed Seed for the random number generator.
 *  \param output Output stream for progress messages.
 *  \param volumeratio If non-null, receives the single-sphere volume when the saturation path completes; it is not written on an early NumCut return.
 *  \param Verbose If true, print voxel-center diagnostics.
 *  \param NumCut Return once this many spheres are inserted; the default runs to saturation.
 *  \return The generated sphere packing. */
SpherePacking GenerateRSAPacking(int dimension, unsigned long nbrsphere, short nbrinsertedlimit, unsigned long trial1, unsigned long trial2, double voxelratio, int seed, std::ostream & output, double * volumeratio=nullptr, bool Verbose=false, size_t NumCut=SIZE_MAX);


/** \brief Convenience wrapper that generates an RSA packing with default trial counts and voxel ratio.
 *  \param dimension  Spatial dimension.
 *  \param nbrsphere  Target number of spheres used to set the sphere radius.
 *  \param parallel   If true, run the insertion routines in parallel.
 *  \param NumCut     Stop and return once this many spheres are inserted (default: run to saturation).
 *  \return           The generated sphere packing. */
inline SpherePacking GenerateRSAPacking(int dimension, unsigned long nbrsphere, bool parallel=false, size_t NumCut=SIZE_MAX)
{
	unsigned long trial = 3 * std::pow(10, dimension);//3000;//temporary changes: 3*std::pow(10, dimension) -> 3000;
	return GenerateRSAPacking((-1)*std::abs(dimension), nbrsphere, parallel?(-3):3, trial, trial, 0.49, GetPreciseClock(), Verbosity>5?std::cout:logfile, nullptr, false, NumCut);
}

/** \brief Saturated RSA packing fraction (the single canonical table shared by PhiC and Rc).
 *  Reference: G. Zhang and S. Torquato, Phys. Rev. E 88, 053312 (2013).
 *  The 1D value is Renyi's car-parking constant 0.7476.
 *  \param d  Spatial dimension.
 *  \return   the saturated packing fraction. */
inline double SaturatedPackingFraction(DimensionType d) {
	switch (d) {
	case 1:	return 0.7476;	//Renyi's constant (1D random sequential adsorption)
	case 2:	return 0.54707;
	case 3:	return 0.38413;
	default:
		std::cout << "SaturatedPackingFraction: out-of-range dimension " << (size_t)d << "\n";
		return 0.2;
	}
}

/** \brief return the radius of a sphere at the saturation density.
 *  \param d           Spatial dimension.
 *  \param numDensity  Number density of the spheres.
 *  \return            the sphere radius at the saturation density. */
double Rc(DimensionType d, double numDensity);

#endif