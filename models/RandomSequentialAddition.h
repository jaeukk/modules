/**
 *	Author	: Ge Zhang
 *	Email	: - 
 *	Created Date	: 2014 
 */

#ifndef RandomSequentialAddition_Included
#define RandomSequentialAddition_Included 

#define MAX_MEM_BOX 40000000000 //Maximal memory space for voxel list 40GB

#include <vector>
#include "PeriodicCellList.h"
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
	/** \brief Build the list of relative offsets used to split a voxel into sub-voxels. (CLAUDE) */
	void getsplitvoxellist(void);

public:
	/** \brief Default constructor: build an empty voxel list. (CLAUDE) */
	VoxelList();
	/** \brief Destructor: release the internal index and split-offset arrays. (CLAUDE) */
	~VoxelList();
	/** \brief Construct a voxel list partitioning the unit box into a regular grid. (CLAUDE)
	 *  \param Dimension     Spatial dimension.
	 *  \param VoxelRank     Number of voxels per side of the unit box.
	 *  \param ExpectedSize  Expected number of voxels, used to reserve storage. */
	VoxelList(int Dimension, voxeltype VoxelRank, size_t ExpectedSize=10000);
	/** \brief Build the initial voxel list, keeping only voxels that do not overlap any sphere. (CLAUDE)
	 *  \param Config  Existing sphere packing the voxels must avoid.
	 *  \param Radius  Radius of the spheres to be inserted. */
	void GetOriginalList(const SpherePacking & Config, double Radius);
	/** \brief Number of voxels currently in the list. (CLAUDE)
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
	/** \brief Recursively test whether a voxel region is fully covered (overlapped) by the spheres. (CLAUDE)
	 *  \param Config      Current sphere packing.
	 *  \param Radius      Radius of the spheres to be inserted.
	 *  \param Center      Center coordinate of the voxel region being tested.
	 *  \param QuaterSize  Quarter of the current voxel edge length.
	 *  \param checklevel  Remaining recursion depth for the deep overlap check.
	 *  \return            true if the whole region is covered (no insertable point), false otherwise. */
	bool CheckOverlapDeep(const SpherePacking & Config, const double & Radius, const double * Center, double QuaterSize, unsigned int checklevel);
	/** \brief Split each voxel and keep insertable sub-voxels, serial version. (CLAUDE)
	 *  \param Config      Current sphere packing.
	 *  \param Radius      Radius of the spheres to be inserted.
	 *  \param checklevel  Recursion depth used in the deep overlap check. */
	void SplitVoxelSerial(const SpherePacking & Config, double Radius, int checklevel=0);
	/** \brief Remove voxels that are fully covered by the spheres without splitting. (CLAUDE)
	 *  \param Config      Current sphere packing.
	 *  \param Radius      Radius of the spheres to be inserted.
	 *  \param checklevel  Recursion depth used in the deep overlap check. */
	void DistillVoxelList(const SpherePacking & Config, double Radius, int checklevel=0);
	/** \brief Split each voxel and keep insertable sub-voxels, parallel version. (CLAUDE)
	 *  \param Config      Current sphere packing.
	 *  \param Radius      Radius of the spheres to be inserted.
	 *  \param checklevel  Recursion depth used in the deep overlap check. */
	void SplitVoxel(const SpherePacking & Config, double Radius, int checklevel=0);
	/** \brief Current voxel edge length. (CLAUDE)
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
	/** \brief Destructor: free the owned sphere packing and voxel list. (CLAUDE) */
	~RSA();

	/** \brief try to insert sphere, if less than NumInsertedLimit spheres inserted in TrialTime trials, stop.
	 *  \param NumInsertedLimit  Minimum number of insertions per round below which the loop stops.
	 *  \param TrialTime         Number of trial insertions per round. */
	void RSA_I(size_t NumInsertedLimit, size_t TrialTime);//try to insert sphere, if less than NumInsertedLimit spheres inserted in TrialTime trials, stop
	/** \brief try to insert sphere until there are NumCut spheres.
	 *  \param NumCut  Target number of spheres at which to stop. */
	void RSA_I(size_t NumCut);//try to insert sphere until there are NumCut spheres
	/** \brief Insert spheres into the voxel list serially until insertions fall below the limit. (CLAUDE)
	 *  \param NumInsertedLimit  Minimum number of insertions per round below which the loop stops.
	 *  \param TrialTime         Number of trial insertions per round. */
	void RSA_II_Serial( size_t NumInsertedLimit, size_t TrialTime);
	/** \brief Insert spheres into the voxel list in parallel until insertions fall below the limit. (CLAUDE)
	 *  \param NumInsertedLimit  Minimum number of insertions per round below which the loop stops.
	 *  \param TrialTime         Number of trial insertions per round. */
	void RSA_II_Parallel( size_t NumInsertedLimit, size_t TrialTime);
	/** \brief Insert spheres into the voxel list, dispatching to the serial or parallel routine. (CLAUDE)
	 *  \param NumInsertedLimit  Minimum number of insertions per round below which the loop stops.
	 *  \param TrialTime         Number of trial insertions per round.
	 *  \param parallel          If true, use the parallel insertion routine. */
	void RSA_II(size_t NumInsertedLimit, size_t TrialTime, bool parallel);
	/** \brief Packing fraction of the current configuration. (CLAUDE)
	 *  \return the volume fraction occupied by the spheres. */
	double Density(void) const;
	/** \brief get original voxel list.
	 *  \param VoxelDiameterRatio  Ratio of the sphere diameter to the initial voxel size. */
	void GetVoxelList(double VoxelDiameterRatio);//get original voxel list
	/** \brief Split the voxel list, dispatching to the serial or parallel routine. (CLAUDE)
	 *  \param parallel    If true, use the parallel splitting routine.
	 *  \param checklevel  Recursion depth used in the deep overlap check. */
	void SplitVoxel(bool parallel, int checklevel=0);
	/** \brief Write the sphere radius and configuration to an output stream. (CLAUDE)
	 *  \param out  Output stream to write to. */
	void OutputSpheres(std::ostream & out);
};

/** \brief return the packing with NumCut particles if it's achieved (implementations of the voxel-list method). (CLAUDE)
 *  \param dimension        Spatial dimension (negative to suppress writing the configuration to file).
 *  \param nbrsphere        Target number of spheres used to set the sphere radius.
 *  \param nbrinsertedlimit Insertions-per-round threshold below which a round stops (negative to run in parallel).
 *  \param trial1           Number of trial insertions per round in the first (voxel-free) stage.
 *  \param trial2           Number of trial insertions per round in the voxel-list stage.
 *  \param voxelratio       Ratio of the sphere diameter to the initial voxel size (negative to re-test at saturation).
 *  \param seed             Seed for the random number generator.
 *  \param output           Output stream for progress messages.
 *  \param volumeratio      If non-null, receives the single-sphere volume.
 *  \param Verbose          If true, print voxel-center diagnostics.
 *  \param NumCut           Stop and return once this many spheres are inserted (default: run to saturation).
 *  \return                 The generated sphere packing. */
SpherePacking GenerateRSAPacking(int dimension, unsigned long nbrsphere, short nbrinsertedlimit, unsigned long trial1, unsigned long trial2, double voxelratio, int seed, std::ostream & output, double * volumeratio=nullptr, bool Verbose=false, size_t NumCut=SIZE_MAX);


/** \brief Convenience wrapper that generates an RSA packing with default trial counts and voxel ratio. (CLAUDE)
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

/** \brief return the radius of a sphere at the saturation density.
 *  \param d           Spatial dimension.
 *  \param numDensity  Number density of the spheres.
 *  \return            the sphere radius at the saturation density. */
double Rc(DimensionType d, double numDensity);

#endif