#ifndef TORQUATOJIAO_INCLUDED
#define TORQUATOJIAO_INCLUDED

//Anything related to linear programming application on jamming hard-sphere packings

#include <etc.h>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <PeriodicCellList.h>
#include <exception>


#ifdef USE_TJ
class LPClass;
class TJJammer
{
public:
	class error : public std::exception
	{
	public:
		std::string message;
		/** \brief Construct the exception with a message string.
		 *  \param a  Null-terminated description of the error. */
		error(const char * a): message(a){}
		/** \brief Return the stored error message.
		 *  \return  C-string with the error description. */
		virtual const char* what() const throw()
		{
			return message.c_str();
		}
	};
	// input file variables - these are the defaults in case input isn't read!
	// See file "HopkinsTJ.standard.input" for details about all of these parameters
	std::string parentDirectory;
	std::string outputFilename;
	std::string tempFolder;
	bool   useNNL ;
	bool autoSave  ;//Whether to save progress in a temporary folder periodically
	int    N  ;                  //The number of spheres
	int    dim  ;                //The dimension of the problem
	bool   maxInitRad  ;     //If true, the radii of the spheres are maximized in the IC (so that the first contact is made)
	double transMax  ;     //The maximum translation of a sphere in a single LP step (lattice vectors?)
	double compMax  ;      //The maximum bulk (compressive) strain in a single step
	double shearMax  ;    //The maximum shear              strain in a single step
	int    strictJam  ;       //1 allows for full box deformations (as allowed by compMax and shearMax), but 0 resticts the compression to be shape-preserving (which will cause collective jamming instead of strict jamming)
	double delta  ;         //The influence sphere radius (?)
	int    maxIters  ;     //The maximum number of times the LP solver may be invoked
	std::string termCriterion;        //Either "maxDisp" or "latticeVol"
	double phiTerm  ;       //The maximum density allowed by the simulation (nice when we want to stop early)
	double termTol  ;    //The termination criterion: if the density changes by less than this, then we are done (?)
	//Alternatively, if the biggest sphere movement (relative to the center of gravity) was smaller than this, then we're done
	double MCStartSweeps  ; //How many MC sweeps to do before the SLP algorithm
	double MCMidSweeps  ;   //How many MC sweeps to do after each LP iteration          
	double mc_dispLimit;         //The displacement limit for a trial move in the MC routine
	double maxDisp;              //The farthest (relative to the total packing displacement) that a sphere has moved in an iteration
	bool   randLatVecs  ;   //When generating a packing, whether or not the lattice vectors should be random (otherwise, cube)
	double maxLatLength  ;  //The maximum length of a lattice vector (at initialization) (?)
	double minLatLength  ; //The minimum length of a lattice vector (at initialization)
	double minAngDeg  ;    //The minimum angle between vectors when starting the lattice.
	double manualMinVol  ; //
	double radLatTol ;
	int    topFail ;
	int    maxVolNum ;
	double inputPhi  ;        //The intial (starting) density
	int    overBoxes  ;          //How many rings (of periodic image cells) to check beyond the actual cell (0 = use L/2 method)
	int    maxBoxNum ;
	bool   randOverlapMove ;
	double resizeSpace ;
	double resizeTol ;
	double NNLextraDist ;
	int    printEvery  ;        //if >0, Write out the configuration to file after every [this many] iterations
	bool   printThisIter_pkg;
	bool   printThisIter_LP;
	int    printPrecision  ;   //Number of decimal places - outFile.precision(printPrecision);
	bool   useTimer  ;          //Whether to time the run (1) or not (0).  Units are microseconds
	bool   scaleUnitRadius  ;   //If yes, scale the packing so that the smalest sphere has radius = 1
	double rescaleIC  ;       //After you make/load an IC, multiply the size of the packing by this (to "rewind" packings)
	int    printLPEvery  ;      //if >0, Write out the LP's solution to file after every [this many] iterations


	//Solver variable:
	std::string solver  ;                //"GLPK" or "Gurobi"

	//GLPK solver params
	double feasibleTol ;
	double pivotTol  ;//Default value used by GLPK, 0<pt<1...increase to improve numerical accuracy...
	int GLPK_PresolveVar ;
	int GLPK_BasisFactType ;

	//Gurobi solver params
	std::string GRB_Method  ;//Which solver to use
	int    GRB_Presolve  ;   //How much presolving to do (-1,0,1,2) = (auto,none,standard,aggresive)
	int    grbThreads ;      //Limit how many threads Gurobi can use at once

	/** \brief Initialize the jammer's configuration defaults, owned pointers, and random-number generator. */
	TJJammer();

	/** \brief Release all owned heap resources (RNG, matrices, coordinate arrays, neighbor lists). */
	~TJJammer();

	/** \brief Copy construction is deleted: TJJammer owns raw heap resources with no deep-copy support. */
	TJJammer(const TJJammer &) = delete;
	/** \brief Copy assignment is deleted: TJJammer owns raw heap resources with no deep-copy support. */
	TJJammer & operator=(const TJJammer &) = delete;

protected:
	// Timer variables:
	class TimerClass 
	{
	public:
		double timePassed;

		/** \brief Record the current time as the timer start. */
		void Start();

		/** \brief Record the current time and store the elapsed interval. */
		void End();

		/** \brief Print the elapsed time to standard output. */
		void Display();
	private:
		time_t startTime;
		time_t endTime;
	};
	class NeighborNode 
	{
	public:
		NeighborNode *next;
		int index;

		/** \brief Iteratively delete the rest of the linked neighbor list (avoids deep recursion / stack overflow on long lists). */
		~NeighborNode()
		{
			NeighborNode *p = next;
			while (p != 0)
			{
				NeighborNode *following = p->next;
				p->next = 0;
				delete p;
				p = following;
			}
		}
	};
	// input file variables - these are the defaults in case input isn't read!
	const gsl_rng_type * RNGtype;
	gsl_rng *RNG;
	TimerClass timer; //See timer.h for structure
	TimerClass saveTimer; //To save the data periodically

	// running variables
	gsl_matrix *lambdas;
	gsl_matrix *inverseLambdas;
	double *localCoords ;
	double *globalCoords ;
	double *radii ;
	double *distTempL ;
	double *distTempG ;
	NeighborNode **neighborListDelta ;
	NeighborNode **neighborListOverlap ;
	int lpIters  ;//The number of LP iterations that have gone by
	double maxSingleMove ;
	double biggestRad ;
	/** \brief Return the Euclidean length of the displacement stored in distTempL after mapping it through lambdas.
	 * \param inputVec Unused by the current implementation.
	 * \return Length of lambdas*distTempL in global Cartesian coordinates. */
	double getGlobalLength(double *inputVec) ;

	/** \brief Volume of the fundamental cell (|det Lambda|).
	 *  \return  Absolute value of the lattice determinant. */
	double getVol();

	/** \brief Rescale the lattice to remove any sphere overlaps, optionally jittering positions.
	 *  \param randOver  If true, leave a small gap and apply random stochastic moves after resizing.
	 *  \return          1 if a resize was performed, 0 otherwise. */
	int resizeIfOverlap(bool randOver) ;

	/** \brief Total volume occupied by all spheres.
	 *  \return  Sum of the d-dimensional volumes of every sphere. */
	double getSphereVol();

	/** \brief Scale lattice and coordinates by rescaleIC to "rewind" a packing. */
	void RescalePacking();

	/** \brief Free and reset both near-neighbor lists. */
	void DeleteOldNNL();

	/** \brief Release the per-packing heap state (matrices, coordinate arrays, neighbor lists) and null the pointers; idempotent. */
	void FreeState();

	/** \brief Rebuild the influence-sphere and overlap near-neighbor lists.
	 *  \return  Total number of neighbor contacts (counted once per pair). */
	int calcNNLs() ;
	double InternalTransMax  ;     //The maximum translation of a sphere in a single LP step (lattice vectors?)
	double InternalCompMax  ;      //The maximum bulk (compressive) strain in a single step
	double InternalShearMax  ;    //The maximum shear              strain in a single step
	double InternalDelta;
	double InternalNNLextraDist;
	/** \brief Uniformly scale all radii by the largest candidate expansion factor computed over the examined sphere/image pairs. */
	void MaximizeRadii();

	/** \brief Adapt the MC trial-move displacement limit toward a target acceptance rate.
	 *  \param accRate  Measured acceptance rate of the last MC sweep. */
	void MC_AdjustDispLimit(double accRate);

	/** \brief Attempt a single random Monte Carlo move of one sphere, rejecting overlaps.
	 *  \param i  Index of the sphere to move.
	 *  \return   Magnitude of the accepted displacement, or 0 if the move was rejected. */
	double MC_Move(int i);

	/** \brief Perform a sweep of Monte Carlo moves and adapt the displacement limit.
	 *  \param numTries  Number of trial moves to attempt.
	 *  \return          Largest accepted displacement during the sweep. */
	double MC_Sweep(int numTries);

	/** \brief Initialize the Monte Carlo displacement limit. */
	void MC_Init();

	/** \brief Write the current packing to a .dat file in Alek's format.
	 *  \param thisOutputFilename  Output path without extension (".dat" is appended).
	 *  \param addIterNum          If true, append the LP iteration number to the filename.
	 *  \return                    1 on success, 0 if the file could not be created. */
	int printPacking(std::string thisOutputFilename, bool addIterNum);

	/** \brief Decide whether the SLP iterations should stop.
	 *  \param lastLastVol  Lattice volume two iterations ago (for the latticeVol criterion).
	 *  \param latVol       Current lattice volume.
	 *  \param phi          Current packing fraction.
	 *  \return             true if a termination criterion is met, false otherwise. */
	bool CheckTermination(double lastLastVol, double latVol, double phi);

	/** \brief Append the timer's elapsed time to TIMEFILE.txt.
	 *  \param timer  Timer whose elapsed time is written. */
	void PrintTime(TimerClass timer);

	/** \brief Delete the temporary auto-save directory and its contents. */
	void ClearSaveData();

	/** \brief Save the current packing and an iteration log to the temporary folder. */
	void SaveProgress();
	SpherePacking result;
public:
	/** \brief Run the sequential-linear-programming jamming iterations until a stop criterion is met. */
	void TJIterations(void);//run iterations until stop criteria is met
	/** \brief Load an initial sphere packing into the jammer's internal state.
	 *  \param p  Sphere packing to use as the initial configuration. */
	void SetPacking(const SpherePacking & p);

	/** \brief Return the jammed packing produced by the iterations.
	 *  \return  The resulting sphere packing. */
	SpherePacking GetPacking(void);
	friend class LPClass;
};

//find if a sphere is a ratteler or not given it has these contacts to non-ratteler spheres
//need to loop through all spheres until no extra rattler are found
//This function used GLPK, which is NOT thread-safe
/** \brief Determine whether a sphere is a rattler given its contacts to non-rattlers.
 *  \param dim       Spatial dimension.
 *  \param contacts  Contact vectors from the sphere to its non-rattler neighbors.
 *  \return          true if the sphere is a rattler (under-constrained), false otherwise. */
bool IsRattler_LP(DimensionType dim, const std::vector<GeometryVector> & contacts);

#else
class TJJammer
{
public:
	// input file variables - these are the defaults in case input isn't read!
	// See file "HopkinsTJ.standard.input" for details about all of these parameters
	std::string parentDirectory;
	std::string outputFilename;
	std::string tempFolder;
	bool   useNNL ;
	bool autoSave  ;//Whether to save progress in a temporary folder periodically
	int    N  ;                  //The number of spheres
	int    dim  ;                //The dimension of the problem
	bool   maxInitRad  ;     //If true, the radii of the spheres are maximized in the IC (so that the first contact is made)
	double transMax  ;     //The maximum translation of a sphere in a single LP step (lattice vectors?)
	double compMax  ;      //The maximum bulk (compressive) strain in a single step
	double shearMax  ;    //The maximum shear              strain in a single step
	int    strictJam  ;       //1 allows for full box deformations (as allowed by compMax and shearMax), but 0 resticts the compression to be shape-preserving (which will cause collective jamming instead of strict jamming)
	double delta  ;         //The influence sphere radius (?)
	int    maxIters  ;     //The maximum number of times the LP solver may be invoked
	std::string termCriterion;        //Either "maxDisp" or "latticeVol"
	double phiTerm  ;       //The maximum density allowed by the simulation (nice when we want to stop early)
	double termTol  ;    //The termination criterion: if the density changes by less than this, then we are done (?)
	//Alternatively, if the biggest sphere movement (relative to the center of gravity) was smaller than this, then we're done
	double MCStartSweeps  ; //How many MC sweeps to do before the SLP algorithm
	double MCMidSweeps  ;   //How many MC sweeps to do after each LP iteration          
	double mc_dispLimit;         //The displacement limit for a trial move in the MC routine
	double maxDisp;              //The farthest (relative to the total packing displacement) that a sphere has moved in an iteration
	bool   randLatVecs  ;   //When generating a packing, whether or not the lattice vectors should be random (otherwise, cube)
	double maxLatLength  ;  //The maximum length of a lattice vector (at initialization) (?)
	double minLatLength  ; //The minimum length of a lattice vector (at initialization)
	double minAngDeg  ;    //The minimum angle between vectors when starting the lattice.
	double manualMinVol  ; //
	double radLatTol ;
	int    topFail ;
	int    maxVolNum ;
	double inputPhi  ;        //The intial (starting) density
	int    overBoxes  ;          //How many rings (of periodic image cells) to check beyond the actual cell (0 = use L/2 method)
	int    maxBoxNum ;
	bool   randOverlapMove ;
	double resizeSpace ;
	double resizeTol ;
	double NNLextraDist ;
	int    printEvery  ;        //if >0, Write out the configuration to file after every [this many] iterations
	bool   printThisIter_pkg;
	bool   printThisIter_LP;
	int    printPrecision  ;   //Number of decimal places - outFile.precision(printPrecision);
	bool   useTimer  ;          //Whether to time the run (1) or not (0).  Units are microseconds
	bool   scaleUnitRadius  ;   //If yes, scale the packing so that the smalest sphere has radius = 1
	double rescaleIC  ;       //After you make/load an IC, multiply the size of the packing by this (to "rewind" packings)
	int    printLPEvery  ;      //if >0, Write out the LP's solution to file after every [this many] iterations


	//Solver variable:
	std::string solver  ;                //"GLPK" or "Gurobi"

	//GLPK solver params
	double feasibleTol ;
	double pivotTol  ;//Default value used by GLPK, 0<pt<1...increase to improve numerical accuracy...
	int GLPK_PresolveVar ;
	int GLPK_BasisFactType ;

	//Gurobi solver params
	std::string GRB_Method  ;//Which solver to use
	int    GRB_Presolve  ;   //How much presolving to do (-1,0,1,2) = (auto,none,standard,aggresive)
	int    grbThreads ;      //Limit how many threads Gurobi can use at once

	/** \brief Stub constructor that reports TJ is disabled. */
	TJJammer()
	{
		std::cerr<<"Error in TJJammer : TJ is not enabled!\n";
	}

protected:
public:
	/** \brief Stub: TJ is disabled; reports an error. */
	void TJIterations(void)//run iterations until stop criteria is met
	{
		std::cerr<<"Error in TJJammer : TJ is not enabled!\n";
	}

	/** \brief Stub: TJ is disabled; reports an error.
	 *  \param p  Ignored sphere packing. */
	void SetPacking(const SpherePacking &)
	{
		std::cerr<<"Error in TJJammer : TJ is not enabled!\n";
	}

	/** \brief Stub: TJ is disabled; reports an error and returns an empty packing.
	 *  \return  An empty SpherePacking. */
	SpherePacking GetPacking(void)
	{
		std::cerr<<"Error in TJJammer : TJ is not enabled!\n";
		return SpherePacking();
	}
};
/** \brief Stub: TJ is disabled; reports an error and returns false.
 *  \param dim       Ignored spatial dimension.
 *  \param contacts  Ignored contact vectors.
 *  \return          Always false. */
inline bool IsRattler_LP(DimensionType, const std::vector<GeometryVector> &)
{
	std::cerr << "Error in IsRattler_LP : TJ is not enabled!\n"; 
		return false;
}
#endif

#endif