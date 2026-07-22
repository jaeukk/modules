/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file GenerateConfigs.h
 * \brief Header file for a ConfigGen class and its derived class.
 * Collections of point configurations and spherepackings that you can simply construct.
 * They can be used as initial configurations.
	//What can we generate?
	----- Lattices/ (Lattice.h)
	** lattice packing
	** Perturbed lattices (Pert/GenPerturb.h)
	** Imperfect lattices with vacancies (Pert/GenPerturb.h)
	** Paracrystals (in progress)

	----- RSA/ (RSA_Gen.h)
	** Saturated RSA (sphere packings)
	----- HSF/ (HSF_Gen.h)
	** equilibrium hard-sphere fluids
	----- DPP/
	** Determinantal point processes (TODO)
	----- Jamming
	** Dense Packing via LS algorithm (TODO)

 * The per-model generator classes were split into subdir headers during the
 * reorganization; this header keeps the ConfigGen base, PoissonPtGen, the
 * Jammed_Gen stub, the predefined-lattice helpers, Echo, and numK_points, and
 * re-includes the split headers at the bottom so that consumers which only
 * #include <GenerateConfigs.h> keep seeing LatticeGen / GenPerturb / GenPtsDefect.
 */

#ifndef GENERATECONFIGS
#define GENERATECONFIGS

#include <omp.h>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <functional>
#include <cstring>
#include <cctype>
#include <cmath>
#include <ctime>
#include <memory>

#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
#include <etc.h>


/** \brief Give basis vectors of well-known lattices
 *  \param[in] latticeType	The name of Bravais lattice, e.g.,
								1D: Integer
								2D: Square, Triangle,
								3D: Cubic, BCC_cubic, FCC, FCC_cubic	 BCC_symmetric,
 *  \param[out] basis	Basis vectors of corresponding Bravais lattice.
*/
void GetPredefinedBasis(const std::string & latticeType, std::vector<GeometryVector> & basis);

/** \brief Give basis vectors and positions of atoms for some predefined lattices (including some Bravais lattices).
 *  \param[in] crystalName It includes some non-Bravais lattices, such as
 *				2D: Honeycomb, Kagome
 *				3D: Diamond.
 *	\param[out] basis	Basis vectors
 *  \param[out] atoms	Relative coordinates of particles (except for one at (0,0,0)). */
void GetPredefinedCrystals(const std::string & crystalName, std::vector<GeometryVector> & basis, std::vector<GeometryVector> & atoms);

/** \brief Compute the number of k-points in spectral densities.
 *  \param FundamentalCell  Fundamental cell whose reciprocal lattice defines the k-points.
 *  \param ks               Output list of k-points.
 *  \return                 the number of k-points. */
size_t numK_points(PeriodicCellList<Empty> & FundamentalCell, std::vector<GeometryVector> & ks);

/** \brief Read a value from an input stream and echo it to an output stream.
 *  \tparam T        Type of the value to read and echo.
 *  \param ifile     Input stream to read from.
 *  \param ofile     Output stream to echo to.
 *  \param variable  Receives the value read from ifile. */
template<typename T> inline void Echo(std::istream & ifile, std::ostream & ofile, T & variable) {
	ifile >> variable;	ofile << variable << "\n";
}

/*
	Callable types to generate configurations
	Instructions to the Constructors
	Ask parameters in the following orders
	dimension -> intensive parameters -> extensive parameters
*/

/** \brief  Abstract class to define a point process.
 */
class ConfigGen
{
protected:
	double rho;		//< number density
	size_t numPrts,	//< number of particles
		numThreads;	//< the number of threads that are used to construct a system.
	DimensionType d;//< space dimension

public:

	/** A default constructor of ConfigGen. */
	ConfigGen() { this->rho = 1.0; this->numPrts = 0; this->d = 1; this->numThreads = 1; };

	/** \brief Fully-parametrized constructor.
	 *  \param dimensions  Spatial dimension.
	 *  \param numPrts     Number of particles.
	 *  \param rho         Number density. */
	ConfigGen(int dimensions, int numPrts, double rho);

	/** \brief Construct a point process interactively from streams.
	 *  Reads dimension, particle number, and number density in that order.
	 *  \param ifile  Input stream supplying the parameters.
	 *  \param ofile  Output stream for prompts and echoes. */
	ConfigGen(std::istream & ifile, std::ostream & ofile);

	/** \brief A pure virtual function to construct a point configuration.
	 *  \param[out] c	A point configuration. */
	virtual void GenerateC(Configuration & c) = 0;

	/** \brief A pure virtual function to construct a sphere packing.
	 *  \param[out] c	A sphere packing. */
	virtual void GenerateP(SpherePacking & c) = 0;

	//Change system parameters
	/** \brief A virtual member function to change number density (must be redefined).
	 *  \param rho	new number density. */
	virtual void SetRho(double rho) { (void)rho; std::cerr << "SetRho() is undefined" << std::endl; };

	/** \brief A virtual member function to change particle number (must be redefined).
	 *  \param num	new particle numbers. */
	virtual void SetNumPrts(int num) { (void)num; std::cerr << "SetNumPrts() is undefined" << std::endl; };

	/** \brief A virtual member function to change the packing fraction (must be redefined).
	 *  \param phi	new packing fraction. */
	virtual void SetPackingFraction(double phi) { (void)phi; std::cerr << "SetPackingFraction() is undefined" << std::endl; };	//must be redefined

	/** \brief Get the number density.
	 *  \return the number density.*/
	virtual double GetRho() { return this->rho; };

	/** \brief Get the particle number.
	 *  \return the particle number.*/
	virtual size_t GetNumPrts() { return this->numPrts; };

	/** \brief Redefin the number of threads that are used to construct realizations.
	 *  \param numThreads  new number of threads. */
	virtual void SetNumThreads(int numThreads) {
		if (numThreads > 0)
			this->numThreads = (size_t)numThreads;
		else
			std::cerr << "SetNumThreads: number of threads must be positive." << std::endl;
	};

	/** \brief A pure virutal member function to input extra parameters.
	 *  \param[in] option	the number of option.
	 *  \param[in] input		An input stream.
	 *  \param[out] output	An output stream.
	 */
	virtual void SetAdditionalOption(const std::string & option, std::istream & input, std::ostream &output) { (void)option; (void)input; output << "Uncrecognized command!\n"; /*by default, accept no additional option*/ }

	/** \brief Get the spatial dimension.
	 *  \return the spatial dimension. */
	DimensionType GetDimension() {
		return this->d;
	}

	/** A default destructor. */
	virtual ~ConfigGen() {}
};

///* Generate an empty periodic fundamental cell */
//struct EmptyBox : public ConfigGen
//{
//protected:
//	Configuration F_Cell;
//	double Cell_Volume;
//	virtual void SetBasis(const std::vector<GeometryVector> & basis);
//public:
//	/*** Constructors **********
//	 * Default Constructor
//	 * Fully parametrized
//	 * Name-parametrized
//	 * CMD interfaces
//	 * Copy
//	 ***************************/
//	EmptyBox();
//	EmptyBox(int dimension, const std::vector<GeometryVector> & BasisVectors);
//	EmptyBox(const std::string & CellName, double CellVolume);
//	EmptyBox(std::istream & ifile, std::ostream & ofile);
//	EmptyBox(const EmptyBox & copy);
//
//	void SetCellVolume(double new_cell_volume) { this->F_Cell.Resize(new_cell_volume); };
//
//	virtual void GenerateC(Configuration & c) { c = Configuration(this->F_Cell); };
//	virtual void GenerateP(SpherePacking & c) { c = SpherePacking(this->F_Cell, 0.0); };
//};


/** \breif A class to construct a simple Poisson point pattern. */
class PoissonPtGen : public ConfigGen {
protected:
	int mode = 0;	//<	0: particle number is fixed
					//< 1: particle number can vary
	RandomGenerator rng;	//< a random number generator

	/** \brief Determine the particle count for a cell of the given volume.
	 *  \param cell_volume Volume of the fundamental cell.
	 *  \return The rounded mean count in fixed-count mode, a Poisson-distributed draw in Poisson mode, or zero for an invalid mode. */
	inline size_t GetExpectedNumPrts(double cell_volume) {
		double mean = cell_volume * this->rho;
		if (mode == 0) {
			return (int)round(mean);
		}
		else if (mode == 1) {
			return rng.PoissonRNG(mean);
		}
		else{
			std::cerr<<"PoissonPtGen::invalid mode!\n";
			return 0;
		}
	}
public:
	/** \brief A default constructor.
	 *  \param rho			The average number density.
	 *  \param mode		0=fixed particle number/ 1=follows Poisson distribution.
	 *  \param seed			A seed for a random numnber generator.	 */
	PoissonPtGen(double rho, int mode = 0, int seed = 0) : ConfigGen() {
		this->rho = rho;
		if (mode == 0 || mode == 1) {
			this->mode = mode;
		}
		else { std::cout << "PoissonPtGen::undefined mode!\n";	return; }
		this->rng.seed(seed);
	}

	/** \brief A constructor via an input stream.
	 *  \param ifile  Input stream supplying the parameters.
	 *  \param ofile  Output stream for prompts and echoes. */
	PoissonPtGen(std::istream & ifile, std::ostream & ofile);

	/** \brief Reset the number density.
	 *  \param rho  new number density. */
	virtual void SetRho(double rho) { this->rho = rho; }

	/** \brief Redefine the mode.
	 *  \param mode  0=fixed particle number / 1=Poisson-distributed particle number. */
	void SetMode(int mode) { this->mode = mode; }

	/** \brief Empty the simulation box and generate a Poisson point pattern of a given number density.
	 *	\param[in,out] c	A periodic simulation box whose volume should be positive!	 */
	virtual void GenerateC(Configuration & c) {
		double volume = c.PeriodicVolume();

		if (volume > 0.0) {
			c.RemoveParticles();
			size_t N = this->GetExpectedNumPrts(volume);
			for (size_t i = 0; i < N; i++) {
				//Randomly place a single point inside the fundamental cell.
				c.Insert("a", this->rng);
			}
		}

	}


	/** \brief Undefined member functions !!
	 *  \param num  Ignored. */
	virtual void SetNumPrts(int num) { (void)num; std::cout << "This is not used\n"; }
	/** \brief Unused override: packing fraction is meaningless for a Poisson point pattern.
	 *  \param phi  Ignored. */
	virtual void SetPackingFraction(double phi) { (void)phi; std::cout << "PoissonPtgen::packing fraction is unnecessary\n"; }
	/** \brief Unused override: sphere-packing generation is not defined for Poisson points.
	 *  \param[out] c  Ignored. */
	virtual void GenerateP(SpherePacking & c) { (void)c; std::cout << "PoissonPtGen::GenerateP is not defined.\n"; }

};



//Jammed packing via the LS algorithms
#ifndef USE_LS
class Jammed_Gen : public ConfigGen {


};
#endif // !USE_LS


// -----------------------------------------------------------------------------
// Backward-compatible facade.  The per-model generator classes that used to live
// in this header now live in per-model subdir headers.  Re-include them here (after
// ConfigGen and the helpers are defined) so existing consumers that only
// #include <GenerateConfigs.h> keep compiling unchanged.  The include guards make
// the reverse include (Lattice.h / GenPerturb.h -> GenerateConfigs.h) a no-op.
// -----------------------------------------------------------------------------
#include "Lattices/Lattice.h"
#include "Pert/GenPerturb.h"


#endif
