/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file GenerateConfigs.cpp
 *	\brief	Implementations of functions declared in GenerateConfigs and derived classes. */


#include "GenerateConfigs.h"

#include <cstring>
#include <cctype>
#include <cmath>


/** \brief Fully-parametrized ConfigGen constructor.
 *  \param dimensions  Spatial dimension.
 *  \param numPrts     Number of particles.
 *  \param rho         Number density. */
ConfigGen::ConfigGen(int dimensions, int numPrts, double rho)
{
	this->d = dimensions;
	this->numPrts = numPrts;
	this->rho = rho;
	this->numThreads = 1;
}

/** \brief Construct a ConfigGen interactively from streams (dimension, particle number, number density).
 *  Validates the stream state and ranges (1 <= dimension <= MaxDimension, positive particle
 *  number, finite positive density); on invalid input it prints a message and leaves the
 *  generator in a safe empty state (dimension 0) rather than using indeterminate values. */
ConfigGen::ConfigGen(std::istream &ifile, std::ostream &ofile)
{
	this->numThreads = 1;
	this->d = 0; this->numPrts = 0; this->rho = 1.0;
	int dim = 0; long long np = 0; double density = 0.0;
	ofile << "Dimensions = ";		if (!(ifile >> dim)) { ofile << "ConfigGen: invalid dimension input.\n"; return; }		ofile << dim << "\n";
	ofile << "Particle Number = ";	if (!(ifile >> np)) { ofile << "ConfigGen: invalid particle-number input.\n"; return; }	ofile << np << "\n";
	ofile << "Number Density = ";	if (!(ifile >> density)) { ofile << "ConfigGen: invalid number-density input.\n"; return; }	ofile << density << "\n";
	if (dim < 1 || dim > (int)MaxDimension) { ofile << "ConfigGen: dimension must be in [1, " << (int)MaxDimension << "].\n"; return; }
	if (np <= 0) { ofile << "ConfigGen: particle number must be positive.\n"; return; }
	if (np > 100000000LL) { ofile << "ConfigGen: particle number must be <= 1e8.\n"; return; }
	if (!(density > 0.0) || !std::isfinite(density)) { ofile << "ConfigGen: number density must be finite and positive.\n"; return; }
	this->d = (DimensionType)dim; this->numPrts = (size_t)np; this->rho = density;
}

/** Give basis vectors of some Bravais lattices. */
void GetPredefinedBasis(const std::string & lattice, std::vector<GeometryVector> & basis) {
	basis.resize(0);

	if (strcmp(lattice.c_str(), "integer") == 0 || strcmp(lattice.c_str(), "Integer") == 0) {
		basis.emplace_back(GeometryVector(static_cast<double>(1.0)));
	}
	else if (strcmp(lattice.c_str(), "square") == 0 || strcmp(lattice.c_str(), "Square") == 0) {
		basis.emplace_back(GeometryVector(1.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 1.0));
	}
	else if (strcmp(lattice.c_str(), "triangle") == 0 || strcmp(lattice.c_str(), "triangular") == 0 || strcmp(lattice.c_str(), "Triangular") == 0) {
		basis.emplace_back(GeometryVector(1.0, 0.0));
		basis.emplace_back(GeometryVector(cos(pi / 3.0), sin(pi / 3.0)));
	}
	else if (strcmp(lattice.c_str(), "cubic") == 0||strcmp(lattice.c_str(), "Cubic") == 0) {
		basis.emplace_back(GeometryVector(1.0, 0.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 1.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 0.0, 1.0));
	}
	else if (strcmp(lattice.c_str(), "fcc") == 0 || strcmp(lattice.c_str(), "FCC") == 0) {//symmetric form
		basis.emplace_back(GeometryVector(0.0, 0.5, 0.5));
		basis.emplace_back(GeometryVector(0.5, 0.0, 0.5));
		basis.emplace_back(GeometryVector(0.5, 0.5, 0.0));
	}
	else if (strcmp(lattice.c_str(), "bcc") == 0 || strcmp(lattice.c_str(), "BCC") == 0) {//symmetric form
		basis.emplace_back(GeometryVector(-0.5, 0.5, 0.5));
		basis.emplace_back(GeometryVector(0.5, -0.5, 0.5));
		basis.emplace_back(GeometryVector(0.5, 0.5, -0.5));
	}
	else if (strcmp(lattice.c_str(), "fcc_cubic") == 0 || strcmp(lattice.c_str(), "FCC_cubic") == 0) {//non-symmetric form
		basis.emplace_back(GeometryVector(1.0, 0.0, 0.0));
		basis.emplace_back(GeometryVector(0.5, sqrt(3.0)*0.5, 0.0));
		basis.emplace_back(GeometryVector(0.5, 0.5 / sqrt(3.0), sqrt(2.0 / 3.0)));
	}
	else if (strcmp(lattice.c_str(), "bcc_cubic") == 0 || strcmp(lattice.c_str(), "BCC_cubic") == 0) {//Non-symmetric form
		basis.emplace_back(GeometryVector(1.0, 0.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 1.0, 0.0));
		basis.emplace_back(GeometryVector(0.5, 0.5, 0.5));
	}

	//else {
	//	std::cout << "Incorrect cell type" << std::endl;
	//}
}

/** Give basis vectors and positions of atomes for predefined periodic point patterns.*/
void GetPredefinedCrystals(const std::string & crystalName, std::vector<GeometryVector> & basis, std::vector<GeometryVector> & atoms) {
	atoms.resize(0);

	GetPredefinedBasis(crystalName, basis);
	//search for nonBravais crystals
	if (basis.size() == 0) {
		if (strcmp(crystalName.c_str(), "honeycomb") == 0 || strcmp(crystalName.c_str(), "HoneyComb") == 0) {
			GetPredefinedBasis("triangular", basis);
			atoms.emplace_back(GeometryVector(1.0/3.0, 1.0/3.0));
		}
		else if (strcmp(crystalName.c_str(), "kagome") == 0 || strcmp(crystalName.c_str(), "Kagome") == 0) {
			GetPredefinedBasis("triangular", basis);
			atoms.emplace_back(GeometryVector(0.5, 0.0));
			atoms.emplace_back(GeometryVector(0.0, 0.5));
		}
		else if (strcmp(crystalName.c_str(), "diamond") == 0 || strcmp(crystalName.c_str(), "Diamond") == 0) {
			GetPredefinedBasis("fcc", basis);
			atoms.emplace_back(GeometryVector(0.25, 0.25, 0.25));
		}
		else if (strcmp(crystalName.c_str(), "hcp") == 0 || strcmp(crystalName.c_str(), "HCP") == 0) {
			basis.clear();
			basis.emplace_back(0.5, -0.5*sqrt(3.0), 0.0);
			basis.emplace_back(0.5, +0.5*sqrt(3.0), 0.0);
			basis.emplace_back(0.0, 0.0, 2.0*sqrt(2.0/3.0));

			atoms.emplace_back(1.0 / 3.0, 1.0 / 3.0, 0.50);
		}
		else {
			std::cout << "Undefined crystal structure "<<std::endl;
		}
	}
}


// LatticeGen implementations live in Lattices/Lattice.cpp; that file is a
// translation-unit fragment #included here so the LatticeGen object code lands in
// this object (models/obj/GenerateConfigs.o).  Consumers that link only
// GenerateConfigs.o keep resolving LatticeGen without any makefile edits.
// IMPORTANT: Lattice.cpp must NOT be compiled as its own translation unit (it is
// deliberately absent from the makefile SRCS); doing so would define LatticeGen's
// members twice and violate the ODR at link time.
#include "Lattices/Lattice.cpp"


PoissonPtGen::PoissonPtGen(std::istream & ifile, std::ostream & ofile) {
	ofile	<< "------------------------------------------\n"
			<< "	Spatially uncorrelated point patterns \n"
			<< "------------------------------------------\n" << std::endl;

	ofile << "Number density = ";	ifile >> this->rho;		ofile << this->rho << "\n";
	ofile << "Modes (0 = cannonical/ 1 = grand cannonical)= ";	ifile >> this->mode;
	if (this->mode == 1 || this->mode == 0)
		ofile << this->mode << "\n";
	else {
		std::cout << "Wrong mode!\n";
		return;
	}

}


#ifdef STRUCTUREFACTOR_INCLUDED
#include "StructureFactor.h"
#endif

size_t numK_points(PeriodicCellList<Empty> & FundamentalCell, std::vector<GeometryVector> & ks) {
#ifdef STRUCTUREFACTOR_INCLUDED

	if (FundamentalCell.PeriodicVolume() > 0) {
		double CircularKMax, LinearKMax;
		std::cout << "Circular Kmax = ";
		std::cin >> CircularKMax;
		std::cout << "Linear Kmax = ";
		std::cin >> LinearKMax;

		ks = GetKs(FundamentalCell, CircularKMax, LinearKMax);
		return ks.size();

	}
	else {
		std::cerr << "FundamentalCell is not valid." << std::endl;
		return 0;
	}
#else
	//Parameters are only used when STRUCTUREFACTOR_INCLUDED is defined.
	(void)FundamentalCell;	(void)ks;
	std::cerr << "Please include \"StructureFactors.h\"\n";
	return 0;
#endif // STRUCTUREFACTOR_INCLUDED

}
