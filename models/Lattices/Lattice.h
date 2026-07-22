/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file Lattice.h
 * \brief A ConfigGen-derived class, LatticeGen, that constructs periodic point
 *        patterns / sphere packings (lattice packings, and the building block for
 *        perturbed lattices and imperfect lattices with vacancies).
 *
 * Split out of GenerateConfigs.h during the per-model reorganization; the class
 * bodies are the code-reviewed (canonical) module versions. ConfigGen and the
 * GetPredefinedBasis / GetPredefinedCrystals helpers live in <GenerateConfigs.h>.
 */

#ifndef MODELS_LATTICE_H
#define MODELS_LATTICE_H

#include <GenerateConfigs.h>

/**	/brief A class to construct a periodic point pattern or sphere packing.
 *  It stores a fundamental cell of a sphere packing. */
class LatticeGen : public ConfigGen
{
	std::vector<size_t> NumInSide;	//< NumInSide[i] = the number of unit cells in the ith direction, for i=0, ..., d-1.

protected:
	double maxPhi;					//< maximal packing fraction.
	SpherePacking F_Cell;			//<	A fundamental cell
//	virtual void SetBasis(const std::vector<GeometryVector> & basis);
public:
	/** A default constructor. */
	LatticeGen() : ConfigGen() { this->maxPhi = 0; this->F_Cell = SpherePacking(); };

	/** A constructor by a SpherePacking.
	 * If particle radii = 0, particles are considered "monodisperse" and the radius are set to be the maximum. */
	LatticeGen(const SpherePacking & unitcell) {
		NumInSide = std::vector<size_t>(unitcell.GetDimension(), 1);
		F_Cell = SpherePacking(unitcell);
		numPrts = F_Cell.NumParticle();
		rho = numPrts / F_Cell.PeriodicVolume();
		d = F_Cell.GetDimension();

		//Particle radii = 0.5 * nearest neighbor distance.
		{
			bool uninitializedPrtRadii = false;		//True = at least one particle has a nonpositive radius.
			for (size_t i = 0; i < this->numPrts; i++)
				uninitializedPrtRadii |= (F_Cell.GetCharacteristics(i) <= 0.0);

			if (F_Cell.CheckOverlap()) {
				std::cout << "Particles are overlapping. Their sizes are adjusted.\n";
				uninitializedPrtRadii = false;
			}

			if (uninitializedPrtRadii) {
				std::cout	<< "Paritlce radii are not initialized."
							<<"Particle radii are set to be nearest neighbor distances.\n";
				/* at least one particle has a nonpositive radius
				=> reset particle radii such that they are half of nearest neighbor distances.*/
				for (size_t i = 0; i < this->numPrts; i++)
					F_Cell.UpdateCharacteristics(0.5*F_Cell.NearestParticleDistance(i), i);
			}
		}

		maxPhi = this->UpdateMaxPhi();
	}

	/** \brief A constructor by name of lattices.
	 *  \param name	The name of predefined periodic packings. */
	LatticeGen(const std::string & name);

	/** A constructor by an input stream. */
	LatticeGen(std::istream & ifile, std::ostream & ofile);

	/** A copy constructor. */
	LatticeGen(const LatticeGen & source) : LatticeGen(source.F_Cell) {};


	/** \brief A member function to update and return maximum packing fraction.
	 *  \return Maximal packing fraction.*/
	double UpdateMaxPhi();

	/** \brief Insert a spherical particle in the fundamental cell.
	 *  \param relativeCoorinates	Relative coordinates of a particle.
	 *  \param radius	Particle radius. The default value is 0.0. */
	void AddPrt(const GeometryVector & relativeCoorinates, double radius = 0.0);

	/** \brief Remove the ith particle in the fundamental cell.
	 *  \param i The index of a particle that is removed. */
	void RemovePrt(size_t i);

	/** \brief Get basis vectors of a fundamental cell.
	 *  \param[out] basis		The basis vectors of the fundamental cell. */
	void GetBasisVectors(std::vector<GeometryVector> & basis) {
		basis.clear();
		for (int i = 0; i < this->d; i++)
			basis.push_back(this->F_Cell.GetBasisVector(i));
	}

	/** \brief Rescale particle size.
	 *  \param newPhi	A new packing fraction. */
	virtual void SetPackingFraction(double newPhi);

	/** \brief Rescale the cell size.
	 *  \param newRho	A new number density. */
	virtual void SetRho(double newRho);

	/** \brief Resize particle radii to attain a maximal packing fraction. */
	void Set2MaxPhi();

	/** \brief Set the number of unit cells.
	 *  \param numInEachSide The number of unit cells in each direction (must be positive). */
	void SetNumCells(size_t numInEachSide) {
		if (numInEachSide == 0) {
			std::cerr << "SetNumCells: number of unit cells must be positive.\n";
			return;
		}
		this->NumInSide = std::vector<size_t>(this->d, numInEachSide);
	};

	/** \brief Set the number of unit cells.
	 *  \param numInSide	The number of unit cells in each direction (entry i for the ith direction); must have exactly d positive entries. */
	void SetNumCells(const std::vector<size_t> & numInSide) {
		if (numInSide.size() != (size_t)this->d) {
			std::cerr << "SetNumCells: expected " << (size_t)this->d << " entries, got " << numInSide.size() << ".\n";
			return;
		}
		for (size_t i = 0; i < numInSide.size(); i++) {
			if (numInSide[i] == 0) {
				std::cerr << "SetNumCells: entry " << i << " must be positive.\n";
				return;
			}
		}
		this->NumInSide = std::vector<size_t>(numInSide);
	};


	/** \brief Generate a periodic point pattern.
 	 *  \param[out] c	A periodic Configuration that contains N[0] x...x N[d-1] unit cells. */
	virtual void GenerateC(Configuration & c) { c = MultiplicateStructure(Configuration(this->F_Cell, "a"), this->NumInSide); }

	/** \brief Generate a periodic sphere packing.
	 *  \param[out] c	A periodic SpherePacking that contains N[0] x...x N[d-1] unit cells. */
	virtual void GenerateP(SpherePacking & c);

};

#endif // MODELS_LATTICE_H
