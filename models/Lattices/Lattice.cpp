/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file Lattice.cpp
 *	\brief Implementations of LatticeGen (declared in Lattice.h).
 *
 * This translation-unit fragment is #included by GenerateConfigs.cpp so that the
 * LatticeGen object code lands in GenerateConfigs.o (consumers that link only
 * models/obj/GenerateConfigs.o keep resolving LatticeGen without makefile edits).
 * GetPredefinedBasis / GetPredefinedCrystals and Echo come from GenerateConfigs.h.
 */

#include "Lattice.h"

/***********************************
	LatticeGen: Constructors
***********************************/

/** \brief Construct a LatticeGen object by lattice name.
 *	\param Name (case-insenstive)	1D:	Integer
									2D: Square, Triangular (triangle), Honeycomb, Kagome
									3D: Cubic, BCC_cubic, FCC, FCC_cubic, BCC_symmetric, Diamond */
LatticeGen::LatticeGen(const std::string & Name) {
	//Change all alphabets to lower cases.
	std::string name(Name);
	std::transform(name.begin(), name.end(), name.begin(), [](unsigned char ch){ return (char)std::tolower(ch); });
	//Obtain basis vectors and positions of atoms
	std::vector<GeometryVector> basis, atoms;
	GetPredefinedCrystals(name, basis, atoms);
	if (basis.empty()) {
		//Unknown lattice: leave a valid empty generator instead of dereferencing &basis[0].
		std::cerr << "LatticeGen: unknown lattice \"" << Name << "\"; constructing an empty generator.\n";
		this->d = 0; this->numPrts = 0; this->rho = 0.0; this->maxPhi = 0.0;
		this->F_Cell = SpherePacking();
		this->SetNumCells(std::vector<size_t>());
		return;
	}
	this->F_Cell = SpherePacking((DimensionType)basis.size(), &basis[0], 1.0);

	//a particle at the origin
	this->F_Cell.Insert(0.0, GeometryVector(static_cast<int>(basis.size())));
	//Insert other particles
	for (auto rel_pos = atoms.begin(); rel_pos != atoms.end(); rel_pos++)
		this->F_Cell.Insert(0.0, *rel_pos);

	this->d	= this->F_Cell.GetDimension();
	this->numPrts = this->F_Cell.NumParticle();
	this->rho = this->numPrts / this->F_Cell.PeriodicVolume();
	this->SetNumCells(1);

	//Max packing fraction
	this->Set2MaxPhi();
}

/** \brief Construct a LatticeGen object by an input stream.
 *	\param[in] ifile	An input stream.
 *	\param[out] ofile	An output stream.
 */
LatticeGen::LatticeGen(std::istream & ifile, std::ostream &ofile) {
	char temp[300] = {};
	std::string cell_type;	double value;	size_t NumInEachSide;
	ofile << "-----------------------------------\n"
		  << "	Choose a lattice to generate     \n"
		  << "-----------------------------------" << std::endl;
	ofile << "Cell Type = ";	Echo<std::string>(ifile, ofile, cell_type);
	{//Define a fundamental cell!
	 //Change all alphabets to the lower cases.
		std::string name(cell_type);
		std::transform(name.begin(), name.end(), name.begin(), [](unsigned char ch){ return (char)std::tolower(ch); });
		//Obtain basis vectors and positions of atoms
		std::vector<GeometryVector> basis, atoms;
		GetPredefinedCrystals(name, basis, atoms);
		if (basis.empty()) {
			//Unknown lattice: leave a valid empty generator instead of dereferencing &basis[0].
			ofile << "LatticeGen: unknown lattice \"" << cell_type << "\"; constructing an empty generator." << std::endl;
			this->d = 0; this->numPrts = 0; this->rho = 0.0; this->maxPhi = 0.0;
			this->F_Cell = SpherePacking();
			this->SetNumCells(std::vector<size_t>());
			return;
		}
		this->F_Cell = SpherePacking(basis.size(), &basis[0], 1.0);

		//a particle at the origin
		this->F_Cell.Insert(0.0, GeometryVector(static_cast<int>(basis.size())));
		//Insert other particles
		for (auto rel_pos = atoms.begin(); rel_pos != atoms.end(); rel_pos++) {
			this->F_Cell.Insert(0.0, *rel_pos);
		}


		this->d = this->F_Cell.GetDimension();

		ofile << "Add additional particles (y/n)?";
		ifile >> temp;
		if (strcmp(temp, "y") == 0) {
			GeometryVector pos(static_cast<int>(this->d));
			size_t NumAdditionalPrts = 0,
				NumInit = this->F_Cell.NumParticle();
			ofile << "How many would you add?";
			ifile >> NumAdditionalPrts;
			for (size_t i = 0; i < NumAdditionalPrts; i++) {
				ofile << "Please input relative coordinates of particle " << i+NumInit <<":\t";
				for (int j = 0; j < this->d; j++)
					ifile >> pos.x[j];

				this->AddPrt(pos);
			}
		}
		this->numPrts = this->F_Cell.NumParticle();

		//Change particle radii to be maximized.
		this->Set2MaxPhi();
	}

	//Resize the fundamental cell and particles.
	{
		ofile << "Rescale the fundamental cell with respect to \n"
			  << "\t D (nearest neighbor distance) \n"
			  << "\t rho (number density) \n"
			  << "\t a (lattice constant).";
		ifile >> temp;
		ofile << temp << " = ";
		ifile >> value;

		double factor;
		if (strcmp(temp, "D") == 0) {
			double D = this->F_Cell.NearestParticleDistance(0);
			for (size_t i = 1; i < this->numPrts; i++)
				D = std::min(D, this->F_Cell.NearestParticleDistance(i));

			factor = value / D;
		}
		else if (strcmp(temp, "rho") == 0) {
			double currRho = this->F_Cell.NumParticle() / this->F_Cell.PeriodicVolume();
			factor = 1.0/pow(value / currRho, 1.0 / this->d);
		}
		else if (strcmp(temp, "a") == 0) {
			double currA = sqrt(this->F_Cell.GetBasisVector(0).Modulus2());
			for (int i = 1; i < this->d; i++)
				currA = std::max(currA, sqrt(this->F_Cell.GetBasisVector(i).Modulus2()));

			factor = value / currA;
		}
		else {
			ofile << "Please choose a correct standard." << std::endl;
			exit(1);
		}

		this->F_Cell.Rescale(factor);
	}




	ofile << "Number of unit cells in each side = ";	Echo<size_t>(ifile, ofile, NumInEachSide);
	assert(NumInEachSide > 0);
	this->SetNumCells(NumInEachSide);
}

/***********************************
	LatticeGen: Member functions
***********************************/

/** \brief Generate a periodic sphere packing.
 *  \param[out] c	A periodic SpherePacking */
void LatticeGen::GenerateP(SpherePacking & c) {
	size_t NumCells = 1;
	for (int i = 0; i < this->d; i++)
		NumCells *= this->NumInSide[i];

	if (NumCells == 1) {
		//There is only one unit cell
		c = SpherePacking(this->F_Cell);
	}
	else {
		//There are multiple unit cells.
		Configuration cc = MultiplicateStructure(Configuration(this->F_Cell, "a"), this->NumInSide);
		c = SpherePacking(cc, 0);

		//Assign particle radii: each fundamental-cell particle has NumCells contiguous replicas.
		size_t idx = 0;
		for (size_t i = 0; i < this->numPrts; i++) {
			for (size_t j = 0; j < NumCells; j++, idx++)
				c.UpdateCharacteristics(this->F_Cell.GetCharacteristics(i), (int)idx);
		}
	}
}

/** Update the maximal packing fraction, but doesn't change particle radii.*/
double LatticeGen::UpdateMaxPhi() {
	double cellLeng = sqrt(this->F_Cell.GetBasisVector(0).Modulus2());
	for (int i = 1; i < this->d; i++)
		cellLeng = std::min(cellLeng, sqrt(this->F_Cell.GetBasisVector(i).Modulus2()));

	double PrtVolume = 0.0;	//Particle volume
	//		NND;		//nearest neighbor distance
	for (size_t i = 0; i < this->numPrts; i++) {
		PrtVolume += HyperSphere_Volume(this->d, 0.5*this->F_Cell.NearestParticleDistance(i));
	}

	this->maxPhi = PrtVolume/this->F_Cell.PeriodicVolume();
	return this->maxPhi;
}

/** Check particle radii to 0.5 of their nearest neighbor distances. */
void LatticeGen::Set2MaxPhi() {
	for (size_t i = 0; i < this->numPrts; i++) {
		this->F_Cell.UpdateCharacteristics(0.5*this->F_Cell.NearestParticleDistance(i), i);
	}
	if (this->F_Cell.CheckOverlap())
		std::cout << "Something wrong\n";
	else {
		this->maxPhi = this->F_Cell.PackingFraction();
		std::cout << "Particle radii are adjusted to 0.5 of nearest neighbor distances.\n"
			<< "New maximal packing fraction is " << this->maxPhi << "\n";
	}


}

/** \brief Add a spherical particle in the fundamental cell.
 *  If the added particle is overlapping with other particles, particle radii are adjusted.
 *  \param relativePosition	Relative coordinates at which a spherical particle is added.
 *  \param radius			Radius of an added particle. */
void LatticeGen::AddPrt(const GeometryVector & relativePosition, double radius) {
	this->numPrts++;
	this->F_Cell.Insert(radius, relativePosition);
	this->rho += 1.0 / this->F_Cell.PeriodicVolume();
	this->UpdateMaxPhi();

	if (this->F_Cell.CheckOverlap()) {
		std::cout << "Particles are overlapping.\n";
		this->Set2MaxPhi();
	}
}

/** Remove the i th particle in the fundamental cell. */
void LatticeGen::RemovePrt(size_t i) {
	this->numPrts--;
	this->F_Cell.DeleteParticle(i);
	this->rho -= 1.0 / this->F_Cell.PeriodicVolume();
	this->UpdateMaxPhi();
}

/** Resize particles. */
void LatticeGen::SetPackingFraction(double targetPackingFraction) {
	if (targetPackingFraction <= this->maxPhi) {
		//Rescale_prtOnly scales radii by `factor`, so packing fraction scales by factor^d.
		double factor = pow(targetPackingFraction / this->F_Cell.PackingFraction(), 1.0 / this->d);
		this->F_Cell.Rescale_prtOnly(factor);
		if (this->F_Cell.CheckOverlap()) {
			this->Set2MaxPhi();
			factor = pow(targetPackingFraction / this->F_Cell.PackingFraction(), 1.0 / this->d);
			this->F_Cell.Rescale_prtOnly(factor);
		}
	}
	else {
		std::cout << "target packing fraction is too high!\n";
	}
}

/** Resize fundamental cell. */
void LatticeGen::SetRho(double numDensity) {
	double factor = 1.0/pow(numDensity / this->rho, 1.0 / this->d);
	this->rho = numDensity;
	this->F_Cell.Rescale(factor);
}
