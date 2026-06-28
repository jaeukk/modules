/**
 *	Author	: Ge Zhang
 *	Email	: 
 *	Created	:	2015? 
 *	Comments: Jaeuk Kim (November 2025).
 */

/** \file PeriodicCellList.h
 *	\brief Header file to define general point patterns under periodic boundary conditions. It and its derived classes (e.g., Configuration and SpherePacking) implement cell-list acceleration techniques. 
 * 			(November 2025)
 *	Some type-related and unintialized warnings are resolved. 
 * */


#ifndef PeriodicCellList_INCLUDED
#define PeriodicCellList_INCLUDED

#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <list>

#include <omp.h>

#include "GeometryVector.h"
#include "RandomGenerator.h"
#include "etc.h"
//#include "etc1.h"

const size_t ParticlePerCell = 1;

//class for point configurations under periodic boundary conditions;
//unit cells can be any Bravais lattice structure.
//Characteristics store additional information on each particle...

//Configuration & SpehrePacking are derived from this.

/** \brief Cell-list accelerator for periodic point configurations parameterized by particle traits. */
template<typename OtherCharacteristics>
class PeriodicCellList
{
	template<typename a> friend class PeriodicCellList;

	//some concepts:
	//there are 3 types of coordinates:
	//1. Cartesian Coordinates
	//2. Relative Coordinates, this is the coordinates relative to the periodic box base vectors.
	//		The periodic box base vectors are PeriodicCellList::BasisVector
	//3. Lattice Coordinates, this is the coordinates relative to the cell lattice(i.e. the base vectors of a cell list cell).
public:
	typedef unsigned long Index;

protected:
	mutable double MaxRadius;//=0 indicates this variable is not valid. call UpdateMaxRadius to update it //moved from SpherePacking
	DimensionType Dimension;
	GeometryVector BasisVector[::MaxDimension];//the base vector of a periodic cell
	signed long CellRank[::MaxDimension];//number of cells along each basis vector
	std::vector<std::vector<size_t > > Cells;//content Particles in the cells
	std::vector<GeometryVector> ParticleRelatives;
	std::vector<GeometryVector> ParticleCartesians;
	std::vector<OtherCharacteristics> ParticleCharacteristics;
	double ApproximateCellSize;// the best choice of cell size obtained in constructor when user input cell size and cell rank. We will adjust CellRank so that their size is roughly equal to this quantity.
	double vol; // the volume of a fundamental cell: add recently

	//reciprocal basis vector, calculated when needed
	mutable bool ReciprocalBasisVectorValid[::MaxDimension];
	mutable GeometryVector ReciprocalBasisVector[::MaxDimension];

	//member data for vicinity lattice
	struct VicinityLatticesMember
	{
		GeometryVector CartesianCoordinates;
		signed long LatticeShift[MaxDimension];
	};
	mutable std::vector<VicinityLatticesMember> VicinityLattices;
	mutable std::map<double, size_t> VicinityLatticePartitions;
	mutable double VicinityLatticesRc2;
	bool SortedVicinityList;//use sorted vicinity lattice list, which takes some time to sort, but makes iteration more efficient

	/** \brief Recompute Cartesian coordinates of particle a from its relative coordinates.
	 *  \param a Index of the particle to update. */
	void UpdateCartesianCoordinates(size_t a)
	{
		if (!(this->ParticleRelatives[a].x[0] < 1 + ::LengthPrecision) || !(this->ParticleRelatives[a].x[0] >= 0 - ::LengthPrecision))
		{
			std::cerr << "Error in PeriodicCellList::UpdateCartesianCoordinates : invalid ParticleRelatives. ParticleRelatives[" << a << "]=" << this->ParticleRelatives[a];
			exit(0);
		}
		this->ParticleCartesians[a] = GeometryVector(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->ParticleCartesians[a].AddFrom((this->BasisVector[i])*(this->ParticleRelatives[a].x[i]));
	}

	/** \brief Wrap particle a back into the primary periodic box.
	 *  \param a Index of the particle to wrap. */
	void PutToPeriodicBox(size_t a)
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->ParticleRelatives[a].x[i] -= std::floor(this->ParticleRelatives[a].x[i]);
	}
	/** \brief Convert a dimensional cell index into a flat 1-dimensional cell index.
	 *  \param index Per-dimension cell indices.
	 *  \return      Flattened cell index. */
	Index ConvertIndex(const signed long * index) const//convert dimensional index to the 1-dimension cell index
	{
		Index result = index[0];
		for (DimensionType i = 1; i<this->Dimension; i++)
		{
			assert(index[i] >= 0);
			assert(index[i]<this->CellRank[i]);
			result *= this->CellRank[i];
			result += index[i];
		}
		return result;
	}

	/** \brief Compute the per-dimension cell index containing a relative coordinate.
	 *  \param[out] index               Filled with the dimensional cell index.
	 *  \param[in]  RelativeCoordinates Point in relative (box) coordinates. */
	void GetFullIndex(signed long * index, const GeometryVector & RelativeCoordinates) const
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			index[i] = (signed long)(std::floor(RelativeCoordinates.x[i] * this->CellRank[i]));
			if (index[i]<0 || index[i] >= this->CellRank[i])
			{
				std::cerr << "Error in Configuration::GetFullIndex : index invalid. index=" << index[i] << ", Relative Coordinate=" << RelativeCoordinates.x[i] << ", Cell Rank=" << this->CellRank[i] << '\n';
				assert(false);
			}
		}
	}


	/** \brief Get the index of the cell that a should be placed in, the LatticeCoordinates of particle must be in box
	 *  \param a   Index of the particle whose containing cell is requested.
	 *  \return    Flattened index of the cell containing particle a. */
	Index GetIndex(const size_t a) const //get the index of the cell that a should be placed in, the LatticeCoordinates of particle must be in box
	{
		signed long index[::MaxDimension] = {0};
		this->GetFullIndex(index, this->ParticleRelatives[a]);
		return this->ConvertIndex(index);
	}


	/** \brief Move a to another cell if necessary, from is the cell that a comes from, return true if moved
	 *  \param a      Iterator to the particle entry within its current cell.
	 *  \param from   The cell that particle a currently comes from.
	 *  \return       True if the particle was moved to a different cell. */
	bool UpdateCellList(std::vector<size_t>::iterator a, std::vector<size_t> & from)//move a to another cell if necessary, from is the cell that a comes from, return true if moved
	{
		//calculate new index
		Index newindex = this->GetIndex((*a));
		if (&(from) == &(this->Cells[newindex]))//don't need to move
			return false;

		//add to new cell
		this->Cells[newindex].push_back(*a);

		//remove from original cell
		*a = from.back();
		from.pop_back();
		return true;
	}


	/** \brief Recursively enumerate cell-lattice points within a spherical shell into the vicinity list.
	 *  \param now         Current cell-lattice point in Cartesian coordinates.
	 *  \param DimIDo      Remaining dimension index to recurse over (0 terminates the recursion).
	 *  \param nowRelative Current lattice shift in relative (cell-lattice) coordinates.
	 *  \param ortho       Gram-Schmidt orthogonalized cell basis vectors.
	 *  \param Rmax2       Squared outer radius of the spherical shell.
	 *  \param Rmin2       Squared inner radius of the spherical shell (already-included lattice points). */
	bool _VicinityLatticeIteration(GeometryVector now, DimensionType DimIDo, GeometryVector nowRelative, GeometryVector * ortho, double Rmax2, double Rmin2) const
	{
		const double VicinityLatticeRmin_Tolerance = 1e-10;
		if (DimIDo == 0)
		{
			double temp = now.Modulus2();
			if (temp<Rmax2)
			{
				bool ToAdd;
				if (temp >(1.0 + VicinityLatticeRmin_Tolerance)*Rmin2)
					ToAdd = true;
				else if (temp < (1.0 - VicinityLatticeRmin_Tolerance)*Rmin2)
					ToAdd = false;
				else
				{
					//cannot tell if the lattice point was already included just from the distance
					//instead, go through the entire VicinityLattices vector to see if it is already added
					VicinityLatticesMember temp;
					for (DimensionType i = 0; i<this->Dimension; i++)
					{
						double & x = nowRelative.x[i];
						temp.LatticeShift[i] = x >= 0 ? static_cast<signed long>(x + 0.5) : static_cast<signed long>(x - 0.5);
					}
					ToAdd = true;
					for (auto iter = this->VicinityLattices.begin(); iter != this->VicinityLattices.end(); iter++)
					{
						bool LatticeShiftIdentical = true;
						for (DimensionType d = 0; d<this->Dimension; d++)
							if (iter->LatticeShift[d] != temp.LatticeShift[d])
							{
								LatticeShiftIdentical = false;
								break;
							}
						if (LatticeShiftIdentical)
						{
							ToAdd = false;
							break;
						}
					}
				}

				if (ToAdd)
				{
					VicinityLatticesMember temp;
					temp.CartesianCoordinates = now;
					for (DimensionType i = 0; i<this->Dimension; i++)
					{
						double & x = nowRelative.x[i];
						temp.LatticeShift[i] = x >= 0 ? static_cast<signed long>(x + 0.5) : static_cast<signed long>(x - 0.5);
					}

					this->VicinityLattices.push_back(temp);
				}
				return true;
			}
			else
				return false;
		}


		GeometryVector k = now;
		for (DimensionType i = 0; i<DimIDo; i++)
		{
			GeometryVector & e = ortho[i];
			k.MinusFrom(e*(k.Dot(e) / e.Dot(e)));
		}
		if (k.Dot(k)>Rmax2)
			return false;
		//right now k is the normal vector of plane (e1, e2,...) in subspace (e1, e2, ... , now)
		GeometryVector a = now - k;//the part of now in that plane
		GeometryVector ethis = this->BasisVector[DimIDo - 1] * (1.0 / this->CellRank[DimIDo - 1]);//base vector that this function will loop through
		for (DimensionType i = 0; i<DimIDo - 1; i++)
		{
			GeometryVector & e = ortho[i];
			a.MinusFrom(e*(a.Dot(e) / e.Dot(e)));
			ethis.MinusFrom(e*(ethis.Dot(e) / e.Dot(e)));
		}
		double coeff = 0;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (std::abs(ethis.x[i])> ::LengthPrecision)
			{
				coeff = a.x[i] / ethis.x[i];
				break;
			}
		}//a equals to ethis*coeff
		assert((a - coeff * ethis).Modulus2()< ::LengthPrecision);

		double CoeffMaxInPlane = std::sqrt((Rmax2 - k.Modulus2()) / ethis.Modulus2());//sqrt(R^2-k^2)/|ethis|
		signed long imin = (signed long)(std::ceil((-1)*CoeffMaxInPlane - coeff));
		signed long imax = (signed long)(std::floor(CoeffMaxInPlane - coeff));

		//int nearest=std::floor(coeff);
		GeometryVector Unit(this->Dimension);
		Unit.x[DimIDo - 1] = 1;
		for (signed long i = imin; i <= imax; i++)
		{
			this->_VicinityLatticeIteration(now + i * this->BasisVector[DimIDo - 1] * (1.0 / this->CellRank[DimIDo - 1]), DimIDo - 1, nowRelative + i * Unit, ortho, Rmax2, Rmin2);
		}
		return true;
	}

public:
	//This class can work in one of the two modes:
	//UseSortedList=true : efficient when the basis vectors almost never change.
	//UseSortedList=false : efficient when the basis vectors change.

	/** \brief Construct an empty list with no dimension set. */
	PeriodicCellList()
	{
		//this->VicinityLatticesRc2=0;
		this->Dimension = 0;
		Index CellCount = 1;//need to calculate overall number of cells
		this->Cells.assign(CellCount, std::vector<size_t >());
		this->MaxRadius = 0.0;
		this->ApproximateCellSize = ::MaxDimension;
		this->SortedVicinityList = false;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		this->vol = -1.0;
	}

	/** \brief Construct a periodic cell list for a given box geometry.
	 *  \param Dimension Spatial dimension.
	 *  \param BasisVectors Basis vectors defining the periodic box.
	 *  \param CellSize Target linear cell size.
	 *  \param UseSortedList Choose sorted vicinity list for faster neighbor checks when geometry is fixed. */
	PeriodicCellList(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList = false)
	{
		//this->VicinityLatticesRc2=0;
		this->Dimension = Dimension;
		Index CellCount = 1;//need to calculate overall number of cells
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			GeometryVector & vec = BasisVectors[i];
			this->CellRank[i] = (signed long)(std::floor(std::sqrt(vec.Modulus2()) / CellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;

			//error check
			if (this->CellRank[i] <= 0)
			{
				std::cerr << "error in Configuration constructor: invalid cell rank. basis vector length squared=" << vec.Modulus2() << ", CellSize=" << CellSize << '\n';
				assert(false);
			}
			this->BasisVector[i] = vec;
			CellCount *= this->CellRank[i];
		}
		this->MaxRadius = 0.0;
		for (DimensionType i = this->Dimension; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = 0;
			this->BasisVector[i] = GeometryVector();
		}
		this->Cells.assign(CellCount, std::vector<size_t >());

		this->ApproximateCellSize = CellSize;
		this->SortedVicinityList = UseSortedList;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		this->vol = ::Volume(BasisVectors, Dimension);
	}

	/** \brief Copy-construct, duplicating geometry and particle data.
	 *  \param source List to copy geometry and particle data from. */
	PeriodicCellList(const PeriodicCellList<OtherCharacteristics> & source)
		: MaxRadius(source.MaxRadius), Dimension(source.Dimension),
		Cells(source.Cells),
		ParticleRelatives(source.ParticleRelatives), ParticleCartesians(source.ParticleCartesians),
		ParticleCharacteristics(source.ParticleCharacteristics),
		ApproximateCellSize(source.ApproximateCellSize), vol(source.vol),
		VicinityLatticesRc2(0.0), SortedVicinityList(source.SortedVicinityList)
	{
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->BasisVector[i] = source.BasisVector[i];
		}
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = source.CellRank[i];
		}
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
	}

	/** \brief Construct from another list while assigning a uniform characteristic.
	 *  \param source Source list providing geometry and particle positions.
	 *  \param cha    Characteristic value assigned to every particle. */
	template<typename add> PeriodicCellList(const PeriodicCellList<add> & source, const OtherCharacteristics & cha)
		: Dimension(source.Dimension),
		Cells(source.Cells),
		ParticleRelatives(source.ParticleRelatives), ParticleCartesians(source.ParticleCartesians),
		ParticleCharacteristics(source.NumParticle(), cha),
		ApproximateCellSize(source.ApproximateCellSize), vol(source.vol),
		VicinityLatticesRc2(0.0), SortedVicinityList(source.SortedVicinityList)
	{
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->BasisVector[i] = source.BasisVector[i];
		}
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = source.CellRank[i];
		}
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
	}


	/** \brief Trivial destructor. */
	virtual ~PeriodicCellList() = default;

	/** \brief Exact comparison of dimension, basis, and particle layout.
	 *  \param right List to compare against.
	 *  \return      True if dimension, basis vectors, and particle positions all match. */
	bool operator == (const PeriodicCellList<OtherCharacteristics> & right) const
	{
		if (right.Dimension != this->Dimension)
			return false;
		if (right.NumParticle() != this->NumParticle())
			return false;
		for (DimensionType d = 0; d<this->Dimension; d++){
			if (right.BasisVector[d] != this->BasisVector[d])
				return false;
		}
		for (size_t i = 0; i<this->NumParticle(); i++){
			if (right.ParticleRelatives[i] != this->ParticleRelatives[i])
				return false;
		}
		return true;
	}

	/** \brief Reset cached vicinity lattice data. */
	void ClearVicinityLatticeList(void)
	{
		this->VicinityLattices.clear();
		this->VicinityLatticesRc2 = 0.0;
		this->VicinityLatticePartitions.clear();
		this->VicinityLatticePartitions.insert(std::make_pair((double)(0.0), (size_t)(0)));
	}


	/** \brief Reduce one basis vector against another to shorten the box, rebuilding the cell list.
	 *  \param ilongest Index of the basis vector to reduce against the others. */
	void TryRefineBasisVectors_inner(int ilongest)
	{
		//find the other basis vector
		int iother;
		signed int movepart = 0;
		for (iother = 0; iother<this->Dimension; iother++)
		{
			if (iother == ilongest) continue;
			movepart = (signed int)(std::floor(0.5 + (this->BasisVector[ilongest].Dot(this->BasisVector[iother])) / (this->BasisVector[iother].Dot(this->BasisVector[iother]))));
			if (movepart != 0)
				break;
		}
		if (movepart == 0)
			return;

		this->BasisVector[ilongest] = this->BasisVector[ilongest] - static_cast<double>(movepart)*this->BasisVector[iother];
		for (size_t i = 0; i<this->ParticleRelatives.size(); i++)
		{
			this->ParticleRelatives[i].x[iother] += this->ParticleRelatives[i].x[ilongest] * movepart;
			this->PutToPeriodicBox(i);
			this->UpdateCartesianCoordinates(i);
		}
		this->ClearVicinityLatticeList();
		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
		size_t NumCells = 1;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double length = std::sqrt(this->BasisVector[i].Modulus2());
			this->CellRank[i] = (signed long)(std::floor(length / this->ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			NumCells *= this->CellRank[i];
		}
		this->Cells.resize(NumCells);
		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();

		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->Cells[this->GetIndex(i)].push_back(i);
		}
	}


	/** \brief Attempt to reduce all basis vectors to a more compact (near-orthogonal) box. */
	void TryRefineBasisVectors(void)
	{
		if (this->Dimension == 1)
			return;
		for (int i = 0; i<this->Dimension; i++)
			this->TryRefineBasisVectors_inner(i);
	}


	/** \brief Return the volume of the fundamental periodic cell.*/
	double PeriodicVolume(void) const
	{
		double result = 0;
		if (this->vol <= 0.0) {
			result = ::Volume(this->BasisVector, this->Dimension);
		}
		else {
			result = this->vol;
		}
		return result;
	}

	/** \brief Return the dimensionality of the configuration box. */
	DimensionType GetDimension(void) const
	{
		return this->Dimension;
	}

	/** \brief Insert a particle with a characteristics at the given relative coordinate (wrapped to the box).
	 *  \param cha               Characteristic value assigned to the new particle.
	 *  \param RelativeCoordinate Position of the new particle in relative (box) coordinates. */
	virtual void Insert(OtherCharacteristics cha, const GeometryVector & RelativeCoordinate)
	{
		//this->Particles.push_back(particle(this->Dimension, cha));
		this->ParticleCartesians.push_back(GeometryVector(this->Dimension));
		this->ParticleRelatives.push_back(GeometryVector(this->Dimension));
		this->ParticleCharacteristics.push_back(cha);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			assert(RelativeCoordinate.x[i] == RelativeCoordinate.x[i]);
			//this->Particles.back().RelativeCoordinates.x[i]=RelativeCoordinate.x[i]-std::floor(RelativeCoordinate.x[i]);
			//this->Particles.back().RelativeCoordinates.x[i]-=std::floor(this->Particles.back().RelativeCoordinates.x[i]);

			//NEED TO subtract the floor TWICE
			this->ParticleRelatives.back().x[i] = RelativeCoordinate.x[i] - std::floor(RelativeCoordinate.x[i]);
			this->ParticleRelatives.back().x[i] -= std::floor(this->ParticleRelatives.back().x[i]);
		}
		this->UpdateCartesianCoordinates(this->NumParticle() - 1);
		//---test code
		int test = this->GetIndex(this->NumParticle() - 1);
		this->Cells[test].push_back(this->NumParticle() - 1);
		//this->Cells[this->GetIndex(this->NumParticle() - 1)].push_back(this->NumParticle() - 1);
	}

	/** \brief Insert a particle at a random position using the supplied RNG.
	 *  \param cha Characteristic value assigned to the new particle.
	 *  \param gen Random generator used to draw the relative coordinates. */
	void Insert(OtherCharacteristics cha, RandomGenerator & gen)
	{
		GeometryVector a(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			a.x[i] = gen.RandomDouble();
		this->Insert(cha, a);
	}

	/** \brief Populate vicinity lattice up to distance Rc (squared cached).
	 *  \param Rc Cutoff distance up to which cell-lattice points are enumerated. */
	void GetVicinityLattice(double Rc) const
	{
		size_t PrevVicinityLatticeNum = this->VicinityLattices.size();
		//prepare orthogonalized basis
		GeometryVector * ortho = new GeometryVector[this->Dimension];
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			ortho[i] = this->BasisVector[i];
			for (DimensionType j = 0; j<i; j++)
			{
				GeometryVector & e = ortho[j];
				ortho[i].MinusFrom(e*(ortho[i].Dot(e) / e.Dot(e)));
			}
		}

		this->_VicinityLatticeIteration(GeometryVector(this->Dimension), this->Dimension, GeometryVector(this->Dimension), ortho, Rc*Rc, this->VicinityLatticesRc2);

		delete[] ortho;

		if (this->SortedVicinityList)
			std::sort(this->VicinityLattices.begin() + PrevVicinityLatticeNum, this->VicinityLattices.end(), [](const VicinityLatticesMember & left, const VicinityLatticesMember & right) -> bool {return left.CartesianCoordinates.Modulus2()<right.CartesianCoordinates.Modulus2(); });

		this->VicinityLatticesRc2 = Rc*Rc;
	}

	/** \brief Debug-print vicinity lattice vectors.
	 *  \param out Output stream to print the lattice vectors to. */
	void _PrintVicinityLattices(std::ostream & out) const
	{
		for (auto iter = this->VicinityLattices.begin(); iter<this->VicinityLattices.end(); iter++)
			out << (*iter).CartesianCoordinates << '\n';
	}

	/** \brief Pick a random particle index; optionally return its cell id.
	 *  \param gen        Random generator used to select the particle.
	 *  \param SourceCell If non-null, receives the cell index containing the chosen particle.
	 *  \return           Index of the randomly chosen particle. */
	size_t GetRandomParticle(RandomGenerator & gen, Index * SourceCell = nullptr) const//pick up a random particle, return it's pointer, if SourceCell is not null, also return the cell that the particle come from
	{
		size_t result;
		result = std::floor(gen.RandomDouble()*this->NumParticle());
		if (SourceCell != nullptr)
			*SourceCell = this->GetIndex(result);
		return result;
	}

	/** \brief Total number of particles stored. */
	size_t NumParticle(void) const
	{
		return this->ParticleRelatives.size();
	}

	/** \brief Access mutable per-particle characteristic.
	 *  \param Num Index of the particle. */
	virtual OtherCharacteristics & GetCharacteristics(size_t Num)
	{
		return this->ParticleCharacteristics[Num];
	}

	/** \brief Access const per-particle characteristic.
	 *  \param Num Index of the particle. */
	virtual const OtherCharacteristics & GetCharacteristics(size_t Num) const
	{
		return this->ParticleCharacteristics[Num];
	}

	/** \brief Relative (box) coordinates of a particle.
	 *  \param Num Index of the particle. */
	const GeometryVector & GetRelativeCoordinates(size_t Num) const
	{
		return this->ParticleRelatives[Num];
	}

	/** \brief Cartesian coordinates of a particle.
	 *  \param Num Index of the particle. */
	const GeometryVector & GetCartesianCoordinates(size_t Num) const
	{
		return this->ParticleCartesians[Num];
	}

	/** \brief Return the length of the longest diagonal of the simulation box */
	//return the length of the longest diagonal of the simulation box
	double MaxDiagonalLength(void) const
	{
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i]);
			else
				maxdist.MinusFrom(this->BasisVector[i]);
		}
		return std::sqrt(maxdist.Modulus2());
	}
	/** \brief Allocate and generate internal cell lists, so that calling IterateThroughNeighbors with the same Rc is thread-safe
	 *  \param Rc Neighbor cutoff radius the internal lists are prepared for. */
	//allocate and generate internal cell lists, so that calling IterateThroughNeighbors with the same Rc is thread-safe
	void PrepareIterateThroughNeighbors(double Rc) const
	{
		//debug temp
		//std::cerr<<"IterateThroughNeighbors for Rc="<<Rc<<'\n';

		Rc = std::abs(Rc);
		//find the true Rc, which consideres the distance from a particle to a lattice point
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
			else
				maxdist.MinusFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
		}
		double newRc = Rc + std::sqrt(maxdist.Modulus2());
		double newRc2 = newRc*newRc;

		auto IterationEnd = this->VicinityLattices.end();
		//Get new VicinityLattice if necessary
		if (this->SortedVicinityList)
		{
			if ((newRc2)>this->VicinityLatticesRc2)
				this->GetVicinityLattice(newRc*1.1);
			IterationEnd = this->VicinityLattices.end();
		}
		else
		{
			auto ub = this->VicinityLatticePartitions.lower_bound(newRc2);
			auto lb = ub;
			if (this->VicinityLatticePartitions.size()>1)
				lb--;
			if ((newRc2)>this->VicinityLatticesRc2)
			{
				this->GetVicinityLattice(newRc);
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, this->VicinityLattices.size()));
				IterationEnd = this->VicinityLattices.end();
			}
			else if (ub->first == newRc2)
			{
				IterationEnd = ub->second + this->VicinityLattices.begin();
			}
			else if (ub->second - lb->second>10)
			{
				auto titer = std::partition(this->VicinityLattices.begin() + lb->second, this->VicinityLattices.begin() + ub->second,
					[&newRc2](const VicinityLatticesMember & a) ->bool
				{
					return a.CartesianCoordinates.Modulus2()<newRc2;
				});
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, titer - this->VicinityLattices.begin()));
				IterationEnd = titer;
			}
			else
				IterationEnd = ub->second + this->VicinityLattices.begin();
		}
	}

	/** \brief Iterate over neighbors within Rc of a point in relative coordinates.
	 *  \param RelativeCoordinate Center point in relative (box) coordinates.
	 *  \param Rc	Neighbor cutoff radius.
	 *  \param callback	A lambda function repeated for each neighbor with shift, lattice shift, periodic index shift, and source index.
	 *  \param pAbort	An optional flag to abort iteration early when set to true. */
	void IterateThroughNeighbors(const GeometryVector & RelativeCoordinate, double Rc, const std::function<void(const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)> callback, bool * pAbort = nullptr) const
	{
		//debug temp
		//std::cerr<<"IterateThroughNeighbors for Rc="<<Rc<<'\n';

		Rc = std::abs(Rc);
		double Rc2 = Rc*Rc;
		GeometryVector StartRelativeCoordinate = RelativeCoordinate;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			StartRelativeCoordinate.x[i] -= std::floor(StartRelativeCoordinate.x[i]);
			StartRelativeCoordinate.x[i] -= std::floor(StartRelativeCoordinate.x[i]);
		}
		//find the true Rc, which consideres the distance from a particle to a lattice point
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
			else
				maxdist.MinusFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
		}
		double newRc = Rc + std::sqrt(maxdist.Modulus2());
		double newRc2 = newRc*newRc;

		auto IterationEnd = this->VicinityLattices.end();
		//Get new VicinityLattice if necessary
		if (this->SortedVicinityList)
		{
			if ((newRc2)>this->VicinityLatticesRc2)
				this->GetVicinityLattice(newRc*1.1);
			IterationEnd = this->VicinityLattices.end();
		}
		else
		{
			auto ub = this->VicinityLatticePartitions.lower_bound(newRc2);
			auto lb = ub;
			if (this->VicinityLatticePartitions.size()>1)
				lb--;
			if ((newRc2)>this->VicinityLatticesRc2)
			{
				this->GetVicinityLattice(newRc);
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, this->VicinityLattices.size()));
				IterationEnd = this->VicinityLattices.end();
			}
			else if (ub->first == newRc2)
			{
				IterationEnd = ub->second + this->VicinityLattices.begin();
			}
			else if (ub->second - lb->second>10)
			{
				auto titer = std::partition(this->VicinityLattices.begin() + lb->second, this->VicinityLattices.begin() + ub->second,
					[&newRc2](const VicinityLatticesMember & a) ->bool
				{
					return a.CartesianCoordinates.Modulus2()<newRc2;
				});
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, titer - this->VicinityLattices.begin()));
				IterationEnd = titer;
			}
			else
				IterationEnd = ub->second + this->VicinityLattices.begin();
		}

		//Get cartesian coordinate of start point
		GeometryVector StartCartesianCoordinate(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			StartCartesianCoordinate.AddFrom(this->BasisVector[i] * StartRelativeCoordinate.x[i]);

		//Get the (dimensional) index of the cell that the start point is in
		signed long index[::MaxDimension] = {0};
		this->GetFullIndex(index, StartRelativeCoordinate);
		//iterate through cells
		for (auto iter = this->VicinityLattices.begin(); iter != IterationEnd; iter++)
		{
			//the vicinity lattice list can be larger than what's needed
			if (this->SortedVicinityList)
			{
				if (iter->CartesianCoordinates.Modulus2()>newRc2)
					break;
			}
			else
			{
				if (iter->CartesianCoordinates.Modulus2()>newRc2)
					continue;
			}


			signed long newindex[::MaxDimension]={0}; //the dimensional index for the cell
			signed long newshift[::MaxDimension]={0}; //the difference of periodic cell of new particle ant current particle, in terms of relative coordinate
			for (DimensionType i = 0; i<this->Dimension; i++)
			{
				signed long n = index[i] + iter->LatticeShift[i];

				if (n >= 0)
				{
					newindex[i] = (n) % this->CellRank[i];
					newshift[i] = (n) / this->CellRank[i];
				}
				else
				{
					newshift[i] = (-1) - (((-1) - n) / this->CellRank[i]);
					newindex[i] = this->CellRank[i] - 1 - (((-1) - n) % this->CellRank[i]);
				}
			}

			const std::vector<size_t > & thiscell = this->Cells[this->ConvertIndex(newindex)];

			if (thiscell.size() != 0)
			{
				GeometryVector LatticeShift(this->Dimension);
				for (int i = 0; i<this->Dimension; i++)
					LatticeShift.AddFrom(static_cast<double>(newshift[i])*this->BasisVector[i]);
				//iterate through Particles
				for (auto iter2 = thiscell.begin(); iter2 != thiscell.end(); iter2++)
				{
					if (pAbort != nullptr)
						if (*pAbort == true)
							return;

					GeometryVector shift = this->ParticleCartesians[*iter2] - StartCartesianCoordinate;

					shift.AddFrom(LatticeShift);

					if (shift.Modulus2()<Rc2)
						callback(shift, LatticeShift, newshift, *iter2);
				}
			}
		}
	}
	/** \brief Move particle Num to a new relative coordinate, updating its cell assignment.
	 *  \param Num                   Index of the particle to move.
	 *  \param newRelativeCoordinate Target position in relative (box) coordinates. */
	void MoveParticle(size_t Num, const GeometryVector & newRelativeCoordinate)
	{
		Index fromIndex;
		fromIndex = this->GetIndex(Num);
		std::vector<size_t > & currentCell = this->Cells[fromIndex];
		std::vector<size_t >::iterator itera = currentCell.begin();
		for (;; itera++)
		{
			assert(itera != currentCell.end());
			if (*itera == Num)
				break;
		}

		this->ParticleRelatives[Num] = newRelativeCoordinate;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			this->ParticleRelatives[Num].x[i] -= std::floor(this->ParticleRelatives[Num].x[i]);
			this->ParticleRelatives[Num].x[i] -= std::floor(this->ParticleRelatives[Num].x[i]);
		}
		this->UpdateCartesianCoordinates(Num);
		this->UpdateCellList(itera, currentCell);
	}

	/** \brief Remove particle Num by swapping it with the last particle and updating cell lists.
	 *  \param Num Index of the particle to remove. */
	void DeleteParticle(size_t Num)
	{
		//particle * a=&this->Particles[Num];
		Index afromIndex = this->GetIndex(Num);
		std::vector<size_t > & acurrentCell = this->Cells[afromIndex];
		std::vector<size_t >::iterator itera = acurrentCell.begin();
		for (;; itera++)
		{
			assert(itera != acurrentCell.end());
			if (*itera == Num)
				break;
		}
		if (Num == this->NumParticle() - 1)
		{
			this->ParticleCartesians.pop_back();
			this->ParticleRelatives.pop_back();
			this->ParticleCharacteristics.pop_back();
			if (itera == acurrentCell.end() - 1)
				acurrentCell.pop_back();
			else
			{
				(*itera) = acurrentCell.back();
				acurrentCell.pop_back();
			}
			return;
		}

		//particle * back=&this->Particles.back();
		size_t backi = this->NumParticle() - 1;
		Index bfromIndex = this->GetIndex(backi);
		std::vector<size_t > & bcurrentCell = this->Cells[bfromIndex];
		std::vector<size_t >::iterator iterb = bcurrentCell.begin();
		for (;; iterb++)
		{
			assert(iterb != bcurrentCell.end());
			if (*iterb == (backi))
				break;
		}

		//(*a)=this->Particles.back();
		//this->Particles.pop_back();
		this->ParticleCartesians[Num] = this->ParticleCartesians.back();
		this->ParticleRelatives[Num] = this->ParticleRelatives.back();
		this->ParticleCharacteristics[Num] = this->ParticleCharacteristics.back();
		this->ParticleCartesians.pop_back();
		this->ParticleRelatives.pop_back();
		this->ParticleCharacteristics.pop_back();

		(*iterb) = Num;
		if (itera == acurrentCell.end() - 1)
			acurrentCell.pop_back();
		else
		{
			(*itera) = acurrentCell.back();
			acurrentCell.pop_back();
		}
	}

	/** \brief Return the i-th basis vector of the periodic box.
	 *  \param i Index of the basis vector. */
	GeometryVector GetBasisVector(::DimensionType i) const
	{
		return this->BasisVector[i];
	}

	/** \brief Return the i-th reciprocal basis vector, computing and caching it if needed.
	 *  \param i Index of the reciprocal basis vector. */
	GeometryVector GetReciprocalBasisVector(::DimensionType i) const
	{
		if (this->ReciprocalBasisVectorValid[i])
			return this->ReciprocalBasisVector[i];
		else
		{

			GeometryVector result = this->BasisVector[i];

			/////////////////////////////////////////////
			std::vector<GeometryVector> ortho;
			ortho.reserve(this->Dimension);
			for (DimensionType j = 0; j<this->Dimension; j++)
			{
				if (j == i) continue;
				GeometryVector temp = this->BasisVector[j];
				for (auto iter = ortho.begin(); iter != ortho.end(); iter++)
					temp.MinusFrom((*iter)*((temp.Dot(*iter)) / (iter->Modulus2())));
				ortho.push_back(temp);
			}
			////////////////////////////////////////////

			for (auto iter = ortho.begin(); iter != ortho.end(); iter++)
				result.MinusFrom((*iter)*((iter->Dot(result)) / (iter->Modulus2())));

			double product = result.Dot(this->BasisVector[i]);
			result.MultiplyFrom(2 * ::pi / product);
			this->ReciprocalBasisVector[i] = result;
			this->ReciprocalBasisVectorValid[i] = true;
			return result;
		}
	}

	/** \brief Obtain indices of cells neighboring the RelativeCoordinates around radius Rc
	 *  \param[in]  RelativeCoordinates Center point in relative (box) coordinates.
	 *  \param[in]  Rc                  Search radius around the center point.
	 *  \param[out] prt_index_list      Filled with the indices of particles in the neighboring cells. */
	//obtain indices of cells neighboring the RelativeCoordinates around radius Rc
	void GetNeighborList(const GeometryVector & RelativeCoordinates, double Rc, std::vector<int> & prt_index_list) const
	{
		prt_index_list.resize(0);
		DimensionType d = this->Dimension;
		std::vector<GeometryVector> reciprocalBasis;
		for (int i = 0; i < d; i++) {
			if (this->ReciprocalBasisVectorValid[i]) {//reciprocal vectors are already computed
				reciprocalBasis.push_back(this->ReciprocalBasisVector[i]);
			}
			else {//reciprocal vectors are yet to computed
				this->GetReciprocalBasisVector(i);
				this->ReciprocalBasisVector[i] = ((double)this->CellRank[i] / (2.0*::pi)) * ReciprocalBasisVector[i];
				reciprocalBasis.push_back(this->ReciprocalBasisVector[i]);
			}
		}
		
		//the cell within which the given point is located
		//the range of cells to span the sphere of radius Rc
		//Cell[i].x[0] = min. of the cell index in i-th relative coordinates
		//Cell[i].x[1] = max. of ...
		std::vector<GeometryVector> Cell;
		double center = 0;
		double change = 0;
		for (DimensionType i = 0; i < d; i++) {
			center = RelativeCoordinates.x[i] * (double)CellRank[i];
			change = Rc*std::sqrt(reciprocalBasis[i].Modulus2());
			double min = floor(center - change);
			double max = floor(center + change) + 1;
			Cell.push_back(GeometryVector(min, max));
		}

		if (this->NumParticle() > 0) {
			signed long *Indices = new signed long[d];
			for (DimensionType i = 0; i < d; i++) {
				Indices[i] = (long)Cell[i].x[0];
			}

			do {
				signed long *Indices_temp = new signed long[d];
				//Apply the periodic boundary condition
				for (int i = 0; i < d; i++) {
					if (Indices[i] >= CellRank[i]) {
						Indices_temp[i] = Indices[i] - CellRank[i];
					}
					else if (Indices[i] < 0) {
						Indices_temp[i] = Indices[i] + CellRank[i];
					}
					else {
						Indices_temp[i] = Indices[i];
					}
				}

				//convert index in the form of a vector into a long index
				signed long long_index = this->ConvertIndex(Indices_temp);
				for (int i = 0; i < this->Cells[long_index].size(); i++) {
					prt_index_list.push_back(this->Cells[long_index].at(i));
				}

				//to the next cell
				Indices[0]++;
				for (int i = 0; i < d - 1; i++) {
					if (Indices[i] == (long)Cell[i].x[1]) {
						Indices[i] = (long)Cell[i].x[0];
						Indices[i + 1] ++;
					}
				}
			} while (Indices[d - 1] < (long)Cell[d - 1].x[1]);
		}
		

				
	}

	/** \brief Recompute cell ranks from the current basis and reassign all particles to cells.*/
	void RefineCellList(void)
	{
		this->ClearVicinityLatticeList();
		size_t NumCells = 1;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double length = std::sqrt(this->BasisVector[i].Modulus2());
			this->CellRank[i] = (signed long)(std::floor(length / this->ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			NumCells *= this->CellRank[i];
		}
		this->Cells.resize(NumCells);
		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();

		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->Cells[this->GetIndex(i)].push_back(i);
		}
	}

	/** \brief Return the target (approximate) linear cell size. */
	double GetCellSize(void)
	{
		return this->ApproximateCellSize;
	}

	/** \brief Enable or disable use of the sorted vicinity lattice list. 
	 *  \param SortedCellList True to use the sorted vicinity lattice list, false otherwise. */
	void SetSortedCellList(bool SortedCellList)
	{
		if (SortedCellList != this->SortedVicinityList)
			this->ClearVicinityLatticeList();
		this->SortedVicinityList = SortedCellList;
	}

	/** \brief Set the target cell size, refining the cell list if cells become too large or small.
	 *  \param CellSize New target linear cell size. */
	void SetCellSize(double CellSize)
	{
		this->ApproximateCellSize = CellSize;
		bool NeedRefine = false;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double celllength = std::sqrt(this->BasisVector[i].Modulus2()) / this->CellRank[i];
			if (celllength > 2 * this->ApproximateCellSize)
				NeedRefine = true;
			else if (celllength < 0.5 * this->ApproximateCellSize && this->CellRank[i] > 1)
				NeedRefine = true;
		}
		if (NeedRefine)
			this->RefineCellList();
	}

	/** \brief Replace the box basis vectors, updating Cartesian coordinates and cell list.
	 *  \param NewBasisVectors Array of new basis vectors (one per dimension). */
	void ChangeBasisVector(const GeometryVector NewBasisVectors[])
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->BasisVector[i] = NewBasisVectors[i];
		this->ClearVicinityLatticeList();
		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->UpdateCartesianCoordinates(i);
		}

		bool NeedRefine = false;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double celllength = std::sqrt(this->BasisVector[i].Modulus2()) / this->CellRank[i];
			if (celllength > 2 * this->ApproximateCellSize)
				NeedRefine = true;
			else if (celllength < 0.5 * this->ApproximateCellSize && this->CellRank[i] > 1)
				NeedRefine = true;
		}
		if (NeedRefine)
			this->RefineCellList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
		this->vol = -1;
		this->vol = this->PeriodicVolume();
	}

	/** \brief Remove all Particles */
	void RemoveParticles(void)//remove all Particles
	{
		this->ParticleCartesians.clear();
		this->ParticleCharacteristics.clear();
		this->ParticleRelatives.clear();

		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();
	}

	/** \brief Scale the box (and cell size) by Factor, isotropically rescaling all particle positions.
	 *  \param Factor Multiplicative scaling factor applied to the box. */
	virtual void Rescale(double Factor)
	{
		this->ApproximateCellSize *= Factor;
		std::vector<GeometryVector> basis;
		for (DimensionType i = 0; i<this->Dimension; i++)
			basis.push_back(this->BasisVector[i] * Factor);
		this->ChangeBasisVector(&basis[0]);
		/*reset volume*/
		this->vol = -1.0;
	}

	/** \brief Rescale the box isotropically so its volume becomes NewVolume.
	 *  \param NewVolume Desired box volume after rescaling. */
	virtual void Resize(double NewVolume)
	{
		double factor = std::pow(NewVolume / this->PeriodicVolume(), 1.0 / this->Dimension);
		this->Rescale(factor);
	}

	/** \brief Convert Cartesian coordinates to relative (box) coordinates.
	 *  \param CartesianCoord Point in Cartesian coordinates.
	 *  \return               The same point in relative (box) coordinates. */
	GeometryVector CartesianCoord2RelativeCoord(const GeometryVector & CartesianCoord) const
	{
		GeometryVector result(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			result.x[i] = CartesianCoord.Dot(this->GetReciprocalBasisVector(i));

		return result*(0.5 / pi);
	}

	/** \brief Convert relative (box) coordinates to Cartesian coordinates. 
	 *  \param RelativeCoord Point in relative (box) coordinates.
	 *  \return              The same point in Cartesian coordinates. */
	GeometryVector RelativeCoord2CartesianCoord(const GeometryVector & RelativeCoord) const
	{
		GeometryVector result(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			result.AddFrom((this->BasisVector[i])*(RelativeCoord.x[i]));

		return result;
	}

	/** \brief Wrap a relative-coordinate vector into its minimum-image representation in [-0.5, 0.5). 
	 *  \param a Relative-coordinate vector wrapped in place to its minimum image. */
	void RelativeCoordToMinimumImage(GeometryVector & a) const
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			a.x[i] -= std::floor(a.x[i]+0.5);
	}

	/** \brief Write dimension, basis, and particle data to a binary stream. 
	 *  \param ofile Output stream to serialize into. */
	void WriteBinary(std::ostream & ofile) const
	{
		ofile.write((char *)(&this->Dimension), sizeof(this->Dimension));
		ofile.write((char *)(&ApproximateCellSize), sizeof(ApproximateCellSize));
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->BasisVector[i].WriteBinary(ofile, this->Dimension);
		int nbrp = this->NumParticle();
		ofile.write((char *)(&nbrp), sizeof(nbrp));
		for (int i = 0; i<nbrp; i++)
		{
			this->ParticleRelatives[i].WriteBinary(ofile, this->Dimension);
			ofile.write((char *)(&this->ParticleCharacteristics[i]), sizeof(OtherCharacteristics));
		}
	}

	/** \brief Construct a PeriodicCellList by reading a binary file (replaces a ReadBinary method). 
	 *  \param ifile        Input stream to deserialize from.
	 *  \param UseSortedList Choose sorted vicinity list for faster neighbor checks when geometry is fixed. */
	// there is no such thing as void ReadBinary(std::istream & ifile)! Use the following constructor instead!
	PeriodicCellList(std::istream & ifile, bool UseSortedList = false)
	{
		ifile.read((char *)(&this->Dimension), sizeof(this->Dimension));
		ifile.read((char *)(&ApproximateCellSize), sizeof(ApproximateCellSize));
		Index CellCount = 1;//need to calculate overall number of cells
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			this->BasisVector[i].ReadBinary(ifile, this->Dimension);
			GeometryVector & vec = this->BasisVector[i];
			this->CellRank[i] = (signed long)(std::floor(std::sqrt(vec.Modulus2()) / ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			//error check
			if (this->CellRank[i] <= 0)
			{
				std::cerr << "error in Configuration constructor: invalid cell rank. basis vector length squared=" << vec.Modulus2() << ", CellSize=" << ApproximateCellSize << '\n';
				assert(false);
			}
			CellCount *= this->CellRank[i];
		}
		for (DimensionType i = this->Dimension; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = 0;
			this->BasisVector[i] = GeometryVector();
		}
		this->Cells.assign(CellCount, std::vector<size_t >());

		this->ApproximateCellSize = ApproximateCellSize;
		this->SortedVicinityList = UseSortedList;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		int nbrp;
		ifile.read((char *)(&nbrp), sizeof(nbrp));
		for (int i = 0; i<nbrp; i++)
		{
			GeometryVector temp;
			temp.ReadBinary(ifile, this->Dimension);
			OtherCharacteristics otemp;
			ifile.read((char *)(&otemp), sizeof(OtherCharacteristics));
			this->Insert(otemp, temp);
		}
		this->vol = ::Volume(this->BasisVector, this->Dimension);
	}

	/** \brief Rotate the box to a standard lower-triangular orientation with descending basis-vector lengths. */
	//rotate the current structure so that
	//the 1st basis vector is (c11, 0, 0, ...)
	//the 2nd basis vector is (c21, c22, 0, ...)
	// ......
	//AND the basis vector length l1>=l2>=l3...
	void StandardizeOrientation(void)
	{
		::StandardizeOrientation(this->BasisVector, this->Dimension);
		this->ChangeBasisVector(this->BasisVector);
	}

	/** \brief PeriodicCellList< of anything > can be converted to PeriocicCellList<Empty> */
	//PeriodicCellList< of anything > can be converted to PeriocicCellList<Empty>
	operator PeriodicCellList<Empty> const () const
	{
		return PeriodicCellList<Empty>(*this, Empty());
	}

	/** \brief Return the distance from fromRelative to the closest particle
	 *  \param fromRelative Query point in relative (box) coordinates.
	 *  \return             Distance from the query point to the nearest particle. */
	//return the distance from fromRelative to the closest particle
	double NearestParticleDistance(const GeometryVector & fromRelative)
	{
		double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension())*1.5;
		std::vector<GeometryVector> neighbors;
		double l = TypicalLength;
		while (neighbors.size() < 1)
		{
			neighbors.clear();
			(*this).IterateThroughNeighbors(fromRelative, l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l *= 1.5;
		}

		std::partial_sort(neighbors.begin(), neighbors.begin() + 1, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2() < right.Modulus2(); });
		return std::sqrt(neighbors[0].Modulus2());
	}

	/** \brief Return the distance from particle No. i to the closest particle
	 *  \param i Index of the particle whose nearest-neighbor distance is requested.
	 *  \return  Distance from particle i to the nearest other particle. */
	//return the distance from particle No. i to the closest particle
	double NearestParticleDistance(size_t i)
	{
		double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension())*1.5;
		std::vector<GeometryVector> neighbors;
		double l = TypicalLength;
		while (neighbors.size() < 2)
		{
			neighbors.clear();
			(*this).IterateThroughNeighbors((*this).GetRelativeCoordinates(i), l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l *= 1.5;
		}

		std::partial_sort(neighbors.begin(), neighbors.begin() + 2, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2() < right.Modulus2(); });
		return std::sqrt(neighbors[1].Modulus2());
	}

	/** \brief Least image displacement : input vector is in relative coordinates
	 *  \param i           Index of the particle.
	 *  \param r0_Relative Reference point in relative (box) coordinates.
	 *  \return            Minimum-image displacement from r0_Relative to particle i, in Cartesian coordinates. */
	//Least image displacement : input vector is in relative coordinates
	GeometryVector displacement_periodic_relative (int i, const GeometryVector & r0_Relative) const {
		if (i < this->NumParticle()) {
			GeometryVector Displacement = ParticleRelatives[i] - r0_Relative;
			for (int k = 0; k < this->Dimension; k++) {
				if (Displacement.x[k] < -0.5) { Displacement.x[k] += 1.0; }
				else if (Displacement.x[k] >0.5) { Displacement.x[k] -= 1.0; }
			}
			Displacement = this->RelativeCoord2CartesianCoord(Displacement);
			return Displacement;
		}
		else {
			std::cout << "PeriodicCellList::displacement_periodic - i is too large \n";
			return GeometryVector(r0_Relative.Dimension);
		}
	}
	
	/** \brief Least image displacement : input vector is in Cartesian coordinates
	 *  \param i  Index of the particle.
	 *  \param r0 Reference point in Cartesian coordinates.
	 *  \return   Minimum-image displacement from r0 to particle i, in Cartesian coordinates. */
	//input vector is in Cartesian coordinates
	GeometryVector displacement_periodic(int i, const GeometryVector & r0) const {
		if (i < this->NumParticle()) {
			GeometryVector Displacement = ParticleCartesians[i] - r0;
			this->CartesianCoord2RelativeCoord(Displacement);
			for (int k = 0; k < this->Dimension; k++) {
				if (Displacement.x[k] < -0.5) { Displacement.x[k] += 1.0; }
				else if (Displacement.x[k] >0.5) { Displacement.x[k] -= 1.0; }
			}
			Displacement = this->RelativeCoord2CartesianCoord(Displacement);
			return Displacement;
		}
		else {
			std::cout << "PeriodicCellList::displacement_periodic - i is too large \n";
			return GeometryVector(r0.Dimension);
		}
	}

	/** \brief Return a uniformly random Cartesian position within the periodic box. */
	GeometryVector GetRandomPosition() {
		RandomGenerator gen;
		GeometryVector result;
		for (int i = 0; i < MaxDimension; i++) {
			result = result + gen.RandomDouble()*this->BasisVector[i];
		}
		return result;
	}

	/** \brief Return squared distance btw particle i and j
	 *  \param i Index of the first particle.
	 *  \param j Index of the second particle.
	 *  \return  Squared minimum-image distance between particles i and j. */
	//return squared distance btw particle i and j
	double distance2(int i, int j) {
		if (i < this->NumParticle() && j < this->NumParticle()) {
			GeometryVector dr = this->ParticleRelatives[i] - this->ParticleRelatives[j];
			for (int k = 0; k < this->Dimension; k++) {
				if (dr.x[k] < -0.5) { dr.x[k] += 1.0; }
				else if (dr.x[k] >0.5) { dr.x[k] -= 1.0; }
			}
			GeometryVector temp = this->RelativeCoord2CartesianCoord(dr);
			return temp.Dot(temp);
		}
		else {
			std::cout << "wrong index\n";
			return 0;
		}
	}
	
	/** \brief The squared distance between i-th particle and test particle
	 *  \param i    Index of the particle.
	 *  \param test Test point in Cartesian coordinates.
	 *  \return     Squared minimum-image distance between particle i and the test point. */
	//the squared distance between i-th particle and test particle
	double distance2(int i, GeometryVector test) {
		if (i < this->NumParticle()) {
			GeometryVector dr = ParticleCartesians[i] - test;
			this->CartesianCoord2RelativeCoord(dr);
			for (int k = 0; k < this->Dimension; k++) {
				if (dr.x[k] < -0.5) { dr.x[k] += 1.0; }
				else if (dr.x[k] >0.5) { dr.x[k] -= 1.0; }
			}
			GeometryVector temp = this->RelativeCoord2CartesianCoord(dr);
			return temp.Modulus2();
		}
		else {
			std::cout << "PeriodicCellList::distance2 - i is too large \n";
			return 0;
		}
	}

	/** \brief Return the distance between particle i and particle j.
	 *  \param i Index of the first particle.
	 *  \param j Index of the second particle.
	 *  \return  Minimum-image distance between particles i and j. */
	double distance(int i, int j) {
		return sqrt(this->distance2(i, j));
	}

};

/** \brief Overlap all the configurations in that std::vector; require that all the configurations in the std::vector have the same basis vectors
 *  \param components Configurations to overlap, all sharing the same basis vectors.
 *  \return           A single list containing the particles of every input configuration. */
//overlap all the configurations in that std::vector
//require that all the configurations in the std::vector have the same basis vectors
template<typename a> PeriodicCellList<a> CombinedList(const std::vector<PeriodicCellList<a> > & components)
{
	if(components.size()==0)
		return PeriodicCellList<a>();
	PeriodicCellList<a> result = components[0];
	for(size_t i=1; i<components.size(); i++)
	{
		const PeriodicCellList<a> & list = components[i];
		for(size_t j=0; j<list.NumParticle(); j++)
			result.Insert(list.GetCharacteristics(j), list.GetRelativeCoordinates(j));
	}
	return result;
}

/** \brief config is a CombinedList of m PeriodicCellLists, extract the nth list
 *  \param config Combined list assembled from m equal-sized sub-lists.
 *  \param m      Total number of sub-lists combined into config.
 *  \param n      Index of the sub-list to extract.
 *  \return       The extracted n-th sub-list. */
//config is a CombinedList of m PeriodicCellLists, extract the nth list
template<typename a>PeriodicCellList<a> SeparateList(const PeriodicCellList<a> & config, size_t m, size_t n)
{
	size_t nt=config.NumParticle();
	if(nt%m!=0)
	{
		std::cerr<<"Error in SeparateList : the total number of particles is not multiple of the number of configurations!\n";
		assert(false);
	}
	size_t np=nt/m;
	size_t start=n*np;
	size_t end=start+np;
	PeriodicCellList<a> result = config;
	result.RemoveParticles();
	for(size_t i=start; i<end; i++)
	{
		result.Insert(config.GetCharacteristics(i), config.GetRelativeCoordinates(i));
	}
	return result;
}

/** \brief Per-particle characteristic holding a short (up to 3-character) atom name.  */
struct AtomName
{
	char name[4];
	/** \brief Construct an empty atom name.  */
	AtomName()
	{
		this->name[0]='\0';
	}

	/** \brief Construct from a C-string, truncating to 3 characters.
	 *  \param n Source C-string atom name. */
	AtomName(const char * n)
	{
		std::strncpy(this->name, n, 3);
		this->name[3]='\0';
	}
};
typedef PeriodicCellList<AtomName> Configuration;

#include <complex>     // PeriodicCellList may be instantiated on std::complex (e.g. weighted S(k))
class SpherePacking;    // fwd decl: ConfigurationIO.h (pulled in by the umbrella below) refers to it

/** \brief Tile a Configuration NumberInEachSide times along every box direction.
 *  \param src             Source configuration to tile.
 *  \param NumberInEachSide Number of replicas along each box direction.
 *  \return                The tiled configuration. */

 Configuration MultiplicateStructure(const Configuration & src, size_t NumberInEachSide);

 /** \brief Tile a Configuration with per-direction replication counts. 
 *  \param src             Source configuration to tile.
 *  \param NumberInEachSide Per-direction number of replicas.
 *  \return                The tiled configuration. */
Configuration MultiplicateStructure(const Configuration & src, const std::vector<size_t> & NumberInEachSide);

/** \brief Return an empty unit cubic box configuration of dimension d. 
 *  \param d        Spatial dimension of the box.
 *  \param CellSize Target linear cell size.
 *  \return         An empty unit cubic box configuration. */
Configuration GetUnitCubicBox(DimensionType d, double CellSize=0.1);

/** \brief A pack of multiple configurations, recorded in files Prefix.ConfigPack and Prefix.ConfigPackIndex */
// a pack of multiple configurations, recorded in files Prefix.ConfigPack and Prefix.ConfigPackIndex
// .ConfigPack file contains configurations in binary format, one by one
// .ConfigPackIndex has the format:
// long long Number of Config, long long location of config0, long long location of config1, ...
class ConfigurationPack
{
private:
	std::string IndexName, PackName;
	long long NConfig;
public:
	/** \brief Construct an empty pack with no associated files. */
	ConfigurationPack()
	{
		this->NConfig=0;
	}

	/** \brief Construct and open a pack using the given file prefix. 
	 *  \param prefix File name prefix identifying the pack files. */
	ConfigurationPack(const std::string & prefix)
	{
		this->Open(prefix);
	}

	/** \brief Open (or attach to) the pack files identified by prefix.
	 *  \param prefix File name prefix identifying the pack files. */
	void Open(const std::string & prefix);

	/** \brief Return the number of configurations currently in the pack. */
	long long NumConfig(void);

	/** \brief Append a configuration to the pack files.
	 *  \param c Configuration to append. */
	void AddConfig(const Configuration & c);

	/** \brief Read and return the i-th configuration from the pack.
	 *  \param i Index of the configuration to read.
	 *  \return  The i-th configuration in the pack. */
	Configuration GetConfig(long long i);

	/** \brief Discard all configurations that are already in the file */
	//discard all configurations that are already in the file
	void Clear(void);
};

//VoxelType should be either char or bool. Char is faster, while bool can work with smaller memory.
typedef char VoxelType;
/** \brief Pixelizing the configuration to MeshSide pixes in each side
 *  \param[in]  c        Configuration to pixelize.
 *  \param[in]  radius   Sphere radius placed on each particle.
 *  \param[in]  MeshSide Number of pixels along each box direction.
 *  \param[out] result   Filled with the occupancy grid (1 = occupied). */
//pixelizing the configuration to MeshSide pixes in each side
void DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> & result);

/** \brief Pixelize the configuration and return the resulting occupancy grid.
 *  \param c        Configuration to pixelize.
 *  \param radius   Sphere radius placed on each particle.
 *  \param MeshSide Number of pixels along each box direction.
 *  \return         The occupancy grid (1 = occupied). */
inline std::vector<VoxelType> DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide)
{
	std::vector<VoxelType> result;
	DigitizeConfiguration(c, radius, MeshSide, result);
	return result;
}

/** \brief Calculate volume fraction if each particle is replaced with a hypersphere of certain radius; do this by pixelizing the configuration to MeshSide pixes in each side
 *  \param c        Configuration whose occupied volume is measured.
 *  \param radius   Sphere radius placed on each particle.
 *  \param MeshSide Number of pixels along each box direction.
 *  \param pbuffer  Optional scratch buffer for the occupancy grid; allocated internally if null.
 *  \return         The occupied (covered) volume. */
//calculate volume fraction if each particle is replaced with a hypersphere of certain radius
//do this by pixelizing the configuration to MeshSide pixes in each side
double Volume(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> * pbuffer = nullptr);


/** \brief Return the minimum D such that if each particle is replaced with a sphere with diameter D, particles a and b will be connected.
 *  \param c Configuration of particles.
 *  \param a Index of the first particle.
 *  \param b Index of the second particle.
 *  \return  Minimum sphere diameter D that connects particles a and b. */
//return the minimum D such that if each particle is replaced with a sphere with diameter D, particles a and b will be connected.
double MinConnectDiameter(const Configuration & c, size_t a, size_t b);
/** \brief Return the minimum D such that if each particle is replaced with a sphere with diameter D, particle n will connect to its periodic image in direction dir.
 *  \param c   Configuration of particles.
 *  \param n   Index of the particle.
 *  \param dir Box direction along which percolation is tested.
 *  \return    Minimum sphere diameter D that connects particle n to its periodic image. */
//return the minimum D such that if each particle is replaced with a sphere with diameter D,
//particle n will connect to its periodic image in direction dir.
double PercolationDiameter(const Configuration & c, size_t n, DimensionType dir);


/** \brief Replace each particle in c with a hypersphere with radius, what's the cluster size InitParticle is in? FoundPeriodicImages: whether or not found more than 1 periodic images of a particle, which means the cluster size is actually infinite rather than the returned value.
 *  \param[in]  c                  Configuration of particles.
 *  \param[in]  InitParticle       Index of the seed particle whose cluster is measured.
 *  \param[in]  radius             Sphere radius placed on each particle.
 *  \param[out] FoundPeriodicImages Set true if more than one periodic image of a particle was found (cluster is effectively infinite).
 *  \return                        Number of particles in the cluster containing InitParticle. */
//Replace each particle in c with a hypersphere with radius,
//what's the cluster size InitParticle is in?
//FoundPeriodicImages: whether or not found more than 1 periodic images of a particle,
//which means the cluster size is actually infinite rather than the returned value.
size_t ClusterSize(const Configuration & c, size_t InitParticle, double radius, bool & FoundPeriodicImages);

/** Find all reciprocal lattice vectors (wavevectors) of a given periodic simulation box.
 *  It is modified to deal with errors.
 *	\param tempList	A periodic simulation box.
 *	\param (CircularKMax, LinearKMax)	Find all reciprocal lattice vectors (excluding parity pairs) whose magnitude is smaller than CircularKMax but greater than CircularKMin, and then fina all integral multiples of the aforementioned vectors whose magnitudes are smaller than LinearKMax.
 *	\param SampleProbability	When <1, randomly choose wavevectors within CircularKMax by the prescribed probability.
 *	\param use_iteratethrough	True = use IterateThrough() function, which may surge memory use. 
 *	\param CircularKMin			0 by default.  
 */

std::vector<GeometryVector> GetKs(const PeriodicCellList<Empty> & tempList, double CircularKMax, double LinearKMax, double SampleProbability = 1.0, bool use_iteratethrough = true, double CircularKMin = 0.0);

// ---------------------------------------------------------------------------
// Backward-compatibility umbrella.  The classes/functions below were split out
// of this file for readability; they are re-included here so existing code
// doing  #include "PeriodicCellList.h"  keeps seeing SpherePacking,
// ClusterMerger and all I/O helpers unchanged.  Include guards make the
// resulting mutual includes safe in any inclusion order.
// ---------------------------------------------------------------------------
#include "SpherePacking.h"
#include "ClusterMerger.h"
#include "ConfigurationIO.h"

#endif
