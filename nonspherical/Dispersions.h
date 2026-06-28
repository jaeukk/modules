/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	Feb. 2019 */


 /** \file Dispersions.h
  * \brief Header file for a Dispersion class 
  * (packings of aribtrarily shaped particles under the periodic boundary conditions).
  *	This is a general extension of PeriodicCellList<T> class, defined in PeriodicCellList.h. 
  * 
  * ------------------------------
  *       To do list
  * ------------------------------
  * Member functions to check whether particles are overlapped or not.
  * Such functions are necessarily to perform the RSA procedure with aribtrarily shaped particles. */


#ifndef DISPERSIONS_H
#define DISPERSIONS_H

#include <memory>	// To use shared_ptr
#include <complex>
#include <algorithm>
#include <bits/stdc++.h>

#include <PeriodicCellList.h>
#include "ParticleShapes.h"
#include <GeometryVector.h>
#include <etc.h>

//< Pointer of the Shape classes = type of particle characteristics of PeriodicCellLsit.	
// Since it uses pointers, users should be very catious when using this class.
//typedef Shape * Shape_;
typedef std::shared_ptr<Shape> Shape_;


/** \class Dispersion class stores information on particle shapes as the shared pointer of Shape classes. 
 *	Particle coordinates describe particle centroids.
 *	Shape_ describes the shape of a certain particle with respect to its "centroid." 
 *  Dispersion class is a friend class of PeriodicCellList, which allows to access its protected members. */
class Dispersion : public PeriodicCellList<Shape_>	
{
	//template<typename a> friend class PeriodicCellList;
public:
	
	/** A default constructor. */
	Dispersion() : PeriodicCellList() {};

	/** A constructor of an empty simulation box.
	 * \param Dimension	Space dimension.
	 * \param BasisVector	An array fo GeometryVectors that describe the basis vectors of the simulation box.
	 * \param cellSize		An approximate size to generate a cell list.
							In each direction, an appropriate cell size is chosen
							such that there are an integral number of cells.
	 * \param UseSortedList	True => Sort Vincinity Cell List to increase efficiency in IterateThroughNeighbors()
							at the expense of time in sorting the list.	 */
	Dispersion(DimensionType Dimension, GeometryVector * BasisVector, double cellSize, bool UseSortedList = false) 
		: PeriodicCellList<Shape_>(Dimension, BasisVector, cellSize, UseSortedList)
	{	this->MaxRadius = 0.0;	}

	/** A constructor from particle centroids in source, and an identical shape in prtShape.
	 * \param source	A realization of a periodic system. Its periodic simulation box and particle centers are taken.
	 * \param prtShape	A pointer of the Shape class that defines a particle shape. */
	template<typename add> Dispersion(const PeriodicCellList<add> & source, const Shape_ & prtShape) : PeriodicCellList(source, prtShape->clone()) {
		//This step is necessary to make sure all pointers of Shape class are independant.
		for (int i = 0; i < source.NumParticle(); i++)
			this->ParticleCharacteristics[i] = prtShape->clone();
	}
	

	/** Directly access to a pointer of a designated Shape object.
	 * However, since return value is a constant pointer, 
	 the user cannot change stored values but it is efficient.
	 * This is an alternative of GetCharacteristics(). Please do not use GetCharacteristics().
	 * \param i	The particle index.
	 * \return A constant shared pointer of the Shape object of particle i. */
	const Shape_ GetShape(size_t i) const {
		if (i >= this->NumParticle())
			throw std::out_of_range("Dispersion::GetShape(size_t)");

		//TODO: add a separate implementation for each derived struct of Shape.
		return this->ParticleCharacteristics[i]->clone();
	}

	/** Give a particle shape.
	 * The returned value is an identical copy of this->ParticleCharacteristics[i], not itself.
	 * Please do not use GetCharacteristics().
	 * \param i	The particle index.
	 * \return A shared pointer of a copied Shape object of particle i. */
	Shape_ GetShape(size_t i) {
		if (i >= this->NumParticle())
			throw std::out_of_range("Dispersion::GetShape(size_t)");
		
		return this->ParticleCharacteristics[i]->clone();
	}

	/** Please do not use. */
	/*virtual const Shape_ & GetCharacteristics(size_t i) const { 
		std::cerr << "Please do not use Dispersion::GetCharacteristics()" << std::endl;
		return this->ParticleCharacteristics[i]->clone();
	}*/
	virtual const Shape_ & GetCharacteristics(size_t i) const {
		std::cerr << "Please do not use Dispersion::GetCharacteristics()" << std::endl;
		return this->ParticleCharacteristics[i];
		//return Shape_(this->ParticleCharacteristics[i]);
	}

	/** Please do not use. */
	virtual Shape_ & GetCharacteristics(size_t i) {
		std::cerr << "Please do not use Dispersion::GetCharacteristics()" << std::endl;
		return this->ParticleCharacteristics[i];
	}
	
	/** Update this->MaxRadius.
	 *  This->MaxRadius = the maximum of R_circum of all particles. */
	void UpdateMaxRadius(void) {
		if (this->MaxRadius == 0.0)	{
			for (size_t i = 0; i < this->NumParticle(); i++)
				this->MaxRadius = std::max(this->MaxRadius, ParticleCharacteristics[i]->GetMaxRadius());
		}
	}
	/** Update this->MaxRadius. */
	void UpdateMaxRadius(void) const {
		if (this->MaxRadius == 0.0) {
			for (size_t i = 0; i < this->NumParticle(); i++)
				this->MaxRadius = std::max(this->MaxRadius, ParticleCharacteristics[i]->GetMaxRadius());
		}
	}

	/** \return The maximal radius of a circumscribing sphere. */
	double GetMaxRadius(void) {
		UpdateMaxRadius();	return this->MaxRadius;
	};

	/** \return The maximal radius of a circumscribing sphere. */
	double GetMaxRadius(void) const {
		UpdateMaxRadius();	return this->MaxRadius;
	};

	/** Insert a particle inside the periodic simulation box.
	 * Position of an added particle is automatically translated to the cell list.
	 * \param ptrShape_	A pointer of particle shape.
	 * \param RelativeCoordinate	A GeometryVector of the particle centroid in relative coordinates.	 */
	virtual void Insert(const Shape_ & prtShape_, const GeometryVector & RelativeCoordinate) {
		size_t idx = this->NumParticle();

		/* Particle position = "prtSHape_::Centroid" + "RelativeCoordinate".  */
		GeometryVector Pos_Relative = RelativeCoordinate + this->CartesianCoord2RelativeCoord(prtShape_->GetCenterOfMass());

		PeriodicCellList::Insert(prtShape_, Pos_Relative);
		this->MaxRadius = std::max(this->MaxRadius, prtShape_->GetMaxRadius());

		/** Place the centroid of the added particle is at (0,0,0).*/
		this->ParticleCharacteristics[idx]->ShiftBy(-1.0 * prtShape_->GetCenterOfMass());
	}

	/** Insert a particle of a given shape at a random position.
	 * \param prtShape_	A shared pointer of Shape object.
	 * \param rng	A random number generator.	 */
	virtual void Insert(const Shape_ & prtShape_, RandomGenerator & rng) {
		GeometryVector pos(this->Dimension);
		for (int i = 0; i < pos.Dimension; i++)
			pos.x[i] = rng.RandomDouble();
		this->Insert(prtShape_, pos);
	}

	/** Check whether a test point is inside a particle of this Dispersion class.
	 * \param[in] RelativeCoordinate 	The position of a test point in the relative coordinates.
	 * \return Particle ID, if the point is inside of a particle. -1, otherwise. */
	int IsInsideParticle(const GeometryVector & RelativeCoordinate) const;

	/** Find the three nearest neighbors.
	 * \param[in] RelativeCoordinate	The position of a test point in the relative coordinates.
	 * \param[out] list	A list of some nearest interfaces. list[i] = GeometryVector(distance, side ID, particle ID)
	 * \param[out] x0_Relative	The point on the interface closest to RelativeCoordinate.
	 * \param[out] n2host	The unit outward normal vector at x0. 	*/
	void FindNearestNeighbors (const GeometryVector & RelativeCoordinate, std::vector<GeometryVector> & list, GeometryVector & x0_Relative, GeometryVector & n2host) const;

	/** A default destructor. */
	~Dispersion() {}
	
	/** Rotate a single particle.
	 * \param i	Particle index.
	 * \param newbasis	A set of orthonormal basis vectors to which particle i will be rotated. */
	virtual void RotateParticle(size_t i, const std::vector<GeometryVector> & newbasis) {
		ParticleCharacteristics[i]->Rotate(newbasis);
	}
	
	/** Rotate a single particle randomly.
	 * \param i Particle index.
	 * \param rng	A random number generator. */
	virtual void RotateParticle(size_t i, RandomGenerator & rng) {
		ParticleCharacteristics[i]->Rotate(rng);
	}
	
	/** Rotate all particle randomly.
	 * \param seed A random seed for a random generator. */
	virtual void RotateAllParticles(int seed = 0) {
		RandomGenerator rng(seed);
		for (size_t i = 0; i < this->NumParticle(); i++)
			this->RotateParticle(i, rng);
	}


	/** Give the packing fraction, assuming particles are not overlapping.
	 * \return Packing fraction. */
	virtual double PackingFraction() const;
	
	/** Rescale all particles and simulation box.
	 * \param LinearFactor	The linear ratio to which everything is scaled. */
	virtual void Rescale(double LinearFactor);
	/** Rescale all particles and simulation box.
	 * \param newVolume	The volume of periodic simulation box to which everything is scaled. */
	virtual void Resize(double NewVolume);

	/** Change packing fraction by scaling only particles.
 	 * \param LinearFactor	The linear ratio to which every particle is scaled. */
	void Rescale_prtOnly(double LinearFactor);	

	/** Change Particle Characteristics.
	 * \param newShape	A shared pointer of Shape class that describes a new particle shape.
	 * \param i	Particle index. */
	void UpdateParticleShape(const Shape_ & newShape, size_t i); 
//	//To implement RSA procedure, CheckVoxelList() function should be defined!!

	/** Compute the spectral density of a particle i at a wavevector k.
	 * \param i	Particle ID.
	 * \param k	Wavevector.
	 * \return Form factor!!
	*/
	inline std::complex<double> GetFormFactor(size_t i, const GeometryVector & k) const {
		return this->ParticleCharacteristics[i]->GetFormFactor(k);	
	}
};


// txt file I/O.

/** Save a Dispersion in a text file.
 * \param conf	A Dipsersion object to save.
 * \param OutputName	File name. */
void WriteDispersion(const Dispersion & conf, const std::string & OutputName);

/** Load a text file into a Dispersion object.
 * \param conf	A Dipsersion object that stores data.
 * \param InputName	File name. */
void ReadDispersion(Dispersion & conf, const std::string & InputName);




///*  cases where all particles are in parallelepiped */
//class Dispersion_parallelepiped : public PeriodicCellList<Parallelepiped> 
//{
//	Dispersion_parallelepiped();
//	/* identical particles */
//	template<typename add> Dispersion_parallelepiped(const PeriodicCellList<add> & source, const Parallelepiped & prtShape);
//	Dispersion_parallelepiped(DimensionType Dimension, GeometryVector * BasisVector, double cellSize, bool UseSortedList = false);
//	virtual Parallelepiped & GetCharacteristics(size_t Num);
//	virtual const Parallelepiped & GetCharacteristics(size_t Num) const;
//
//	virtual void Insert(const Parallelepiped & prtShape, const GeometryVector & RelativeCoordinate);
//	virtual void Insert(const GeometryVector & sideLengths, const GeometryVector & RelativeCoordinate);
//
//	virtual void RotateParticle(size_t i, const std::vector<GeometryVector> & newbasis);	// rotate particle i with a new orientation
//	virtual void RotateParticle(size_t i);													// rotate particle i with a random orientation
//	virtual void RotateAllParticles();														// rotate all particles in random oriengtations
//};
//
////typedef shared_ptr<Cube> Cube_;
//class Dispersion_cube : public PeriodicCellList<Cube> 
//{
//	//template<typename a> friend class PeriodicCellList;
//
//public:
//	Dispersion_cube() : PeriodicCellList<Cube> (){};
//	/* identical particles */
//	template<typename add> 
//	Dispersion_cube(const PeriodicCellList<add> & source, Cube & prtShape) : PeriodicCellList(source, prtShape) {
//		this->MaxRadius = prtShape.GetMaxRadius();
//	}
//	Dispersion_cube(DimensionType Dimension, GeometryVector * BasisVector, double cellSize, bool UseSortedList = false);
//	
//	virtual Cube & GetCharacteristics(size_t Num);
//	virtual const Cube & GetCharacteristics(size_t Num) const;
//
//	virtual void Insert(const Cube & prtShape, const GeometryVector & RelativeCoordinate);
//	virtual void Insert(const double sideLength, const GeometryVector & RelativeCoordinate) {
//		this->Insert(Cube(RelativeCoordinate.Dimension, sideLength), RelativeCoordinate);
//	}
//
//	void UpdateMaxRadius(void) const;
//	double GetMaxRadius(void) const;
//
//	virtual void RotateParticle(size_t i, const std::vector<GeometryVector> & newbasis);	// rotate particle i with a new orientation
//	virtual void RotateParticle(size_t i, RandomGenerator & rng);													// rotate particle i with a random orientation
//	virtual void RotateAllParticles();														// rotate all particles in random oriengtations
//
//	virtual double PackingFraction() const;
//	virtual void Rescale(double LinearFactor);			//Change simulation box and particle at a certain linear ratio
//	virtual void Resize(double NewVolume);				//Change simulation box to a NewVolume, scale particle size accordingly
//	virtual void Rescale_prtOnly(double LinearFactor);	//Change particle sizes by a factor of ScaleFactor
//	virtual void UpdateCharacteristics(const Cube & newShape, size_t i);	//Replace particle i with newShape 
//																			//To implement RSA procedure, CheckVoxelList() function should be defined!!
//
//};

#endif // !DISPERSIONS_H
