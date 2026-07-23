#ifndef SpherePacking_INCLUDED
#define SpherePacking_INCLUDED

/** \file SpherePacking.h
 *  \brief Polydisperse sphere packing (SpherePacking, derived from
 *         PeriodicCellList<double>) and its on-disk pack container
 *         SpherePackingPack.  Split out of PeriodicCellList.h for readability. */

#include "PeriodicCellList.h"   // base PeriodicCellList<double>, Configuration

/** \brief Polydisperse sphere packing : a periodic point pattern whose per-particle characteristic is a sphere radius */
//polydisperse sphere packing
class SpherePacking : public PeriodicCellList<double>
{
public:
	/** \brief Construct an empty sphere packing. */
	SpherePacking() : PeriodicCellList()//, MaxRadius(0.0)
	{
	}

	/** \brief Construct an empty sphere packing for a given box geometry.
	 *  \param Dimension    Spatial dimension.
	 *  \param BasisVectors Basis vectors defining the periodic box.
	 *  \param CellSize     Target linear cell size.
	 *  \param UseSortedList Choose sorted vicinity list for faster neighbor checks when geometry is fixed. */
	SpherePacking(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList=false);

	/** \brief Construct from another list, assigning a uniform sphere radius to every particle.
	 *  \param source Source list providing geometry and particle positions.
	 *  \param radius Sphere radius assigned to every particle. */
	template<typename add> SpherePacking(const PeriodicCellList<add> & source, const double & radius) : PeriodicCellList(source, radius)
	{
		this->MaxRadius=radius;
	}

	/** \brief Access the mutable radius of sphere Num.
	 *  \param Num Index of the sphere. */
	virtual double & GetCharacteristics(size_t Num)
	{
		//mutable access may enlarge a radius; invalidate the cached MaxRadius so
		//UpdateMaxRadius() recomputes it before the next neighbor-cutoff use.
		this->MaxRadius = 0.0;
		return this->ParticleCharacteristics[Num];
	}

	/** \brief Access the const radius of sphere Num.
	 *  \param Num Index of the sphere. */
	virtual const double & GetCharacteristics(size_t Num) const
	{
		return this->ParticleCharacteristics[Num];
	}

	/** \brief Construct a sphere packing through a binary stream.
	 *  \param ifile        Input stream to deserialize from.
	 *  \param UseSortedList Choose sorted vicinity list for faster neighbor checks when geometry is fixed. */
	SpherePacking(std::istream & ifile, bool UseSortedList=false) : PeriodicCellList<double>(ifile, UseSortedList)
	{
		this->MaxRadius=0.0;
		for(size_t i=0; i<this->NumParticle(); i++)
			this->MaxRadius=std::max(this->MaxRadius, this->GetCharacteristics(i));
	}

	/** \brief Recompute the cached maximum sphere radius if it is currently invalid. */
	void UpdateMaxRadius(void) const
	{
		if (this->MaxRadius == 0.0)
		{
			for (size_t i = 0; i < this->NumParticle(); i++)
				this->MaxRadius = std::max(this->MaxRadius, this->ParticleCharacteristics[i]);
		}
	}

	/** \brief Return the largest sphere radius in the packing. */
	double GetMaxRadius(void) const
	{
		UpdateMaxRadius();
		return this->MaxRadius;
	}

	/** \brief Insert a sphere of the given radius at a relative coordinate.
	 *  \param radius            Radius of the sphere to insert.
	 *  \param RelativeCoordinate Position of the sphere center in relative (box) coordinates. */
	virtual void Insert(double radius, const GeometryVector & RelativeCoordinate);

	/** \brief Check if the voxel is COMPLETELY occupied
	 *  \param cr       Voxel center in relative (box) coordinates.
	 *  \param halfsize Half the voxel side length (in relative units).
	 *  \param Radius   Cartesian radius of the upcoming sphere.
	 *  \return         True if the voxel is completely covered by existing spheres. */
	//check if the voxel is COMPLETELY occupied
	//no need to distinguish partial availability and full availability
	//relative coord and halfsize
	//Radius : cartesian radius of the upcoming sphere
	bool CheckVoxelOverlap(const GeometryVector & cr, double halfsize, double Radius) const;

	/** \brief Returns 0 if voxel is completely empty, 2 if voxel is fully occupied, 1 otherwise
	 *  \param cr       Voxel center in relative (box) coordinates.
	 *  \param halfsize Half the voxel side length (in relative units).
	 *  \param Radius   Cartesian radius of the upcoming sphere.
	 *  \return         0 if empty, 2 if fully occupied, 1 if partially occupied. */
	//returns 0 if voxel is completely empty
	//2 if voxel is fully occupied
	//1 otherwise
	//relative coord and halfsize
	//Radius : cartesian radius of the upcoming sphere
	int CheckVoxelOverlap2(const GeometryVector & cr, double halfsize, double Radius) const;


	/** \brief Check whether a sphere of the given radius placed at cr overlaps any existing sphere.
	 *  \param cr     Sphere center in relative (box) coordinates.
	 *  \param Radius Cartesian radius of the upcoming sphere.
	 *  \return       True if the new sphere overlaps any existing sphere. */
	//relative coord
	//Radius : cartesian radius of the upcoming sphere
	bool CheckOverlap(const GeometryVector & cr, double Radius) const;

	/** \brief Check if there are overlapping spheres.
	 *  \return false = no overlapping. */
	//check if there are overlapping spheres
	/* \return false = no overlapping. */
	bool CheckOverlap(void) const
	{
		UpdateMaxRadius();
		for (size_t i = 0; i < this->NumParticle(); i++)
		{
			bool Overlap = false;
			const std::vector<double> * pSphereRadii = &this->ParticleCharacteristics;
			double Radius = this->GetCharacteristics(i);
			this->IterateThroughNeighbors(this->GetRelativeCoordinates(i), this->GetCharacteristics(i) + MaxRadius, [&Overlap, &Radius, &pSphereRadii, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long *, const size_t Sourceparticle) -> void
			{
				if (Sourceparticle != i || LatticeShift.Modulus2() != 0)
				{
					double OverlapRadius = Radius + pSphereRadii->at(Sourceparticle);
					if (shift.Modulus2() < OverlapRadius*OverlapRadius)
						Overlap = true;
				}
			}, &Overlap);
			if (Overlap)
				return true;
		}
		return false;
	}

	/** \brief Return the packing fraction (total sphere volume divided by box volume). */
	double PackingFraction(void) const
	{
		double v1=0.0;
		DimensionType d=(DimensionType) this->GetDimension();
		for(size_t i=0; i<this->NumParticle(); i++)
			v1+= ::HyperSphere_Volume(d, this->GetCharacteristics(i));
		return v1/this->PeriodicVolume();
	}

	/** \brief Rescale system box and particle sizes....
	 *  \param Factor Multiplicative scaling factor applied to the box and sphere radii. */
	//-----rescale system box and particle sizes....
	virtual void Rescale(double Factor)
	{
		this->ApproximateCellSize *= Factor;
		this->MaxRadius *= Factor;
		std::vector<GeometryVector> basis;
		for (DimensionType i = 0; i<this->Dimension; i++)
			basis.push_back(this->BasisVector[i] * Factor);
		this->ChangeBasisVector(&basis[0]);
		
		for (size_t i = 0; i < this->NumParticle(); i++) {
			this->ParticleCharacteristics[i] *= Factor;
		}
		/* inform that the volume is updated. */
		this->vol = -1.0;
	}

	/** \brief Rescale the box and spheres isotropically so the box volume becomes NewVolume.
	 *  \param NewVolume Desired box volume after rescaling. */
	virtual void Resize(double NewVolume)
	{
		double factor = std::pow(NewVolume / this->PeriodicVolume(), 1.0 / this->Dimension);
		this->Rescale(factor);
	}

	/** \brief Rescale particle sizes only
	 *  \param factor Multiplicative scaling factor applied to sphere radii only. */
	//Rescale particle sizes only
	void Rescale_prtOnly(double factor) {
		for (size_t i = 0; i < this->NumParticle(); i++) {
			this->ParticleCharacteristics[i] *= factor;
		}
		this->MaxRadius *= factor;
	}

	/** \brief Set the radius of particle index to newRadius, updating the cached maximum radius.
	 *  \param newRadius New radius assigned to the sphere.
	 *  \param index     Index of the sphere to update. */
	void UpdateCharacteristics(double newRadius, int index) {
		this->ParticleCharacteristics[index] = newRadius;
		this->MaxRadius = std::max(this->MaxRadius, newRadius);
	}
	
	// Newly defined for First-passage time simulations

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

};

/** \brief A SpherePacking counterpart of ConfigurationPack, using .SpherePack and .SpherePackIndex files */
/*	 a SpherePacking counterpart of ConfigurationPack
*.SpherePack files and .SpherePackIndex
*/
class SpherePackingPack
{
private:
	std::string IndexName, PackName;
	long long NConfig;
public:
	/** \brief Construct an empty pack with no associated files.*/
	SpherePackingPack() { this->NConfig = 0; }

	/** \brief Construct and open a pack using the given file prefix.
	 *  \param prefix File name prefix identifying the pack files. */
	SpherePackingPack(const std::string & prefix) {
		this->Open(prefix);
	}

	/** \brief Open (or attach to) the pack files identified by prefix.
	 *  \param prefix File name prefix identifying the pack files. */

	 void Open(const std::string &prefix);

	 /** \brief Return the number of sphere packings currently in the pack. */
	long long NumConfig(void);

	/** \brief Append a sphere packing to the pack files.
	 *  \param c Sphere packing to append. */
	void AddConfig(const SpherePacking &c);

	/** \brief Read and return the i-th sphere packing from the pack. 
	 *  \param i Index of the sphere packing to read.
	 *  \return  The i-th sphere packing in the pack. */
	SpherePacking GetConfig(long long i);

	/** \brief Discard all sphere packings currently stored in the files. */
	void Clear(void);
};

#endif
