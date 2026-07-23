/** \file SpherePacking.cpp
 *  \brief Implementations of SpherePacking and SpherePackingPack member
 *         functions declared in SpherePacking.h.  Split out of
 *         PeriodicCellList.cpp for readability. */
#include "SpherePacking.h"
#include <exception>

// Implementations of member functions of the SpherePacking class
// for First-Passage-Time simulations
/** \brief Return the ID of the particle containing the given relative coordinate, or -1 if none.
 *  \param RelativeCoordinate Test point in relative (box) coordinates.
 *  \return                   ID of the containing particle, or -1 if the point is outside every sphere. */
int SpherePacking::IsInsideParticle(const GeometryVector & RelativeCoordinate) const
{
	int idx = -1;
	bool inside = false;
	const std::vector<double> * pSphereRadii = &this->ParticleCharacteristics;

	this->IterateThroughNeighbors(RelativeCoordinate, this->GetMaxRadius(), 
	[&pSphereRadii, &inside, &idx, this] (const GeometryVector & shift, const GeometryVector &, const signed long *, const size_t Sourceparticle) -> void{
		//GeometryVector pos_Cartesian = this->RelativeCoord2CartesianCoord(shift);
		
//		double r = pSphereRadii->at(Sourceparticle)
//		double R = std::sqrt(shift.Modulus2());

		if(pSphereRadii->at(Sourceparticle)*pSphereRadii->at(Sourceparticle) > shift.Modulus2() ){
			idx = Sourceparticle;
			inside = true;
		}
	});
	return idx;
}

/** \brief Find the nearest sphere interfaces to a test point and the closest interface point and normal. 
 *  \param[in]  RelativeCoordinate Test point in relative (box) coordinates.
 *  \param[out] list               Nearest interfaces, each as GeometryVector(distance, side ID, particle ID).
 *  \param[out] x0_Relative        Interface point closest to the test point, in relative coordinates.
 *  \param[out] n2host             Unit outward normal vector at x0. */
void SpherePacking::FindNearestNeighbors (const GeometryVector & RelativeCoordinate, std::vector<GeometryVector> & list, GeometryVector & x0_Relative, GeometryVector & n2host) const{
	//TODO: fix this!
	size_t num = 1;
	for (size_t i=0; i < this->Dimension; i++){
		num *= 2;
	}
	//std::min((size_t)3, this->NumParticle());
	if (this->NumParticle() == 0) {
		//empty packing: nothing to search (avoids /0 and an unbounded search loop)
		list.clear();
		return;
	}
	double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension());
	std::vector<size_t> neighbors;
	std::vector<GeometryVector> shifts;
	// search neighbors
	while (neighbors.size() < num)
	{
		neighbors.clear();
		shifts.clear();
		(*this).IterateThroughNeighbors(RelativeCoordinate, TypicalLength, [&shifts, &neighbors](const GeometryVector & shift, const GeometryVector &, const signed long *, const size_t Sourceparticle) ->void
		{	
			neighbors.push_back(Sourceparticle);
			shifts.push_back(shift);
		});
		TypicalLength *= 1.5;
	}

	list.resize(0);
	//std::vector<GeometryVector> origins, ns;
	double smallest_distance = 0.0;

	//search for the closest points to each neighbor
	for (size_t i = 0; i < neighbors.size(); i++){
		size_t prt_idx = neighbors[i];
		GeometryVector pos_wrt_center = -1.0* shifts[i];	

		std::vector<GeometryVector> list_i;
		GeometryVector x0, n;
		double r = std::sqrt(pos_wrt_center.Modulus2());
		list_i.emplace_back(std::abs(r-this->GetCharacteristics(prt_idx)), 0.0, prt_idx);
		list_i.emplace_back(std::abs(r+this->GetCharacteristics(prt_idx)), 0.0, prt_idx);

		// find the smallest distance to neighboring particles.
		if (i==0)
			smallest_distance = 1.0 + list_i[0].x[0];

		//skip the 1/r direction/interface-point update when the test point
		//coincides with a particle center (r==0 would produce NaN).
		if (r > 0.0)
		{
			if (
				(list_i[0].x[0]<list_i[1].x[0])
			&& (list_i[0].x[0] < smallest_distance)
			){
				smallest_distance = list_i[0].x[0];
				n2host = 1.0/r*pos_wrt_center ;
				x0_Relative = (1.0 - this->GetCharacteristics(prt_idx)/r) * shifts[i];	//In Cartesian coordinates.
			}
			else if(
				(list_i[0].x[0] > list_i[1].x[0])
			&& (list_i[1].x[0] < smallest_distance)
			){
				smallest_distance = list_i[1].x[0];
				n2host = -1.0/r*pos_wrt_center ;
				x0_Relative = (1.0 + this->GetCharacteristics(prt_idx)/r) * shifts[i];
				// In Cartesian coordinates.
			}
		}
		list.insert(list.end(), list_i.begin(), list_i.end());
	}

	// sort the list of distances...
	std::sort(list.begin(), list.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.x[0] < right.x[0];} );
	list.resize(3);

	// up to here, x0_Relative = the Cartesian coordinates of interface point closest to a test point with respect to this test point. 
	x0_Relative = this->CartesianCoord2RelativeCoord(x0_Relative) + RelativeCoordinate;
	this->RelativeCoordToMinimumImage(x0_Relative);
}


// SpherePackingPack member functions.  Mirror of ConfigurationPack (see
// PeriodicCellList.cpp), but keyed on the .SpherePack / .SpherePackIndex files
// and storing SpherePacking objects (radius-bearing) instead of Configurations.

/** \brief Open (or attach to) the .SpherePack / .SpherePackIndex files identified by prefix.
 *  \param prefix File name prefix identifying the pack files. */
void SpherePackingPack::Open(const std::string &prefix) {
	this->IndexName = prefix + ".SpherePackIndex";
	this->PackName = prefix + ".SpherePack";
	std::fstream ifile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	if (ifile.good())
	{
		ifile.read((char*)(&NConfig), sizeof(NConfig));
		if (ifile.eof())
			NConfig = 0;
	}
	else
		NConfig = 0;
}

/** \brief Return the number of sphere packings currently in the pack. */
long long SpherePackingPack::NumConfig(void)
{
	return NConfig;
}

/** \brief Append a sphere packing to the pack files.
 *  \param c Sphere packing to append. */
void SpherePackingPack::AddConfig(const SpherePacking & c)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary);
	if (indexFile.good() == false)
	{
		indexFile.clear();
		indexFile.close();
		indexFile.clear();
		indexFile.open(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	}
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::ate);
	if (packFile.good() == false)
	{
		packFile.clear();
		packFile.close();
		packFile.clear();
		packFile.open(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	}
	this->NConfig++;
	indexFile.seekp(0, std::fstream::beg);
	indexFile.write((char*)(&NConfig), sizeof(NConfig));
	long long loc = packFile.tellp();
	indexFile.seekp(0, std::fstream::end);
	indexFile.write((char*)(&loc), sizeof(loc));
	c.WriteBinary(packFile);
	indexFile.close();
	if (indexFile.fail())
		std::cerr << "Warning in SpherePackingPack::AddConfig : Index file failing!\n";
	packFile.close();
	if (packFile.fail())
		std::cerr << "Warning in SpherePackingPack::AddConfig : Pack file failing!\n";
}

/** \brief Read and return the i-th sphere packing from the pack.
 *  \param i Index of the sphere packing to read.
 *  \return  The i-th sphere packing in the pack. */
SpherePacking SpherePackingPack::GetConfig(long long i)
{
	if (i < 0 || i >= this->NConfig)
	{
		std::cerr << "Error in SpherePackingPack::GetConfig : Index is out of range!\n";
		exit(1);
	}
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::binary);
	long long loc;
	indexFile.seekg(sizeof(long long)*(i + 1));
	indexFile.read((char*)(&loc), sizeof(loc));
	packFile.seekg(std::streampos(loc));
	return packFile.good() ? SpherePacking(packFile) : SpherePacking();
}

/** \brief Discard all sphere packings currently stored in the files. */
void SpherePackingPack::Clear(void)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	this->NConfig = 0;
}

/** \brief Construct an empty sphere packing for a given box geometry. 
 *  \param Dimension    Spatial dimension.
 *  \param BasisVectors Basis vectors defining the periodic box.
 *  \param CellSize     Target linear cell size.
 *  \param UseSortedList Choose sorted vicinity list for faster neighbor checks when geometry is fixed. */
//polydisperse sphere packing
SpherePacking::SpherePacking(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList) : PeriodicCellList<double>(Dimension, BasisVectors, CellSize, UseSortedList)
{
	//std::cerr<<"SpherePacking CellSize="<<CellSize<<'\n';
	MaxRadius=0;
}

/** \brief Insert a sphere of the given radius at a relative coordinate.
 *  \param radius            Radius of the sphere to insert.
 *  \param RelativeCoordinate Position of the sphere center in relative (box) coordinates. */
void SpherePacking::Insert(double radius, const GeometryVector & RelativeCoordinate)
{
	UpdateMaxRadius();
	this->ParticleCartesians.push_back(GeometryVector(this->Dimension));
	this->ParticleRelatives.push_back(GeometryVector(this->Dimension));
	this->ParticleCharacteristics.push_back(radius);
	for(DimensionType i=0; i<this->Dimension; i++)
	{
		assert( RelativeCoordinate.x[i]==RelativeCoordinate.x[i] );

		//NEED TO subtract the floor TWICE
		this->ParticleRelatives.back().x[i]=RelativeCoordinate.x[i]-std::floor(RelativeCoordinate.x[i]);
		this->ParticleRelatives.back().x[i]-=std::floor(this->ParticleRelatives.back().x[i]);
	}
	this->UpdateCartesianCoordinates(this->NumParticle()-1);
	this->Cells[this->GetIndex(this->NumParticle()-1)].push_back(this->NumParticle()-1);

	if(this->MaxRadius<radius)
		this->MaxRadius=radius;
}

/** \brief Check whether the voxel at cr (with given halfsize) is completely covered by existing spheres.  (mainly used in RandomSequentialAddition.{h,cpp})
 *  \param cr       Voxel center in relative (box) coordinates.
 *  \param halfsize Half the voxel side length (in relative units).
 *  \param Radius   Cartesian radius of the upcoming sphere.
 *  \return         True if the voxel is completely covered by existing spheres. */
//relative coord and halfsize
//Radius : cartesian radius of the upcoming sphere
bool SpherePacking::CheckVoxelOverlap(const GeometryVector & cr, double halfsize, double Radius) const
{
	UpdateMaxRadius();
	bool Overlap=false;
	const GeometryVector * ba=this->BasisVector;
	DimensionType dim=this->Dimension;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius, [&Overlap, &halfsize, &ba, dim, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector &, const signed long *, const size_t Sourceparticle) -> void
	{
		//additional processing of shift
		GeometryVector ss(shift);
		for(DimensionType d=0; d<dim; d++)
			if(ss.Dot(ba[d])>0.0)
				ss.AddFrom( halfsize*ba[d]);
			else
				ss.MinusFrom( halfsize*ba[d]);

		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);
		if(ss.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	return Overlap;
}

/** \brief Returns 0 if voxel is completely empty, 2 if voxel is fully occupied, 1 otherwise. (mainly used in RandomSequentialAddition.{h,cpp})
 *  \param cr       Voxel center in relative (box) coordinates.
 *  \param halfsize Half the voxel side length (in relative units).
 *  \param Radius   Cartesian radius of the upcoming sphere.
 *  \return         0 if empty, 2 if fully occupied, 1 if partially occupied. */
//returns 0 if voxel is completely empty
//2 if voxel is fully occupied
//1 otherwise
int SpherePacking::CheckVoxelOverlap2(const GeometryVector & cr, double halfsize, double Radius) const
{
	UpdateMaxRadius();
	bool Overlap=false;
	bool PartialOverlap=false;
	const GeometryVector * ba=this->BasisVector;
	DimensionType dim=this->Dimension;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	double maxVoxelDiagonal=this->MaxDiagonalLength()*halfsize;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius+maxVoxelDiagonal, [&maxVoxelDiagonal, &PartialOverlap, &Overlap, &halfsize, &ba, dim, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector &, const signed long *, const size_t Sourceparticle) -> void
	{
		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);

		double temp=OverlapRadius+maxVoxelDiagonal;
		if(shift.Modulus2()<temp*temp)
			PartialOverlap=true;
		//additional processing of shift
		GeometryVector ss(shift);
		for(DimensionType d=0; d<dim; d++)
			if(ss.Dot(ba[d])>0.0)
				ss.AddFrom( halfsize*ba[d]);
			else
				ss.MinusFrom( halfsize*ba[d]);

		if(ss.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	if(Overlap)
		return 2;
	else if(PartialOverlap)
		return 1;
	else
		return 0;
}


/** \brief Check whether a sphere of the given radius placed at cr overlaps any existing sphere. 
 *  \param cr     Sphere center in relative (box) coordinates.
 *  \param Radius Cartesian radius of the upcoming sphere.
 *  \return       True if the new sphere overlaps any existing sphere. */
//relative coord and halfsize
//Radius : cartesian radius of the upcoming sphere
bool SpherePacking::CheckOverlap(const GeometryVector & cr, double Radius) const
{
	UpdateMaxRadius();
	bool Overlap=false;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius, [&Overlap, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector &, const signed long *, const size_t Sourceparticle) -> void
	{
		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);
		if(shift.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	return Overlap;
}


