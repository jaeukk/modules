/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: April 14th 2019 */

/** \file Dispersion.cpp 
 * \brief Function implementations for the Dispersion classes. */

#include "Dispersions.h"

 /** Give a particle shape.
  * @param i	The particle index.
  * @return A constant shared pointer of the particle shape of particle i. */
 /*	virtual Shape_ & GetCharacteristics(size_t i) const {
 if (i >= this->NumParticle())
 throw std::out_of_range("Dispersion::GetCharacteristics(size_t)");

 Shape_ prt = this->ParticleCharacteristics[i]->clone();
 std::cout << "this->ParticleCharacteristics[i] " << this->ParticleCharacteristics[i].use_count() << "\n";
 std::cout << "prt " << prt.use_count() << "\n";
 std::cout << *prt;
 return prt;
 }*/
 /** Give a particle shape.
 * @param i The particle index.
 * @return A shared pointer of the particle shape of particle i. */
 /*	virtual Shape_ & GetCharacteristics(size_t i) {
 if (i >= this->NumParticle())
 throw std::out_of_range("Dispersion::GetCharacteristics(size_t)");

 Shape_ prt = this->ParticleCharacteristics[i]->clone();
 return prt;
 }*/


/** @return Ratio of all particle volumes to the volume of simulation box. */
double Dispersion::PackingFraction() const {
	double prtVol = 0.0;
	for (int i = 0; i < this->NumParticle(); i++) {
		prtVol += this->ParticleCharacteristics[i]->GetVolume();
	}
	return prtVol / this->PeriodicVolume();
}

//-----rescale system box and particle sizes....

/** Resize all particles and the simulation box. */
void Dispersion::Rescale(double LinearFactor)
{	//Resize the simulation box.
	PeriodicCellList::Rescale(LinearFactor);
	//Resize Particles.
	this->MaxRadius *= LinearFactor;
	for (int i = 0; i < this->NumParticle(); i++) {
		this->ParticleCharacteristics[i]->Scale(LinearFactor);
	}
}

/** Resize all particles and the simulation box. */
void Dispersion::Resize(double NewVolume)
{
	double factor = std::pow(NewVolume / this->PeriodicVolume(), 1.0 / this->Dimension);
	this->Rescale(factor);
}

//Rescale particle sizes only
void Dispersion::Rescale_prtOnly(double LinearFactor) {
	for (int i = 0; i < this->NumParticle(); i++) {
		this->ParticleCharacteristics[i]->Scale(LinearFactor);
	}
	this->MaxRadius *= LinearFactor;
}

void Dispersion::UpdateParticleShape(const Shape_ & newShape, size_t index) {
	this->ParticleCharacteristics[index] = newShape->clone();
	this->MaxRadius = std::max(this->MaxRadius, newShape->GetMaxRadius());
}

void WriteDispersion(const Dispersion & conf, const std::string & OutputName){

	std::string name = OutputName + std::string(".Dispersion");
	std::fstream output(name.c_str(), std::fstream::out);
	
	output << std::setprecision(10);
	output << "##Dispersion##\n";
	output << "#Particle\t"<<conf.NumParticle()<<"\n";
	output << conf.GetDimension() << "\n";
	for (int i = 0; i < conf.GetDimension(); i++) {
		output <<conf.GetBasisVector(i);
	}

	for (size_t i = 0; i < conf.NumParticle(); i++) {
		Shape_ ptr_ = conf.GetShape(i);
		ptr_->ShiftBy(conf.GetCartesianCoordinates(i));
		output << * ptr_;
	}
	output.close();
}

void ReadDispersion(Dispersion & conf, const std::string & InputName){
	std::string name = InputName + std::string(".Dispersion");
	std::fstream ifile (name.c_str(), std::fstream::in);

	if (ifile.good() == false)
	{
		std::cerr << "Error in ReadFunction : input stream is not good!\n" << InputName << "\n";
		return;
	}

	std::string dump;
	std::getline(ifile, dump);
	if (strcmp(dump.c_str(), "##Dispersion##") == 0){
		// set the number of particles
		int num_particles = 0;
		ifile >> dump;
		if(strcmp(dump.c_str(), "#Particle")==0){
			ifile >> num_particles;
		}

		DimensionType d = 0;
		GeometryVector * basis;
		ifile >> d;
		basis = new GeometryVector[d];

		// set dimension
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				ifile >> basis[i].x[j];
			}
			basis[i].SetDimension(d);
		}

		conf = Dispersion(d, basis, 1.0);
		for(int i=0; i<num_particles; i++)
		{
			char dump[255] = {};

			ifile >> dump;
			if(strcmp(dump, "[Polygon]")==0){
				conf.Insert(std::shared_ptr<Shape> ( new Polygon(ifile)), GeometryVector(0,0));
			}
			else if (strcmp(dump, "[Polyhedron]")==0){
				conf.Insert(std::shared_ptr<Shape> ( new Polyhedron(ifile)), GeometryVector(0,0,0));
			}
			else{
				std::cerr << "Undefined derived class of Shape\n";
				return;
			}

		}

	}
	else{
		std::cerr << "Error in ReadDispersion : File format is inappropriate!\n" << InputName << "\n";
		return;		
	}
}

int Dispersion::IsInsideParticle(const GeometryVector & RelativeCoordinate) const
{
	int idx = -1;
	this->IterateThroughNeighbors(RelativeCoordinate, this->GetMaxRadius(), 
	[&idx, this] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void{
		//GeometryVector pos_Cartesian = this->RelativeCoord2CartesianCoord(shift);
		if(this->ParticleCharacteristics[Sourceparticle]->IsPtInside(-1.0*shift)){
			idx = Sourceparticle;
		}
	});
	return idx;
}

void Dispersion::FindNearestNeighbors (const GeometryVector & RelativeCoordinate, std::vector<GeometryVector> & list, GeometryVector & x0_Relative, GeometryVector & n2host) const{
	//TODO: fix this!
	size_t num = 1;
	for (size_t i=0; i < this->Dimension; i++){
		num *= 2;
	}
	//std::min((size_t)3, this->NumParticle());
	double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension());
	std::vector<size_t> neighbors;
	std::vector<GeometryVector> shifts;
	// search neighbors
	while (neighbors.size() < num)
	{
		neighbors.clear();
		shifts.clear();
		(*this).IterateThroughNeighbors(RelativeCoordinate, TypicalLength, [&shifts, &neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
		{	
//			bool repeated = false;
//			for (size_t i=0; i<neighbors.size() && !repeated; i++){
//				repeated = (neighbors[i] == Sourceparticle);
//			}
//			if(!repeated){
				neighbors.push_back(Sourceparticle);
				shifts.push_back(shift);
//			}
		});
		TypicalLength *= 1.5;
	}

	list.resize(0);
	//std::vector<GeometryVector> origins, ns;
	double smallest_distance = 0.0;

	for (size_t i = 0; i < neighbors.size(); i++){
		size_t prt_idx = neighbors[i];
//		GeometryVector pos_wrt_center = RelativeCoordinate - this->ParticleRelatives[prt_idx];
//		this->RelativeCoordToMinimumImage(pos_wrt_center); 
//		pos_wrt_center = this->RelativeCoord2CartesianCoord(pos_wrt_center);
		GeometryVector pos_wrt_center = -1.0* shifts[i];		

		std::vector<GeometryVector> list_i;
		GeometryVector x0, n;
		double distance_i = this->ParticleCharacteristics[prt_idx]->GetDistance2Boundary(pos_wrt_center, list_i, x0, n);

		for(auto v = list_i.begin(); v != list_i.end(); v++){
			v->SetDimension(3);
			v->x[2] = (double)prt_idx;
		}
		list.insert(list.end(), list_i.begin(), list_i.end());

		//origins.push_back(		this->CartesianCoord2RelativeCoord(x0 + this->ParticleCartesians[prt_idx]));
		//ns.push_back(n);
		if( i == 0 || smallest_distance > list_i[0].x[0]){
			smallest_distance = list_i[0].x[0];
			x0_Relative = x0 + shifts[i];// now in the Cartesian coordinates.
			n2host = n;
		}
	}

	std::sort(list.begin(), list.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.x[0] < right.x[0];} );
	list.resize(3);

	//size_t idx = std::find(neighbors.begin(), neighbors.end(), (size_t)round(list[0].x[2])) - neighbors.begin();
	//x0_Relative = origins[idx];
	//this->RelativeCoordToMinimumImage(x0_Relative);
	//n2host = ns[idx];

	// up to here, x0_Relative = the Cartesian coordinates of interface point closest to a test point with respect to this test point. 
	x0_Relative = this->CartesianCoord2RelativeCoord(x0_Relative) + RelativeCoordinate;
	this->RelativeCoordToMinimumImage(x0_Relative);
}


//Dispersion_cube::Dispersion_cube(DimensionType Dimension, GeometryVector * BasisVector, double cellSize, bool UseSortedList)
//	: PeriodicCellList<Cube>(Dimension, BasisVector, cellSize, UseSortedList)
//{
//	this->MaxRadius = 0.0;
//};
//
//Cube & Dispersion_cube::GetCharacteristics(size_t i) {
//	return this->ParticleCharacteristics[i];
//	//return Cube_(new Cube(this->ParticleCharacteristics[i]));
//}
//
//const Cube & Dispersion_cube::GetCharacteristics(size_t i) const {
//	return this->ParticleCharacteristics[i];
//	//return Cube_(new Cube(this->ParticleCharacteristics[i]));
//}
//
//void Dispersion_cube::UpdateMaxRadius(void) const {
//	if (this->MaxRadius == 0.0)
//	{
//		for (int i = 0; i < this->NumParticle(); i++)
//			this->MaxRadius = std::max(this->MaxRadius, this->ParticleCharacteristics[i].GetMaxRadius());
//	}
//}
//
//double Dispersion_cube::GetMaxRadius(void) const {
//	this->UpdateMaxRadius();
//	return this->MaxRadius;
//}
//
//void Dispersion_cube::Insert(const Cube & prtShape, const GeometryVector & RelativeCoordinate) {
//	UpdateMaxRadius();
//	this->ParticleCartesians.push_back(GeometryVector(this->Dimension));
//	this->ParticleRelatives.push_back(GeometryVector(this->Dimension));
//	this->ParticleCharacteristics.push_back(Cube(prtShape));
//	for (DimensionType i = 0; i<this->Dimension; i++)
//	{
//		assert(RelativeCoordinate.x[i] == RelativeCoordinate.x[i]);
//
//		//NEED TO subtract the floor TWICE
//		this->ParticleRelatives.back().x[i] = RelativeCoordinate.x[i] - std::floor(RelativeCoordinate.x[i]);
//		this->ParticleRelatives.back().x[i] -= std::floor(this->ParticleRelatives.back().x[i]);
//	}
//	this->UpdateCartesianCoordinates(this->NumParticle() - 1);
//	this->Cells[this->GetIndex(this->NumParticle() - 1)].push_back(this->NumParticle() - 1);
//
//	if (this->MaxRadius<prtShape.GetMaxRadius())
//		this->MaxRadius = prtShape.GetMaxRadius();
//}
//
//
//void Dispersion_cube::RotateAllParticles() {
//	RandomGenerator rng;
//	for (size_t i = 0; i < this->NumParticle(); i++) {
//		this->RotateParticle(i, rng);
//	}
//}
//
//void Dispersion_cube::RotateParticle(size_t i, RandomGenerator & rng) {
//	if (i < this->NumParticle()) {
//		this->ParticleCharacteristics[i].Rotate(rng);
//	}
//	else
//		std::cerr << "Dispersion::RotateParticle(.), index is out of range" << std::endl;
//}
//
//void Dispersion_cube::RotateParticle(size_t i, const std::vector<GeometryVector> & newbasis) {
//	if (i < this->NumParticle()) {
//		this->ParticleCharacteristics[i].Rotate(newbasis);
//	}
//	else
//		std::cerr << "Dispersion::RotateParticle(.), index is out of range" << std::endl;
//}
//
//double Dispersion_cube::PackingFraction() const {
//	double prtVol = 0.0;
//	for (int i = 0; i < this->NumParticle(); i++) {
//		prtVol += this->ParticleCharacteristics[i].GetVolume();
//	}
//	return prtVol / this->PeriodicVolume();
//}
//
////-----rescale system box and particle sizes....
//void Dispersion_cube::Rescale(double LinearFactor)
//{
//	this->ApproximateCellSize *= LinearFactor;
//	this->MaxRadius *= LinearFactor;
//	std::vector<GeometryVector> basis;
//	for (DimensionType i = 0; i<this->Dimension; i++)
//		basis.push_back(this->BasisVector[i] * LinearFactor);
//	this->ChangeBasisVector(&basis[0]);
//
//	for (int i = 0; i < this->NumParticle(); i++) {
//		this->ParticleCharacteristics[i].Scale(LinearFactor);
//	}
//}
//void Dispersion_cube::Resize(double NewVolume)
//{
//	double factor = std::pow(NewVolume / this->PeriodicVolume(), 1.0 / this->Dimension);
//	this->Rescale(factor);
//}
////Rescale particle sizes only
//void Dispersion_cube::Rescale_prtOnly(double LinearFactor) {
//	for (int i = 0; i < this->NumParticle(); i++) {
//		this->ParticleCharacteristics[i].Scale(LinearFactor);
//	}
//	this->MaxRadius *= LinearFactor;
//}
//
//void Dispersion_cube::UpdateCharacteristics(const Cube & newShape, size_t index) {
//	this->ParticleCharacteristics[index] = Cube(newShape);
//	this->MaxRadius = std::max(this->MaxRadius, newShape.GetMaxRadius());
//}

