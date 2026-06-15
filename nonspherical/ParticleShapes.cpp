/**
*	Author	: Jaeuk Kim
*	Email	: phy000.kim@gmail.com
*	Date	:	Feb. 2019 */

/** \file ParticleShapes.cpp
 * \brief Function implementations for "Shape" and its children structs */

#include "ParticleShapes.h"
#include <algorithm>



void Shape::Rotate(RandomGenerator & rng) {
	if (this->Dimension == 2 || this->Dimension == 3) {
		std::vector<GeometryVector> newcoordbasis;
		if (this->Dimension == 2) {
			double phi = 2 * pi*rng.RandomDouble();
			newcoordbasis.emplace_back(GeometryVector(cos(phi), sin(phi)));
			newcoordbasis.emplace_back(GeometryVector(-sin(phi), cos(phi)));
		}
		else {
			//See X1 Z2 X3 rotation in https://en.wikipedia.org/wiki/Euler_angles
			double alpha = 2.0*pi*rng.RandomDouble(), gamma = 2.0*pi*rng.RandomDouble();
			double z = 2.0*(rng.RandomDouble() - 0.5);
			double c1 = cos(alpha), c2 = z,						c3 = cos(gamma);
			double s1 = sin(alpha), s2 = sqrt(fabs(1 - z * z)), s3 = sin(gamma);

			newcoordbasis.emplace_back(GeometryVector(c2, -c3 * s2, s2*s3));
			newcoordbasis.emplace_back(GeometryVector(c1*s2, c1*c2*c3 - s1 * s3, -c3 * s1 - c1 * c2*s3));
			newcoordbasis.emplace_back(GeometryVector(s1*s2, c1*s3 + c2 * c3 * s1, c1*c3 - c2 * s1*s3));
			/*
			//Gram-Schmidt process
			newcoordbasis[0] = 1.0 / sqrt(newcoordbasis[0].Modulus2()) * newcoordbasis[0];
			newcoordbasis[1] = newcoordbasis[1] - newcoordbasis[1].Dot(newcoordbasis[0]) * newcoordbasis[0];
			newcoordbasis[1] = 1.0 / sqrt(newcoordbasis[1].Modulus2()) * newcoordbasis[1];
			newcoordbasis[2] = newcoordbasis[2] - newcoordbasis[2].Dot(newcoordbasis[0]) * newcoordbasis[0] - newcoordbasis[2].Dot(newcoordbasis[1])*newcoordbasis[1];
			newcoordbasis[2] = 1.0 / sqrt(newcoordbasis[2].Modulus2())*newcoordbasis[2];
			*/
		}

		this->Rotate(newcoordbasis);
	}
	else {
		std::cerr << "Rotate(.) is not defined for d>3 cases\n";
	}
};

void Shape::Scale(double linearFactor) {
	this->R_circum *= linearFactor;
	this->R_inscr *= linearFactor;
	double multiplier = linearFactor;
	for (int i = 2; i < this->Dimension; i++)	
		multiplier *= linearFactor;
	
	this->area *= multiplier;
	this->volume *= multiplier * linearFactor;
}

void Shape::Scale2Vol(double newVolume) {
	double factor = newVolume / this->GetVolume();
	switch (this->GetDimension()) {
	case 1:
		break;
	case 2:
		factor = std::sqrt(factor);	break;
	case 3:
		factor = std::cbrt(factor); break;
	default:
		factor = std::pow(factor, 1.0 / this->GetDimension());	break;
	}
	this->Scale(factor);
}

bool CheckNewCoordinates(int Dimension, const std::vector<GeometryVector> & newbasis) {
	const double error = 1e-13;	// 
	/* new basis vector should be orthonormal */
	if (Dimension == newbasis.size()) {
		bool check = true;
		if (Dimension == 2) {
			check &= fabs(newbasis[0].Modulus2() - 1.0) < error;
			check = check && ((fabs(newbasis[0].x[0] - newbasis[1].x[1])< error) && (fabs(newbasis[0].x[1] + newbasis[1].x[0])< error));
			//return check;
		}
		else if (Dimension == 3) {
			/* check normal */
			for (int i = 0; i < Dimension; i++) {
				check &= fabs(newbasis[i].Modulus2() - 1.0) < error;
			//	std::cout << newbasis[i];
			//	std::cout << newbasis[i].Modulus2() - 1.0 << std::endl;
			}
				
			
			if (check) {
				/* check orthogonality */
				check &= fabs(newbasis[0].Dot(newbasis[1])) < error;
				check &= fabs(newbasis[1].Dot(newbasis[2])) < error;
				check &= fabs(newbasis[2].Dot(newbasis[0])) < error;
				if (check) {
					/* check orientation */
					GeometryVector crossXY(newbasis[0].x[1] * newbasis[1].x[2] - newbasis[0].x[2] * newbasis[1].x[1], newbasis[0].x[2] * newbasis[1].x[0] - newbasis[0].x[0] * newbasis[1].x[2], newbasis[0].x[0] * newbasis[1].x[1] - newbasis[0].x[1] * newbasis[1].x[0]);
					check &= (crossXY - newbasis[2]).Modulus2() < error;
				}
				//else {
				//	return check;
				//}
			}
			//else
			//	return check;
		}
		return check;
	}
	else {
		std::cerr << "Dimensions are inconsistent!\n";
		return false;
	}
}

/* -------------------------------------------------
 * PlanarFigure protected member functions
 * -------------------------------------------------*/

// Check	(1) all vertices are on the same plane and
//			(2) vertices are aligned in the same orientations
void PlanarFigure::CheckVertices() {
	GeometryVector temp((int)3);
	this->Trust_V = true;
	size_t NUM = this->num_Vs;
	//compute the normal vector
	if (this->num_Vs > 2) {
		if (this->Dimension == 2) {
			double sum = 0;
			for (size_t i=0; i<this->num_Vs; i++){
				size_t i_1 = (i == 0)? this->num_Vs - 1 : i - 1;

				sum += (this->V[i_1]->x[0] * this->V[i]->x[1]) - (this->V[i]->x[0] * this->V[i_1]->x[1]);
			}
			this->Trust_V = (std::abs(sum) >1e-20);
			this->normal = GeometryVector(0.0, 0.0, (sum > 0.0)?1.0:-1.0);

			// bool firstNonzeroProduct = true;
			// for (size_t i = 0; i < NUM; i++) {
			// 	size_t j = (i+1 == NUM) ? 0 : i + 1;	//next vertex
			// 	/*double DotProduct = V[i]->Dot(*V[j]);
			// 	if (DotProduct != 0)	//when the inner product of two neighboring vertices is nonzero.
			// 	*/
			// 	double cross = (V[i]->x[0])*(V[j]->x[1]) - (V[j]->x[0])*(V[i]->x[1]);
			// 	if(cross!=0.0)
			// 	{   /*
			// 		temp.x[2] = (V[i]->x[0])*(V[j]->x[1]) - (V[j]->x[0])*(V[i]->x[1]);
			// 		temp = (1.0 / sqrt(temp.Modulus2()))*temp;
			// 		*/
			// 		temp.x[2] = cross / abs(cross);
			// 		if (firstNonzeroProduct) {
			// 			this->normal = temp;	firstNonzeroProduct = false;
			// 		}
			// 		else {
			// 			this->Trust_V &= this->normal.IsEqual(temp);
			// 		}

			// 	}
			// }
			this->r_perp = 0.0;
		}
		else if (this->Dimension == 3) {
			/*temp = (*V[NUM - 1]).Cross(*V[0]);
			this->normal = (1.0 / sqrt(temp.Modulus2()))*temp;

			for (int i = 1; i < NUM; i++) {
				temp = (*V[i - 1]).Cross(*V[i]);
				temp = (1.0/ sqrt(temp.Modulus2()))*temp;
				this->Trust_V &= this->normal.IsEqual(temp);
			}*/
			temp = (this->E[NUM - 1]).Cross(this->E[0]);

			this->normal = (1.0 / sqrt(temp.Modulus2()))*temp;
 
			for (int i = 1; i < NUM; i++) {
				temp = (this->E[i - 1]).Cross(this->E[i]);
				temp = (1.0 / sqrt(temp.Modulus2()))*temp;
				this->Trust_V &= this->normal.IsEqual(temp);
			}
			this->r_perp = this->normal.Dot(*V[0]);
		}
		else {
			std::cerr << "PlanarFigure is not defined in " << this->Dimension << "dimensions!"<<std::endl;
		}
	}
	else
		std::cerr << "PlanarFigure cannot be defined with less than 3 vertices!"<<std::endl;

	//if (!this->Trust_V)
	//	std::cerr << "PlanarFigure:: orientations of vertices are inconsistent." << std::endl;
}



// Update two types of vectors: E_i = (V_i - V_{i-1})/2 and R_i = (V_i+V_{i-1})/2
void PlanarFigure::Update_ER_vectors() {
	this->E.resize(this->num_Vs);
	this->R.resize(this->num_Vs);
	size_t j = 0;

	this->E[0] = 0.5*(*(this->V[0]) - *(this->V[this->num_Vs - 1]));
	this->R[0] = 0.5*(*(this->V[0]) + *(this->V[this->num_Vs - 1]));

	for (int i = 1; i < this->num_Vs; i++) {
//		j = (i == 0) ? this->num_Vs - 1 : i - 1;
		this->E[i] = 0.5*(*(this->V[i]) - *(this->V[i-1]));
		this->R[i] = 0.5*(*(this->V[i]) + *(this->V[i-1]));
	}
	this->update_ER = true;
}

/**********************************
 * Inherited member functions
 **********************************/
// Update the largest distnace from vertices to the center.
void PlanarFigure::UpdateR() {
	double temp_max = 0, temp_min=0, tempV2 = 0;
	for (int i = 0; i < this->num_Vs; i++) {
		tempV2 = (this->V[i])->Modulus2();
		temp_max = std::max(tempV2, temp_max);		
	}
	this->R_circum = sqrt(temp_max);
	this->R_inscr = 0.0; // sqrt(temp_min);
}

void PlanarFigure::UpdateVolume() {
	GeometryVector sum((int)3);

	if(this->Dimension == 2){
		//Shoelace formula
		for (size_t i=0; i<this->num_Vs; i++){
			size_t i_1 = (i == 0)? this->num_Vs - 1 : i - 1;

			sum.x[2] += ((this->V[i_1]->x[0] * this->V[i]->x[1]) - (this->V[i]->x[0] * this->V[i_1]->x[1]));
		}
		this->volume = 0.5*(this->normal.Dot(sum));	
	}
	else if(this->Dimension == 3){
		for (size_t i = 0; i < this->num_Vs; i++) {
		size_t j = (i == 0) ? this->num_Vs - 1 : i - 1;
		
			sum = sum + (this->V[j])->Cross(*(this->V[i]));

		}
		this->volume = 0.5*(this->normal.Dot(sum));	
	}
	// for (size_t i = 0; i < this->num_Vs; i++) {
	// 	size_t j = (i == 0) ? this->num_Vs - 1 : i - 1;
	// 	if (this->Dimension == 2)
	// 		sum.x[2] += (this->V[j]->x[0] * this->V[i]->x[1]) - (this->V[i]->x[0] * this->V[j]->x[1]);

	// 	else if (this->Dimension == 3)
	// 		sum = sum + (this->V[j])->Cross(*(this->V[i]));

	// }
	// this->volume = 0.5*(this->normal.Dot(sum));
}

void PlanarFigure::UpdateSurfaceArea() {
	double sum = 0.0;
	if (!this->update_ER)
		Update_ER_vectors();

	for (int i = 0; i < this->num_Vs; i++) {
		sum += sqrt(this->E[i].Modulus2());
	}
	this->area = 2.0 * sum;
}


/* -------------------------------------------------
 * PlanarFigure constructors
 * -------------------------------------------------*/

PlanarFigure::PlanarFigure(GeometryVector * v, size_t num_vertices) : Shape::Shape() {
	this->num_Vs = num_vertices;
	this->Dimension = v[0].Dimension;	
	for (int i = 0; i < this->num_Vs; i++)
		this->V.push_back(&v[i]);
	
	this->Update_ER_vectors();
	this->CheckVertices();
}

PlanarFigure::PlanarFigure(GeometryVector * v, std::vector<size_t> vert_idx) : Shape::Shape(){
	this->num_Vs = vert_idx.size();
	this->Dimension = v[0].Dimension;
	for (int i = 0; i < this->num_Vs; i++)
		this->V.push_back(&v[vert_idx[i]]);

	this->Update_ER_vectors();
	this->CheckVertices();
}

/* -------------------------------------------------
* PlanarFigure member functions
* -------------------------------------------------*/
GeometryVector PlanarFigure::GetCenterOfMass() {
	GeometryVector x0;
	if (this->Dimension == 2) {			x0 = GeometryVector(0.0, 0.0);}
	else if (this->Dimension == 3) {	x0 = this->r_perp * this->normal;}
	else {								std::cerr << "PlanarFigure:: wrong dimension\n";}

	GeometryVector center((int)this->Dimension);
	//Divide a planar figure into triangles
	//Compute the center of mass of each triangle, and compute the center of mass of the entire shape.
	for (int i = 0; i < this->num_Vs; i++) {
		size_t j = (i == 0) ? this->num_Vs - 1 : i - 1;
		GeometryVector v = *this->V[i] - x0, u = *this->V[j] - x0;
		GeometryVector C = (1.0 / 3.0 *(v+u))+x0;
		double area;
		{
			GeometryVector sum((int)3);		
			if (this->Dimension == 2)
				sum.x[2] += (u.x[0] * v.x[1]) - (v.x[0] * u.x[1]);

			else if (this->Dimension == 3)
				sum = sum + u.Cross(v);
			area = 0.5*this->normal.Dot(sum);
		}
		center = center + area * C;
	}
	center = 1.0 / this->GetVolume() * center;
	Centroid = center;
	return center;
}
void PlanarFigure::Print(std::ostream & ofile) const {
	ofile << "#Vertices\t" << this->num_Vs << "\n";
	for (int i = 0; i < this->num_Vs; i++)
		ofile << *(this->V[i]);
}

void PlanarFigure::GetVertices(std::vector<GeometryVector> & Verts) const {
	Verts.resize(this->num_Vs);
	for (int i = 0; i < this->num_Vs; i++) {
		Verts[i] = *(this->V[i]);
	}
}

void PlanarFigure::Rotate(const std::vector<GeometryVector> & newcoordinatebasis) {
	//if (CheckNewCoordinates(this->dimension, newcoordinatebasis)) {
		//vertices
		//for (int i = 0; i < this->num_Vs; i++) {
		//	GeometryVector temp((int)this->Dimension);
		//	for (int j = 0; j < this->Dimension; j++) {
		//		temp.x[j] += newcoordinatebasis[j].Dot(* this->V[i]);
		//	}
		//	*(this->V[i]) = temp;
		//}
		this->update_ER = false;
		// normal vector
		if(this->Dimension==3){
			GeometryVector temp((int)this->Dimension);
			for (int j = 0; j < this->Dimension; j++) {
				temp.x[j] += newcoordinatebasis[j].Dot(this->normal);
			}
			this->normal = temp;
		}

	//}
	//else {
	//	std::cerr << "input basis vectors are incorrect!\n";
	//}
};

void PlanarFigure::Scale(double linearfactor) {
	Shape::Scale(linearfactor);
//	this->R_circum *= linearfactor;
//	this->volume *= linearfactor * linearfactor;
//	this->area *= linearfactor;
	this->r_perp *= linearfactor;
	for (int i = 0; i < this->num_Vs; i++) {
//		*(this->V[i]) = linearfactor * (*(this->V[i]));
		if (this->update_ER) {
			this->E[i] = linearfactor * this->E[i];
			this->R[i] = linearfactor * this->R[i];
		}
	}
}


/* -------------------------------------------------
 * Polygon constructors 
 * -------------------------------------------------*/
Polygon::Polygon(const Polygon & source) : Polygon::Polygon(&(source.Vertices[0]), source.Vertices.size()){
	this->R_circum = source.R_circum;
	this->R_inscr = source.R_inscr;
	this->volume = source.volume;
	this->area = source.area;
}

Polygon::Polygon(const std::vector<GeometryVector> & vertices) : Polygon::Polygon(&vertices[0], vertices.size()) {
}

Polygon::Polygon(const GeometryVector * vertices, size_t num_vertices) : Shape::Shape() {
	init(vertices, num_vertices);
}

void Polygon::init(const GeometryVector * vertices, size_t num_vertices){
	this->name = "Polygon";
	bool TwoDimensions = true;
	//confirm dimensionality of vertices
	this->Vertices.resize(num_vertices);
	for (size_t i = 0; i < num_vertices; i++) {
		TwoDimensions &= (vertices[i].Dimension == 2);
		GeometryVector vec = vertices[i];
		this->Vertices[i] = vec;
	}

	if (TwoDimensions) {
		this->Dimension = 2;
		// shift vertices by the first vertex to reduce the round-off errors.
		GeometryVector temp_shift(this->Vertices[0]);
		for(size_t i=0; i<this->Vertices.size(); i++){
			this->Vertices[i] = this->Vertices[i] - temp_shift;
		}
		this->Gamma = PlanarFigure(&this->Vertices[0], this->Vertices.size());

		this->UpdateVolume();
		this->UpdateSurfaceArea();
		//Translate the center of mass to the origin
		GeometryVector CM = this->GetCenterOfMass(); 
		this->ShiftBy(-1.0*CM);
		this->Centroid = CM + temp_shift;
		this->UpdateR();

		this->Gamma.Update_ER_vectors();
		this->Gamma.CheckVertices();
		if (!this->Gamma.Trust_V)
			std::cerr << "Polygon::Polygon is ill-defined\n";
	}
	else
		std::cerr << "Vertices in a Polygon should be in two dimensions\n";
}

Polygon::Polygon(std::istream & ifile) : Shape::Shape(){

	std::string phrase("#Vertices");
	size_t num_vertices = 0;
	std::vector<GeometryVector> vertices;

	std::string buffer; std::istringstream iss;

	while (num_vertices == 0 && std::getline(ifile, buffer)){
		if(buffer[0] == '#'){
			num_vertices = (size_t) std::stoi(buffer.substr(phrase.size()));	
		}
	}

	if (num_vertices > 0){
		for (size_t i = 0; i < num_vertices; i++){
			std::getline(ifile, buffer);
			GeometryVector temp((int)2);
			iss.str(buffer);
			iss.clear();

			for (int j=0; j<2; j++)
				iss >> temp.x[j];
			
			if (ifile.eof() == false)
				vertices.push_back(temp);
		}
	}
	else{
		std::cerr << "Polygon::Polygon(ifile): Wrong input\n";
		return ;
	}

	/** The following is the same as the implementation of 
	 * Polygon::Polygon(const GeometryVector * vertices, size_t num_vertices).*/
	init(&vertices[0], num_vertices);
}

void Polygon::SetPositiveOrientation() 
{
	if(this->Gamma.normal.x[2] < 0){
		//Change the ordering of vertices
		std::reverse(this->Vertices.begin(), this->Vertices.end());

		this->Gamma.normal = (-1.0)*this->Gamma.normal;
		this->Gamma.Update_ER_vectors();
	}
}

/* -------------------------------------------------
* Polygon member functions
* -------------------------------------------------*/
void Polygon::UpdateR() {
	this->Gamma.UpdateR();	
	this->R_circum = this->Gamma.GetMaxRadius();
	GeometryVector normal; normal.SetDimension(2);

	double r = 0;
	for (size_t i = 0; i < this->Vertices.size(); i++) {
		//normal vector of the side connecting V[i] and V[i+1]
		normal.x[0] = this->Gamma.E[i].x[1];
		normal.x[1] = -1.0* this->Gamma.E[i].x[0];
		normal = 1.0 / sqrt(normal.Modulus2())*normal;

		if (i == 0)
			r = std::abs(normal.Dot(this->Vertices[0]));
		else
			r = std::min(std::abs(normal.Dot(this->Vertices[i])), r);
	}
	this->R_inscr = r;
}
void Polygon::UpdateVolume() {
	this->Gamma.UpdateVolume(); this->volume = this->Gamma.GetVolume();
}
void Polygon::UpdateSurfaceArea() {
	this->Gamma.UpdateSurfaceArea(); this->area = this->Gamma.GetSurfaceArea();
}
void Polygon::ShiftBy(const GeometryVector & x) {
	//translate vertices
	if(this->Centroid.Dimension == 0){
		//centroid is undefined yet.
		for (auto v = this->Vertices.begin(); v != this->Vertices.end(); v++)
			*v = *v + x;
	
		//translate R vectors (E vectors are invariant under translations)
		this->Gamma.Update_ER_vectors();
		//this->Gamma.CheckVertices();
	}
	else	//When the centroid is defined.
		this->Centroid = this->Centroid + x;
}
GeometryVector Polygon::GetCenterOfMass() {
	if(this->Centroid.Dimension == 0)
		return this->Gamma.GetCenterOfMass();
	else
		return this->Centroid;
}

void Polygon::GetVertices(std::vector<GeometryVector> & vertices) const {
	vertices.resize(0); vertices = std::vector<GeometryVector>(this->Vertices);
	for (auto v = vertices.begin(); v!=vertices.end(); v++)
		*v = *v + this->Centroid;
}
void Polygon::Rotate(const std::vector<GeometryVector> & newbasis) {
	//rotate vertices
	for (auto v = this->Vertices.begin(); v != this->Vertices.end(); v++) {
		GeometryVector temp((int)this->Dimension);
		for (int j = 0; j < this->Dimension; j++) {
			temp.x[j] += newbasis[j].Dot(*v);
		}
		*v = temp;
	}
	//rotate normal vectors
	this->Gamma.Rotate(newbasis);
}
//void Polygon::Rotate(RandomGenerator & rng) : Shape::Rotate(rng) {};
void Polygon::Scale(double factor) {
	Shape::Scale(factor);
	this->Gamma.Scale(factor);
	//this->volume *= factor * factor;
	//this->area *= factor;
	//this->R_circum *= factor;
	//scale vertices
	for (auto v = this->Vertices.begin(); v != this->Vertices.end(); v++)
		*v = factor * (*v);
}

Polygon ReadPolygon(std::istream & ifile) 
{
	std::string phrase("#Vertices");
	size_t num_vertices = 0;
	std::vector<GeometryVector> vertices;

	std::string buffer; std::istringstream iss;

	while (num_vertices == 0 && std::getline(ifile, buffer)){
		if(buffer[0] == '#'){
			num_vertices = (size_t) std::stoi(buffer.substr(phrase.size()));	
		}
	}

	if (num_vertices > 0){
		for (size_t i = 0; i < num_vertices; i++){
			std::getline(ifile, buffer);
			GeometryVector temp((int)2);
			iss.str(buffer);
			iss.clear();

			for (int j=0; j<2; j++)
				iss >> temp.x[j];
			
			if (ifile.eof() == false)
				vertices.push_back(temp);
		}
		return Polygon(vertices);
	}
	else{
		std::cerr << "Polygon::ReadPolygon: Wrong input\n";
		return Polygon();
	}


}

/** Check whether a given point is inside this Polygon object. 
 * @return true = pos is inside this Polygon. */
bool Polygon::IsPtInside(const GeometryVector & pos) {
	//Make sure the position vector is in 2D.
	assert(pos.Dimension == 2);	

	/* The winding number algorithm is implemented.
		* winding number (wc) of a point = how many times the vertices wrap around the point.
		* if wc = 0, the point is outside of the polygon.	 */

	int wc = 0;
	/* Consider a horizontal line segment from pos to infinity. */

	/* Check whether sides are crossing the horizontal line above. */

	double x = pos.x[0];
	double y = pos.x[1];
	for (size_t i = 0; i < this->Vertices.size(); i++) {

		size_t i_1 = (i==0)? this->Vertices.size()-1: i-1;
		bool is_Ei_crossing = false;
		if ( (this->Vertices[i].x[1] - y) * (this->Vertices[i_1].x[1] - y ) < 0.0){
			is_Ei_crossing = (x < std::min(this->Vertices[i].x[0], this->Vertices[i_1].x[0]));
			
			if( (!is_Ei_crossing) && 
				(x < std::max(this->Vertices[i].x[0], this->Vertices[i_1].x[0])) ) {
				//check the intersection of the horizontal line and the line segment of V_i V_{i-1}
				double t = (y-this->Vertices[i].x[1]) / (this->Vertices[i_1].x[1] - this->Vertices[i].x[1]);

				double x_intersect = this->Vertices[i].x[0] + t*(this->Vertices[i_1].x[0]-this->Vertices[i].x[0]);
				is_Ei_crossing = (x < x_intersect);
			}
			else if((this->Vertices[i].x[1] == y)  && (y == this->Vertices[i_1].x[1])){
				is_Ei_crossing = x < std::min(this->Vertices[i].x[0], this->Vertices[i_1].x[0]);
			}
		}

		if (is_Ei_crossing){
			if (this->Gamma.E[i].x[1] >= 0.0)
				wc++;
			else
			{
				wc--;
			}
		}
	} 


	if (wc == 0){
		return false;
	}
	else{
		return true;
	}
}

bool compare_temp(GeometryVector a, GeometryVector b) {
	assert(a.Dimension == b.Dimension);
	return a.x[0] < b.x[0];
}

/**	Compute the closest distance from a given point to the side of this polygon.
 * @param pos	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
 * @param[out] list_dist	A list of distances in the ascending order.
 * @param[out] n_outward	The unit outward normal vector at x0. 
 * @return distance.	 */
double Polygon::GetDistance2Boundary(const GeometryVector & pos, std::vector<GeometryVector> & list_dist, GeometryVector & x0, GeometryVector & n_outward) const{
	double distance = this->GetMaxRadius();
	list_dist.resize(0);
	list_dist.reserve(this->Vertices.size());
	GeometryVector side;
	for (size_t i=0; i<this->Vertices.size(); i++){
		double dist_Ei = 0.0;
		GeometryVector E = 2.0*this->Gamma.E[i];
		/* "t E_i + V_i" is the position of the projection of "pos" onto the infinitely long line connecting V_i and V_{i-1}.*/
		double t = - E.Dot(this->Vertices[i] - pos) / E.Modulus2();
		GeometryVector pt_polygon;
		
		/* if the projected point is between V_i and V_{i-1}, */
		if((-1.0<=t) && (t<= 0.0) ){
			pt_polygon = t * E + this->Vertices[i];
			dist_Ei = sqrt((pt_polygon - pos).Modulus2());
		}
		/* if the projected point is outside the line segment bewteen V_i and V_{i-1}, */
		else{
			double l1 = (pos - this->Vertices[i]).Modulus2();
			double l2 = (pos - (this->Vertices[i]- E)).Modulus2();

			if (l1 > l2){
				dist_Ei = sqrt(l2);
				pt_polygon = this->Vertices[i]- E;
			}
			else{
				dist_Ei = sqrt(l1);
				pt_polygon = this->Vertices[i];
			}
		}
		list_dist.emplace_back(dist_Ei, (double)i);

		if (i==0 || distance > dist_Ei){
			distance = dist_Ei;
			x0 = pt_polygon;
			side = E;
		}
		// else if (distance > dist_Ei){
		// 	distance = dist_Ei;
		// 	x0 = pt_polygon;
		// }
	}
	
	/* sort distances in the increasing order. */
	std::sort(list_dist.begin(), list_dist.end(), compare_temp);

	x0 = x0 + this->Centroid;

	if(this->Gamma.normal.x[2] > 0){
		/*A positive orientation*/
		n_outward = GeometryVector(side.x[1], -side.x[0]);
	}
	else{
		/*A negative orientation*/
		n_outward = GeometryVector(-side.x[1], side.x[0]);
	}
	n_outward = 1.0/(sqrt(n_outward.Modulus2()))*n_outward;

	return distance;
}

/**	Compute the closest distance from a given point to the side of this polygon.
 * @param pos	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
 * @param[out] list_dist	A list of distances in the ascending order.
 * @return distance.	 */
double Polygon::GetDistance2Boundary(const GeometryVector & pos) const{
	std::vector<GeometryVector> dump;
	GeometryVector x0, n;
	return Polygon::GetDistance2Boundary(pos, dump, x0, n);
}


/* -------------------------------------------------
 * Polyhedron constructors 
 * -------------------------------------------------*/

Polyhedron::Polyhedron(const std::vector<GeometryVector>& vertices, const std::vector<std::vector<size_t>>& idx) : Shape::Shape()
{
	this->Dimension = 3;
	this->name = "Polyhedron";
	//Check dimensions of vertices, and then input them
/* 	GeometryVector guessedInteriorPt((int)3);
	this->Vertices.resize(0);	this->Vertices.reserve(vertices.size());
	for (auto v = vertices.begin(); v != vertices.end(); v++) {
		if (v->Dimension == 3) {
			this->Vertices.push_back(*v);
			guessedInteriorPt = guessedInteriorPt + *v;
		}
		else
			std::cerr << "Polyhedron::Polyhedron\tInput vertices should be in 3 dimensions\n";
	}
	guessedInteriorPt = 1.0 / (double)this->Vertices.size() * guessedInteriorPt;
	//Translate the origin to the guessed interior point
	for (auto v = this->Vertices.begin(); v != this->Vertices.end(); v++)
		*v = *v - guessedInteriorPt; */
	
	this->Vertices.resize(0);	this->Vertices.reserve(vertices.size());
	for (auto v = vertices.begin(); v != vertices.end(); v++) {
		if (v->Dimension == 3) {
			this->Vertices.push_back(*v);
		}
		else
			std::cerr << "Polyhedron::Polyhedron\tInput vertices should be in 3 dimensions\n";
	}

	this->Vertex_idx = std::vector<std::vector<size_t>>(idx);
	if (this->Vertex_idx.size() > 3) {
		//Construct each face
		this->Gamma.resize(0);
		this->Gamma.reserve(this->Vertex_idx.size());
		for (auto v_idx = this->Vertex_idx.begin(); v_idx != this->Vertex_idx.end(); v_idx++) {
			PlanarFigure temp(&this->Vertices[0], *v_idx);
			this->Gamma.emplace_back(temp);
		}

		//Translate the center of mass to the origin
		GeometryVector CM = this->GetCenterOfMass();
		this->ShiftBy((-1.0)*CM);
		this->Centroid = CM;
		//Update the radii of circumscribing and inscribing spheres.
		this->UpdateR();
		//Update the surface area.
		this->UpdateSurfaceArea();

		//Check Polyhedron
		if (this->CheckEulerIdentity()) {
			this->CheckOrientations();
			if (!this->CheckClosedSuface())
				std::cerr<<"A Polyhedron must be surrounded by faces.\n";
		}		
	}
	else
		std::cerr << "Polyhedron::Polyhedron\tThere should be at least 4 faces" << std::endl;

}

Polyhedron::Polyhedron(std::istream & ifile) : Shape::Shape(){
	
	std::string phrase("#Vertices");
	size_t num_vertices = 0, num_faces = 0;
	
	/** Parameters for a basic constructor.*/
	std::vector<GeometryVector> vertices;
	std::vector<std::vector<size_t>> idx;
	
	
	std::string buffer; std::istringstream iss;
	/** Read vertices */
	while (num_vertices == 0 && std::getline(ifile, buffer)){
		if(buffer[0] == '#'){
			num_vertices = (size_t) std::stoi(buffer.substr(phrase.size()));
		}
	}

	if (num_vertices > 0){
		for (size_t i = 0; i < num_vertices; i++){
			std::getline(ifile, buffer);
			GeometryVector temp((int)3);
			iss.str(buffer);
			iss.clear();

			for (int j=0; j<3; j++)
				iss >> temp.x[j];
			
			if (ifile.eof() == false)
				vertices.push_back(temp);
		}
	}
	else{
		std::cerr << "Polygon::Polygon(ifile): Wrong input\n";
		return ;
	}

	/** Read idx. */

	while (num_faces == 0 && std::getline(ifile, buffer)){
		if(buffer[0] == '#'){
			num_faces = (size_t) std::stoi(buffer.substr(std::string("#Faces").size()));
		}
	}

	if (num_faces > 0){
		for (size_t i = 0; i < num_faces; i++){
			std::getline(ifile, buffer);
			std::vector<size_t> idx_face;
			iss.str(buffer);
			iss.clear();
			size_t temp = 0;

			iss >> temp;
			while(iss){
				idx_face.push_back(temp);
				iss >> temp;
			}

			idx.push_back(idx_face);

//			if (ifile.eof() == false)
//				vertices.push_back(temp);
		}
	}
	else{
		std::cerr << "Polyhedron::Polyhedron(ifile): Wrong input\n";
		return ;
	}


	/** The following is the same as the implementation of 
	 * Polyhedron::Polyhedron(const std::vector<GeometryVector>& vertices, const std::vector<std::vector<size_t>>& idx).*/
	this->Dimension = 3;
	this->name = "Polyhedron";
	
	this->Vertices=std::vector<GeometryVector>(vertices);
	this->Vertex_idx = std::vector<std::vector<size_t>>(idx);
	if (this->Vertex_idx.size() > 3) {
		//Construct each face
		this->Gamma.resize(0);
		this->Gamma.reserve(this->Vertex_idx.size());
		for (auto v_idx = this->Vertex_idx.begin(); v_idx != this->Vertex_idx.end(); v_idx++) {
			PlanarFigure temp(&this->Vertices[0], *v_idx);
			this->Gamma.emplace_back(temp);
		}

		//Translate the center of mass to the origin
		GeometryVector CM = this->GetCenterOfMass();
		this->ShiftBy((-1.0)*CM);
		this->Centroid = CM;
		//Update the radii of circumscribing and inscribing spheres.
		this->UpdateR();
		//Update the surface area.
		this->UpdateSurfaceArea();

		//Check Polyhedron
		if (this->CheckEulerIdentity()) {
			this->CheckOrientations();
			if (!this->CheckClosedSuface())
				std::cerr<<"A Polyhedron must be surrounded by faces.\n";
		}		
	}
	else
		std::cerr << "Polyhedron::Polyhedron\tThere should be at least 4 faces" << std::endl;


}

 /* -------------------------------------------------
  * Polyhedron member functions
  * -------------------------------------------------*/

/** Update this->R_circum and this->R_inscr with the largest distance from the centroid to a vertex,
 * and the smallest r_perp for faces, respectively. */
void Polyhedron::UpdateR() {
	double R2 = 0.0;
	for (auto iter = this->Vertices.begin(); iter != this->Vertices.end(); iter++)
		R2 = std::max(R2, iter->Modulus2());
	this->R_circum = sqrt(R2);

	double r = 0.0;
	for (auto face = this->Gamma.begin(); face != this->Gamma.end(); face++) {
		r = (face == this->Gamma.begin()) ? face->r_perp : std::min(face->r_perp, r);
	}
	this->R_inscr = r;
}
void Polyhedron::UpdateVolume() {
	double vol = 0.0;
	for (auto i = this->Gamma.begin(); i != this->Gamma.end(); i++)
		vol += i->GetVolume()*i->r_perp;
	this->volume = vol / 3.0;
}
void Polyhedron::UpdateSurfaceArea() {
	double surface_area = 0.0;
	for (auto i = this->Gamma.begin(); i != this->Gamma.end(); i++)
		surface_area += i->GetVolume();
	this->area = surface_area;
}
bool Polyhedron::CheckEulerIdentity() const {
	bool EulerIdentity, ClosedSurface;
	//Euler identity
	size_t V = this->Vertices.size(), F = this->Gamma.size(), E = 0;
	for (auto i = this->Gamma.begin(); i != this->Gamma.end(); i++)
		E += i->num_Vs;

	if (E % 2 != 0) {EulerIdentity = false;	}
	else {
		E /= 2;
		EulerIdentity = (V + F - E == 2);
	}

	return EulerIdentity;
}
void Polyhedron::ShiftBy(const GeometryVector & x) {
	//translate all vertices
	if(this->Centroid.Dimension == 0){
		for (auto v = this->Vertices.begin(); v != this->Vertices.end(); v++)
			*v = *v + x;
		//translate all R vectors of each face (if they are required)
		for (auto gamma = this->Gamma.begin(); gamma != this->Gamma.end(); gamma++) {
			gamma->Update_ER_vectors();
			//update r_perp
			gamma->CheckVertices();
		}
	}
	else
	{
		this->Centroid = this->Centroid + x;
	}
	
}
GeometryVector Polyhedron::GetCenterOfMass() {
	if(this->Centroid.Dimension == 0){
		GeometryVector CM((int)this->Dimension);
		double volume = 0.0;
		//A polygon is divided into pyramids whose bases are the faces of the polygon
		for (auto gamma = this->Gamma.begin(); gamma != this->Gamma.end(); gamma++) {
			// The center of mass of the asscoiate pyramid whose base is gamma
			GeometryVector cm_pyramid = 0.75*gamma->GetCenterOfMass();
			// The volume of this pyramid
			double vol_pyramid = 1.0 / 3.0* gamma->r_perp * gamma->GetVolume();
			CM = CM + vol_pyramid * cm_pyramid;		
			volume += vol_pyramid;
		}
		return 1.0 / volume * CM;

	}
	else{
		return this->Centroid;
	}
}

bool Polyhedron::CheckClosedSuface() {
	bool ClosedSurface = false;
	if (this->Gamma.size() > 3) {
		GeometryVector sum((int)3);
		for (auto iter = this->Gamma.begin(); iter != this->Gamma.end(); iter++) {
			sum = sum + iter->GetVolume()*iter->normal;
		}
		//ClosedSurface = (sum.IsEqual(GeometryVector((int)3)));
		ClosedSurface = (sqrt(sum.Modulus2()) < 1e-17);
	}
	else {
		std::cerr << "Polyhedron:: at least 3 faces are required " << std::endl;
	}
	
	return ClosedSurface;
}

void Polyhedron::CheckOrientations() {
	for (size_t idx = 0; idx < this->Gamma.size(); idx++) {
		if (this->Gamma[idx].r_perp < 0) {//should invert the orientation
			std::reverse(this->Vertex_idx[idx].begin(), this->Vertex_idx[idx].end());
			this->Gamma[idx] = PlanarFigure(&(this->Vertices[0]), this->Vertex_idx[idx]);
			this->Gamma[idx].Update_ER_vectors();
			this->Gamma[idx].CheckVertices();
		}
	}
}

void Polyhedron::Print(std::ostream & ofile) const {
	ofile << "[" << this->name << "]\n";

	ofile << "#Vertices\t" << this->Vertices.size() << "\n";
	for (auto iter = this->Vertices.begin(); iter != this->Vertices.end(); iter++)
		ofile << (*iter) + this->Centroid;

	ofile << "#Faces\t" << this->Vertex_idx.size() << "\n";
	for (auto iter = this->Vertex_idx.begin(); iter != this->Vertex_idx.end(); iter++) {
		for (auto elem = iter->begin(); elem != iter->end(); elem++)
			ofile << *elem << "\t";
		ofile << "\n";
	}
}
void Polyhedron::GetVertices(std::vector<GeometryVector> & vertices, std::vector<std::vector<size_t>> & idx) const {
	vertices.resize(0);	idx.resize(0);
	vertices = std::vector<GeometryVector>(this->Vertices);
	for (auto v = vertices.begin(); v != vertices.end(); v++)
		*v = *v + this->Centroid;

	idx = std::vector<std::vector<size_t>>(this->Vertex_idx);
}

void Polyhedron::Rotate(const std::vector<GeometryVector>& newbasis){
	if (newbasis.size() == 3) {
		//Rotate this->Vertices
		for (auto v = this->Vertices.begin(); v != this->Vertices.end(); v++) {
			GeometryVector temp((int)this->Dimension);
			for (int j = 0; j < this->Dimension; j++) {
				temp.x[j] += newbasis[j].Dot(*v);
			}
			*v = temp;
		}
		//Rotate normal vectors and E and R vectors of each face
		for (auto gamma = this->Gamma.begin(); gamma != this->Gamma.end(); gamma++) {
			gamma->Rotate(newbasis);
			gamma->Update_ER_vectors();
		}
	}
	else
		std::cerr << "Polyhedron::Rotate\tWrong dimensions of newbasis\n";
}

void Polyhedron::Scale(double factor){
	Shape::Scale(factor);
//	this->volume *= factor * factor*factor;
//	this->area *= factor * factor;
//	this->R_circum *= factor;
//	this->R_inscr *= factor;
	for (auto v = this->Vertices.begin(); v != this->Vertices.end(); v++)
		*v = factor * (*v);

	//scale each face
	for (auto gamma = this->Gamma.begin(); gamma != this->Gamma.end(); gamma++)
		gamma->Scale(factor);
}

/** Check whether a point is inside this Polyhedron or not.
 * Here, we assume that all normal vectors of faces are outward. */
bool Polyhedron::IsPtInside(const GeometryVector & pos) {
	assert(pos.Dimension == 3);
	double R = this->GetMaxRadius(),
		r = this->GetInscrRadius(),
		dist2pos = sqrt(pos.Modulus2());

	if (dist2pos > R) {
		return false;
	}
	else {
		// pos is inside the circumscribing sphere.
		
		if (dist2pos < r) {
			return true;
		}
		else {
			//We need to scrutinize.
			bool isInside = true;

			for (auto face = this->Gamma.begin(); face != this->Gamma.end() && isInside; face++) {
				isInside = face->r_perp >= face->normal.Dot(pos);
			}
			return isInside;
		}
	}

}

/**	Compute the closest distance from a given point to the facets of this polygon.
 * TODO: implement it!
 * @param pos	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
 * @return distance.	 */
double Polyhedron::GetDistance2Boundary(const GeometryVector & pos, std::vector<GeometryVector> & list_dist, GeometryVector & x0, GeometryVector & n_outward) const
{
	std::cout<<"Please implement!!\n";
	return 0;
}

double Polyhedron::GetDistance2Boundary(const GeometryVector & pos) const
{
	std::vector<GeometryVector> dump;
	GeometryVector x0, n;
	return Polyhedron::GetDistance2Boundary(pos, dump, x0, n);
}

///* -------------------------------------------------
// * Parallelepiped constructor 
// * -------------------------------------------------*/
//
//Parallelepiped::Parallelepiped(const Parallelepiped & source) : Shape::Shape(source), a(source.a), b(source.b) {
//}
//
//Parallelepiped::Parallelepiped(DimensionType Dimension, const std::vector<GeometryVector> & sides) {
//	if (Dimension == sides.size()) {
//		this->Dimension = Dimension;
//		this->a = std::vector<GeometryVector>(sides);
//		switch (this->Dimension)
//		{
//		case 1:
//			this->num_edges = 1;
//			this->num_vertices = 2;
//			this->num_faces = 0;	break;
//		case 2:
//			this->num_edges = 4;
//			this->num_vertices = 4;
//			this->num_faces = 1;	break;
//		case 3:
//			this->num_edges = 12;
//			this->num_vertices = 8;
//			this->num_faces = 6;	break;
//		default:
//			std::cout << "Please implement for d<=3 cases\n";
//			this->num_edges = 1;
//			this->num_vertices = 1;
//			this->num_faces = 1;
//			break;
//		}
//		this->UpdateR();
//		this->area = -1.0;
//		this->vol = -1.0;
//	}
//	else {
//		std::cerr << "Dimensions and the number of sides are inconsistent\n";
//	}
//}
//
//Parallelepiped::Parallelepiped(DimensionType Dimension, double L) : Parallelepiped::Parallelepiped(Dimension, std::vector<GeometryVector> (Dimension, GeometryVector(static_cast<int> (Dimension)))){
//	for (int i = 0; i < this->Dimension; i++) {
//		GeometryVector e((int)this->Dimension);	e.x[i] = L;
//		this->a[i] = e;
//	}
//
//	/* update volume */
//	this->vol = 1.0;
//	for (int i = 0; i < Dimension; i++)
//		this->vol *= L;
//	/* update surface area*/
//	if (Dimension == 2)
//		this->area = 4.0*L;
//	else if (Dimension == 3)
//		this->area = 6.0*L*L;
//	else if (Dimension == 1)
//		this->area = 0.0;
//	else {
//		std::cerr << "Dimension is too high\n";
//	}
//	/* update exclusion radius*/
//	this->R_circum = 0.5 * sqrt(Dimension) * L;
//}
//
//void Parallelepiped::GetVertices(std::vector<GeometryVector> & vertices) const {
//	vertices.resize(0);
//	vertices.reserve(this->num_vertices);
//	int count = 0;
//
//	n_naryNumber pos(this->Dimension, 2);
//	for (int i = 0; i < this->num_vertices; i++) {
//		GeometryVector x((int)this->Dimension);
//		
//		for (int j = 0; j < this->Dimension; j++)
//			x = x + (pos.GetNthDigit(j) - 0.5)*this->a[j];
//
//		vertices.push_back(x);
//		if(i < this->num_vertices-1)
//			pos++;
//	}
//}
//
//void Parallelepiped::Print() const {
//	for (DimensionType i = 0; i < this->Dimension; i++)
//		std::cout << this->a[i];
//}
////Dimension > 4 is undefined
//double Parallelepiped::GetVolume() const {
//	if (this->vol < 0)
//		this->vol = ::Volume(&this->a[0], (DimensionType)this->Dimension);
//
//	return this->vol;
//}
////Dimension > 3 is undefined
//double Parallelepiped::GetSurfaceArea() const {
//	if (this->area < 0) {
//		if (this->Dimension == 1)
//			this->area = 0;
//		else if (this->Dimension == 2)
//			this->area = 2.0*(sqrt(this->a[0].Modulus2()) + sqrt(this->a[1].Modulus2()));
//		else if (this->Dimension == 3) {
//			this->area = sqrt(this->a[0].Modulus2() * this->a[1].Modulus2() - this->a[0].Dot(this->a[1])*this->a[0].Dot(this->a[1]));
//			this->area += sqrt(this->a[1].Modulus2() * this->a[2].Modulus2() - this->a[1].Dot(this->a[2])*this->a[1].Dot(this->a[2]));
//			this->area += sqrt(this->a[2].Modulus2() * this->a[0].Modulus2() - this->a[2].Dot(this->a[0])*this->a[2].Dot(this->a[0]));
//			this->area *= 2.0;
//		}
//		else {
//			std::cerr << "GetSurfaceArea() is undefined \n";
//		}
//	}
//	return this->area;
//}
//
//void Parallelepiped::Rotate(const std::vector<GeometryVector> & newcoordinatebasis) {
//	if (CheckNewCoordinates(this->Dimension, newcoordinatebasis)) {
//		for (int i = 0; i < this->Dimension; i++) {
//			GeometryVector temp((int)this->Dimension);
//			for (int j = 0; j < this->Dimension; j++) {
//				temp.x[j] += newcoordinatebasis[j].Dot(this->a[i]);
//			}
//			this->a[i] = temp;
//		}
//	}
//	else {
//		std::cerr << "input basis vectors are incorrect!\n";
//	}
//};
//
//void Parallelepiped::Scale(double factor){
//	Shape::Scale(factor);
//	for (int i = 0; i < this->Dimension; i++) {
//		this->a[i] = factor * this->a[i];
//	}
//}
//
//void Parallelepiped::UpdateR() const{
//	GeometryVector maxdist((int)this->Dimension);
//	for (DimensionType i = 0; i<this->Dimension; i++)
//	{
//		if (maxdist.Dot(this->a[i])>0)
//			maxdist.AddFrom(this->a[i]);
//		else
//			maxdist.MinusFrom(-1.0*this->a[i]);
//	}
//	this->R_circum = 0.5*std::sqrt(maxdist.Modulus2());
//}
//
//
///* -------------------------------------
// * Cube  member functions 
// * -------------------------------------*/
//Cube::Cube(const Cube & source) : Parallelepiped::Parallelepiped(source), L(source.L), halfL(source.halfL){
//}
//
//Cube::Cube(DimensionType Dimension, double L) : Parallelepiped::Parallelepiped(Dimension, L), L(L), halfL(0.5*L) {
//}
//
////double Cube::GetVolume() const {return Parallelepiped::GetVolume();}
////void Cube::Rotate(std::vector<GeometryVector> & newbasis) { Parallelepiped::Rotate(newbasis); }
////void Cube::Rotate(RandomGenerator & rng) { Shape::Rotate(rng); }
