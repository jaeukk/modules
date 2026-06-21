#ifndef DISCREPANCY_H
#define DISCREPANCY_H

#include <vector>
#include <functional>
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <etc.h>

void MergeSort(std::vector<double> &stack);
//sorts data in stack in the ascending order according to std_col-th component.
void MergeSort(std::vector<GeometryVector> & stack, int std_col);

//sorts index of particles in the ascending order according to std_col-th components of relative coordinates.
template<typename T> void Merge_Sort(PeriodicCellList<T> & stack, int std_col) {
	if (stack.GetDimension() < std_col + 1) {
		std::cout << "Discrepancy::MergeSort - wrong standard column \n";	return;
	}

	std::vector<GeometryVector> tempConfig;	//relative coordinates
	std::vector<T> tempCharacteristics;
	tempConfig.reserve(stack.NumParticle());
	tempCharacteristics.reserve(stack.NumParticle());

	if (stack.NumParticle() < 1e+6 && stack.GetDimension()>1) {
		//Make an index
		GeometryVector index_vec(stack.GetDimension());
		const int indexing_dimension = (0 == std_col) ? 1 : 0;
		for (int i = 0; i < stack.NumParticle(); i++) {
			tempConfig.push_back(stack.GetRelativeCoordinates(i) + index_vec); //<= the relative coordinates indicate the original ordering
			tempCharacteristics.push_back(stack.GetCharacteristics(i));
			index_vec.x[indexing_dimension]++;
		}

		MergeSort(tempConfig, std_col);
		stack.RemoveParticles();
		int t = 0;
		for (int i = 0; i < tempConfig.size(); i++) {
			t = (int)floor(tempConfig[i].x[indexing_dimension]);	//here t indicates the index of the particle
			tempConfig[i].x[indexing_dimension] -= (double)t;
			stack.Insert(tempCharacteristics[t], tempConfig[i]);
		}
		tempConfig.clear();
		tempCharacteristics.clear();
	}
	else if (stack.GetDimension() == 1 && std_col == 0) {//1D system
		//Make an index
		GeometryVector index_vec(2);
		const int indexing_dimension = 1;
		for (int i = 0; i < stack.NumParticle(); i++) {
			tempConfig.push_back(GeometryVector(stack.GetRelativeCoordinates(i).x[0],0) + index_vec); //<= the relative coordinates indicate the original ordering
			tempCharacteristics.push_back(stack.GetCharacteristics(i));
			index_vec.x[indexing_dimension]++;
		}

		MergeSort(tempConfig, std_col);
		stack.RemoveParticles();
		int t = 0;
		for (int i = 0; i < tempConfig.size(); i++) {
			t = (int)round(tempConfig[i].x[indexing_dimension]);	//here t indicates the index of the particle
			stack.Insert(tempCharacteristics[t], GeometryVector(tempConfig[i].x[0]));
		}
		tempConfig.clear();
		tempCharacteristics.clear();
	}
	else {//Need to think other method...
		std::cout << "Discrepancy::MergeSort - consider the case where 1D configuration and the number of particles is more than 1 million.";
	}
}


//For a given point configuration, indices of particles satifying w(x-x0)<r are returned.
//w(x) : argument x is in Cartesian coordinates (not relative coordinates)
//x0 : Cartesian coordinates
template<typename T> 
std::vector<size_t> Crop(const PeriodicCellList<T> &config, const std::function<double(const GeometryVector &)> w, const GeometryVector &x0, double r) {
	double prefactor = config.GetDimension();
	std::vector<size_t> result;
		result.reserve(config.NumParticle());
	
		GeometryVector r0(config.CartesianCoord2RelativeCoord(x0));
	
		//redefine the window indicator function considering the periodic boundary condition
		//auto W = [&w,&config](const GeometryVector &X)->double {
			//GeometryVector temp = config.CartesianCoord2RelativeCoord(X);
			//config.RelativeCoordToMinimumImage(temp);
			//return w(config.RelativeCoord2CartesianCoord(temp));
		//};

		config.IterateThroughNeighbors(r0, prefactor*r, [&result, &r, &w](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)->void {
			if (w(shift) < r) {//sourceparticle belongs to the window
				result.push_back(Sourceparticle);
			}}
		);
		return result;
}

//Using the window rotated by R (rotation matrix),
//return Cartesian coordinates of points satisfying w(R(x-x0))<r with respect to the rotated coordinates systems.
//x0 : in Cartesian coordinates
template<typename T>
std::vector<GeometryVector> Crop(const PeriodicCellList<T> &config, const std::function<double(const GeometryVector &)> w, const GeometryVector & x0, double r, std::vector<GeometryVector> R) {
	if (R.size() == 0) { //check empty R
		std::vector<GeometryVector> t;
		std::cout << "Discrepancy::Crop : R is empty \n";
		return  t;
	}
	else {//Check the dimension of vectors of R
		int i = 0;	bool coincident = true;
		while (i < R.size() && coincident) {
			if (i == 0) { coincident = (config.GetDimension() == R[0].Dimension);}
			else {	coincident = (R[i].Dimension == R[i - 1].Dimension);}
			i++;
		}
		if (coincident) {
			
			auto w_ = [&w, &R](const GeometryVector X) -> double {
				GeometryVector x_(X); for (int i = 0; i < x_.Dimension; i++) { x_.x[i] = X.Dot(R[i]); } return w(x_);};
			std::vector<size_t> indices= Crop(config, w_, x0, r);
			std::vector<GeometryVector> cartesianCoord; cartesianCoord.reserve(indices.size());
			
			for (size_t j = 0; j < indices.size(); j++) {
				//GeometryVector temp(config.GetCartesianCoordinates(indices[j]) - x0); //before the transform
				//GeometryVector temp_(temp);
				GeometryVector temp = config.CartesianCoord2RelativeCoord(config.GetCartesianCoordinates(indices[j]) - x0); //before the transform
				config.RelativeCoordToMinimumImage(temp);
				GeometryVector temp_(config.RelativeCoord2CartesianCoord(temp));

				for (int k = 0; k < temp.Dimension; k++) {
					temp.x[k] = temp_.Dot(R[k]);
				}
				cartesianCoord.push_back(temp);
			}
			return cartesianCoord;
		}
		else {
			std::vector<GeometryVector> t;
			std::cout << "Discrepancy::Crop : sth wrong in R \n";
			return t;
		}




	}



}

double vanDerCorputSequence(int n, int basis);

//Compute d-dimensional van der Corput sequence
//n = numbering
//basis[0] = Max;
//basis[i] = basis for i-th coordinates
GeometryVector vanDerCorputSequence(int n, const std::vector<int> &basis);

//compute a set of van der Corput points (n points)
std::vector<GeometryVector> vanDerCorputPts(int n, const std::vector<int> &basis);


//Compute the star discrepancy
template<typename T>
double Star_discrepancy(const PeriodicCellList<T> & config, int std_col = 0, bool isSorted = true);

//measure the maximum discrepancy with a cubic window
template<typename T>
double max_discrepnacy(const PeriodicCellList<T> & config, int std_col = 0, bool isSorted = true);

//measure the maximum discrepancy with a arbitrary shaped window
template<typename T>
double max_discrepancy(const PeriodicCellList<T> & config, const std::function<double(const GeometryVector &)> w);


#endif