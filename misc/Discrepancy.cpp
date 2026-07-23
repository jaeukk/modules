#include <Discrepancy.h>

//My dirty version of merge sort [in the increasing order]
void MergeSort(std::vector<double> &stack) {
	if (stack.size() > 2) {
		int size = stack.size();

		int bin = 1;
		int pre_bin = 1;
		while (pre_bin < size) {
			if (bin > 1) {
				std::vector<double> tempList;
				tempList.reserve(size);

				int j1 = 0; //index for sublist1
				int j2 = pre_bin; //index for sublist2
				int j1_max = pre_bin;
				int j2_max = (bin>size) ? size : bin;
				while (tempList.size()<(size_t)size && j1 < size) {
					double temp;

					if (j1 < j1_max && j2<j2_max) {
						if (stack[j1] > stack[j2]) {
							//In ascending order (true), tempList[j1] > tempList[j2]. In descending order(false), tempList[j1] < tempList[j2]
							temp = stack[j2];
							j2++;
						}
						else {
							temp = stack[j1];
							j1++;
						}
					}
					else if (j2 == j2_max) {//elements left in sublist 1 will be larger than the last element in the sublist 2
						temp = stack[j1]; j1++;
					}
					else {//elements left in sublist 2 will be larger than the last element in the sublist 1
						temp = stack[j2]; j2++;
					}
					tempList.push_back(temp);

					if (j1 == j1_max && j2 == j2_max) { //complete to sort both sublists
						j1_max = ((j1_max + bin)>size) ? size : j1_max + bin;
						j1 = ((j1 + pre_bin) > size) ? size : j1 + pre_bin;
						j2_max = ((j2_max + bin) > size) ? size : j2_max + bin;
						j2 = ((j2 + pre_bin) > size) ? size : j2 + pre_bin;
					}
				}

				//copy the contents of tempList into stack
				for (int k = 0; k < size; k++) {
					stack.at(k) = tempList[k];
				}
			}

			pre_bin = bin;
			bin *= 2;
		}
	}
	else if (stack.size()==2){// stack size <=2
		std::vector<double> temp(stack);
		stack[0] = std::min(temp[0], temp[temp.size() - 1]);
		stack[stack.size()-1] = std::max(temp[0], temp[temp.size() - 1]);
		//std::cout << "MergeSort::stack is too small";
	}
	else {
		//do nothing
		//stack is too small
	}

}


//My dirty version of merge sort (non-recursive)
void MergeSort(std::vector<GeometryVector> & stack, int std_col) {
	if (stack.empty()) {//don't dereference stack[0] on an empty vector
		return;
	}
	if (std_col < 0 || std_col+1 > stack[0].Dimension) {
		std::cout << "Discrepancy::MergeSort - wrong standard column \n";	return;
	}

	int size = stack.size();

	int bin = 1;
	int pre_bin = 1;
	while (pre_bin < size) {
		if (bin > 1) {
			std::vector<GeometryVector> tempList;
			tempList.reserve(size);

			int j1 = 0; //index for sublist1
			int j2 = pre_bin; //index for sublist2
			int j1_max = pre_bin;
			int j2_max = (bin>size) ? size : bin;
			while (tempList.size()<(size_t)size && j1 < size) {
				GeometryVector temp;

				if (j1 < j1_max && j2<j2_max) {
					if (stack[j1].x[std_col] > stack[j2].x[std_col]) {
						//In ascending order (true), tempList[j1] > tempList[j2]. In descending order(false), tempList[j1] < tempList[j2]
						temp = stack[j2];
						j2++;
					}
					else {
						temp = stack[j1];
						j1++;
					}
				}
				else if (j2 == j2_max) {//elements left in sublist 1 will be larger than the last element in the sublist 2
					temp = stack[j1]; j1++;
				}
				else {//elements left in sublist 2 will be larger than the last element in the sublist 1
					temp = stack[j2]; j2++;
				}
				tempList.push_back(temp);

				if (j1 == j1_max && j2 == j2_max) { //complete to sort both sublists
					j1_max = ((j1_max + bin)>size) ? size : j1_max + bin;
					j1 = ((j1 + pre_bin) > size) ? size : j1 + pre_bin;
					j2_max = ((j2_max + bin) > size) ? size : j2_max + bin;
					j2 = ((j2 + pre_bin) > size) ? size : j2 + pre_bin;
				}
			}

			//copy the contents of tempList into stack
			for (int k = 0; k < size; k++) {
				stack.at(k) = tempList[k];
			}
		}

		pre_bin = bin;
		bin *= 2;
	}
}


//template<typename T>
//std::vector<size_t> Crop(const PeriodicCellList<T> &config, const std::function<int(const GeometryVector &)> w, const GeometryVector &CartesianCoordinates, double R) {

//}


//template<typename t> 
//std::vector<size_t> Crop(const PeriodicCellList<t> & config, const GeometryVector &L, const GeometryVector &U, int std_col, bool isSorted) {
//	if (L.Dimension != U.Dimension) {
//		std::cout << "Discrepancy::Crop - wrong dimensions in L and U \n";
//		return ;
//	}
//	else if (L.Dimension != config.GetDimension()) {
//		std::cout << "Discrepancy::Crop - wrong dimensions in LU and config \n";
//		return ;
//	}
//	else {
//		const DimensionType d = config.GetDimension();
//		const size_t Num = config.NumParticle();
//		bool doesCorrectBounds = true;
//		for (DimensionType i = 0; i < d && doesCorrectBounds; i++) {
//			doesCorrectBounds = (L.x[i] < U.x[i]) && (U.x[i] < 1.0);
//		}
//		if (!doesCorrectBounds) {
//			std::cout << "Discrepancy::Crop - wrong bounds L and U";
//			return;
//		}
//
//		PeriodicCellList<t> temp(config);
//		if (!isSorted) {
//			MergeSort(temp, std_col);
//		}
//
//		int lowerIndex = 0;
//		int upperIndex = Num;
//
//		//Determine indices of particles whose std_col - th component is within the range.
//		for (lowerIndex = 0; (lowerIndex < Num)  && (temp.GetRelativeCoordinates(lowerIndex)<L.x[std_col]) ;lowerIndex++) {}
//		for (upperIndex = Num - 1; (upperIndex>=0) && (temp.GetRelativeCoordinates(upperIndex)>U.x[std_col]) ;upperIndex--){}
//		
//		//Template for the result configuration
//		GeometryVector *basis;	basis = new GeometryVector[d];
//		//Think about the way using pointer...
//		for (auto i = 0; i < d; i++) {
//			basis[i] = GeometryVector(temp.GetBasisVector(i));
//		}
//		PeriodicCellList<t> result(d, basis, temp.GetCellSize());
//
//		//can be improved by sorting particle indices according to the other components.
//		//But here, use a bit brutal force method.
//		for (int i = lowerIndex; i <= upperIndex; i++) {
//			bool isWithinTheCube = true;
//			for (int j = 0; (j < d) && isWithinTheCube; j++) {
//				if (j != std_col) 
//					isWithinTheCube = (L.x[j]<temp.GetRelativeCoordinates(j)) && (temp.GetRelativeCoordinates(j)<U.x[j]);
//			}
//
//			if (isWithinTheCube) {
//				result.Insert(temp.GetCharacteristics(i), temp.GetRelativeCoordinates(i));
//			}
//		}
//		return result;
//	}
//}
double vanDerCorputSequence(int n, int basis) {
	if (basis < 2)//radix must be >= 2; otherwise temp%basis / temp/basis divide by zero
		return 0.0;
	double n_basis = 0;
	double digit = 1.0;
	int temp = n;
	while (temp > 0) {
		digit /= (double)basis;
		n_basis += (temp%basis)*digit;
		temp = (int)floor(temp / basis);
	}
	return n_basis;
}

GeometryVector vanDerCorputSequence(int n, const std::vector<int> &basis) {
	GeometryVector result(static_cast<DimensionType>(basis.size()));
	if (basis.size() < 2) {
		std::cout << "Discrepancy::vanDerCorputSequence : wrong input\n";
		return result;
	}
	else if (basis[0] <= 0) {
		std::cout << "Discrepancy::vanDerCorputSequence : basis[0] must be positive\n";
		return result;
	}
	else {
		result.x[0] = (double)n / basis[0];
		for (size_t i = 1; i < basis.size(); i++) {
			result.x[i] = vanDerCorputSequence(n, basis[i]);
		}

		return result;

	}


}

std::vector<GeometryVector> vanDerCorputPts(int n, const std::vector<int> &basis) {
	std::vector<GeometryVector> result;
	if (n < 0) {//negative n would reserve a huge size_t and the loop would produce nothing
		std::cout << "Discrepancy::vanDerCorputPts : n must be nonnegative\n";
		return result;
	}
	result.reserve(n);

	for (int i = 1; i < n+1; i++) {
		result.push_back(vanDerCorputSequence(i, basis));
	}
	return result;
}
