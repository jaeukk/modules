#include "RandomOrganization.h"


//Berthier's model



//----------------For square systems
void IsActive(std::vector<GeometryVector> & List, std::vector<int> & result, double exclusionR, double SystSize) {
	//without parallel computing
	std::vector<int> TempList(List.size(), 0);
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 1; i < (int)List.size(); i++) {
		for (int j = 0; j < i; j++) {
			GeometryVector temp = List[i] - List[j];
			for (int k = 0; k < temp.Dimension; k++) {//Implement periodic boundary conditions
				if (fabs(temp.x[k] > SystSize*0.5)) {
					temp.x[k] = SystSize - fabs(temp.x[k]);
				}
			}

			if (sqrt(temp.Modulus2()) < exclusionR) {
				TempList.at(i)++; TempList.at(j)++;
			}
		}
	}

	//Output 
	result.resize(0);
	result.reserve(List.size());
	for (int i = 0; i < (int)TempList.size(); i++) {
		if (TempList.size() > 0)	result.push_back(i);
	}

	
}

bool B_Model(std::vector<GeometryVector> & List, double exclusionR, double delta, double SystSize) {
	std::vector<int> index;
	IsActive(List, index, exclusionR, SystSize);

	if (index.size() > 0) {
		DimensionType dim = List[0].Dimension;
		for (int i = 0; i < (int)index.size(); i++) {
			int j = index[i];
			GeometryVector dr = delta*RandomUnitVector(dim);
			List.at(j) = List[j] + dr;
		}
		return true;
	}
	else {		return false;	}

}



//------------------For general periodic systems




void IsActive(const Configuration & List, std::vector<size_t> & ActivePrt_index, double exclusionR) {
	double diameter2 = 4.0*exclusionR*exclusionR;
	size_t Num = List.NumParticle();
	ActivePrt_index.resize(0);
	ActivePrt_index.reserve(Num);
	
	//Count(i) is the number that i-th particle is verified as active : Count(i) == 0 => particle i is passive
	std::vector<int> Count(Num, 0);
	for (size_t i = 0; i < Num; i++) {
		//#pragma omp for parallel
		if (Count[i] == 0) {//Check whether i-th particle is active or not, when i-th particle may not be active particle
			List.IterateThroughNeighbors(List.GetRelativeCoordinates(i), 2.0*exclusionR, 
				[&Count, &i, &diameter2, &Num](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)->void 
			{
				if (Sourceparticle != i && shift.Modulus2() < diameter2 && Sourceparticle < Num) {
					Count.at(Sourceparticle)++;
					Count.at(i)++;
				}
					
			});
		}
	}
	
	for (size_t i = 0; i < Count.size(); i++) {
		if (Count[i] > 0) {
			ActivePrt_index.push_back(i);
		}
	}
}

void IsActive(const PeriodicCellList<double> &List, std::vector<size_t> & ActivePrt_index, double exclusionR) {
	//List.ParticleCharacteristics = particle exclusion raodis
	size_t Num = List.NumParticle();
	ActivePrt_index.resize(0);
	ActivePrt_index.reserve(Num);
	
	std::vector<int> Count(Num, 0);
	for (size_t i = 0; i < Num; i++) {
		//#pragma omp for parallel
		if (Count[i] == 0) {//Check whether i-th particle is active or not, when i-th particle may not be active particle
			List.IterateThroughNeighbors(List.GetRelativeCoordinates(i), 2.0*exclusionR,
				[&Count, &i, &List, &Num](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)->void
			{
				if (Sourceparticle != i && sqrt(shift.Modulus2()) < List.GetCharacteristics(i)+List.GetCharacteristics(Sourceparticle) && Sourceparticle < Num) {
					Count.at(Sourceparticle)++;
					Count.at(i)++;
				}

			});
		}
	}

	for (size_t i = 0; i < Count.size(); i++) {
		if (Count[i] > 0) {
			ActivePrt_index.push_back(i);
		}
	}


}

void B_model(Configuration &List, double delta, const std::vector<size_t> & ActivePrt_index, RandomGenerator & g) {
	GeometryVector new_r;
	DimensionType d = List.GetDimension();
	//assume periodic cell is a simple cubic..
//#pragma omp parallel for
	for (int i = 0; i < (int)ActivePrt_index.size(); i++) {
		int t = ActivePrt_index[i];
		//new_r = List.GetRelativeCoordinates(i) + delta*RandomUnitVector(d,g);
		List.MoveParticle(t, List.GetRelativeCoordinates(t) + delta*RandomUnitVector(d, g));
	}
}

void B_model(PeriodicCellList<double> &List, double delta, const std::vector<size_t> & activePrt_index, RandomGenerator &g) {
	GeometryVector new_r;
	DimensionType d = List.GetDimension();
	//assume periodic cell is a simple cubic..
	//#pragma omp parallel for
	for (int i = 0; i < (int)activePrt_index.size(); i++) {
		int t = activePrt_index[i];
		//new_r = List.GetRelativeCoordinates(i) + delta*particle radius*RandomUnitVector(d,g);
		List.MoveParticle(t, List.GetRelativeCoordinates(t) + delta*List.GetCharacteristics(t)* RandomUnitVector(d, g));
	}
}




void B_model(Configuration & List, double Rc, double delta, std::vector<size_t> & ActivePrt_index, RandomGenerator & g) {

	IsActive(List, ActivePrt_index, Rc);
	
	B_model(List, delta, ActivePrt_index,g);
}