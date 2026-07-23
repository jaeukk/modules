/** Author	:	Jaeuk Kim
*	Email	:	phy000.kim@gmail.com
*	Date	:	March	2021 */

/** \file TwoPhaseDispersion.cpp
*	\brief Function implementations of TwoPhaseDispersion class. */

#include "TwoPhaseDispersion.h"

void TwoPhaseDispersion::SetResolutions(const std::vector<size_t> & num_voxels) {
	if (num_voxels.size() == this->packing.GetDimension()) {
		this->num_side_voxels = std::vector<size_t>(num_voxels);
		this->num_total_voxels = 1;
		for (auto iter = num_voxels.begin(); iter != num_voxels.end(); iter++)
			this->num_total_voxels *= (*iter);
	}

}

void TwoPhaseDispersion::OptimizeSize()
{
	std::cerr << "OptimizeSize() is Not implemented yet\n";
	//TODO:: implement to adjust cell sizes to be omptimal for FFTW.
}

//void TwoPhaseDispersion::ComputeIndicator() {
//	/* generate grid */
//	std::vector<GeometryVector> basis;	/* basis vector for the grid*/
//	for (int i = 0; i < this->packing.GetDimension(); i++) {
//		basis.emplace_back(this->packing.GetBasisVector(i));
//		basis[i] = 1.0 / (double)this->num_side_voxels[i] * basis[i];
//	}
//	Configuration grid;
//	{
//		Configuration single_cell(this->packing.GetDimension(), &basis[0], 1.0);
//		GeometryVector x0(static_cast<int> (this->packing.GetDimension()));
//		for (int i = 0; i < x0.Dimension; i++)
//			x0.x[i] = 0.5;
//		single_cell.Insert("a", x0);
//		grid = MultiplicateStructure(single_cell, this->num_side_voxels);
//	}
//	grid.PrepareIterateThroughNeighbors(this->packing.GetMaxRadius());
//
//	/* prt_indicator[i] = grid coordinates of the particle phase in a raw-major format. */
//	this->prt_indicator.resize(grid.NumParticle(), 0); 
//	size_t index = 0;
//
//	/* check whether the center of each grid is inside the particle phase */
//	bool DoOverlap = this->packing.CheckOverlap();
//	double cov_fraction = 0;
//	if (DoOverlap) {
//		/* this estimates how many grid numbers will be counted. */
//		cov_fraction = this->packing.PackingFraction(); 
//	}
//
//	if (DoOverlap && (cov_fraction > 0.9)) {
//		/* we need extra space for prt_indicator temporarily. */
//		double extra_frac = (cov_fraction > 1.1) ? cov_fraction : 1.1;
//		int est_num_grids = (int)round(extra_frac*grid.NumParticle());
//		this->prt_indicator.resize(est_num_grids, 0);
//#pragma omp parallel for
//		for (int i = 0; i < this->packing.NumParticle(); i++) {
//			GeometryVector center = this->packing.GetRelativeCoordinates(i);
//			double r = this->packing.GetCharacteristics(i);
//
//			grid.IterateThroughNeighbors(center, r, [&basis, &grid, &r, &index, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long *periodicshift, size_t SourceParticle)
//				-> void {
//				GeometryVector prt = grid.GetRelativeCoordinates(SourceParticle);
//				std::vector<size_t> component(prt.Dimension, 0);	/*coordinates of grid*/
//				for (int j = 0; j < prt.Dimension; j++) {
//					component[j] = (size_t)floor(prt.x[j] * this->num_side_voxels[j]);
//				}
//				size_t temp = this->Components2Index(component);
//
//#pragma omp critical
//				{
//					if (index < this->prt_indicator.size())
//						this->prt_indicator[index] = temp;
//					else
//						this->prt_indicator.push_back(temp);
//					index++;
//				}
//			});
//		}
//	}
//	else {
//		/* we don't need to consider overflow. */
//#pragma omp parallel for
//		for (int i = 0; i < this->packing.NumParticle(); i++) {
//			GeometryVector center = this->packing.GetRelativeCoordinates(i);
//			double r = this->packing.GetCharacteristics(i);
//
//			grid.IterateThroughNeighbors(center, r, [&basis, &grid, &r, &index, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long *periodicshift, size_t SourceParticle)
//				-> void {
//					GeometryVector prt = grid.GetRelativeCoordinates(SourceParticle);
//					std::vector<size_t> component(prt.Dimension, 0);	/*coordinates of grid*/
//					for (int j = 0; j < prt.Dimension; j++) {
//						component[j] = (size_t)floor(prt.x[j] * this->num_side_voxels[j]);
//					}
//					size_t temp = this->Components2Index(component);
//
//#pragma omp critical /* if openmp 3.1 is avaiable, "#pragma omp atomic capture" may be used instead. */
//					{
//						this->prt_indicator[index] = temp;
//						index++;
//					}
//				}
//			);
//		}
//	}
//
//	/* erase all grid centers inside the matrix (reference) phase. */
//	this->prt_indicator.resize(index);
//	std::cout << "----------- generated the grid ---------\n";
//	std::cout << "Total number of voxels = " << grid.NumParticle() << "\n";
//	std::cout << "The number of voxels in the particle phase = " << this->prt_indicator.size() << "\n";
//	/* if there are overlapping particles, check whehter there are duplicate voxels. */
//	if (DoOverlap) {
//		std::cout << "Some duplicate voxels are removed \n";
//		std::cout << "The number of voxels in the particle phase = " << this->prt_indicator.size() << "\n";
//
//		std::sort(this->prt_indicator.begin(), this->prt_indicator.end());
//		for (auto iter = this->prt_indicator.begin(); iter != this->prt_indicator.end() - 1; iter++) {
//			if (*iter == *(iter + 1))
//				this->prt_indicator.erase(iter + 1);
//		}
//	}
//	std::cout << "The volume fraction of the particle phase is \n"
//		<< "\t" << this->packing.PackingFraction() << " (from SpherePacking object)\n"
//		<< "\t" << this->prt_indicator.size() / (double)grid.NumParticle() << " (from prt_indicator)\n";
//	update_indicator = true;
//}

void TwoPhaseDispersion::ComputeIndicator() {
	int dim = this->packing.GetDimension();
	/* basis vectors for the grid points (in relative coordinates)*/
	std::vector<double> basis(dim, 0);	
	size_t tot_num_grids = 1;
	for (int i = 0; i < this->packing.GetDimension(); i++) {
		basis[i] = 1.0 / (double)this->num_side_voxels[i];
		tot_num_grids *= this->num_side_voxels[i];
	}

	/* prt_indicator[i] = grid coordinates of the particle phase in a raw-major format. */
	this->prt_indicator.resize(tot_num_grids, 0); 	

	/* make IterateThroughNeighbors thread-safe */
	double Rmax = this->packing.GetMaxRadius();
	this->packing.PrepareIterateThroughNeighbors(Rmax);

	/* other preparations */
	size_t num_grids_in_part = 0;

	/* iterate over all grid points */
#pragma omp parallel for schedule(dynamic)
	for (size_t idx = 0; idx < tot_num_grids; idx++) {
		/* convert index of a grid point to the position of the grid */
		std::vector<size_t> comp_grid;
		this->Index2Components(idx, comp_grid);
		GeometryVector coord_grid(dim);
		for (int i = 0; i < dim; i++)
			coord_grid.x[i] = (0.5+(double)comp_grid[i]) * basis[i];

		/*If the grid point is inside a particle */
		if (this->packing.CheckOverlap(coord_grid, 0)) {
#pragma omp critical
			{
				this->prt_indicator[num_grids_in_part] = idx;
				num_grids_in_part++;
			}
		}
	}

	this->prt_indicator.resize(num_grids_in_part);
	std::cout << "The volume fraction of the particle phase is \n"
		<< "\t" << this->packing.PackingFraction() << " (from SpherePacking object)\n"
		<< "\t" << this->prt_indicator.size() / (double)tot_num_grids << " (from prt_indicator)\n";
	update_indicator = true;
}

double TwoPhaseDispersion::WeightedPhaseIndicator(size_t index) const
{	
	//TODO:: refine this function to account for the fraction of particles inside the cell!!

	if(use_weights){
		return this->weights[index];
	}
	else{
		return 1.0L;
	}
}

void TwoPhaseDispersion::Print(const std::string & prefix) const {
	std::fstream output(prefix + ".DisArray", std::fstream::out);
	int d = this->packing.GetDimension();
	output << "##DiscretizedArray##\n";
	output << "#dimension\n";
	output << packing.GetDimension() << "\n";
	output << "#resolutions\n";
	for (int i = 0; i < d; i++)
		output << this->num_side_voxels[i]<<"\t";
	output << "\n";
	output << "#num_prt_pixels\n";
	output << this->num_total_voxels <<"\n";
	output << "#basis";
	for (int i = 0; i < d; i++) {
		output << this->packing.GetBasisVector(i);
	}

	std::vector<size_t> comp;
	for (size_t i = 0; i < prt_indicator.size(); i++) {
		this->Index2Components(prt_indicator[i], comp);
		for (int j = 0; j < d; j++) {
			output << comp[j] << "\t";
		}
		output << this->WeightedPhaseIndicator(i)<<"\n";
	}
	output.close();
}

void TwoPhaseDispersion::Read(const std::string & prefix) {
	std::string name = prefix + std::string(".DisArray");
	std::fstream ifile (name.c_str(), std::fstream::in);

	if (ifile.good() == false){
		std::cerr << "TwoPhaseDispersion::Read : input stream is not good!\n" << name << "\n";
		return;
	}

	std::string dump;
	std::getline(ifile, dump);
	this->prt_indicator.clear();
	this->weights.clear();

	size_t num_prt_pixels = 0;
	DimensionType d=0;

	/* read some basic info. */
	if (strcmp(dump.c_str(), "##DiscretizedArray##") == 0){
		ifile >> dump;
		if(strcmp(dump.c_str(), "#dimension")==0){
			ifile >> d;
		}
		ifile >> dump;

		std::vector<size_t> num_pixels(d, 0);
		if(strcmp(dump.c_str(), "#resolutions")==0){
			for (int i=0; i<d; i++){
				ifile >> num_pixels[i];
			}
		}
		ifile >> dump;

		if(strcmp(dump.c_str(), "#num_prt_pixels")==0){
			ifile >> num_prt_pixels;
			this->prt_indicator.resize(num_prt_pixels);
			this->weights.resize(num_prt_pixels);
		}
		ifile >> dump;		

		if(strcmp(dump.c_str(), "#basis")==0){
			std::vector<GeometryVector> basis(d);
			for (int i = 0; i < d; i++) {
				for (int j = 0; j < d; j++) {
					ifile >> basis[i].x[j];
				}
				basis[i].SetDimension(d);
			}

			this->packing = SpherePacking(d, basis.data(), 1.0);
		}		

		this->SetResolutions(num_pixels);

	}

	{
		std::vector<size_t> size = this->GetResolutions();
		for (size_t i=0; i<d; i++)
			std::cout<<size[i]<<"\t";
		std::cout<<"\n";
	}

	/* read pixel coordinates and the corresponding weight. */
	std::istringstream iss;
	size_t id = 0;
	while (std::getline(ifile, dump)){
		std::vector<size_t> components(d, 0);
		double val = 0;

		iss.str(dump);
		iss.clear();

		for(int i=0; i<d; i++){
			iss >> components[i];
		}
		iss >> val;

		if(ifile.eof() == false && id < this->prt_indicator.size()){
			this->prt_indicator[id] = this->Components2Index(components);
			this->weights[id] = val;
			id++;
		}
	}
	this->use_weights = true;
	this->update_indicator = true;
};

TwoPhaseDispersion::TwoPhaseDispersion(const SpherePacking & packing, 
	size_t num_pixels): packing(packing){
	DimensionType d = this->packing.GetDimension();
	this->SetResolutions(std::vector<size_t>(d, num_pixels));
}

TensorField TwoPhaseDispersion::GetKs() const
{
	TensorField result(packing.GetDimension(), 1, GetResolutions());
	std::vector<double> L;
	for (int i = 0; i < packing.GetDimension(); i++) {
		L.push_back(sqrt(packing.GetBasisVector(i).Modulus2()));
	}

	for (size_t k_idx = 0; k_idx < result.num_total_pixels; k_idx++) {
		std::vector<size_t> comp;
		Index2Components(k_idx, comp);
		GeometryVector k(static_cast<int>(comp.size()));

		for (size_t i = 0; i < comp.size(); i++) {
			k.x[i] = (double)comp[i];
			if (2*comp[i] > result.num_points[i])
				k.x[i] -= (double)result.num_points[i];
			k.x[i] /= L[i];
			
			result[i][k_idx][0] = k.x[i];
		}
	}
	result *= 2.0*pi;
	//result.Print("testks");
	return result;
}