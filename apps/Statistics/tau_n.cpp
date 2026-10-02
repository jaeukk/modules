#include "tau_n.h"

/* if  */
double tau_n::gn(const Configuration &c, const std::vector<GeometryVector> &xs, size_t prt_idx, std::vector<int> overlap){
	assert(xs.size() == this->n-1);
	assert(overlap.size() == this->n-1);

	/* the effective number of particles around x_{i-1} */
	std::vector<double> ns(this->n,0);

	/* at prt_idx */
	std::vector<int> prts;
	GeometryVector x1 = c.GetRelativeCoordinates(prt_idx);
	//c.GetNeighborList(x1, this->dR, prts);
	c.IterateThroughNeighbors(x1, this->dR, [&prts](const GeometryVector &, const GeometryVector &, const signed long *, const size_t prt_index) -> void { prts.emplace_back((int) prt_index); });
	ns[0] = (double)prts.size();

	for (size_t i=0; i < xs.size(); i++){
		GeometryVector xi = x1 + xs[i];
		if (overlap[i] < 0){
			// no overlap
			//c.GetNeighborList(xi, this->dR, prts);
			prts.resize(0);
			c.IterateThroughNeighbors(xi, this->dR, [&prts](const GeometryVector &, const GeometryVector &, const signed long *, const size_t prt_index) -> void { prts.emplace_back((int) prt_index); });
			
			ns[i+1] = (double) prts.size();
		}
		else {
			// if x_i is overlapped with one of previous points,
			// the effective number of points will be (number of particles near the overlapped point) - 1.
			ns[i+1] = ns[overlap[i]] - 1.;
		}

		if (ns[i+1] < 0.5){
			// the number of particles around points are too small.
			return 0.0;
		}
	}

	/* compute */
	double val = 1.0;
	for (size_t i=1; i<ns.size(); i++){
		val *= ns[i];
	} 
	return val;
}

double tau_n::gn_from_nlist(const nlist &nlist_, const std::vector<int> &idx_relative_pos, const std::vector<int> overlap)
{
	assert(idx_relative_pos.size() == this->n-1);
	assert(overlap.size() == this->n-1);

	double result = 0;
	double N = 0;
	for (auto iter = nlist_.bins.begin(); iter != nlist_.bins.end(); ++iter){
		std::vector<double> ns(this->n, 0.0);

		ns[0] = iter->second;
		for (size_t i = 0; i<idx_relative_pos.size(); i++){
			int idx_pos = idx_relative_pos[i] + iter->first;

			if (overlap[i] < 0){
				// no overlap
				ns[i+1] = (double)nlist_.get(idx_pos);
			}
			else{
				// if idx_pos is overlapped with one of previous points,
				ns[i+1] = ns[overlap[i]]-1.0;
			}
		}
		double val = ns[0];
		for (size_t i = 1; i < ns.size(); i++){
			val *= ns[i];
		}
		result += val;
		N += ns[0];
	}

	return result/N;
}

tau_n::tau_n(std::ostream &ofile, std::istream &ifile)
{
	ofile << "order n = ";
	ifile >> this->n;
	if (this->n < 2 or this->n > 4){
		ofile << "n is defined only at 2, 3, or 4";
		return ;
	}

	ofile << "ConfigPack name = ";
	ifile >> this->cp_name;
	this->CP = ConfigurationPack(this->cp_name);
	ofile << "Nc = ";
	ifile >> this->Nc;
	if (this->Nc > (size_t)this->CP.NumConfig()){
		this->Nc = this->CP.NumConfig();
		ofile << "$Nc is too large, so it is replaced with ";  
	}
	ofile << this->Nc << "\n";
	
	ofile << "num. threads = ";
	ifile >> this->num_threads;

	ofile << "output prefix = ";
	ifile >> this->savename;

	ofile << "prefix of binary file = ";
	ifile >> this->bin_name;
	bin_name += ".tau_"+std::to_string(this->n);
	ofile << bin_name << std::endl;

	/** If the binary file exists, then load the parameters. (next two inputs are ignored.) */
	if (this->ReadBinary()){
		double temp;
		ifile >> temp;	ofile << temp <<"Rmax is ignored\n";
		ifile >> temp;	ofile << temp <<"dR is ignored\n";
	}
	/* If the binary file doesn't exist, then input parameters here. */
	else{
		ofile << "Start new calculations\n";
		this->curr_idx = 0;
		ofile << "Rmax = ";
		ifile >> this->Lmax;	ofile << this->Lmax << std::endl;

		ofile << "dx = ";		
		ifile >> this->dx;	ofile << this->dx << std::endl;
		this->dR = 0.5*this->dx;

		this->normalization_factor = 1.0/pow(this->dx, (double)(this->n - 1));
		this->num_half = (int) floor(this->Lmax/this->dx);
		this->len_row = 2*this->num_half + 1 ;
		this->tot_num_tuples = (this->num_half+1)*pow(this->len_row, n-2);

		this->result.resize(this->num_half + 1, GeometryVector(0.,0.,0.,0.));
	}
	{
		ofile <<"\n\n\n\n\n\n";
		ofile << "num_half = "<< this->num_half <<"\n";
		ofile << "len_row = "<<this->len_row <<"\n";
		ofile << "normalization factor = "<<this->normalization_factor << "\n";
		ofile << "Lmax = " << this->Lmax << "\n";
		ofile << "dx = " << this->dx << "\n";
		ofile << "len. val = "<< this->result.size() <<"\n";


	}


}

void tau_n::WriteBinary() const
{
	std::fstream ofile(this->bin_name.c_str(), std::fstream::out);

	ofile.write((char *)(&this->curr_idx), sizeof(this->curr_idx));
	ofile.write((char *)(&this->num_half), sizeof(this->num_half));
	ofile.write((char *)(&this->n), sizeof(this->n));
	ofile.write((char *)(&this->dx), sizeof(this->dx));
	ofile.write((char *)(&this->Lmax), sizeof(this->Lmax));
	ofile.write((char *)(&this->Nc), sizeof(this->Nc));
	
	
	for (size_t i = 0; i < this->result.size(); i++) {
		for (size_t j=0; j < 2; j++)
			ofile.write((char *)(&this->result[i].x[j]), sizeof(this->result[i].x[j]));
	}
	ofile.close();
}

bool tau_n::ReadBinary()
{
	std::fstream ifile(this->bin_name.c_str(), std::fstream::in);
	if (ifile.good()){
		std::cout << "Load progress from "<< bin_name << "\n";

		ifile.read((char *)(&this->curr_idx), sizeof(this->curr_idx));
		ifile.read((char *)(&this->num_half), sizeof(this->num_half));
		ifile.read((char *)(&this->n), sizeof(this->n));
		ifile.read((char *)(&this->dx), sizeof(this->dx));
		ifile.read((char *)(&this->Lmax), sizeof(this->Lmax));
		ifile.read((char *)(&this->Nc), sizeof(this->Nc));
		
		this->len_row = 2*this->num_half+1;
		this->normalization_factor = 1.0/pow(this->dx, (double)(this->n - 1));
		this->tot_num_tuples = (this->num_half+1) * std::pow(this->len_row, this->n-2);
		this->dR = 0.5*this->dx;

		this->result.resize(this->num_half + 1, GeometryVector(0,0,0,0));
		
		for (size_t i = 0; i < this->result.size(); i++) {
			for (size_t j=0; j < 2; j++)
				ifile.read((char *)(&this->result[i].x[j]), sizeof(this->result[i].x[j]));
		}
		ifile.close();
		return true;
	}
	else{
		std::cout << "Binary file doesn't exist.\n";
		ifile.close();
		return false;
	}
}

void tau_n::Write_tau_n()
{
	char filename[300];
	sprintf(filename, "%s-tau%d", this->savename.c_str(), (int)this->n);

	std::vector<std::string> comments;
	comments.push_back("direct calculations of tau_n from configurations");
	comments.push_back("source = " + this->cp_name + "\n");
	comments.push_back("Nc = " + std::to_string(this->Nc) + "\n");
	comments.push_back("x = L, y = tau_n, dy = err. in y");
	comments.push_back("x\t y\t dy");

	std::vector<GeometryVector> table(this->result.size(), GeometryVector(0,0,0));
	double dv = pow(this->dx, this->n-1);
	/* cumulative sums */
	table[0] = GeometryVector(this->result[0]);
	for (size_t i = 1; i < table.size(); i++){
		table[i].x[1] = 2.*this->result[i].x[1] + table[i-1].x[1];
		table[i].x[2] = 2.*this->result[i].x[2] + table[i-1].x[2];
	}
	/* multiply volume element */
	for (size_t i = 0; i < table.size(); i++){
		table[i].x[0] = this->result[i].x[0];
		table[i].x[1] *= dv;
		table[i].x[2] *= dv;
	}

	WriteFunction(table, filename, comments);
}

/** \brief Order by absolute value (sorts bin offsets so the last one is the largest |offset|).
 *  \param i  First offset.
 *  \param j  Second offset.
 *  \return |i| < |j|. */
inline bool comp(int i, int j){ return std::abs(i) < std::abs(j); }	


void tau_n::Run(std::ostream &ofile, time_t TimeLimit, bool output_gn){
	omp_set_num_threads(this->num_threads);
	/* Generate hash tables mapping the index of bin to the number of particles therein. */
	int num_bins = (int) round( std::sqrt(this->CP.GetConfig(0).GetBasisVector(0).Modulus2()) / this->dx );
	std::vector<nlist> arr_nlist(this->Nc, nlist(num_bins));
	double factor = 0.0;
	//ofile << "tau_n::Run(...) is initializing nlist\n";
	#pragma omp parallel for 
	for (size_t i=0; i < this->Nc; i++){
		Configuration c = this->CP.GetConfig(i);
		if (i == 0){
			factor = this->normalization_factor / this->Nc;
		}
		/* Get nlist */
		std::vector<int> xs (c.NumParticle());
		for (size_t j = 0; j < c.NumParticle(); j++){
			int idx = (int) floor(c.GetCartesianCoordinates(j).x[0]/this->dx);
			arr_nlist[i].update(idx, arr_nlist[i].get(idx)+1);
		}
	} 

	std::vector<std::set<int>> pair_dists(this->Nc);
	if (this->n > 2){
		/* if n = 3 or higher, we compute a set of unique pair distances in each configuration, which is used to skip some calculations of g_n later.*/
	#pragma omp parallel for 
		for (size_t idx = 0; idx < arr_nlist.size(); idx ++){
			pair_dists[idx].insert(0);
			for (auto iter1 = arr_nlist[idx].bins.begin(); iter1 != arr_nlist[idx].bins.end(); ++ iter1){
				for (auto iter2 = arr_nlist[idx].bins.begin(); iter2 != arr_nlist[idx].bins.end(); ++ iter2){
					int dist = iter1->first - iter2->first;
					if (dist > 0){
						dist = std::min(dist, arr_nlist[idx].period - dist);
						if (dist < this->num_half + 2){
							pair_dists[idx].insert(dist);
						}
					}
				}
			}
		}
	}
	ofile << "tau_n::Run(...) completed initialization of nlist\n";

	/* output stream for saving gn function.*/
	std::fstream out_gn((this->savename + "_g" + std::to_string(this->n) +".txt").c_str(), std::fstream::out);
	if (output_gn){
		for (size_t i =0; i < this->n-1; i++){
			out_gn << "#x_1"<< std::to_string(i+2) <<"\t";
		}
		out_gn << "#gn\t#err_gn\n";
	}
	ofile << "Start from " << this->curr_idx <<"th tuples of points" <<std::endl;
	ofile << "Computing tau_"<<this->n<<std::endl;
	progress_display pd (this->tot_num_tuples);
	pd += this->curr_idx;


	size_t idx_positions = this->curr_idx;
	if(idx_positions == 0){
		/* Initialize results */
		for (size_t i=0; i < this->result.size(); i++){
			this->result[i].x[0] = this->dx * (double) (i+1); 
		}
	}

	/* variables for time estimations */
	auto start_loop = std::time(nullptr);
	time_t exp_time = 0;
	std::vector<int> prev_relative_pos;

	/* iterate over positions */
	for (; idx_positions < this->tot_num_tuples; idx_positions ++ ){
		// convert index to tuples of bin indices
		std::vector<int> curr_relative_pos;
		this->idx2pos((int)idx_positions, curr_relative_pos);
		// sort "absolute values" of relative positions
		std::sort(curr_relative_pos.begin(), curr_relative_pos.end(), comp);
		int idx_L = std::abs(curr_relative_pos[curr_relative_pos.size()-1]);	// get L

		// compare wether the same 


		// get L 
		//int idx_L = std::abs(* std::max_element(curr_relative_pos.begin(), curr_relative_pos.end(), comp));
		
		/** check whether points are overlapped. 
		 * if (idx_overlapped[i] == -1), no overlap
		 * elif (idx_overlapped[i] <= i), xs[i] is overlapped with point (idx_overlapped[i]).
		 * 
		 * If x0, x1 (=xs[0]), and x3 (=xs[2]) are overlapped,    
		 * idx_overlapped = {0,-1, 1}.
		 * 
		 * If x0 and x3 (=xs[2]) are overlapped, and x1 and x2 are overlapped, then 
		 * idx_overlapped = {-1,1, 0}.
		*/
		std::vector<int> idx_overlapped(this->n-1,-1);
		{
			for (size_t i=0; i < curr_relative_pos.size(); i++){
				if (curr_relative_pos[i] == 0){
					idx_overlapped[i] = 0;
				}			
				for (int j = (int)i-1; j >= 0; j--){
					if (curr_relative_pos[j] == curr_relative_pos[i]){
						idx_overlapped[i] = j+1;
						break;
					} 
				}
			}
		}

		double gn = 0.0, err_gn = 0.0;
#pragma omp parallel for reduction(+: gn, err_gn) 
		/* In each configuration, */
		for (size_t i = 0; i < this->Nc; i++){
			bool NeedToCompute = true;
			if ( this->n > 2){
			/* check whethter all elements in curr_relative_pos are possible pair distances. */
				for (size_t j = 0; j < curr_relative_pos.size(); j++){
					NeedToCompute &= (pair_dists[i].find(abs(curr_relative_pos[j])) != pair_dists[i].end());
				}
			}

			if (NeedToCompute){
				/* gn of a single configuration */
				double num_tuples = this->gn_from_nlist(arr_nlist[i], curr_relative_pos, idx_overlapped);
				gn += num_tuples;
				err_gn += num_tuples*num_tuples;
			}
		}

		gn *= factor;
		err_gn *= factor*factor;
		err_gn = std::sqrt(abs(err_gn / (double)this->Nc - gn*gn) / (this->Nc - 1.0));

		if (output_gn){
			for (size_t i =0; i < this->n-1; i++){
				out_gn << this->dx*curr_relative_pos[i] << "\t";
			}
			out_gn << std::setprecision(8);
			out_gn << gn << "\t" <<err_gn <<"\n";
		}

		double y = 1.-gn;
		this->result.at(idx_L).x[1] += y*y;
		this->result.at(idx_L).x[2] += 2.*y*err_gn;

		/* Time estimations */
		{
			auto curr_time = std::time(nullptr);
			if (exp_time == 0){
				exp_time = std::max(curr_time - start_loop, (long int)1);
	//			exp_time = (curr_time - start_loop > 1)? curr_time - start_loop:1;
			}

			/** If the program doesn't have time to add three more configurations, it save current progress and quite the program */
			if (TimeLimit - curr_time < 3*exp_time){
				ofile << "time over"<<std::endl;
				break;
			}
			pd ++ ;		
		}
	}

	out_gn.close();
	ofile << "\nComputation Done!\n";
	ofile << "It took " << std::time(nullptr) - start_loop << " seconds\n";
	ofile << "Saving progress ...\n";
	{
		this->curr_idx = idx_positions;
		this->WriteBinary();
	}
	ofile << "Progress is safely saved. \n\n";

	ofile << "Writing tau_n ... \n";
	this->Write_tau_n();
	ofile << "Done\n";

}

/* obsolete version of tau_n::Run(...) */
// void tau_n::Run(std::ostream &ofile, time_t TimeLimit, bool output_gn)
// {
// 	/* Load Configurations */
// 	std::vector<Configuration> c_list(this->Nc, GetUnitCubicBox(1));
// 	for (int i=0; i < this->Nc; i++){
// 		c_list[i] = Configuration(this->CP.GetConfig(i));
// 	}
// 	double factor = this->normalization_factor / ((double)c_list[0].NumParticle() *this->Nc);
	
// 	/* output stream for saving gn function.*/
// 	std::fstream out_gn((this->savename + "_g" + std::to_string(this->n) +".txt").c_str(), std::fstream::out);
// 	if (output_gn){
// 		for (int i =0; i < this->n-1; i++){
// 			out_gn << "#x_1"<< std::to_string(i+2) <<"\t";
// 		}
// 		out_gn << "#gn\t#err_gn\n";
// 	}
	
// 	ofile << "Start from " << this->curr_idx <<"th tuples of points" <<std::endl;
// 	ofile << "Computing tau_"<<this->n<<std::endl;
// 	progress_display pd (this->tot_num_tuples);
// 	pd += this->curr_idx;


// 	size_t idx_positions = this->curr_idx;
// 	if(idx_positions == 0){
// 		/* Initialize results */
// 		for (int i=0; i < this->result.size(); i++){
// 			this->result[i].x[0] = this->dx * (double) (i+1); 
// 		}
// 	}
	
// 	/* variables for time estimations */
// 	auto start_loop = std::time(nullptr);
// 	time_t exp_time = 0;

// 	omp_set_num_threads(this->num_threads);

// 	/* iterate over positions */
// 	for (; idx_positions < this->tot_num_tuples; idx_positions ++ ){
// 		// convert index to tuples 
// 		std::vector<int> relative_pos;
// 		this->idx2pos((int)idx_positions, relative_pos);
			
// 		// get L 
// 		int idx_L = std::abs(* std::max_element(relative_pos.begin(), pos.end(), comp));

// 		/* Positions relative to one particle are written in the relative coordinates. */
// 		std::vector<GeometryVector> xs(this->n-1, GeometryVector(0.0000));
// 		/** check whether points are overlapped. 
// 		 * if (idx_overlapped[i] == -1), no overlap
// 		 * elif (idx_overlapped[i] <= i), xs[i] is overlapped with point (idx_overlapped[i]).
// 		 * 
// 		 * If x0, x1 (=xs[0]), and x3 (=xs[2]) are overlapped,    
// 		 * idx_overlapped = {0,-1, 1}.
// 		 * 
// 		 * If x0 and x3 (=xs[2]) are overlapped, and x1 and x2 are overlapped, then 
// 		 * idx_overlapped = {-1,1, 0}.
// 		*/
// 		std::vector<int> idx_overlapped(this->n-1,-1);
// 		{
// 			for (int i=0; i<pos.size(); i++){
// 				xs[i] = c_list[0].CartesianCoord2RelativeCoord(
// 					GeometryVector(this->dx * pos[i]) );

				
// 				if (pos[i] == 0){
// 				idx_overlapped[i] = 0;
// 				}			
// 				for (int j = (int)i-1; j >= 0; j--){
// 					if (pos[j] == pos[i]){
// 						idx_overlapped[i] = j+1;
// 						break;
// 					} 
// 				}
// 			}
// 		}
	
// 		double gn = 0.0, err_gn = 0.0;
// 	#pragma omp parallel for reduction(+: gn, err_gn) 
// 		/* In each configuration, */
// 		for (int i=0; i<this->Nc; i++){
// 			/* gn of a single configuration */
// 			double num_tuples = 0.0;
// 			/* iterate over centers */
// 			for (size_t j = 0; j < c_list[i].NumParticle(); j++){
// 				double y = this->gn(c_list[i], xs, j, idx_overlapped);
// 				num_tuples += y; 
// 			}
// 			//gn_temp *= 
// 			gn += num_tuples;
// 			err_gn += num_tuples*num_tuples;
// 		}

// 		gn *= factor;
// 		err_gn *= factor*factor;
// 		err_gn = std::sqrt(abs(err_gn / (double)this->Nc - gn*gn) / (this->Nc - 1.0));

// 		if (output_gn){
// 			for (int i =0; i < this->n-1; i++){
// 				out_gn << this->dx*relative_pos[i] << "\t";
// 			}
// 			out_gn << std::setprecision(8);
// 			out_gn << gn << "\t" <<err_gn <<"\n";
// 		}

// 		double y = 1.-gn;
// 		this->result.at(idx_L).x[1] += y*y;
// 		this->result.at(idx_L).x[2] += 2.*y*err_gn;

// 		/* Time estimations */
// 		{
// 			auto curr_time = std::time(nullptr);
// 			if (exp_time == 0){
// 				exp_time = std::max(curr_time - start_loop, (long int)1);
// 	//			exp_time = (curr_time - start_loop > 1)? curr_time - start_loop:1;
// 			}

// 			/** If the program doesn't have time to add three more configurations, it save current progress and quite the program */
// 			if (TimeLimit - curr_time < 3*exp_time){
// 				ofile << "time over"<<std::endl;
// 				break;
// 			}
// 			pd ++ ;		
// 		}
// 	}

// 	out_gn.close();
// 	ofile << "\nComputation Done!\n";

// 	ofile << "Saving progress ...\n";
// 	{
// 		this->curr_idx = idx_positions;
// 		this->WriteBinary();
// 	}
// 	ofile << "Progress is safely saved. \n\n";

// 	ofile << "Writing tau_n ... \n";
// 	this->Write_tau_n();
// 	ofile << "Done\n";

// }

void tau_n::Test1(std::ostream &ofile)
{

	/* For debugging of idx <-> positions*/
	ofile << "Sanity check of idx <-> positions \n";
	RandomGenerator rng(std::time(nullptr));
	int count = 0;
	for (size_t i = 0; i < this->tot_num_tuples; i ++){
		std::vector<int> pos;
		int idx1 = i; //(int) std::round(rng.RandomDouble() * model.tot_num_tuples);

		this->idx2pos(idx1, pos);

		int idx2 = this->pos2idx(pos);

		bool out_of_range = false;
		out_of_range |= pos[0] < 0; 
		for (auto iter = pos.begin(); iter != pos.end(); iter++)
			out_of_range |= abs(*iter) > this->num_half;

		if ( out_of_range ||( idx1 != idx2)){
			count++;
			ofile << idx1 << "\t" << idx2 << ":\t";
			for (auto iter = pos.begin(); iter != pos.end(); iter++) ofile << *iter <<",\t";

			ofile << std::endl;
		}

	}
	ofile << "trial = "<< this->tot_num_tuples <<", error = " << count << "\n";
}

void tau_n::Test2(std::ostream &ofile)
{

	for (size_t idx_positions = 0; idx_positions < tot_num_tuples; idx_positions ++ ){
		// convert index to tuples 
		std::vector<int> pos;
		this->idx2pos((int)idx_positions, pos);
			
		/** check whether points are overlapped. 
		 * if (idx_overlapped[i] == -1), no overlap
		 * elif (idx_overlapped[i] <= i), xs[i] is overlapped with point (idx_overlapped[i]).
		 * 
		 * If x0, x1 (=xs[0]), and x3 (=xs[2]) are overlapped,    
		 * idx_overlapped = {0,-1, 1}.
		 * 
		 * If x0 and x3 (=xs[2]) are overlapped, and x1 and x2 are overlapped, then 
		 * idx_overlapped = {-1,1, 0}.
		*/
		
		std::vector<int> idx_overlapped(this->n-1,-1);
		for (size_t i = 0; i< idx_overlapped.size(); i++){
			if (pos[i] == 0){
				idx_overlapped[i] = 0;
			}			
			for (int j = (int)i-1; j >= 0; j--){
				if (pos[j] == pos[i]){
					idx_overlapped[i] = j+1;
					break;
				} 
			}
		}


		/* Debug overlapps */
		ofile << 0 <<"\t";
		for (size_t i = 0; i < pos.size(); i++){
			ofile << pos[i] << "\t";	
		}
		ofile << std::endl;
		
		ofile << "\t";
		for (size_t i = 0; i < idx_overlapped.size(); i++){
			ofile << idx_overlapped[i] << "\t";	
		}
		ofile << std::endl << std::endl << std::endl;
	}
}
