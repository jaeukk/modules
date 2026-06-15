#include "LinealPathFunction.h"

void ChordLengthDistribution(const std::function<const SpherePacking(size_t i)> & GetC, size_t num_c, std::vector<GeometryVector> & p0, std::vector<GeometryVector> & p1, double Rmax, double dR, size_t num_sample_lines, int seed)
{
	p0 = std::vector<GeometryVector>((int)std::round(Rmax/dR), GeometryVector(0,0,0,0));
	for (size_t i = 0; i < p0.size() ; i++){
		p0[i].x[0] = dR*i;
	}
	p1 = std::vector<GeometryVector> (p0);
	// In p0 and p1, x[3] is for storing distribution of each configuration temporarily

	RandomGenerator rng(seed);
	std::cout<<"Computing Chord-Length Distribution";
	progress_display pd(num_c);
	for (size_t i = 0; i < num_c; i++){
		SpherePacking c = GetC(i);
		c.SetCellSize(c.GetMaxRadius()*2.0);
		c.PrepareIterateThroughNeighbors(c.GetMaxRadius());
		std::vector<size_t> l_matrix, // l_matrix = line legnths in a matrix phase
			l_particle;	//l_particle = line lengths in a particle phase
		GetLinealSegments(c, l_matrix, l_particle, 2*Rmax, dR, num_sample_lines, rng);

		for (auto iter = l_matrix.begin(); iter != l_matrix.end(); iter++){
			if (*iter < p0.size())
				p0[*iter-1].x[3] ++; 
		}
		for (auto iter = l_particle.begin(); iter != l_particle.end(); iter++){
			if (*iter < p1.size())
				p1[*iter-1].x[3] ++; 
		}

		/* Processing to extract the chord-length distribution of each configuration*/
		double tot_num_l_matrix = 0.0, tot_num_l_particle = 0.0;
		#pragma omp parallel for reduction(+: tot_num_l_matrix, tot_num_l_particle)
		for (int j=0; j< (int)p0.size(); j++){
			tot_num_l_matrix += (double)p0[j].x[3];
			tot_num_l_particle += (double)p1[j].x[3];
		}
		
		double p0_ = 0.0, p1_ = 0.0;
		#pragma omp parallel for private(p0_, p1_)
		for (int j=0; j< (int)p0.size(); j++){
			p0_ = p0[j].x[3] / tot_num_l_matrix / dR;
			p1_ = p1[j].x[3] / tot_num_l_particle / dR;
			p0[j].x[1] += p0_;
			p0[j].x[2] += p0_*p0_;
			p1[j].x[1] += p1_;
			p1[j].x[2] += p1_*p1_;

			//initialization for the next configuration
			p0[j].x[3] = 0.0;	
			p1[j].x[3] = 0.0;			
		}
		pd++;
	}

	// summarize the data
	for (size_t i = 0; i < p0.size(); i++){
		p0[i].x[1] /= (double) num_c;
		p0[i].x[2] = std::sqrt( (p0[i].x[2] / num_c - p0[i].x[1] * p0[i].x[1]) / (num_c-1));

		p1[i].x[1] /= (double) num_c;
		p1[i].x[2] = std::sqrt( (p1[i].x[2] / num_c - p1[i].x[1] * p1[i].x[1]) / (num_c-1));
	}
	std::cout<<"done!\n";
}

void GetLinealSegments(const SpherePacking &c, std::vector<size_t> &lengths0, std::vector<size_t> &lengths1, double Lmax, double dL, size_t num_sample_lines, RandomGenerator &rng)
{
	lengths0.resize(0);
	lengths1.resize(0);
	lengths0.reserve(num_sample_lines * static_cast<size_t>( std::round(Lmax/dL/20)));
	lengths1.reserve(num_sample_lines * static_cast<size_t>( std::round(Lmax/dL/20)));

	int d = (int)c.GetDimension();
	int size = (int)std::round(Lmax/dL);

	for (size_t i = 0; i < num_sample_lines; i++){
		// points 1 and 2; u = p2-p1
		// p1 and p2 are in relative coordinates;
		// u is a vector of length dL, in relative coordinates.
		GeometryVector p1, p2, u;
		{
			GeometryVector temp1(static_cast<int>(d)), temp2(static_cast<int>(d));

			for (int j=0; j<d; j++){
				temp1.x[j] = rng.RandomDouble();
				temp2.x[j] = rng.RandomDouble_normal(1.0);
			}
			p1 = temp1;

			temp2 = 1./std::sqrt(temp2.Modulus2()) * temp2;
			u = c.CartesianCoord2RelativeCoord(temp2);
			p2 = p1 + Lmax * u;
			u = dL * u;

			// std::cout << d <<"\n";
			// std::cout << temp1;
			// std::cout << temp2;
		}

		// along the line segment betwee p1 and p2, check whether a point is inside the particle phase or not.
		std::vector<int> phases(size, 0);

		#pragma omp parallel for shared(phases)
		for(int j=0; j< size; j++){
			GeometryVector x = p1 + j*u;
			if(c.CheckOverlap(x,0.0)){
				/* Inside the particle phase */
				phases[j] = 1;
			}
			/* In a matrix phase */
			//phases[j] = 0;
		}

		// Find the continuous segmants
		int curr_phase = 0, prev_phase = phases[0];
		int leng = 1;
		for (size_t j = 1; j < phases.size(); j++){
			curr_phase = phases[j];
			bool is_not_at_end = (j < (phases.size()-1));
			if ((curr_phase == prev_phase) && is_not_at_end){
				// same phase and before end of the line
				leng ++ ;
			}
			else{
				// end of the line, or phase boundary
				if (!is_not_at_end){
					leng++;
				}

				if (prev_phase == 0){
					lengths0.push_back(leng);
				}
				else{
					lengths1.push_back(leng);
				}
				leng = 1;
			}

			if ((j==phases.size()-1) || (curr_phase != prev_phase)){
			}
			else{
				// same phase
			}
			prev_phase = curr_phase;
		}



		/* sanity check */
		if(false){
			int l1 = 0, l2 = 0;
			for (int j = 0; j<size; j++){
				std::cout << phases[j];
			}
			std::cout << "\n";

			std::cout << "phase 0:\n";
			for (auto iter = lengths0.begin(); iter != lengths0.end(); iter++){
				std::cout << *iter <<"\t";
				l1 += *(iter);
			}
			std::cout << "\n";
			std::cout << "phase 1:\n";
			for (auto iter = lengths1.begin(); iter != lengths1.end(); iter++){
				l2 += *(iter);
				std::cout << *iter <<"\t";
			}
			std::cout << "\n";

			std::cout << "p1=\t" << p1;
			std::cout << "u=\t" << u;


			std::cout << "total length = " << (int)phases.size() <<"\t";
			std::cout << "\t (lengths in phases 1 and 2) = " << (l1+l2) << "\n";
		}
	}
}

void LinealPathFunction(const std::function<const SpherePacking(size_t i)> &GetC, size_t num_c, std::vector<GeometryVector> &L0, std::vector<GeometryVector> &L1, double Rmax, double dR, size_t num_sample_lines, int seed)
{
	L0 = std::vector<GeometryVector>((int)std::round(Rmax/dR), GeometryVector(0,0,0,0));
	for (size_t i = 0; i < L0.size() ; i++){
		L0[i].x[0] = dR*i;
	}
	L1 = std::vector<GeometryVector> (L0);
	// In L0 and L1, x[3] is for storing distribution of each configuration temporarily

	double phi0 = 0.0, phi1 = 0.0;

	RandomGenerator rng(seed);
	std::cout<<"Computing Lineal-Path Function";
	progress_display pd(num_c);
	for (size_t i = 0; i < num_c; i++){
		SpherePacking c = GetC(i);
		c.SetCellSize(c.GetMaxRadius()*2.0);
		c.PrepareIterateThroughNeighbors(c.GetMaxRadius());
		std::vector<size_t> l_matrix, // l_matrix = line legnths in a matrix phase
			l_particle;	//l_particle = line lengths in a particle phase
		GetLinealSegments(c, l_matrix, l_particle, 2*Rmax, dR, num_sample_lines, rng);

		/* some basic statistics */
		double phi0_ = 0.0, phi1_ = 0.0;
		double lc0_ = 0.0, lc1_ = 0.0;
		for (size_t j = 0; j < l_matrix.size(); j++){
			phi0_ += l_matrix[j];
		}
		for (size_t j = 0; j < l_particle.size(); j++){
			phi1_ += l_particle[j];
		}
		lc0_ = phi0_ ; /// (double) l_matrix.size();
		lc1_ = phi1_ ; /// (double) l_particle.size();
		phi0_ = phi0_ / (phi0_+phi1_);
		phi1_ = 1.0 - phi0_;

		phi0 += phi0_;
		phi1 += phi1_;

		/* Extract the lineal-path function of each configuration */
		// in matrix phase 
		std::vector<double> sum(L0.size(),0.0), sum1(L0.size(),0.0); 
		for (auto iter = l_matrix.begin(); iter != l_matrix.end(); iter++){
			if (*iter < sum.size()){
				sum[*iter-1] += 1;
				sum1[*iter-1] += (*iter)-0.5;
			}
		}
		for (size_t j = 0; j < L0.size(); j++){
			/* sum in the reverse order */
			size_t idx = L0.size() - j - 1;
			if (j > 0){
				sum1[idx] += sum1[idx+1];
				sum[idx] += sum[idx+1];
			}

			double l0 = phi0_ *(sum1[idx] - (double)idx * sum[idx]) / lc0_;
			L0[idx].x[1] += l0;
			L0[idx].x[2] += l0*l0;
		}

		// in particle phase
		sum = std::vector<double>(L1.size(),0.0);
		sum1 = std::vector<double>(L1.size(),0.0); 
		for (auto iter = l_particle.begin(); iter != l_particle.end(); iter++){
			if (*iter < sum.size()){
				sum[*iter-1] += 1;
				sum1[*iter-1] += (*iter-0.5);
			}
		}
		for (size_t j = 0; j < L1.size(); j++){
			/* sum in the reverse order */
			size_t idx = L1.size() - j - 1;
			if (j > 0){
				sum1[idx] += sum1[idx+1];
				sum[idx] += sum[idx+1];
			}

			double l1 = phi1_ * (sum1[idx] - (double)idx * sum[idx]) / lc1_;
			L1[idx].x[1] += l1;
			L1[idx].x[2] += l1*l1;
		}
		pd++;
	}

	// summarize the data
	for (size_t i = 0; i < L0.size(); i++){
		L0[i].x[1] /= (double) num_c;
		L0[i].x[2] = std::sqrt( (L0[i].x[2] / num_c - L0[i].x[1] * L0[i].x[1]) / (num_c-1));

		L1[i].x[1] /= (double) num_c;
		L1[i].x[2] = std::sqrt( (L1[i].x[2] / num_c - L1[i].x[1] * L1[i].x[1]) / (num_c-1));
	}
	std::cout<<"done!\n";
	phi0 /= num_c;
	phi1 /= num_c;
	std::cout << "phi of the matrix phase = "<< phi0 << "\n"; 
	std::cout << "phi of the particle phase = "<< phi1 << "\n";

}
