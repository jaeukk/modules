#include "gn.h"

#include <algorithm>
#include <stdexcept>


gn::gn(std::ostream & ofile, std::istream & ifile){
	ofile << "order n = ";
	ifile >> this->n;
	if (this->n < 2 or this->n > 4){
		throw std::invalid_argument("gn: n is defined only at 2, 3, or 4");
	}

	ofile << "binary file name = ";
	ifile >> this->bin_name;
	bin_name += ".g"+std::to_string(this->n);
	ofile << bin_name << std::endl;
	std::fstream infile(bin_name.c_str(), std::fstream::in | std::fstream::binary);
	/** If the binary file exists, then load the parameters. */
	if (infile.good()){
		this->ReadBinary(infile);
		infile.close();
		double temp;
		ifile >> temp;	ofile << temp <<" is ignored\n";
		ifile >> temp;	ofile << temp <<" is ignored\n";
	}
	/* If the binary file doesn't exist, then input parameters here. */
	else{
		this->curr_c_idx = 0;

		ofile << "Rmax = ";
		ifile >> this->Rmax;	ofile << this->Rmax << std::endl;
		ofile << "dR = ";
		ifile >> this->dx;	ofile << this->dx << std::endl;

		this->num_half = (int) floor(this->Rmax/this->dx);
		this->len_row = 2*this->num_half + 1 ;

		size_t total_length = (this->num_half+1)*pow(this->len_row, n-2);
		this->sum1.assign(total_length, 0.0);
		this->sum2.assign(total_length, 0.0);
	}
	ofile << "\n\n\n\n\n\n";
	ofile << "num_half = "<< this->num_half <<"\n";
	ofile << "len_row = "<< this->len_row <<"\n";
	ofile << "Rmax = " << this->Rmax << "\n";
	ofile << "dR = " << this->dx << "\n";
	ofile << "num. bins = "<< this->sum1.size() <<"\n";
	ofile << "configurations so far = "<< this->curr_c_idx <<"\n";

	/** Set the number of threads */
	size_t num_threads;
	ofile << "Num Threads = ";
	ifile >> num_threads;	ofile << num_threads << "\n";
	this->SetNumThreads(num_threads);

	std::string name;
	ofile << "Input ConfigPack file = ";
	ifile >> name;
	this->CP = ConfigurationPack(name);
	ofile << "\nNum Configs = ";
	ifile >> this->targ_c_idx;
	if(this->targ_c_idx > (size_t)this->CP.NumConfig()){
		this->targ_c_idx = this->CP.NumConfig();
		ofile <<"This ConfigurationPack contains only " << this->targ_c_idx << " configurations, which is set to be the num_config\n";
	}


	ofile << "Output file name = ";
	ifile >> this->savename;	ofile << this->savename << "\n";
}


void gn::WriteBinary(std::ostream & ofile ) const {
	ofile.write((char *)(&this->curr_c_idx), sizeof(this->curr_c_idx));
	ofile.write((char *)(&this->num_half), sizeof(this->num_half));
	ofile.write((char *)(&this->n), sizeof(this->n));
	ofile.write((char *)(&this->dx), sizeof(this->dx));
	ofile.write((char *)(&this->Rmax), sizeof(this->Rmax));

	// Legacy layout: per bin, the ensemble means of g_n and g_n^2.
	const double c = (this->curr_c_idx > 0) ? (double)this->curr_c_idx : 1.0;
	for (size_t i = 0; i < this->sum1.size(); i++) {
		double m1 = this->sum1[i] / c, m2 = this->sum2[i] / c;
		ofile.write((char *)(&m1), sizeof(m1));
		ofile.write((char *)(&m2), sizeof(m2));
	}
};

void gn::ReadBinary(std::istream & ifile ){
	ifile.read((char *)(&this->curr_c_idx), sizeof(this->curr_c_idx));
	ifile.read((char *)(&this->num_half), sizeof(this->num_half));
	ifile.read((char *)(&this->n), sizeof(this->n));
	ifile.read((char *)(&this->dx), sizeof(this->dx));
	ifile.read((char *)(&this->Rmax), sizeof(this->Rmax));
	this->len_row = 2*this->num_half+1;

	size_t total_length = (this->num_half+1) * std::pow(this->len_row, n-2);
	this->sum1.assign(total_length, 0.);
	this->sum2.assign(total_length, 0.);

	for (size_t i = 0; i < this->sum1.size(); i++) {
		double m1 = 0.0, m2 = 0.0;
		ifile.read((char *)(&m1), sizeof(m1));
		ifile.read((char *)(&m2), sizeof(m2));
		this->sum1[i] = m1 * this->curr_c_idx;
		this->sum2[i] = m2 * this->curr_c_idx;
	}
}

void gn::Moments(size_t i, double & mean, double & err) const {
	const double c = (double)this->curr_c_idx;
	mean = (c > 0) ? this->sum1[i] / c : 0.0;
	err = (c > 1) ? std::sqrt(std::max(0.0, this->sum2[i] / c - mean * mean) / (c - 1.0)) : 0.0;
}

void gn::Run(std::ostream &ofile, time_t TimeLimit)
{
	ofile << "Start from " << this->curr_c_idx <<"th configuration " <<std::endl;
	ofile << "Computing g"<<this->n<<std::endl;
	progress_display pd (this->targ_c_idx);
	for (size_t i = 0 ; i<this->curr_c_idx; i++)
		pd ++ ;

	auto start_loop = std::time(nullptr);
	time_t exp_time = 0;

	for (size_t i = this->curr_c_idx; i < this->targ_c_idx; i++){
		this->AddConfig(this->CP.GetConfig(i));
		auto curr_time = std::time(nullptr);
		if (exp_time == 0){
			exp_time = (curr_time - start_loop > 1)? curr_time - start_loop:1;
		}

		/** If the program doesn't have time to add three more configurations, it saves the current progress and quits. */
		if (TimeLimit - curr_time < 3*exp_time){
			ofile << "time over"<<std::endl;
			break;
		}
		pd ++ ;
	}
	ofile << "expected time / config = "<< exp_time <<std::endl;
	ofile << "\n\n\nsaving progress" <<std::endl;
	std::fstream outfile(this->bin_name.c_str(), std::fstream::out | std::fstream::binary);
	this->WriteBinary(outfile);
	outfile.close();

	ofile << "\nProgress is saved safely"<<std::endl;
	// Legacy tested the loop counter, which reported an unfinished run when the last configuration hit the time limit.
	if (this->curr_c_idx >= this->targ_c_idx){
		if (this->print_gn){
			ofile << "Write gn " << std::endl;
			this->Write_gn(this->savename);
		}
		ofile << "Computation is done\n";

	}
	else{
		ofile << "Computation requires more time\n";
		ofile << "Continue ... " <<std::endl;
	}
	ofile << "\n\n";
	ofile << "Compute tau_n\n";
	std::vector<GeometryVector> result;
	this->Compute_tau_n(result);
	{
		std::string filename = this->savename + "_tau" + std::to_string(this->n);
		std::vector<std::string> comm = {
			"Up to " + std::to_string(this->curr_c_idx) + "th configuration",
			"R\t#tau_n\t#dtau_n"
			};
		WriteFunction(result, filename, comm);
	}


	ofile << "\n\n";
	ofile << "Compute tau_n for long-range correlations \n";
	this->Compute_tau_n_excluded(result);
	{
		std::string filename = this->savename + "_tau" + std::to_string(this->n) + "_exc";
		std::vector<std::string> comm = {
			"Up to " + std::to_string(this->curr_c_idx) + "th configuration",
			"Rmin = " + std::to_string(this->Rmin),
			"R\t#tau_n_exc\t#dtau_n_exc"
			};
		WriteFunction(result, filename, comm);
	}

}

void gn::Compute_tau_n(std::vector<GeometryVector> & result)
{
	result.clear();
	result.resize(this->num_half, GeometryVector(0.,0.,0.));
	double dv = pow(this->dx, this->n-1) ;
	for (int i=0; i<this->num_half; i++){
		result[i].x[0] = this->dx * (i+1);
	}
	int nh = this->num_half, lr = (int)this->len_row;
	if (this->n == 2){
#pragma omp parallel for
		for (int i2 = 0; i2 < this->num_half; i2++){
			double m, dy;
			this->Moments(i2, m, dy);
			double y = std::abs(m - 1.0);
			result[i2].x[1] += 2*y*y;
			result[i2].x[2] += 4.*y*dy;
		}
	}
	else if (this->n==3){

# pragma omp parallel for shared(result)
		for (int i2 = 0; i2 < nh + 1; i2 ++){
			for (int i3 = -nh; i3 < nh+1 ; i3 ++){
				int ri = std::max(i2, std::abs(i3));
				if (ri < nh){
					double m, dy;
					this->Moments((size_t)i2*lr + (i3+nh), m, dy);
					double y = std::abs(m - 1.0);
# pragma omp atomic
					result[ri].x[1] += y*y;
# pragma omp atomic
					result[ri].x[2] += 2.*y*dy;
				}
			}
		}
	}
	else if (this->n==4){

# pragma omp parallel for shared(result)
		for (int i2 = 0; i2 < nh + 1; i2 ++){
			for (int i3 = -nh; i3 < nh+1 ; i3 ++){
				int rj = std::max(i2, std::abs(i3));
				for (int i4 = -nh; i4 < nh+1; i4 ++ ){
					int ri = std::max(rj, std::abs(i4));

					if (ri < nh){
						double m, dy;
						this->Moments(((size_t)i2*lr + (i3+nh))*lr + (i4+nh), m, dy);
						double y = std::abs(m - 1.0);
	# pragma omp atomic
						result[ri].x[1] += y*y;
	# pragma omp atomic
						result[ri].x[2] += 2.*y*dy;
					}
				}
			}
		}
	}

	/** Cumulative integration */
	result[0].x[1] *= dv;
	result[0].x[2] *= dv;
	for (int i = 1; i < this -> num_half; i++){
		result[i].x[1] *= dv;
		result[i].x[2] *= dv;

		result[i].x[1] += result[i-1].x[1];
		result[i].x[2] += result[i-1].x[2];
	}
}

void gn::Compute_tau_n_excluded(std::vector<GeometryVector> & result)
{
	int nh = this->num_half, lr = (int)this->len_row;
	result.clear();
	result.resize(nh, GeometryVector(0.,0.,0.));
	double dv = pow(this->dx, this->n-1) ;
	int i_ = 0;	// indicates the lowest index to be added.
	for (int i=0; i<nh; i++){
		result[i].x[0] = this->dx * (i+1);
		if (result[i].x[0] < this->Rmin){
			i_ = i;
		}
	}

	if (this->n == 2){
#pragma omp parallel for
		for (int i2 = i_; i2 < nh; i2++){
			double m, dy;
			this->Moments(i2, m, dy);
			double y = std::abs(m - 1.0);
			result[i2].x[1] += 2*y*y;
			result[i2].x[2] += 4.*y*dy;
		}
	}
	else if (this->n==3){

# pragma omp parallel for shared(result)
		for (int i2 = i_; i2 < nh + 1; i2 ++){
			for (int i3 = -nh; i3 < nh+1 ; i3 ++){
				int ri = std::max(i2, std::abs(i3));
				bool isLongRange3 = (std::abs(i3) > i_) && ( std::abs(i2-i3) > i_ );
				if ( isLongRange3  && (ri < nh)){
					double m, dy;
					this->Moments((size_t)i2*lr + (i3+nh), m, dy);
					double y = std::abs(m - 1.0);
# pragma omp atomic
					result[ri].x[1] += y*y;
# pragma omp atomic
					result[ri].x[2] += 2.*y*dy;
				}
			}
		}
	}
	else if (this->n==4){

# pragma omp parallel for shared(result)
		for (int i2 = i_; i2 < nh + 1; i2 ++){
			for (int i3 = -nh; i3 < nh+1 ; i3 ++){
				bool isLongRange3 = (std::abs(i3) > i_) && ( std::abs(i2-i3) > i_ );
				if (isLongRange3){
					int rj = std::max(i2, std::abs(i3));
					for (int i4 = -nh; i4 < nh+1; i4 ++ ){
						int ri = std::max(rj, std::abs(i4));
						bool isLongRange4 = (std::abs(i4) > i_) && (std::abs(i2-i4) > i_) && (std::abs(i3-i4) > i_);
						if (isLongRange4 && (ri < nh)){
							double m, dy;
							this->Moments(((size_t)i2*lr + (i3+nh))*lr + (i4+nh), m, dy);
							double y = std::abs(m - 1.0);
		# pragma omp atomic
							result[ri].x[1] += y*y;
		# pragma omp atomic
							result[ri].x[2] += 2.*y*dy;
						}
					}

				}
			}
		}
	}

	/** Cumulative integration (bins below i_ are empty; start at 1 so result[i-1] exists) */
	result[0].x[1] *= dv;
	result[0].x[2] *= dv;
	for (int i = std::max(i_, 1); i < nh; i++){
		result[i].x[1] *= dv;
		result[i].x[2] *= dv;

		result[i].x[1] += result[i-1].x[1];
		result[i].x[2] += result[i-1].x[2];
	}
	result.erase(result.begin(), result.begin()+i_);
}

namespace {

/** \brief An occupied bin around one particle. */
struct BinCount {
	int bin;	///< bin index i + num_half in [0, len_row)
	int count;	///< number of particles in the bin (x_1 itself included in bin num_half)
};

/** \brief Floor of q/N for a possibly negative q.
 *  \param q  Numerator.
 *  \param N  Positive denominator.
 *  \return floor(q/N). */
inline long FloorDiv(long q, long N){
	return (q >= 0) ? q / N : -((-q + N - 1) / N);
}

}

void gn::AddConfig(const Configuration & c){
	if (c.GetDimension() != 1)
		throw std::invalid_argument("gn::AddConfig: only 1D configurations are supported");
	const long N = (long)c.NumParticle();
	if (N == 0)
		throw std::invalid_argument("gn::AddConfig: empty configuration");
	const double L = std::sqrt(c.GetBasisVector(0).Modulus2());
	const double dx = this->dx;
	const int nh = this->num_half;
	const size_t lr = this->len_row;
	// Search up to the outer edge of the outermost bins, so every bin is fully sampled.
	const double rcut = (nh + 0.5) * dx;
	const double rho = N / L;
	const double C = 1.0 / (N * std::pow(rho * dx, (double)(this->n - 1)));

	/* Sorted positions; pos(q) is the unwrapped position of the q-th point of the periodic sequence. */
	std::vector<double> x(N);
	for (long i = 0; i < N; i++){
		double r = c.GetRelativeCoordinates(i).x[0];
		x[i] = (r - std::floor(r)) * L;
	}
	std::sort(x.begin(), x.end());
	auto pos = [&x, N, L](long q) -> double { long w = FloorDiv(q, N); return x[q - w*N] + w*L; };

	/* Run-length binned neighbor list of every particle (CSR: entries[off[p]..off[p+1])).
	 * Unwrapped positions increase with q, so bins arrive in ascending order. */
	std::vector<BinCount> entries;
	entries.reserve((size_t)N * std::min(lr, (size_t)(2.0 * rho * rcut + 3.0)));
	std::vector<size_t> off(N + 1);
	for (long p = 0; p < N; p++){
		off[p] = entries.size();
		long q = p;
		while (x[p] - pos(q - 1) <= rcut) --q;
		for (; pos(q) - x[p] <= rcut; ++q){
			long b = std::lround((pos(q) - x[p]) / dx);
			if (b < -nh || b > nh) continue;
			int bin = (int)b + nh;
			if (entries.size() > off[p] && entries.back().bin == bin)
				entries.back().count++;
			else
				entries.push_back({bin, 1});
		}
	}
	off[N] = entries.size();

	if (this->n == 2){
		std::vector<double> tot(lr, 0.0);
		for (const BinCount & e : entries) tot[e.bin] += e.count;
		for (int i = 0; i <= nh; i++){
			// Symmetrized over +-x_21; bin 0 excludes x_1 itself.
			double v = C * (0.5 * (tot[nh + i] + tot[nh - i]) - (i == 0 ? (double)N : 0.0));
			this->sum1[i] += v;
			this->sum2[i] += v * v;
		}
		this->curr_c_idx++;
		return;
	}

	/* Inverted index over x_21 >= 0: bucket[boff[i2-nh] .. boff[i2-nh+1]) = (particle, occupancy). */
	std::vector<size_t> boff(nh + 2, 0);
	for (const BinCount & e : entries)
		if (e.bin >= nh) boff[e.bin - nh + 1]++;
	for (int k = 0; k <= nh; k++) boff[k + 1] += boff[k];
	std::vector<std::pair<long, int>> bucket(boff[nh + 1]);
	{
		std::vector<size_t> fill(boff.begin(), boff.end() - 1);
		for (long p = 0; p < N; p++)
			for (size_t k = off[p]; k < off[p + 1]; k++)
				if (entries[k].bin >= nh) bucket[fill[entries[k].bin - nh]++] = {p, entries[k].count};
	}

	/* Each x_21 slab is owned by one thread: no atomics, and a fixed summation order, so the
	 * result is independent of the thread count.  Below ~1e5 operations a parallel region
	 * costs more than the kernel (ponytail: fixed threshold; tune if the core count changes). */
	const size_t slab_len = (this->n == 3) ? lr : lr * lr;
	const double mean_len = (double)entries.size() / N;
	const double work = (double)entries.size() * std::pow(mean_len, (double)(this->n - 2));
	if (this->slab.size() < (size_t)omp_get_max_threads())
		this->slab.resize(omp_get_max_threads());
	const bool n4 = (this->n == 4);

#pragma omp parallel if(work > 1e5)
	{
		std::vector<double> & buf = this->slab[omp_get_thread_num()];
		if (buf.size() != slab_len) buf.assign(slab_len, 0.0);
		std::vector<size_t> touched;
		// Contributions are positive integers, so buf[t] == 0 marks a bin not yet touched in this slab.
		auto add = [&buf, &touched](size_t t, double v){
			if (buf[t] == 0.0) touched.push_back(t);
			buf[t] += v;
		};

#pragma omp for schedule(dynamic)
		for (int i2 = nh; i2 <= 2 * nh; i2++){
			for (size_t k = boff[i2 - nh]; k < boff[i2 - nh + 1]; k++){
				const long p = bucket[k].first;
				const double n2 = bucket[k].second - (i2 == nh);
				if (n2 <= 0) continue;
				const BinCount * e = entries.data() + off[p];
				const size_t m = off[p + 1] - off[p];
				for (size_t a = 0; a < m; a++){
					const int i3 = e[a].bin;
					const double n3 = e[a].count - (i3 == nh) - (i3 == i2);
					if (n3 <= 0) continue;
					if (!n4){
						add(i3, n2 * n3);
						continue;
					}
					const double f = n2 * n3;
					const size_t row = (size_t)i3 * lr;
					for (size_t b = 0; b < m; b++){
						const int i4 = e[b].bin;
						const double n4_ = e[b].count - (i4 == nh) - (i4 == i2) - (i4 == i3);
						if (n4_ > 0) add(row + i4, f * n4_);
					}
				}
			}
			/* Flush this slab into the ensemble sums and reset the touched bins. */
			const size_t base = (size_t)(i2 - nh) * slab_len;
			for (size_t t : touched){
				const double v = buf[t] * C;
				this->sum1[base + t] += v;
				this->sum2[base + t] += v * v;
				buf[t] = 0.0;
			}
			touched.clear();
		}
	}
	this->curr_c_idx++;
}

void gn::Write_gn(const std::string & prefix){
	auto write_row = [this](std::ostream & ofile, size_t idx){
		double m, dy;
		this->Moments(idx, m, dy);
		ofile << std::scientific << std::setprecision(8) << m;
		if (this->curr_c_idx > 1)
			ofile << "\t" << dy;
		ofile << "\n";
	};
	if (this->n == 2){
		/* g2 */
		std::fstream ofile (prefix + "_g2.txt", std::fstream::out);
		ofile << "#x = x_21, y = g_2(x), dy = dg_2(x)\n";
		ofile << "#x\t#y\t#dy\n";
		for (int i2 = 0; i2 <= this->num_half; i2++){
			double x21 = (double)i2 * this->dx;
			ofile << std::fixed << std::setprecision(3) << x21 <<"\t" ;
			write_row(ofile, i2);
		}
		ofile.close();
	}
	else if (this->n == 3){
		/* g3 */
		std::fstream ofile (prefix + "_g3.txt", std::fstream::out);
		ofile << "#x2 = x_21, x3=x_31 y = g_3(x2,x3), dy = dg_3(x2,x3)\n";
		ofile << "#x2\t#x3\t#y\t#dy\n";

		size_t idx = 0;
		for (int i2 = 0; i2 <= this->num_half; i2++){
			double x21 = (double)i2 * this->dx;
			for (int i3 = -1*this->num_half; i3 <= this->num_half; i3 ++ ){
				double x31 = (double)i3 * this->dx;
				ofile << std::fixed << std::setprecision(3) << x21 << "\t" << x31 <<"\t" ;
				write_row(ofile, idx++);
			}
		}
		ofile.close();
	}
	else if (this->n == 4){

#pragma omp parallel for
		for (int i2 = 0; i2 <= this->num_half; i2++){
			double x21 = (double)i2 * this->dx;
			char x21_str[64];
			snprintf(x21_str, sizeof(x21_str), "%0.3f", x21);
			std::fstream ofile (prefix + "_g4_x21_" + x21_str + ".txt", std::fstream::out);

			ofile << "#x2 = x_21, x3=x_31, x4=x_41, y = g_4(x2,x3,x4), dy = dg_4(x2,x3,x4)\n";
			ofile << "#x3\t#x4\t#y\t#dy\n";
			size_t idx = this->len_row*this->len_row*i2;
			for (int i3 = -1*this->num_half; i3 <= this->num_half; i3 ++ ){
				double x31 = (double)i3 * this->dx;
				for (int i4 = -1*this->num_half; i4 <= this->num_half; i4 ++ ){
					double x41 = (double)i4 * this->dx;
					ofile << std::fixed << std::setprecision(3) << x31 <<"\t"<<x41 <<"\t";
					write_row(ofile, idx++);
				}
			}
			ofile.close();
		}

	}
}
