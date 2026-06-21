/**	Author	:	Jaeuk Kim 
 * 	Email	:	phy000.kim@gmail.com
 * 	Date	:	March	2021
*/

/**	\file main.cpp
 * \brief	Function implementations of performing 
 * the fast-Fourier-transform-based homogenization (electrostatic) for two-phase media.
*/


#include "FFH.h"
#include "TwoPhaseDispersion.h"
#include <GeometryVector.h>
//#include <StructureFactor.h>
#include <vector>
#include <functional>


size_t FFH_verbosity = 2;

/** \brief Read a list of DisArray files, perform homogenization, and compute statistics.
 *  \param argc  Argument count.
 *  \param argv  [1] DisArray prefix (printf format), [2] number of files, [3] epsilon2,
 *               [4] output prefix, [5] number of threads.
 *  \return      0 on success. */
int main(int argc, char ** argv) {
	char tempstring[1000] = {};
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[200] = {};
	char dir[] = "./results/n1/";
	


	RandomGenerator rng(0);
	fftw_init_threads();
	std::string dump;
	
	std::string load_prefix(argv[1]);
	int num_files = atoi(argv[2]);
	double epsilon1 = 1.0L, epsilon2 = atof(argv[3]);
	std::string save_prefix(argv[4]);
	int num_threads = atoi(argv[5]);
	omp_set_num_threads(num_threads);

	int num_pixels;
	std::vector<GeometryVector> epsilon_eff;
	epsilon_eff.reserve(num_files);
	
	for (int i = 0; i < num_files; i++){
		sprintf(name, load_prefix.c_str(), i);
		
		/* read a configuration */
		TwoPhaseDispersion config(name);
		DimensionType d = config.packing.GetDimension();
		
		/* set up the problem */
		StaticDielectricProblem stat_prob(config, epsilon1, epsilon2);
		stat_prob.SetMethod(true);

		{
			ConstantTensor E0(d, 1);
			E0[0] = 1.0L;
			stat_prob.DefineProperties(epsilon1, epsilon2);
			stat_prob.Prepare(E0);
		}

		/* compute ...  */
		stat_prob.Run();
		std::complex<double> epsilone = stat_prob.ComputeEffectiveEpsilon();

		epsilon_eff.emplace_back((double)i, epsilone.real(), epsilone.imag());
	}

	std::vector<std::string> comments;
	comments.emplace_back("#ID\tRe[eps_eff]\tIm[eps_eff]");
	WriteFunction(epsilon_eff, save_prefix, comments);

	if(num_files > 1){
		/* Compute average */
		double av_real = 0.0, av_image = 0.0;
		double sq_real = 0.0, sq_image = 0.0;

		for(int i = 0; i < num_files; i++){
			av_real += epsilon_eff[i].x[1];
			av_image += epsilon_eff[i].x[2];

			sq_real += epsilon_eff[i].x[1] * epsilon_eff[i].x[1];
			sq_image += epsilon_eff[i].x[2] * epsilon_eff[i].x[2];
		}
		av_real /= (double)num_files;
		av_image /= (double)num_files;
		sq_real = sqrt(sq_real/(double)num_files - av_real * av_real );
		sq_image = sqrt(sq_image/(double)num_files - av_image * av_image );
		sq_real /= sqrt(num_files - 1);
		sq_image /= sqrt(num_files - 1);

		std::ofstream outfile;
		outfile.open(save_prefix+".txt", std::ios_base::app);
		outfile << "Av. real part = "<< av_real <<" pm " << sq_real <<"\n"; 
		outfile << "Av. imag part = "<< av_image <<" pm " << sq_image <<"\n"; 
		outfile.close();
	}


	return 0;
}
