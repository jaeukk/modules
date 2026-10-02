/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: Feburary 2023 */

/** \file main.cpp
 *	\brief Implementations for computing g_3 and g_4 functions for 1D configurations 
 */


#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

//-------Core part--------
#include <etc.h>
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
//#include "../StructureFactor.h"

#include "tau_n.h"
#include "gn.h"

size_t Verbosity = 3;

/** \brief Hole probability H_V(R) and void exclusion probability E_V(R) of an ensemble,
 *  sampled at random test points.
 *  \param getc              getc(i) returns the i-th configuration.
 *  \param numconfigs        Number of configurations.
 *  \param dR                Bin width in R [length].
 *  \param[out] Ev           (R, E_V(R), error) per bin.
 *  \param[out] Hv           (R, H_V(R), error) per bin.
 *  \param SamplesPerConfig  Random test points per configuration. */
void Measure_Hole_Prob(std::function<const Configuration(size_t i)> getc, size_t numconfigs,
	double dR, std::vector<GeometryVector> & Ev, std::vector<GeometryVector> & Hv, size_t SamplesPerConfig = 1000);

/** \brief CLI entry point; argv[1] selects the mode (test, gn_1D, taun_1D, hole, overlap),
 *  and the remaining parameters are read from stdin.
 *  gn_1D [hours] [n]: time limit in hours (default 10); a third argument "n" skips writing g_n.
 *  taun_1D [hours]: time limit in hours (default 10).
 *  \param argc  Argument count.
 *  \param argv  Argument vector.
 *  \return 0 on success, 1 on invalid input. */
int main (int argc, char ** argv){

	std::string tempstring;
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[200] = {};
	omp_set_num_threads(4);

	if (argc==1){
		ofile << "Quit because of no input\n";
		return 0;
	}

	if (strcmp(argv[1], "test") == 0){
		/* test */
		int nh;
		size_t n =4;
		ofile << "nh = ";
		ifile >> nh;
		int lr = 2*nh+1;
		size_t total_length = (nh+1)*pow(lr, n-2);
		ofile << nh <<std::endl;
		ofile << lr <<std::endl;
		ofile << total_length <<std::endl;
	//# pragma omp parallel for 


		int idx = 0;
		for (int i2 = nh; i2 < 2*nh + 1; i2 ++){
			std::cout << "------" << std::endl <<std::endl;
			for (int i3 = 0; i3 < 2*nh+1 ; i3 ++){
				std::cout << std::endl;
				for (int i4 = 0; i4 < 2*nh+1; i4 ++ ){
					std::vector<size_t> indices1 = {(size_t)(i2-nh),(size_t)i3,(size_t)i4};
					std::vector<size_t> indices2;
					size_t idx1 = vector2idx(indices1, lr);
					idx2vector(idx1, indices2, 4, lr);
					size_t idx2 = vector2idx(indices2, lr);
					
					// #pragma omp critical
					// 	{
					// 		std::cout << i2 <<"\t"<<i3<<"\t"<<i4<<"\t"<<idx1 << "\t" << idx2 <<std::endl;
					// 		std::cout << indices1[0] <<"\t"<< indices1[1] <<"\t"<< indices1[2] << std::endl;
					// 	}	
					int tt = (int)idx1;
					ofile << indices1[0] <<"\t"<< indices1[1] <<"\t" << indices1[2] <<"\t" << tt <<"\t" ;

					ofile << (int)idx2 << "\t" << idx <<std::endl;

					idx ++;


					if (idx1 >= total_length){
						#pragma omp critical
						{
							ofile << "overflow: " << i2 <<"\t"<<i3<<"\t"<<i4<<"\t"<<idx1 <<" / "<<total_length <<std::endl;
						}
					}

					bool correct = true;
					for (size_t i=0; i<3; i++){
						correct &= (indices1[i] == indices2[i]);
					}
					if (! correct){	
						#pragma omp critical
						{
							ofile << "conversion error: " << i2 <<"\t"<<i3<<"\t"<<i4<<"\t"<<idx1  <<std::endl;
						}					
					}
					if (! (idx1 == idx2)){	
						#pragma omp critical
						{
							ofile << "vector2idx: " << i2 <<"\t"<<i3<<"\t"<<i4<<"\t"<<idx2  <<std::endl;
						}					
					}
						
				}
			}
		}
		return 0;
	}

	else if (strcmp(argv[1], "gn_1D") == 0){
		/* Time limit */
		double timelimit_in_hour = 10;
		if (argc > 2){
			timelimit_in_hour = (double) std::atof(argv[2]);
		}

		auto ProgramStart = std::time(nullptr);
		auto TimeLimit = ProgramStart + timelimit_in_hour*3600;
		ofile << "Time limit = "<< timelimit_in_hour << " hours or, "<< timelimit_in_hour*3600 <<" sec"<<std::endl;
		
		gn model;
		try{
			model = gn(ofile, ifile);
		}
		catch (const std::invalid_argument & e){
			ofile << e.what() << std::endl;
			return 1;
		}
		// Legacy tested argv[2] (the time limit), so "n" was never honored.
		if (argc > 3){
			if (strcmp(argv[3], "n") == 0){
				model.DisablePrint_gn();		
			}
		}

		try{
			model.Run(ofile, TimeLimit);
		}
		catch (const std::invalid_argument & e){
			ofile << e.what() << std::endl;
			return 1;
		}
	}
	else if (strcmp(argv[1], "Kirkwood")==0){
		/* Compute the Kirkwood approximation of g_n from g_n */

		/* Load g2 */

		/* Load gn */

		/* Compute gn approximation */

		/* Compare gn approximation */

	}
	else if (strcmp(argv[1], "taun_1D") == 0){
		/* Time limit */
		double timelimit_in_hour = 10;
		if (argc == 3){
			timelimit_in_hour = (double) std::atof(argv[2]);
		}
		auto ProgramStart = std::time(nullptr);
		auto TimeLimit = ProgramStart + timelimit_in_hour*3600;
		ofile << "Time limit = "<< timelimit_in_hour << " hours or, "<< timelimit_in_hour*3600 <<" sec"<<std::endl;
		
		tau_n model(ofile, ifile);
		model.Run(ofile, TimeLimit, false);
		//model.Test2(ofile);
	}
	else if (strcmp(argv[1], "hole") == 0){
		ofile << "-------------------------------------- \n";
		ofile << "Compute void-related statistics, Ev Hv \n";
		ofile << "-------------------------------------- \n";
		
		double dR;
		long long idx1, idx2;
		size_t num_threads, samples, Nc;
		ofile << "Input ConfigPack file = ";
		ifile >> name;
		ConfigurationPack cp(name);
		ofile << "\nIndex to begin = ";
		ifile >> idx1;
		ofile << "\nIndex to end = ";
		ifile >> idx2;
		if(idx2 > cp.NumConfig()){
			idx2 = cp.NumConfig();
			ofile <<"This ConfigurationPack contains only " << idx2 << " configurations, which is set to be the num_config\n";
		}
		Nc = idx2 - idx1;

		ofile << "\nNum Threads = ";
		ifile >> num_threads;	omp_set_num_threads(num_threads);
		ofile << "\ndR = ";
		ifile >> dR;
		ofile << "\nSamples per configuration = ";
		ifile >> samples;
		ofile << "\nsave prefix = ";
		ifile >> tempstring;

		auto getc = [&cp, &idx1](size_t i)->Configuration{
			return cp.GetConfig(i+idx1);
		};		
		std::vector<GeometryVector> Hv;
		std::vector<GeometryVector> Ev;
		Measure_Hole_Prob(getc, Nc, dR, Ev, Hv, samples);

		std::string filename = tempstring + "__Hv";
		{
			std::vector<std::string> comm = {
				"ConfigPack: "+std::string(name)+"\n",
				"Index range: "+std::to_string(idx1)+"\t"+std::to_string(idx2)+"\n",
				"R\t#Hv#dHv"
			};
			WriteFunction(Hv, filename, comm);
		}

		filename = tempstring + "__Ev";
		{
			std::vector<std::string> comm = {
				"ConfigPack: "+std::string(name)+"\n",
				"Index range: "+std::to_string(idx1)+"\t"+std::to_string(idx2)+"\n",
				"R\t#Ev#dEv"
			};
			WriteFunction(Ev, filename, comm);
		}
	}
	else if (strcmp(argv[1], "overlap") == 0){
		ofile << "-------------------------------------- \n";
		ofile << "Compute overlap function Q(t;a) \n";
		ofile << "	in NV ensembles. \n";
		ofile << "-------------------------------------- \n";

		double sigma_max, dsigma;
		long long idx1, idx2;
		size_t num_average = 10;
		size_t num_threads, Nc;
		ofile << "Input ConfigPack file = ";
		ifile >> name;
		ConfigurationPack cp(name);
		ofile << "\nIndex to begin = ";
		ifile >> idx1;
		ofile << "\nIndex to end = ";
		ifile >> idx2;
		if(idx2 > cp.NumConfig()){
			idx2 = cp.NumConfig();
			ofile <<"This ConfigurationPack contains only " << idx2 << " configurations, which is set to be the num_config\n";
		}
		Nc = idx2 - idx1;

		ofile << "\nNum Threads = ";
		ifile >> num_threads;	omp_set_num_threads(num_threads);
		ofile << "\nMaximum overlap radius = ";
		ifile >> sigma_max;
		ofile << "\nIncrement in overlap radius = ";
		ifile >> dsigma;
		ofile << "\nsave prefix = ";
		ifile >> tempstring;

		auto getc = [&cp, &idx1](size_t i)->Configuration{
			return cp.GetConfig(i+idx1);
		};		
		
		/* initial configuration */
		double rho = 0.0, length = 0.0, dlength = 0.0, N = 0.0;
		std::vector<Configuration> list_c0;
		for (size_t i = 0; i< num_average; i++){
			list_c0.emplace_back(getc(i));
			
			if (i == 0){
				Configuration c0 = getc(i);
				rho = c0.NumParticle() / c0.PeriodicVolume();
				length = std::pow(rho, -1.0/c0.GetDimension());
				dlength = dsigma * length;
				N = c0.NumParticle();
			}			
		}

		/* initialize the list */
		size_t num_bins = (size_t)std::floor(sigma_max/dsigma);
		std::vector<std::vector<double>> time_series(Nc, std::vector<double> (2*num_bins + 1, 0.0L));
		
		for (size_t t = 0; t<Nc; t++){
			Configuration ci = getc(t);	/* target system */
			time_series [t][0] = (double)t;
			
			/* load a reference system from list_c0 */
			for (size_t ref = 0; (ref < t+1) && (ref < list_c0.size()) ; ref++){
				
				/* compute the overlap function for a give configuration.
					Q[j] = Q(t_i, a=dsigma*(j+1)) */
				std::vector<double> Q(num_bins,0.0L);
				Configuration * c0 = &list_c0.at(ref);

			#pragma omp parallel for 
				for (size_t prt_idx = 0; prt_idx < c0->NumParticle(); prt_idx++){
					GeometryVector diff = c0->GetRelativeCoordinates(prt_idx) - ci.GetRelativeCoordinates(prt_idx);
					c0->RelativeCoordToMinimumImage(diff);
					double distance = std::sqrt( c0->RelativeCoord2CartesianCoord(diff).Modulus2());
					size_t idx = (size_t) std::floor(distance / dlength );

					if (idx < Q.size()){
			#pragma omp atomic
						Q.at(idx) ++;
					}				
				}

				for (size_t j = 1; j<Q.size(); j++){
					Q.at(j) += Q.at(j-1);
				}
				size_t delta_t = t;
				if (t < ref){
					ofile << "something is wrong\n";
				}
				else{
					delta_t -= ref;
				}
				for (size_t j = 0; j<Q.size(); j++){
					Q[j] /= N;

					time_series[delta_t][2*j+1] += Q[j];
					time_series[delta_t][2*j+2] += Q[j]*Q[j];
				}
			}
		}

	#pragma omp parallel for
		for(size_t i = 0; i < time_series.size(); i++){
			double num = ((time_series.size() - num_average) > i)? (double)num_average : (double) (time_series.size() - i + 1); 
			for (size_t j = 0; j < num_bins; j++){
				time_series[i][2*j+1] /= num ;
				time_series[i][2*j+2] /= num ;
				if (num > 0.5){
					time_series[i][2*j+2] = std::sqrt((time_series[i][2*j+2] - time_series[i][2*j+1]*time_series[i][2*j+1]) / (num - 1.0));
				}
			}
		}


		/* output */
		{
			std::string filename = tempstring + "_OverlapFunction.txt";
			std::fstream ofile (filename, std::fstream::out);
			
			ofile << "# Overlap function to check coherence between configurations\n";
			ofile << "# dr = " << dlength << "\n";
			ofile << "# diff. idx";
			for (size_t i =0; i< num_bins; i++){
				ofile << "\t dr *"<<i;
				ofile << "\t err";
			}
			ofile << "\n";
			for (size_t i =0; i<time_series.size(); i++){
				ofile <<std::fixed<< std::setprecision(1) << time_series[i][0];
				for (size_t j=1; j<time_series[i].size(); j++){
					ofile <<"\t"<<std::scientific << std::setprecision(8)<< time_series[i][j];
				}
				ofile <<"\n";
			}

			ofile.close();
		}	
	}
	return 0;
}



void Measure_Hole_Prob(std::function<const Configuration(size_t i)> getc, size_t numconfigs, double dR, std::vector<GeometryVector> & Ev, std::vector<GeometryVector> & Hv, size_t SamplePerConfig){
	/* Resize Ev and Hv */
	DimensionType dim;
	double Rc = 0.0, Rc_mean = 0.0, Rc_error = 0.0;
	{
		Configuration c = getc(0);
		dim = c.GetDimension();
		GeometryVector x(dim);
		for (int i=0; i< dim; i++){
			x = x + c.GetBasisVector(i);
		}
		double L = 0.5*std::sqrt(x.Modulus2());
		Ev.resize((size_t)round(L/dR), GeometryVector(0,0,0,0));
		Hv.resize((size_t)round(L/dR), GeometryVector(0,0,0,0));
		std::cout<< L <<" : " <<Ev.size() <<std::endl;
	}
	double rho = 0.0;
#pragma omp parallel for
	for (size_t i = 0; i< numconfigs; i++){
		Configuration c = getc(i);
#pragma omp atomic
		rho += c.NumParticle() / c.PeriodicVolume();
	}
	rho /= (double)(numconfigs);

#pragma omp parallel for
	for (size_t i=0; i<Hv.size(); i++){
		double R = dR*(0.5+i);
		//Hv[i].x[3] = 1.0/dR;
		Hv[i].x[0] = R;
		Ev[i].x[0] = R;
	}

	size_t UPPER = 0;
	double acc = 1.0/SamplePerConfig;

	std::cout<< "Computing Hole Probabilities:"<<std::endl;
	progress_display pd(numconfigs);
#pragma omp parallel for
	for (size_t i = 0; i< numconfigs; i++){
		double rc = 0.0;
		std::vector<double> NR(Ev.size(), 0.0);
		RandomGenerator rng(i);
		Configuration c = getc(i);
		GeometryVector x0(dim);
		size_t upper = 0; // define the support!
		for (size_t j = 0; j < SamplePerConfig; j++){
			for (DimensionType k = 0; k < dim ; k++)
				x0.x[k] = rng.RandomDouble();

			double dist = c.NearestParticleDistance(x0);
			size_t id = (size_t)floor(dist/dR);
			NR.at(id) += acc ;
			upper = (upper < id)? id : upper;
			rc = (rc < dist) ?  dist : rc;
		}
#pragma omp critical
		{
			UPPER = (upper > UPPER)? upper: UPPER;		
		}
#pragma omp critical
		{
			Rc = (Rc < rc)? rc: Rc;			
		}
#pragma omp atomic
			Rc_mean += rc;
#pragma omp atomic 
			Rc_error += rc*rc;

		/** Add Hv */
		for (size_t j = 0; j<=upper; j++){
			double val_Gv = NR[j]; 
#pragma omp atomic
			Hv[j].x[1] += val_Gv;
#pragma omp atomic 
			Hv[j].x[2] += val_Gv*val_Gv;
		}
		/** Add Ev */
		double val_Ev = 0.0;
		for (size_t j = 0; j <= upper; j++){
			size_t jj = upper - j;
//			std::cout<<jj<<std::endl;
			val_Ev += NR[jj];
#pragma omp atomic
			Ev[jj].x[1] += val_Ev;
#pragma omp atomic 
			Ev[jj].x[2] += val_Ev*val_Ev;
		}

#pragma omp critical
		{
			pd ++ ;
		}

	}

	/** Summarize data */
	Ev.resize(UPPER+2);
	Hv.resize(UPPER+2);
	double Nc = (double)(numconfigs);
#pragma omp parallel for
	for (size_t i = 0; i<Hv.size(); i++){
		Hv[i].x[1]/=Nc*dR;
		Hv[i].x[2]/=Nc*dR*dR;
		Hv[i].x[2] = std::sqrt((Hv[i].x[2] - Hv[i].x[1]*Hv[i].x[1])/(Nc-1.0));
		Hv[i].x[3] = 0.0;
	}

#pragma omp parallel for
	for (size_t i = 0; i<Ev.size(); i++){
		Ev[i].x[1]/=Nc;
		Ev[i].x[2]/=Nc;
		Ev[i].x[2] = std::sqrt((Ev[i].x[2] - Ev[i].x[1]*Ev[i].x[1])/(Nc-1.0));
	}	

	/** Largest holes */
	Rc_mean /= Nc;
	Rc_error = std::sqrt(Rc_error / Nc - Rc_mean*Rc_mean);
	Rc_error /= std::sqrt(Nc-1.0);
	std::cout << "R_c = " << Rc <<std::endl;
	std::cout << "<R_c> = " << Rc_mean << "+-" << Rc_error << std::endl;
}




