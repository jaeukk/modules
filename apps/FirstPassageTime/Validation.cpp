/** \file Validation.cpp
  * \brief Off-line debug/validation harness (not called by main). */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>



#include <etc.h>
#include <GeometryVector.h>
#include <RandomGenerator.h>

#include <ParticleShapes.h>

#include <PeriodicCellList.h>
#include <Dispersions.h>

#include "FirstPassageTime.h"
#include <chrono>

void Test_interface(double x, double y, size_t num){

	double count_total_trials = (double) num, count2host = 0.0;
	double sigma1 = 1.0, sigma2 = 4.0;
	
	FirstPassageTimeParameters param(sigma1, sigma2, 2);
	Dispersion test;
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/Z2-tri");
	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/checkerboard/regular2");
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/honeycomb-0.50");
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/2D_HSF_0.20_100-f_0.5-0");
	Trajectory one;
	RandomGenerator rng(0);
	one.sigma_host = sigma1;	one.sigma_particle = sigma2;
	one.dimension = 2;
	one.pos_curr = GeometryVector(x,y);
	std::cout << "initial position = " << one.pos_curr;
	for(size_t i=0; i<num; i++){
		param.PerformSimulation2D_test(test, 1, one, rng);
		if(test.IsInsideParticle(one.pos_curr) == -1){
			count2host ++;
		}
		one.pos_curr = GeometryVector(x,y);
	}
	std::cout<<"p_1 = "<< count2host / count_total_trials <<"\n";

}



void Test_simulation(double x, double y, size_t num_steps){
	double sigma1 = 1.0, sigma2 = 0.0;
	FirstPassageTimeParameters param(sigma1, sigma2, 2);
	Dispersion test;
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/Z2-tri");
	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/checkerboard/regular2");
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/honeycomb-0.50");
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/2D_HSF_0.20_100-f_0.5-0");
	Trajectory one;
	RandomGenerator rng(0);
	one.sigma_host = sigma1;	one.sigma_particle = sigma2;
	one.dimension = 2;
	one.pos_curr = GeometryVector(x,y);
	one.delta_pos = GeometryVector(0,0);
	std::cout << "initial position = " << one.pos_curr;

	param.PerformSimulation2D_test(test, num_steps, one, rng);

	std::cout<<"R2 total, R2_host, R2_prt, R2_s, tau_s\n";
	std::cout << one.R2_host / one.R2_particle <<"\n";
//	if (i==0){
//		sprintf(temp, "data = np.array([[%1.7e, %1.7e, %1.7e, %1.7e, %1.7e, %1.7e,]", one.GetR2Total(), one.R2_host, one.R2_particle, one.R2_s, one.tau_s, one.GetEffectiveConductivity());
//	}
//	else{
//		sprintf(temp, ",\n[%1.7e, %1.7e, %1.7e, %1.7e, %1.7e, %1.7e,]", one.GetR2Total(), one.R2_host, one.R2_particle, one.R2_s, one.tau_s, one.GetEffectiveConductivity());
//	}

//	ofile<<temp;
//	ofile<<"])\n";

	 
}

void Test_p1(){
	FirstPassageTimeParameters param(1.0, 10.0, 2);
	double R = 1.0;
	
	std::vector<GeometryVector> p1_1, p1_2, taus_1, taus_2;

	for (double r = 0.01; r<R; r+=0.01){
		p1_1.emplace_back(r, param.p1_2D(r, true));
		p1_2.emplace_back(r, param.p1_2D(r, false));

		taus_1.emplace_back(r, param.tau_s_2D(r, R, true));
		taus_2.emplace_back(r, param.tau_s_2D(r, R, false));
	}
	WriteFunction(p1_1, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/p1_1");
	WriteFunction(p1_2, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/p1_2");
	WriteFunction(taus_1, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/taus_1");
	WriteFunction(taus_2, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/taus_2");
}

void Test_Dispersion(){
	/* read an example dispersion. */
	Dispersion test;
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/honeycomb-0.50");
	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/2D_HSF_0.20_100-f_0.5-0");

	std::cout << * test.GetShape(86);
	std::cout << "\n";
	std::cout << test.GetCartesianCoordinates(86);
	RandomGenerator rng(1);
	std::vector<GeometryVector> inside, outside;
	inside.reserve(1000);
	outside.reserve(1000);

	for (int i=0; i<1000; i++){
		GeometryVector p_test(rng.RandomDouble()/10., rng.RandomDouble()/10.);
		GeometryVector p_Cartesian = test.RelativeCoord2CartesianCoord(p_test);
		p_Cartesian.SetDimension(3);
		p_Cartesian.x[2] = i;
		if(test.IsInsideParticle(p_test)==-1){
			outside.push_back(p_Cartesian);
		}
		else{
			inside.push_back(p_Cartesian);
		}
	}

	WriteFunction(inside, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/pts_inside");
	WriteFunction(outside, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/pts_outside");

	std::vector<GeometryVector> centers, origins, radii;


	for (int i=0; i<30; i++){
		GeometryVector p_test(rng.RandomDouble()/10., rng.RandomDouble()/10.);
		GeometryVector p_Cartesian = test.RelativeCoord2CartesianCoord(p_test);
		p_Cartesian.SetDimension(3);
		p_Cartesian.x[2] = i;

		std::vector<GeometryVector> list;
		GeometryVector x0, n;

		test.FindNearestNeighbors(p_test, list, x0, n);
		
		centers.push_back(p_Cartesian);
		origins.push_back(test.RelativeCoord2CartesianCoord(x0));
		radii.emplace_back(list[0].x[0], list[1].x[0], list[2].x[0]);

	}
	WriteFunction(centers, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/test_particles");
	WriteFunction(origins, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/test_x0");
	WriteFunction(radii, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/test_radii");
}

void Test_IsInside(){
	std::vector<GeometryVector> vs;
	char temp_str[300]="";
	vs.emplace_back(-0.5, -0.5);
	vs.emplace_back(0.5, -0.5);
	vs.emplace_back(0.5, 0.5);
	vs.emplace_back(0.4, 0.3);
	vs.emplace_back(-0.5, 0.5);

	Polygon temp(vs);
	std::fstream ofile("test_isinside.txt", std::fstream::out);



	ofile<< "vs = np.array([";
	for (size_t i=0; i < vs.size(); i++){
		snprintf(temp_str, sizeof(temp_str), "[%0.6f,%0.6f],\n", vs[i].x[0],vs[i].x[1]);

		ofile<< temp_str;
	}
	snprintf(temp_str, sizeof(temp_str), "[%0.6f,%0.6f]])\n\n", vs[0].x[0],vs[0].x[1]);
	ofile << temp_str;


	int num = 1000;
	RandomGenerator rng(312);
	std::vector<GeometryVector> inside, outside;
	inside.reserve(num);
	outside.reserve(num);

	GeometryVector cm = temp.GetCenterOfMass();
	for (int i=0; i <num ; i++){
		GeometryVector ptr(2.0*(rng.RandomDouble()-0.5), 2.0*(rng.RandomDouble()-0.5));

		if(temp.IsPtInside(ptr-cm)){
			inside.emplace_back(ptr.x[0],ptr.x[1],(double)i);
		}
		else{
			outside.emplace_back(ptr.x[0],ptr.x[1],(double)i);
		}
	}

	ofile<<"ptr1=np.array([";
	for (size_t i=0; i<inside.size(); i++){
		if (i == 0)
			snprintf(temp_str, sizeof(temp_str), "[%0.6f,%0.6f,%0.1f]", inside[i].x[0], inside[i].x[1], inside[i].x[2]);
		else{
			snprintf(temp_str, sizeof(temp_str), ",\n[%0.6f,%0.6f,%0.1f]", inside[i].x[0], inside[i].x[1], inside[i].x[2]);}
		ofile<< temp_str;
	}
	ofile<<"])\n\n";

	ofile<<"ptr2=np.array([";
	for (size_t i=0; i<outside.size(); i++){
		if (i == 0)
			snprintf(temp_str, sizeof(temp_str), "[%0.6f,%0.6f,%0.1f]", outside[i].x[0], outside[i].x[1],outside[i].x[2]);
		else{
			snprintf(temp_str, sizeof(temp_str), ",\n[%0.6f,%0.6f,%0.1f]", outside[i].x[0], outside[i].x[1],outside[i].x[2]);}
		ofile<< temp_str;
	}
	ofile<<"])\n\n";
	ofile.close();
}


void Test_distance(){
	std::vector<GeometryVector> vs;
	char temp_str[300]="";
	vs.emplace_back(-0.5, -0.5);
	vs.emplace_back(0.5, -0.5);
	vs.emplace_back(0.5, 0.5);
	vs.emplace_back(0.4, 0.7);
	vs.emplace_back(-0.5, 0.5);

	Polygon temp(vs);
	std::fstream ofile("test_distance.txt", std::fstream::out);



	ofile<< "vs = np.array([";
	for (size_t i=0; i < vs.size(); i++){
		snprintf(temp_str, sizeof(temp_str), "[%0.6f,%0.6f],\n", vs[i].x[0],vs[i].x[1]);

		ofile<< temp_str;
	}
	snprintf(temp_str, sizeof(temp_str), "[%0.6f,%0.6f]])\n\n", vs[0].x[0],vs[0].x[1]);
	ofile << temp_str;


	int num = 1000;
	RandomGenerator rng(312);
	std::vector<GeometryVector> inside, outside;
	inside.reserve(num);
	outside.reserve(num);

	GeometryVector cm = temp.GetCenterOfMass();
	for (int i=0; i <num ; i++){
		GeometryVector ptr(2.0*(rng.RandomDouble()-0.5), 2.0*(rng.RandomDouble()-0.5));
		GeometryVector ptr_center = ptr-cm;
		std::vector<GeometryVector> list;
		GeometryVector x0, n;
		double distance = temp.GetDistance2Boundary(ptr_center, list, x0, n);
		if(temp.IsPtInside(ptr_center)){
			inside.emplace_back(ptr.x[0],ptr.x[1], -1.0*distance, (double)i);
		}
		else{
			outside.emplace_back(ptr.x[0],ptr.x[1], distance, (double)i);
		}

		if (i < 10){
			std::cout<< "----\t"<<i <<"\t----\n";
			std::cout<< "x = " << ptr;
			std::cout<< "x0 = "<< x0 + cm;
			std::cout<< "\n";
			for (auto l = list.begin(); l != list.end(); l++)
				std::cout << *l;
			
			
			std::cout<<"\n\n";
		}
	}

	ofile<<"ptr1=np.array([";
	for (size_t i=0; i<inside.size(); i++){
		if (i == 0)
			snprintf(temp_str, sizeof(temp_str), "[%0.6f,%0.6f,%0.6f,%0.1f]", inside[i].x[0], inside[i].x[1], inside[i].x[2], inside[i].x[3]);
		else{
			snprintf(temp_str, sizeof(temp_str), ",\n[%0.6f,%0.6f,%0.6f,%0.1f]", inside[i].x[0], inside[i].x[1], inside[i].x[2], inside[i].x[3]);}
		ofile<< temp_str;
	}
	ofile<<"])\n\n";

	ofile<<"ptr2=np.array([";
	for (size_t i=0; i<outside.size(); i++){
		if (i == 0)
			snprintf(temp_str, sizeof(temp_str), "[%0.6f,%0.6f,%0.6f,%0.1f]", outside[i].x[0], outside[i].x[1],outside[i].x[2],outside[i].x[3]);
		else{
			snprintf(temp_str, sizeof(temp_str), ",\n[%0.6f,%0.6f,%0.6f,%0.1f]", outside[i].x[0], outside[i].x[1],outside[i].x[2],outside[i].x[3]);}
		ofile<< temp_str;
	}
	ofile<<"])\n\n";
	ofile.close();
}
