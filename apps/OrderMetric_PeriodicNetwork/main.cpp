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

#include <SpectralDensity_Dispersion.h>

//#include "../StructureFactor.h"
//#include "../NumberVariance.h"
//#include "../PairCorrelation.h"


void WriteFunction(const std::vector<GeometryVector> & result, std::ostream & ofile, const std::vector<std::string> & fields);
void WriteFunction(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::vector<std::string> & fields);
void ReadFunction(std::vector<GeometryVector> & result, std::istream & ifile, size_t NumColumns);
void ReadFunction(std::vector<GeometryVector> & result, const std::string & InFilePrefix, size_t NumColumns);

void WriteToContourf(const std::vector<GeometryVector> & vec_k, const std::vector<GeometryVector> &s_k, const std::string &prefix, int length);
void PrintContour(Shape & shape, std::string name, double max, double increment);
void Show(Shape & shape);
void PrintAlongLine(Shape & shape, const GeometryVector & axes, std::string name, double max, double increment);

bool compare(GeometryVector a, GeometryVector b) {
	assert(a.Dimension == b.Dimension);
	return a.x[a.Dimension - 1] < b.x[b.Dimension - 1];
}


double Sinc(double x) {
	if (fabs(x) < 1e-10)
		return 1.0;
	else
		return sin(x) / x;
}

double mt(const GeometryVector & k, double L) {
	double result = 1.0L;
	for (int i = 0; i < k.Dimension; i++)
		result *= (L*Sinc(0.5*k.x[i] * L));

	return result;
}

void test_triangles(std::ostream &, std::istream &);
void test_squares(std::ostream &, std::istream &);
void test_cubes(std::ostream &, std::istream &);
void test_formfactor(std::ostream &, std::istream &);

Polyhedron GetCube(double L);

typedef std::shared_ptr<Shape> Shape_;

Shape_ GenSquare(double L);
Shape_ GenCube(double L);



size_t Verbosity = 2;  int main(int argc, char ** argv) {
	std::string tempStr;
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	RandomGenerator rng(0);
	int num_threads = 6;
	omp_set_num_threads(num_threads);

	if(argc > 6){
		/* =====================================
		 *		Compute spectral density of dispersions.
		 * =====================================*/

		/*	argv[1]	load prefix
			argv[2]	number of configs
			argv[3]	save prefix
			argv[4] circular kmax
			argv[5] linear kmax
			argv[6] dk		*/
		
		std::string load_prefix(argv[1]);
		int num_configs = std::atoi(argv[2]);
		std::string save_prefix(argv[3]);
		double circular_k = std::atof(argv[4]);
		double linear_k = std::atof(argv[5]);
		double dk = std::atof(argv[6]);
		
		std::vector<std::string> comments;
		{
			char temp[300];
			sprintf(temp, "Load files: %s\t from 0 to %d\n", load_prefix.c_str(), num_configs);
			comments.emplace_back(temp);
			
			sprintf(temp, "Save the computed spectral density at %s\n", save_prefix.c_str());
			comments.emplace_back(temp);

			sprintf(temp, "Circular K_max = %0.4f\n", circular_k);
			comments.emplace_back(temp);

			sprintf(temp, "Linear K_max = %0.4f\n", linear_k);
			comments.emplace_back(temp);

			sprintf(temp, "dk = %0.4f\n", dk);
			comments.emplace_back(temp);
		}

		for(size_t i=0; i<comments.size(); i++)
			ofile << comments[i];

		comments.emplace_back("k\tspectral density\terror_x\terror_y");
		auto getc = [& load_prefix](size_t i) -> Dispersion{
			char name[200] = {};
			// load_prefix is a user printf pattern; it must use an int conversion (e.g. %d/%04d).
			snprintf(name, sizeof(name), load_prefix.c_str(), static_cast<int>(i));
			Dispersion c;
			ReadDispersion(c, name);
			return c;
		};
		
		std::vector<GeometryVector> result;
		IsotropicSpectralDensity_Dispersion(getc, num_configs, circular_k, linear_k, result, dk);
		WriteFunction(result, save_prefix, comments);
	}
	else if(argc == 2){
		std::string load_prefix(argv[1]);
		Dispersion c;
		ReadDispersion(c, load_prefix);
	}
	else{
		//Polyhedron cube = GetCube(1.0);
		//ofile << cube;
		//test_cubes(ofile, ifile);
		//test_triangles(ofile, ifile);
		if(true)
			test_formfactor(ofile, ifile);

		if(false){
			//a packing with a cube
			Shape_ temp = GenCube(1.0);
			ofile << *temp;

			temp->ShiftBy(GeometryVector(0.1,0.0,0.0));

			std::vector<GeometryVector> basis;
			basis.emplace_back(1.0,0.0,0.0);
			basis.emplace_back(0.0,1.0,0.0);
			basis.emplace_back(0.0,0.0,1.0);
			Dispersion packing_test(3, &basis[0], 1.0);
			packing_test.Insert(temp, GeometryVector(0,0.1,0));

			Shape_ temp2 = packing_test.GetShape(0); 
			ofile << * temp2;
		}

		if(false){
			ofile << "2D version!!\n";
			Shape_ temp = GenSquare(0.5);
			ofile << *temp;
			temp->ShiftBy(GeometryVector(0.1,0.0));

			std::vector<GeometryVector> basis;
			basis.emplace_back(1.0,0.0);
			basis.emplace_back(0.0,1.0);
			Dispersion packing_test(2, &basis[0], 1.0);
			packing_test.Insert(temp, GeometryVector(0,0.1));

			Shape_ temp2 = packing_test.GetShape(0); 
			ofile << "Before transform:\n";
			ofile << * temp2;

			temp2->ShiftBy(GeometryVector(1.0,1.0));
			ofile << "After transform 1:\n";
			ofile << * temp2;
			
			temp2 = packing_test.GetShape(0); 
			ofile << "Before transform 2:\n";
			ofile << * temp2;

			std::string name = "square_output" + std::string(".txt");
			std::fstream output(name.c_str(), std::fstream::out);

			output << * temp2;
			output.close();

		}

		if(false)
		{
			ofile << "=========================\n";
			ofile << "	Test 2D Disperion I/O test \n";
			ofile << "=========================\n";
			Shape_ temp = GenSquare(0.1);
			
			std::vector<GeometryVector> basis;	
			basis.emplace_back(1.0,0.0);
			basis.emplace_back(0.0,1.0);
			Dispersion packing_test(2, &basis[0], 1.0);
			packing_test.Insert(temp, GeometryVector(0,0.1));
			packing_test.Insert(temp, GeometryVector(0.1,0.3));

			//WriteDispersion(packing_test, "test_2D");


	//		Dispersion test;
	//		ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/test3");
	//		WriteDispersion(test, "/home/jaeukk/Projects/OrderMetrics/test4");
	//		ofile << * (test.GetShape(0));
	//		ofile << * (test.GetShape(0));
		}

		if(false)
		{
			ofile << "=========================\n";
			ofile << "	Test 3D Disperion I/O test \n";
			ofile << "=========================\n";
			Shape_ temp = GenCube(0.1);
			
			std::vector<GeometryVector> basis;	
			basis.emplace_back(1.0,0.0,0.0);
			basis.emplace_back(0.0,1.0,0.0);
			basis.emplace_back(0.0,0.0,1.0);
			Dispersion packing_test(3, &basis[0], 1.0);

			packing_test.Insert(temp, GeometryVector(0.0, 0.1, 0.2));
			packing_test.Insert(temp, GeometryVector(0.1, 0.3, 0.0));

			WriteDispersion(packing_test, "test_3D");

			Dispersion test;

			{
				std::string name("test_3D.Dispersion");
				std::fstream Ifile (name.c_str(), std::fstream::in);

				std::string dump;

				while (Ifile.eof() == false){
					std::getline(Ifile, dump);
					ofile << dump <<std::endl;			
				}
			}

			ReadDispersion(test, "test_3D");

			ofile << test.GetCartesianCoordinates(0);
			ofile << test.GetCartesianCoordinates(1);

		}

		//test_squares(ofile, ifile);
		//test_cubes(ofile, ifile);

	}



	return 0;
}

void test_formfactor(std::ostream & ofile, std::istream & ifile){
	
	/* compute formfactors of a random triangle. */
	std::vector<GeometryVector> v_tr;
	v_tr.emplace_back(-1.0,1.0);
	v_tr.emplace_back(3.0,-1.0);
	v_tr.emplace_back(0.1,3.0);

	Polygon triangle(v_tr);
	
	ofile << triangle;
	
	/* a list of 2D random vectors.*/
	std::vector<GeometryVector> ks;
	ks.emplace_back(4.6235983,2.2084670);
	ks.emplace_back(0.8336992,1.1380447);
	ks.emplace_back(12.6005627,7.7597935);
	ks.emplace_back(3.2477628,8.1011854);
	ks.emplace_back(8.7447961,9.1326413);
	ks.emplace_back(1.6967074,3.4584658);
	ks.emplace_back(9.6394202,2.7467691);
	ks.emplace_back(16.5577468,18.2850515);
	ks.emplace_back(16.0122479,5.1673137);
	ks.emplace_back(4.1876083,12.2343984);
	ks.emplace_back(5.3404423,8.9151078);
	ks.emplace_back(11.6527179,5.0938620);
	ks.emplace_back(12.4160427,5.0924797);
	ks.emplace_back(2.0118595,5.1645708);
	ks.emplace_back(14.5241260,17.1944929);
	ks.emplace_back(7.8611688,9.2116279);
	ks.emplace_back(13.7056868,8.3468536);
	ks.emplace_back(8.1586370,4.3567845);
	ks.emplace_back(3.7621995,15.6497999);
	ks.emplace_back(12.9048235,11.3182128);
	ks.emplace_back(1.6243016,10.4284878);
	ks.emplace_back(4.9809728,19.7295559);
	ks.emplace_back(5.6472651,16.7331007);
	ks.emplace_back(14.3817718,15.6148567);
	ks.emplace_back(15.0179111,13.6630471);
	ks.emplace_back(6.1969314,19.7132957);
	ks.emplace_back(16.6575406,19.8772827);
	ks.emplace_back(12.8507761,9.2404402);
	ks.emplace_back(11.8366760,12.1753568);
	ks.emplace_back(16.6853618,10.9591459);
	ks.emplace_back(19.6507115,3.0527435);
	ks.emplace_back(8.1704238,14.3970467);
	ks.emplace_back(13.0642585,19.5851434);
	ks.emplace_back(6.5218854,6.0758404);
	ks.emplace_back(8.9315610,2.5109316);
	ks.emplace_back(18.4243708,12.6174785);
	ks.emplace_back(0.6140416,7.5409491);
	ks.emplace_back(6.7500367,17.0803252);
	ks.emplace_back(0.9615481,14.4410247);
	ks.emplace_back(8.1461685,8.4513576);
	ks.emplace_back(7.4416962,1.0872819);
	ks.emplace_back(7.9612719,1.4114950);
	ks.emplace_back(17.5286884,2.1217478);
	ks.emplace_back(19.3040106,6.1587468);
	ks.emplace_back(16.7018220,19.5475775);
	ks.emplace_back(12.3105698,14.2021371);
	ks.emplace_back(9.7252915,5.3648251);
	ks.emplace_back(17.7741883,9.6281073);
	ks.emplace_back(16.9552726,7.4955174);
	ks.emplace_back(12.3372856,6.7480653);

	ofile << "np.array([\n";
	char temp[100]={};
	for (size_t i=0; i<ks.size(); i++){
		std::complex<double> value = triangle.GetFormFactor(ks[i]);
		sprintf(temp, "[%0.7f, %0.7f],\n", value.real(), value.imag());
		ofile << temp;
	}
	ofile << "]";
}

Shape_ GenSquare(double L){
	std::vector<GeometryVector> vs;
/* 	vs.emplace_back(0.5, -0.5);
	vs.emplace_back(0.5, 0.5);
	vs.emplace_back(-0.5, 0.5);
	vs.emplace_back(-0.5, -0.5);
 */
	vs.emplace_back(-0.5, -0.5);
	vs.emplace_back(-0.5, 0.5);
	vs.emplace_back(0.5, 0.5);
	vs.emplace_back(0.5, -0.5);

	Polygon temp(vs);
	temp.Scale(L);
	return std::shared_ptr<Shape> (new Polygon (temp));
}

Shape_ GenCube(double L){
	GeometryVector v1(1.0,0.0,0.0);
	GeometryVector v2(0.0,1.0,0.0);
	GeometryVector v3(0.0,0.0,1.0);
	v1 = L*v1; v2 = L*v2; v3 = L*v3;
	GeometryVector x0 = 0.5*(v1 + v2 + v3);

	std::vector<GeometryVector> vertices;	
	std::vector<std::vector<size_t>> index;
	/*0*/	vertices.emplace_back(0.0*v1 - x0);
	/*1*/	vertices.emplace_back(v2 - x0);
	/*2*/	vertices.emplace_back(v1 + v2 - x0);
	/*3*/	vertices.emplace_back(v1 - x0);
	/*4*/	vertices.emplace_back(v3 - x0);
	/*5*/	vertices.emplace_back(v1 + v3 - x0);
	/*6*/	vertices.emplace_back(v2 + v3 - x0);
	/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
	{
		index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
		index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
		index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
		index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
		index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
		index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
	}
	Shape_ temp = std::shared_ptr<Shape> (new Polyhedron(vertices, index));
	return temp;
}

void test_triangles(std::ostream & ofile, std::istream & ifile){
	
	ofile << "triangles...\n";

	std::vector<GeometryVector> vs;
	vs.emplace_back(-0.5, -0.5/sqrt(3.0));
	vs.emplace_back(+0.5, -0.5/sqrt(3.0));
	vs.emplace_back(0.0, 1.0/sqrt(3.0));

	Polygon triangle(vs);
	ofile << "first triangle \n";
	ofile << triangle;

	RandomGenerator rng;
	for (int i=0; i< 10; i++){
		GeometryVector displacement (rng.RandomDouble(), rng.RandomDouble());
		std::vector<GeometryVector> test_vs = vs;

		ofile <<"-----------"<< i<<"th triangle -----------\n";
		for (auto x = test_vs.begin(); x!= test_vs.end(); x++){
			*x = *x + displacement;
			ofile << *x;
		}
		ofile << "\n";
		Polygon test_triangle (test_vs);
		//ofile << test_triangle;
		ofile << "centroid = " << test_triangle.GetCenterOfMass();
		ofile << "shift = " << displacement;
		ofile << "\n\n";

	}

	{
		GeometryVector displacement (0.0, 0.5/sqrt(3.0));
		std::vector<GeometryVector> test_vs = vs;

		for (auto x = test_vs.begin(); x!= test_vs.end(); x++){
			*x = *x + displacement;
		}
		Polygon test_triangle (test_vs);

		ofile << "11th triangle generated\n";
		ofile << test_triangle;
	}
}

void test_squares(std::ostream & ofile, std::istream & ifile){
	
	ofile << "squares...\n";

	std::vector<GeometryVector> vs;
/* 	vs.emplace_back(0.5, -0.5);
	vs.emplace_back(0.5, 0.5);
	vs.emplace_back(-0.5, 0.5);
	vs.emplace_back(-0.5, -0.5);
 */
	vs.emplace_back(-0.5, -0.5);
	vs.emplace_back(-0.5, 0.5);
	vs.emplace_back(0.5, 0.5);
	vs.emplace_back(0.5, -0.5);

	Polygon square(vs);
	ofile << "first square \n";
	ofile << square;

	{
		std::fstream of ("square1.txt", std::fstream::out);
		of << square;
	}

	{
		ofile << "Check filestream\n";
		std::fstream iF ("square1.txt", std::fstream::in);
		std::string dump;

		while (std::getline(iF, dump)){
			ofile << dump;
			ofile << "\n";
		}
		ofile << "\n\n\n";
	}

	{
		ofile << "Read square1.txt file\n";
		std::fstream iF ("square1.txt", std::fstream::in);
		char dump[255];
		//std::getline(iF, dump);
		iF >> dump;
		Polygon what;
		if (strcmp(dump, "[Polygon]")==0){
			what = ReadPolygon(iF);
			ofile <<"Succeed to read\n";
			ofile << what;
		}
		else{
			ofile << "falied to read\n";
		}
	}


	ofile << "set orientation positive\n";
	square.SetPositiveOrientation();
	ofile << square;


	for (double kx =1e-9; kx<10; kx+=0.1){
		GeometryVector k (kx, 1.0);
		double analytic = mt(k, 1.0); 
		double numeric_r = square.GetFormFactor(k).real();
		double numeric_i = square.GetFormFactor(k).imag();
		ofile << kx <<"\t"<< analytic <<"\t("<< numeric_r<<", " << numeric_i  << ")\t"<< abs(analytic - numeric_r) <<"\n";
	} 
}



Polyhedron GetCube(double L){
	GeometryVector v1(1.0,0.0,0.0);
	GeometryVector v2(0.0,1.0,0.0);
	GeometryVector v3(0.0,0.0,1.0);
	v1 = L*v1; v2 = L*v2; v3 = L*v3;
	GeometryVector x0 = 0.5*(v1 + v2 + v3);

	std::vector<GeometryVector> vertices;	
	std::vector<std::vector<size_t>> index;
	/*0*/	vertices.emplace_back(0.0*v1 - x0);
	/*1*/	vertices.emplace_back(v2 - x0);
	/*2*/	vertices.emplace_back(v1 + v2 - x0);
	/*3*/	vertices.emplace_back(v1 - x0);
	/*4*/	vertices.emplace_back(v3 - x0);
	/*5*/	vertices.emplace_back(v1 + v3 - x0);
	/*6*/	vertices.emplace_back(v2 + v3 - x0);
	/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
	{
		index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
		index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
		index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
		index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
		index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
		index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
	}

	return Polyhedron(vertices, index);
}

void test_cubes(std::ostream & ofile, std::istream & ifile){
	
	ofile << "cubes...\n";

	GeometryVector v1(1.0,0.0,0.0);
	GeometryVector v2(0.0,1.0,0.0);
	GeometryVector v3(0.0,0.0,1.0);
	GeometryVector x0 = 0.5*(v1 + v2 + v3);

	std::vector<GeometryVector> vertices;	
	std::vector<std::vector<size_t>> index;
	/*0*/	vertices.emplace_back(0.0*v1 - x0);
	/*1*/	vertices.emplace_back(v2 - x0);
	/*2*/	vertices.emplace_back(v1 + v2 - x0);
	/*3*/	vertices.emplace_back(v1 - x0);
	/*4*/	vertices.emplace_back(v3 - x0);
	/*5*/	vertices.emplace_back(v1 + v3 - x0);
	/*6*/	vertices.emplace_back(v2 + v3 - x0);
	/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
	{
		index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
		index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
		index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
		index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
		index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
		index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
	}

	ofile <<"Input vertices \n";
	for (size_t i=0; i<vertices.size(); i++)
		ofile << vertices[i];

	ofile <<"Generate a cube \n";
	ofile << "The centroid of each face is \n";
	for (int i=0; i< 6; i++){
		GeometryVector x0 = vertices[index[i][0]];
		for (int j=1; j<4; j++)
			x0 = x0 + vertices[index[i][j]];

		x0 = 0.25 * x0;
		ofile << "face " << i << ": "<< x0;
	}

	Polyhedron cube1 (vertices, index);

	ofile << cube1;


	ofile << "centroid = "<<cube1.GetCenterOfMass() <<"\n";
	ofile << "volume = "<< cube1.GetVolume() <<"\n"; 

	ofile << "Test form factor functions!!\n";

	auto mtilde = [](const GeometryVector & k) -> double{
	 	double result = 1.0;
	 	for (int i=0; i<k.Dimension; i++){
	 	
	 		double x = k.x[i] * 0.5;
	 		result *= Sinc(x);
	 	}

	 	return result;		
	};

	for (double k=0.0; k<5.0; k+=0.01){
	 	GeometryVector K(k, 0.0, 0.0);
	 	double val = cube1.GetFormFactor(K).real();
	 	double val2 = mtilde(K); 

	 	GeometryVector temp(k, val, val2, val-val2);
	 	ofile << temp;
	}



}

void WriteToContourf(const std::vector<GeometryVector> & vec_k, const std::vector<GeometryVector> &s_k, const std::string &prefix, int length) {

	if (vec_k.size() == s_k.size()) {
		if (vec_k[0].Dimension == 2) {
			std::string file = prefix + "X_.txt";
			std::ofstream X_(file);
			file = prefix + "Y_.txt";
			std::ofstream Y_(file);
			file = prefix + "S_k.txt";
			std::ofstream S(file);

			for (size_t i = 0; i < vec_k.size(); i++) {
				X_ << vec_k[i].x[0] << "\t";
				Y_ << vec_k[i].x[1] << "\t";
				S << s_k[i].x[0] << "\t";
				if ((i + 1) % length == 0) {
					X_ << "\n";
					Y_ << "\n";
					S << "\n";
				}
			}
			X_.close();	Y_.close(); S.close();
		}
		else if (vec_k[0].Dimension == 3) {
			std::string file = prefix + "X_.txt";
			std::ofstream X_(file);
			file = prefix + "Y_.txt";
			std::ofstream Y_(file);
			file = prefix + "Z_.txt";
			std::ofstream Z_(file);
			file = prefix + "S_k.txt";
			std::ofstream S(file);

			for (size_t i = 0; i < vec_k.size(); i++) {
				X_ << vec_k[i].x[0] << "\t";
				Y_ << vec_k[i].x[1] << "\t";
				Z_ << vec_k[i].x[2] << "\t";
				S << s_k[i].x[0] << "\t";
				if ((i + 1) % length == 0) {
					X_ << "\n";
					Y_ << "\n";
					Z_ << "\n";
					S << "\n";
				}
			}
			X_.close();	Y_.close(); Z_.close();  S.close();
		}
	}
	else {
		std::cout << " wrong input \n";
	}



}

void PrintContour(Shape & shape, std::string name, double max, double increment) {
	std::vector<GeometryVector> ks, formfactor_Re, formfactor_Im, formfactor_Abs;
	double dk = increment;
	double kMax = max;
	int d = shape.GetDimension();
	int size = (int)ceil(pow(kMax / dk, d));
	int length = 0;
	ks.reserve(size);
	formfactor_Re.reserve(size);
	formfactor_Im.reserve(size);
	formfactor_Abs.reserve(size);
	if (d == 2) {
		for (double kx = -kMax; kx < kMax; kx += dk) {
			length++;
			for (double ky = -kMax; ky < kMax; ky += dk) {
				ks.emplace_back(GeometryVector(kx, ky));
			}
		}
	}
	else if (d == 3) {
		for (double kz = 0.0; kz < kMax; kz += dk) {
			length++;
			for (double ky = 0.0; ky < kMax; ky += dk)
				for (double kx = 0.0; kx < kMax; kx += dk)
					ks.emplace_back(GeometryVector(kx, ky, kz));

		}
	}


	for (auto v = ks.begin(); v != ks.end(); v++) {
		std::complex<double> Ik = shape.GetFormFactor(*v);
		formfactor_Re.emplace_back(GeometryVector(std::real(Ik)));
		formfactor_Im.emplace_back(GeometryVector(std::imag(Ik)));
		formfactor_Abs.emplace_back(GeometryVector(std::abs(Ik)));
	}


	WriteToContourf(ks, formfactor_Re, (name + std::string("_FormFactor_real")).c_str(), length);
	WriteToContourf(ks, formfactor_Im, (name + std::string("_FormFactor_imag")).c_str(), length);
	WriteToContourf(ks, formfactor_Abs, (name + std::string("_FormFactor_abs")).c_str(), length);
}

void Show(Shape & shape) {
	std::cout << "Vertices: \n";
	shape.Print();
	std::cout << "Volume: \n";
	std::cout << shape.GetVolume() << "\n";
	std::cout << "Surface Area: \n";
	std::cout << shape.GetSurfaceArea() << "\n";
	std::cout << "Rmax: \n";
	std::cout << shape.GetMaxRadius() << "\n\n\n";
}

void PrintAlongLine(Shape & shape, const GeometryVector & axes, std::string name, double max, double increment) {
	GeometryVector dk = 1.0 / sqrt(axes.Modulus2())*axes;
	std::vector<GeometryVector> result;
	for (double K = 0.0; K < max; K += increment) {
		GeometryVector k = K * dk;
		std::complex<double> sk = shape.GetFormFactor(k);
		result.emplace_back(GeometryVector(K, std::real(sk), std::imag(sk), std::abs(sk)));
	}
	WriteFunction(result, (name + std::to_string(dk.x[0]) + "_" + std::to_string(dk.x[1]) + "_" + std::to_string(dk.x[2])).c_str());
}