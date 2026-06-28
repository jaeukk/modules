#include "RandomGenerator.h"
#include "GeometryVector.h"
#include <cmath>
#include <ctime>
void RandomGenerator::initialize(void)
	{
		this->min=gen.min();
		this->max=gen.max();
	}
void RandomGenerator::seed(int s)
	{
		this->gen.seed(s);
	}
RandomGenerator::RandomGenerator()
	{
		this->initialize();
		this->seed(std::time(NULL));
	}
RandomGenerator::RandomGenerator(int Seed)
	{
		this->initialize();
		this->seed(Seed);
	}
double RandomGenerator::RandomDouble(void)//[0,1)
	{
		return (double)(gen()-min)/((double)(max-min)+1);
	}
double RandomGenerator::RandomDouble(const std::function<double(double)> &inverseCp)
	{
		return inverseCp(RandomDouble());
	}

int RandomGenerator::PoissonRNG(double mean){
	if (mean <= 0.0) {
		return 0;
	}
	std::poisson_distribution<int> dist(mean);
	return dist(this->gen);
}

// RandomUnitVector: documented on the declaration in RandomGenerator.h.
GeometryVector RandomUnitVector(DimensionType d, RandomGenerator & gen) {
	GeometryVector result(static_cast<int>(d));
	double theta = 0;
	double z = 0;
	switch (d)
	{
	case 1:
		result.x[0] = (gen.RandomDouble() < 0.5)? -1:1;
		break;
	case 2:
		theta = 2.0*pi*gen.RandomDouble();
		result.x[0] = cos(theta);
		result.x[1] = sin(theta);
		break;
	case 3:
		//Equal-area projection was implemented.
		theta = 2.0*pi*gen.RandomDouble();
		z = 2.0*(gen.RandomDouble() - 0.5);
		result.x[0] = sqrt(1.0 - z*z)*cos(theta);
		result.x[1] = sqrt(1.0 - z*z)*sin(theta);
		result.x[2] = z;
		break;
	default:
		for (int i = 0; i < d; i++)
			result.x[i] = gen.RandomDouble_normal(1.0);// 2.0*(gen.RandomDouble() - 0.5);

		double mod = 1.0 / sqrt(result.Dot(result));
		for (int i = 0; i < d; i++)
			result = mod * result;

		break;
	}
	return result;
}

GeometryVector RandomUnitVector(DimensionType d) {
	RandomGenerator gen;
	return RandomUnitVector(d, gen);

}

//random unit vector written into a raw buffer (rejection sampling).
//differs from RandomUnitVector (raw double* output + rejection method).
void GetRandomVector(int Dimension, RandomGenerator & gen, double * result)
{
	double r;
	do
	{
		r=0;
		for(int i=0; i<Dimension; i++)
		{
			result[i]=2*gen.RandomDouble()-1;
			r+=result[i]*result[i];
		}
	}
	while(r>1);
	for(int i=0; i<Dimension; i++)
		result[i]/=std::sqrt(r);

	return;
}

//random vector uniformly distributed in the cube [-1,1]^d.
void RandomVectorInCube(DimensionType d, RandomGenerator & gen, GeometryVector & result)
{
	result.SetDimension(d);
	for (DimensionType i = 0; i < d; i++)
		result.x[i] = 2.0*(gen.RandomDouble() - 0.5);
}
