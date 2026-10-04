#ifndef RANDOMGENERATOR_INCLUDED
#define RANDOMGENERATOR_INCLUDED
#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <random>
#include <functional>
#include "etc.h"          // DimensionType

class GeometryVector;     // forward decl: returned by value by RandomUnitVector below

/** \brief Random-number utilities for uniform, Gaussian, and Poisson draws. */
class RandomGenerator
{
private:
	boost::mt19937 gen;
	boost::mt19937::result_type min, max;
	/** \brief Cache the engine's min and max output values. */
	void initialize(void);
public:
	/** \brief Seed the underlying generator with a specific integer.
	 *  \param s  Seed value for the Mersenne-Twister engine. */
	void seed(int s);

	/** \brief Construct a generator seeded from current time. */
	RandomGenerator();

	/** \brief Construct a generator with an explicit seed value.
	 *  \param Seed  Seed value for the underlying engine. */
	RandomGenerator(int Seed);

	/** \brief Draw a uniform deviate in [0,1). */
	double RandomDouble(void);

	/** \brief Draw a deviate using an inverse cumulative distribution. 
	 *  \param inverseCp Inverse CDF mapping from uniform [0,1) to the target distribution. */
	double RandomDouble(const std::function<double(double)> &inverseCp);

	/** \brief Draw a Gaussian deviate with mean 0 and supplied standard deviation.
	 *  \param std  Standard deviation of the normal distribution. */
	double RandomDouble_normal(double std) {
		boost::normal_distribution<> nd(0.0, std);
		boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(this->gen, nd);
		return var_nor();
	}

	/** \brief Draw a Poisson deviate with given mean.
	 *  \param mean  Mean (lambda) of the Poisson distribution. */
	int PoissonRNG(double mean); 

	/** \brief Draw an integer from the full underlying engine range. */
	inline int RandomInt(void) {return this->gen();}
};

/** \brief Random unit vector whose direction is statistically isotropic.
 *  Exact per-dimension sampling (closed-form for d=1,2,3; Gaussian-normalized
 *  otherwise).
 *  \param d Space dimension.  \param gen Random number generator. */
GeometryVector RandomUnitVector(DimensionType d, RandomGenerator & gen);
GeometryVector RandomUnitVector(DimensionType d);

/** \brief Random unit vector written into a raw buffer (rejection sampling).
 *  Unlike RandomUnitVector (which returns a GeometryVector via exact per-
 *  dimension sampling), this fills result[0..Dimension-1] in place by drawing
 *  uniformly in the cube, rejecting points outside the unit ball, then
 *  normalizing.
 *  \param Dimension Dimension of the target space.
 *  \param gen Random number generator to sample from.
 *  \param result Output buffer (length Dimension) to receive the components. */
void GetRandomVector(int Dimension, RandomGenerator & gen, double * result);

/** \brief Fill result with a random vector uniformly distributed in the cube [-1,1]^d.
 *  \param d       Dimension of the vector to generate.
 *  \param gen     Random number generator to sample from.
 *  \param result  Output vector to receive the random components. */
void RandomVectorInCube(DimensionType d, RandomGenerator & gen, GeometryVector & result);


#endif
