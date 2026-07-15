/**
 *	Editor	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: Nov. 2023 */

/** \file LinealPathFunction.h 
 * \brief Header file for computing the lineal path function from two-phase media consisting of (penetrable) spheres.
*/

#ifndef __LINEALPATHFUNCTION_H__
#define __LINEALPATHFUNCTION_H__

#include <etc.h>
#include <PeriodicCellList.h>
#include <GeometryVector.h>


/** \brief Compute the chord-length distribution from a single distribution of Spheres.
 *  \param GetC	A function to load SpherePacking objects for analysis.
 *  \param num_c	Number of configurations.
 *  \param p0 	chord-length distribution in matrix phase
 *  \param p1 	chord-length distribution in particle phase
 *  \param Rmax maximum line length
 *  \param dR 	increment in line length
 *  \param num_sample_lines 	The number of random lines
 *  \param seed random seeds for generating lines. */
void ChordLengthDistribution(const std::function<const SpherePacking(size_t i)> & GetC, size_t num_c, std::vector<GeometryVector> & p0, std::vector<GeometryVector> & p1, double Rmax, double dR, size_t num_sample_lines, int seed = 0);

/** \brief A function to compute the lineal path lengths from a distribution of spheres
 *  \param c 	A SpherePacking object for analysis (must be optimized before input)
 *  \param lengths0	A vector of lineal segment lengths in a matrix phase (in the unit of dL)
 *  \param lengths1	A vector of lineal segment lengths in a particle phase (in the unit of dL)
 *  \param Lmax 	Sizes of the maximum lineal path
 *  \param dL 		Increment in lengths
 *  \param num_sample_lines 	The number of sampling lines
 *  \param rng 	Random number generator for sampling lines. */
void GetLinealSegments(const SpherePacking & c, std::vector<size_t> & lengths0, std::vector<size_t> & lengths1, double Lmax, double dL, size_t num_sample_lines, RandomGenerator & rng);


/** \brief Compute the lineal-path function via the MC method.
 *  \param GetC	A function to load SpherePacking objects for analysis.
 *  \param num_c	Number of configurations.
 *  \param L0 	lineal-path function in matrix phase
 *  \param L1 	lineal-path function in particle phase
 *  \param Rmax maximum line length
 *  \param dR 	increment in line length
 *  \param num_sample_lines 	The number of random lines
 *  \param seed random seeds for generating lines. */
void LinealPathFunction(const std::function<const SpherePacking(size_t i)> & GetC, size_t num_c, std::vector<GeometryVector> & L0, std::vector<GeometryVector> & L1, double Rmax, double dR, size_t num_sample_lines, int seed = 0);


#endif // __LINEALPATHFUNCTION_H__