/**
 *	Editor	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: June 2024 */


/**
 * \file ContactStatistics.h
 * \brief Header file for a
 */

#ifndef Z_INCLUDED
#define Z_INCLUDED

#include "PeriodicCellList.h"
#include <vector>
#include <omp.h>
#include <math.h>


/**
 * \brief Measure a distribution of coordination numbers in a range between MinDistance and MaxDistance.
 * \param GetConfigsFunction  A lambda function to generate the point configurations.
 * \param NumConfigs          The number of point configurations.
 * \param MinDistance         Neighbors closer than this distance are not counted as contacts.
 * \param MaxDistance         Neighbors farther than this distance are not counted as contacts.
 * \param Result              Output distribution P(Z) of coordination numbers as (Z, P(Z), error in P(Z), 0).
 */
void ContactNumbers(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs,
	double MinDistance, double MaxDistance, std::vector<GeometryVector> & Result);

/**
 * \brief Measure cumulative coordination number as a function of distance.
 * \param GetConfigsFunction  A lambda function to generate the point configurations.
 * \param NumConfigs          The number of point configurations.
 * \param r_min               Smallest distance at which Z(r) is reported.
 * \param r_max               Largest distance up to which Z(r) is computed.
 * \param dr                  Bin width in distance (initial bin width when logscale is used).
 * \param Result              Output cumulative coordination number as (r, Z(r), error in Z(r)).
 * \param logscale            If true, use logarithmically growing (geometric) bin widths. (CLAUDE)
 * \return                    The actual minimum distance used as r_min.
 */
double CumulativeCoordinationNumber(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double r_min, double r_max, double dr, std::vector<GeometryVector> &Result, bool logscale=false);

#endif