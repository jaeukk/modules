/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file VolumeFractionVariance.h
 *	\brief Local volume-fraction variance sigma_V^2(R) of sphere packings.
 *
 *	Split out of NumberVariance.{h,cpp}: while the local *number* variance counts
 *	points in a window, the local *volume-fraction* variance measures the
 *	fluctuation of the covered (particle-phase) volume fraction inside a spherical
 *	window. Ported from the full top-level NumberVariance; the per-window
 *	local-volume helpers (I_1 / I_config) are file-local in the .cpp. */


#ifndef __VOLUMEFRACTIONVARIANCE_H__
#define __VOLUMEFRACTIONVARIANCE_H__

#include <vector>
#include <functional>

/* header files in $(cores) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>


/** \brief Monte-Carlo local volume-fraction variance sigma_V^2(R) of sphere packings,
 *  using prescribed sampling-window centers.
 *  \param GetConfig       Callable returning sphere packing i.
 *  \param numConfig       Number of configurations in the ensemble.
 *  \param[out] result     result[i].x = { R, sigma_V^2(R), error, running average }.
 *  \param Rmax            Largest window radius.
 *  \param dR              Resolution in R.
 *  \param samplingPoints  Window centers in relative coordinates.
 *  \param num_threads     Number of threads used (default 1). */
void MC_vfv(const std::function<SpherePacking(size_t i)> & GetConfig, size_t numConfig, std::vector<GeometryVector> & result, double Rmax, double dR, std::vector<GeometryVector> & samplingPoints, size_t num_threads = 1);

/** \brief Monte-Carlo local volume-fraction variance sigma_V^2(R) of sphere packings,
 *  using num randomly generated sampling-window centers.
 *  \param GetConfig       Callable returning sphere packing i.
 *  \param numConfig       Number of configurations in the ensemble.
 *  \param[out] result     result[i].x = { R, sigma_V^2(R), error, running average }.
 *  \param Rmax            Largest window radius.
 *  \param dR              Resolution in R.
 *  \param numSampling     Number of random window centers per configuration.
 *  \param num_threads     Number of threads used (default 1). */
void MC_vfv(const std::function<SpherePacking(size_t i)> & GetConfig, size_t numConfig, std::vector<GeometryVector> & result, double Rmax, double dR, size_t numSampling, size_t num_threads = 1);


#endif
