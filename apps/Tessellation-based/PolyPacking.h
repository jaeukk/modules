#ifndef TESS_POLYPACKING_H
#define TESS_POLYPACKING_H

/** \file PolyPacking.h
 *  \brief Turn a monodisperse hard-sphere packing into a polydisperse one by
 *  growing each sphere to a target fraction of its Voronoi-cell volume.
 *
 *  The `MakeMono2Poly*` family, relocated out of the legacy grab-bag
 *  etc1.{h,cpp} (its only tiling consumer is TilingProblemCLI, which calls
 *  `MakeMono2Poly_parallel`).  The cell-volume computations rely on CGAL /
 *  voro++ (Voronoi.h), so the Voronoi-based overloads are only *defined* when
 *  `USE_CGAL` is set; the RSA-substrate overload is always available but
 *  delegates to a CGAL-guarded one.
 */

#include <vector>
#include <PeriodicCellList.h>

/** \brief adjust radii of monodisperse HS packings s.t. sphere volume/Voronoi cell = targetPhi
 *  \param config        Input monodisperse hard-sphere packing.
 *  \param resultConfig  Resulting configuration, data will be reset.
 *  \param targetPhi     Target ratio of sphere volume to Voronoi-cell volume.
 *  \param[out] Volumes  Volumes[i]: volume of a Voronoi cell of ith particle in config.
 *  \param numVertices   Optional output: number of Voronoi-cell vertices per particle.
 *  \return              the maximum packing fraction achievable without overlapping. */
double MakeMono2Poly(const SpherePacking & config, SpherePacking & resultConfig, double targetPhi, std::vector<double> &Volumes, std::vector<int> * numVertices = nullptr);

/** \brief build an RSA substrate of the requested size, then MakeMono2Poly it.
 *  \param d             Spatial dimension.
 *  \param numParticles  Number of particles in the RSA substrate.
 *  \param boxSize       Linear size to which the substrate is rescaled.
 *  \param phi0          Packing fraction of the RSA substrate.
 *  \param resultConfig  Resulting configuration, data will be reset.
 *  \param targetPhi     Target ratio of sphere volume to Voronoi-cell volume.
 *  \param[out] Volumes  Per-particle Voronoi-cell volumes.
 *  \param seed          Offset added to the random-generator seed.
 *  \return              the maximum packing fraction achievable without overlapping. */
double MakeMono2Poly(DimensionType d, size_t numParticles, double boxSize, double phi0, SpherePacking &resultConfig, double targetPhi, std::vector<double> &Volumes, int seed = 0);

/** \brief parallel variant (uses voro++ when USE_VORO, otherwise per-cell GetVoronoiCell).
 *  \param config        Input monodisperse hard-sphere packing.
 *  \param resultConfig  Resulting configuration, data will be reset.
 *  \param targetPhi     Target ratio of sphere volume to Voronoi-cell volume.
 *  \param[out] Volumes  Per-particle Voronoi-cell volumes.
 *  \param numThreads    Number of threads used for the cell computation.
 *  \param numVertices   Optional output: number of Voronoi-cell vertices per particle.
 *  \return              the maximum packing fraction achievable without overlapping. */
double MakeMono2Poly_parallel(const SpherePacking & config, SpherePacking & resultConfig, double targetPhi, std::vector<double> &Volumes, int numThreads, std::vector<int> * numVertices = nullptr);

#endif
