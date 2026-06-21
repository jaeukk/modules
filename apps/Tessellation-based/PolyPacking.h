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

//adjust radii of monodisperse HS packings s.t. sphere volume/Voronoi cell = targetPhi
//Volumes[i]: volume of a Voronoi cell of ith particle in config.
//resultConfig: resulting configuration, data will be reset
//return the maximum packing fraction achievable without overlapping.
double MakeMono2Poly(const SpherePacking & config, SpherePacking & resultConfig, double targetPhi, std::vector<double> &Volumes, std::vector<int> * numVertices = nullptr);

//build an RSA substrate of the requested size, then MakeMono2Poly it.
double MakeMono2Poly(DimensionType d, size_t numParticles, double boxSize, double phi0, SpherePacking &resultConfig, double targetPhi, std::vector<double> &Volumes, int seed = 0);

//parallel variant (uses voro++ when USE_VORO, otherwise per-cell GetVoronoiCell).
double MakeMono2Poly_parallel(const SpherePacking & config, SpherePacking & resultConfig, double targetPhi, std::vector<double> &Volumes, int numThreads, std::vector<int> * numVertices = nullptr);

#endif
