/**
 * \file smoke.cpp
 * \brief Minimal self-contained smoke test / usage example for the `core` module.
 *
 * Demonstrates the three foundational types every other module builds on:
 *   - GeometryVector       : a point/vector in up to MaxDimension Euclidean dims
 *   - PeriodicCellList<T>  : a periodic point pattern with O(1) neighbour queries
 *   - RandomGenerator      : reproducible Mersenne-Twister RNG
 *
 * Build & run:  make && ./smoke.out
 */

#include <cstdio>
#include <cmath>
#include "etc.h"
#include "GeometryVector.h"
#include "PeriodicCellList.h"
#include "RandomGenerator.h"

/** \brief Smoke-test entry point: build a periodic point pattern and count neighbours.
 *  \return  0 on success. */
int main()
{
    const DimensionType dim = 3;
    const double L = 10.0;

    // Cubic periodic box: basis vectors are L * e_i.
    GeometryVector basis[::MaxDimension];
    for (DimensionType i = 0; i < dim; ++i) {
        basis[i].SetDimension(dim);
        for (DimensionType j = 0; j < dim; ++j)
            basis[i].x[j] = (i == j) ? L : 0.0;
    }
    PeriodicCellList<Empty> cfg(dim, basis, /*CellSize*/ 1.0);

    // Insert N random points with a seeded RNG (reproducible run-to-run).
    RandomGenerator rng(12345);
    const size_t N = 200;
    for (size_t n = 0; n < N; ++n) {
        GeometryVector rel;
        rel.SetDimension(dim);
        for (DimensionType i = 0; i < dim; ++i)
            rel.x[i] = rng.RandomDouble();        // fractional coords in [0,1)
        cfg.Insert(Empty(), rel);
    }

    // Count neighbours of particle 0 within radius R -- exercises the cell list.
    const double R = 2.0;
    size_t neighbours = 0;
    cfg.IterateThroughNeighbors(cfg.GetRelativeCoordinates(0), R,
        [&](const GeometryVector & /*shift*/, const GeometryVector & /*latShift*/,
            const signed long * /*periodicShift*/, const size_t source) {
            if (source != 0) ++neighbours;        // exclude the particle itself
        });

    std::printf("core smoke test\n");
    std::printf("  dimension      = %d\n", (int)dim);
    std::printf("  particles      = %zu\n", cfg.NumParticle());
    std::printf("  number density = %.4f\n", cfg.NumParticle() / std::pow(L, (int)dim));
    std::printf("  neighbours of particle 0 within R=%.1f : %zu\n", R, neighbours);
    std::printf("OK\n");
    return 0;
}
