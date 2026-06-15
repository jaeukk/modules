/**
 * @file g2_sk_demo.cpp
 * @brief End-to-end proof consumer: build a point pattern with `core`, then
 *        compute g(r) and S(k) with `pairstat`.
 *
 * This file lives OUTSIDE any single module and depends only on the public
 * module headers (`<...>` resolved via -I$(core) -I$(pairstat) in the makefile).
 * It demonstrates the recommended way to consume the modular library.
 *
 * Build & run:  make && ./g2_sk_demo.out
 */

#include <cstdio>
#include <vector>
#include <functional>
#include <etc.h>
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
#include <PairCorrelation.h>
#include <StructureFactor.h>

// Build one Poisson configuration: N random points in a cubic box of side L.
static Configuration MakePoisson(DimensionType dim, double L, size_t N,
                                 int seed)
{
    GeometryVector basis[::MaxDimension];
    for (DimensionType i = 0; i < dim; ++i) {
        basis[i].SetDimension(dim);
        for (DimensionType j = 0; j < dim; ++j)
            basis[i].x[j] = (i == j) ? L : 0.0;
    }
    Configuration c(dim, basis, /*CellSize*/ 1.0);

    RandomGenerator rng(seed);
    for (size_t n = 0; n < N; ++n) {
        GeometryVector rel;
        rel.SetDimension(dim);
        for (DimensionType i = 0; i < dim; ++i)
            rel.x[i] = rng.RandomDouble();
        c.Insert(AtomName("A"), rel);
    }
    return c;
}

int main()
{
    const DimensionType dim = 3;
    const double L = 10.0;
    const size_t N = 1000;          // number density rho = 1.0
    const size_t NumConfigs = 5;

    // A generator that hands pairstat a fresh (but reproducible) config each call.
    std::function<const Configuration(size_t)> GetConfig =
        [=](size_t i) -> const Configuration {
            return MakePoisson(dim, L, N, 1000 + (int)i);
        };

    // ---- g(r) ----  Result[i] = (r, g2, err_r, err_g2)
    std::vector<GeometryVector> g2;
    IsotropicTwoPairCorrelation(GetConfig, NumConfigs, /*MaxDistance*/ 4.0, g2);

    // ---- S(k) ----  Result[i] = (k, S(k), err_k, err_S)
    std::vector<GeometryVector> sk;
    IsotropicStructureFactor(GetConfig, NumConfigs,
                             /*CircularKMax*/ 10.0, /*LinearKMax*/ 20.0, sk);

    std::printf("g2/Sk demo : %zu configs of N=%zu in L=%.1f^%d (rho=%.2f)\n",
                NumConfigs, N, L, (int)dim, N / std::pow(L, (int)dim));
    std::printf("  g2 bins  = %zu   (first: r=%.3f g2=%.3f)\n",
                g2.size(), g2.empty() ? 0 : g2[0].x[0], g2.empty() ? 0 : g2[0].x[1]);
    std::printf("  S(k) pts = %zu   (first: k=%.3f S=%.3f)\n",
                sk.size(), sk.empty() ? 0 : sk[0].x[0], sk.empty() ? 0 : sk[0].x[1]);
    // For an uncorrelated (Poisson) pattern, g2 ~ 1 and S(k) ~ 1 away from k=0.
    std::printf("OK\n");
    return 0;
}
