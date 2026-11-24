/** \file test_contact_grid.cpp
 * \brief The contact grid of MultiStealthyPotential against the generic neighbor walk.
 * The grid must give the bits of the walk in energy and forces (1-3D, orthogonal and skew cells, unequal radii,
 * coincident particles, coordinates on cell faces, cells holding several periodic images, cells filled in another
 * order than the particle index), an optimization must end at the same bits with either, and with the grid the
 * three entry points must stay consistent and independent of the thread count. */
#include "fixtures.h"
#include <StructureOptimization.h>
#include <cstring>
#include <ctime>
#include <limits>

size_t Verbosity = 0;
using namespace threading_test;

namespace {
int failures = 0;

/** \brief Report one check.
 *  \param ok     Whether the check passed.
 *  \param label  What was checked. */
void Check(bool ok, const std::string & label)
{
    failures += !ok;
    std::cout << (ok ? "ok   " : "FAIL ") << label << '\n';
}

/** \brief Whether two results are bitwise equal in energy and every force component.
 *  \param a  First result.
 *  \param b  Second result.
 *  \return True if identical. */
bool Identical(const Result & a, const Result & b)
{
    if (a.forces.size() != b.forces.size() || std::memcmp(&a.energy, &b.energy, sizeof(double)) != 0) return false;
    for (size_t i = 0; i < a.forces.size(); ++i)
        if (std::memcmp(a.forces[i].x, b.forces[i].x, sizeof(a.forces[i].x)) != 0) return false;
    return true;
}

/** \brief A packing of the fixture with unequal radii, one coincident pair and particles on cell faces.
 *  \param dim   Dimension.
 *  \param n     Particle number.
 *  \param skew  Skew cell.
 *  \return The packing. */
SpherePacking Awkward(int dim, size_t n, bool skew)
{
    SpherePacking base = Packing(dim, {n}, 40 + dim, skew);
    GeometryVector basis[3];
    for (int d = 0; d < dim; ++d) basis[d] = base.GetBasisVector(d);
    SpherePacking packing(dim, basis, 1.0);
    for (size_t i = 0; i < n; ++i) {
        GeometryVector position = base.GetRelativeCoordinates(i);
        if (i == 1) position = base.GetRelativeCoordinates(0);	// coincident with particle 0
        if (i == 2) position.x[0] = 0.0;	// on the cell face
        if (i == 3) position.x[0] = -1e-18;	// wraps to 1.0 before it wraps to 0
        if (i == 4) for (int d = 0; d < dim; ++d) position.x[d] = 1.0 - 1e-16;
        packing.Insert(base.GetCharacteristics(i) * (0.6 + 0.4 * double(i % 5) / 4.0), position);
    }
    return packing;
}
}

/** \brief Run every check.
 *  \return Number of failed checks. */
int main()
{
    try {
        for (int dim = 1; dim <= 3; ++dim)
            for (bool skew : {false, true})
                for (size_t n : {size_t(60), size_t(700), size_t(3000)}) {
                    const std::vector<size_t> sizes = {n / 3, n - n / 3};
                    const auto packing = Awkward(dim, n, skew);
                    Configuration config(packing, "a");
                    const std::string label = "dim " + std::to_string(dim) + (skew ? " skew" : " orthogonal") + " N " + std::to_string(n);
                    auto grid = MakePotential(packing, sizes, 4, false, Backend::Direct);
                    auto walk = MakePotential(packing, sizes, 4, false, Backend::Direct);
                    walk->SetContactGridEnabled(false);
                    const Result a = Evaluate(*grid, config), b = Evaluate(*walk, config);
                    Check(grid->ContactGridInUse() && !walk->ContactGridInUse(), label + ": grid in use, and off when disabled");
                    Check(Identical(a, b), label + ": grid gives the bits of the neighbor walk");
                    // The three entry points share the kernel: fused = Energy then AllForce, bitwise.
                    grid->SetConfiguration(config);
                    Result separate;
                    separate.energy = grid->Energy();
                    grid->AllForce(separate.forces);
                    Check(Identical(a, separate), label + ": EnergyAndAllForce equals Energy then AllForce bitwise");
                    GeometryVector single;
                    grid->Force(single, n / 2);
                    Check(std::memcmp(single.x, a.forces[n / 2].x, sizeof(single.x)) == 0, label + ": Force(i) equals AllForce bitwise");
                    // Contact alone (one species without constraints): the thread count does not change a bit. With
                    // constraints the force gate may choose another path for another team, so there only a repetition
                    // at the same thread count must agree bitwise (FINUFFT backend: contact chunks between transforms).
                    {
                        Result first;
                        bool same = true;
                        for (int threads : {1, 4, 8}) {
                            omp_set_num_threads(threads);
                            QuietOutput quiet;
                            MultiStealthyPotential contact(dim, 100.0);
                            contact.ParallelNumThread = threads;
                            contact.SetPacking(packing);
                            contact.SetNumSpecies(1);
                            for (size_t i = 0; i < n; ++i) contact.SetSpeciesOfParticle(i, 0);
                            contact.SetPacking(packing);
                            const Result r = Evaluate(contact, config);
                            if (threads == 1) first = r;
                            same = same && Identical(r, first) && contact.ContactGridInUse() && r.energy > 0.0;
                        }
                        omp_set_num_threads(4);
                        Check(same, label + ": contact alone is bitwise the same for 1, 4 and 8 threads");
                        auto one = MakePotential(packing, sizes, 4, false), two = MakePotential(packing, sizes, 4, false);
                        Check(Identical(Evaluate(*one, config), Evaluate(*two, config)) && one->ContactGridInUse(),
                            label + ": FINUFFT backend with the grid repeats bitwise");
                    }
                    // Switching the kernel on an existing potential rebuilds the contact data.
                    walk->SetContactGridEnabled(true);
                    Check(Identical(Evaluate(*walk, config), a), label + ": enabling the grid later gives the grid result");
                }
        // Small cells: one or two cells per side, so a neighbor cell appears as several periodic images.
        for (int dim = 1; dim <= 3; ++dim)
            for (size_t n : {size_t(2), size_t(8), size_t(20)}) {
                const std::vector<size_t> sizes = {n / 2, n - n / 2};
                auto packing = Packing(dim, sizes, 9);
                for (size_t i = 0; i < packing.NumParticle(); ++i) packing.GetCharacteristics(i) *= 1.4;	// contacts across images
                Configuration config(packing, "a");
                auto small = MakePotential(packing, sizes, 2, false, Backend::Direct);
                auto walk = MakePotential(packing, sizes, 2, false, Backend::Direct);
                walk->SetContactGridEnabled(false);
                const Result a = Evaluate(*small, config), b = Evaluate(*walk, config);
                Check(small->ContactGridInUse() && Identical(a, b) && b.energy > 0.0,
                    "dim " + std::to_string(dim) + " N " + std::to_string(n) + ": several periodic images, bits of the neighbor walk");
            }
        // Cells filled in another order than the particle index (particles moved one by one), and a second state
        // through the same potentials (stored configuration reused).
        {
            const std::vector<size_t> sizes = {300, 300};
            const auto packing = Packing(2, sizes, 11);
            Configuration config(packing, "a"), moved(config);
            RandomGenerator random(5);
            for (size_t step = 0; step < 400; ++step) {
                GeometryVector target(2);
                for (int d = 0; d < 2; ++d) target.x[d] = random.RandomDouble();
                moved.MoveParticle((step * 7) % moved.NumParticle(), target);
            }
            // The source configurations have answered neighbor queries at other radii, which reorders their cached
            // lattice lists; a copy starts from a clean cache, and so must the reused stored configuration.
            for (Configuration * state : {&config, &moved})
                for (double radius : {6.0, 2.3, 4.1})
                    state->IterateThroughNeighbors(state->GetRelativeCoordinates(0), radius,
                        [](const GeometryVector &, const GeometryVector &, const signed long *, size_t) {});
            auto grid = MakePotential(packing, sizes, 4, false, Backend::Direct);
            auto walk = MakePotential(packing, sizes, 4, false, Backend::Direct);
            auto reusedWalk = MakePotential(packing, sizes, 4, false, Backend::Direct);
            walk->SetContactGridEnabled(false);
            walk->ReuseConfigurationStorage = false;
            reusedWalk->SetContactGridEnabled(false);
            bool same = true;
            for (const Configuration * state : {&moved, &config, &moved}) {
                const Result reference = Evaluate(*walk, *state);
                same = same && Identical(Evaluate(*grid, *state), reference) && Identical(Evaluate(*reusedWalk, *state), reference);
            }
            Check(same && grid->ContactGridInUse(), "cells in arrival order, queried sources, three states in turn: bits of the neighbor walk on fresh copies");
        }
        // The property that matters for a run: an optimization ends at the same bits with the grid and with the walk.
        {
            ProgramStart = std::time(nullptr);
            TimeLimit = std::numeric_limits<time_t>::max();
            const std::vector<size_t> sizes = {400, 400, 400};
            const auto packing = Packing(2, sizes, 17);
            Configuration end[2] = {Configuration(packing, "a"), Configuration(packing, "a")};
            for (int g = 0; g < 2; ++g) {
                auto potential = MakePotential(packing, sizes, 4, false);
                potential->SetContactGridEnabled(g == 1);
                QuietOutput quiet;
                RelaxStructure_NLOPT(end[g], *potential, 0.0, 0, 0.0, 120);
            }
            bool same = end[0].NumParticle() == end[1].NumParticle(), movedAway = false;
            const Configuration start(packing, "a");
            for (size_t i = 0; same && i < end[0].NumParticle(); ++i) {
                same = std::memcmp(end[0].GetRelativeCoordinates(i).x, end[1].GetRelativeCoordinates(i).x, sizeof(end[0].GetRelativeCoordinates(i).x)) == 0;
                movedAway = movedAway || std::memcmp(end[0].GetRelativeCoordinates(i).x, start.GetRelativeCoordinates(i).x, sizeof(start.GetRelativeCoordinates(i).x)) != 0;
            }
            Check(same && movedAway, "optimization of 360 evaluations: same end point, bit for bit, with the grid and with the walk");
        }
        // No repulsion: no grid, zero contact.
        {
            const std::vector<size_t> sizes = {50, 50};
            const auto packing = Packing(2, sizes, 3);
            Configuration config(packing, "a");
            auto off = MakePotential(packing, sizes, 2, false, Backend::Direct, 0.0);
            Evaluate(*off, config);
            Check(!off->ContactGridInUse(), "no repulsion: grid not built");
        }
    }
    catch (const std::exception & error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << (failures == 0 ? "ALL PASSED\n" : "FAILURES\n");
    return failures;
}
