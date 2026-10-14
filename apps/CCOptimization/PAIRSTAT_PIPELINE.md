# In-memory multicomponent PairStat pipeline

The adapter `MultiComponentPairStat.{h,cpp}` connects the accepted optimization
configuration to the existing `pairstat/` library. **No PairStat source changes**
or new statistics estimators are required. It does not launch PairStatCLI.

```cpp
#include "MultiComponentPairStat.h"

potential.SetConfiguration(final_configuration);
MultiComponentPairStatOptions options;
options.kMax = 2.0;
options.kBinWidth = 0.05;
options.rMax = 2.0;
options.rStep = 0.1;
options.windowSamples = 1000;
options.seed = 17;
auto results = ComputeMultiComponentPairStatistics(potential, options);
// Consume results in memory; writing a final report is optional.
WriteMultiComponentPairStatistics(std::cout, results, options);
```

The potential supplies a `SpherePacking` snapshot containing the actual periodic
basis, relative positions, and per-particle radii. The adapter validates that species
index lists form a complete, disjoint partition and creates species subsets in RAM.
There is one result for the combined packing (`species=-1`) followed by each declared
species, including empty species (zero count, empty tables). Labels and unequal
radii are preserved; this is not restricted to five equally populated components.
The current adapter supports 2D and 3D.

## Reused calculations

| Result | Existing PairStat function | Interpretation |
|---|---|---|
| Isotropic structure factor | `IsotropicStructureFactor` | Reciprocal-point average, all eligible modes up to Kmax, normalization by each subset's particle count |
| Local number variance | `MC_NV_single` | Spherical observation windows, prescribed relative centers |
| Local volume-fraction variance | `MC_vfv` | Sphere-window intersections with actual particle radii |

`GetSamplingPoints` supplies seeded random centers reused across the combined packing,
all species, and both variance estimators. Results retain PairStat's native four-column
tables. The combined structure factor is calculated from all particles; it is not the
average of individual-component structure factors. Likewise, total variances are
computed directly, preserving cross-component effects.

## CLI

Run `soft_core_stealthy.out multi`, choose `ground`, and enter, for example:

```text
combined
export binary
pairstat 2.0 0.05 2.0 0.1 1000 17
run
```

The `pairstat` arguments are `Kmax dK Rmax dR samples seed`; analysis is disabled
unless this command is supplied. The existing `# threads` input controls the OpenMP
runtime used by PairStat. Each newly accepted configuration generates final statistics
in `<prefix>_PairStat_<success-index>.txt`, alongside the optional configuration export.
Statistics also work with `export off`. The index is zero-based and follows the
success pack when resuming. Existing successes are not reanalyzed automatically.

The adapter never writes or rereads an intermediate configuration. Legacy restart
packs are still saved by the optimizer. `combined` additionally disables the old
post-run species-pack splitting path, which otherwise rereads success packs; it does
not disable per-species analysis in this adapter. MD and the preexisting incomplete
filename-input branch are not extended here.

## Interpretation and limitations

- K and R use Cartesian inverse-length and length units, respectively. These are
  per-configuration estimates, not ensemble averages or new acceptance criteria.
- PairStat samples radii `dR, 2*dR, ..., ceil(Rmax/dR)*dR`; the last radius can
  slightly exceed Rmax. Choose window ranges with finite-cell effects in mind.
- `MC_vfv` adds sphere-window intersection volumes. Its phase-volume interpretation
  assumes nonoverlapping spheres. Each result includes `overlaps` from the existing
  overlap check; any overlap, including roundoff-sized contact, sets this flag.
- The volume-variance ensemble error is unavailable for one configuration; the
  library's zero placeholder must not be interpreted as zero uncertainty.
- Seeded centers are reproducible, but threaded floating-point reductions need not
  be bitwise repeatable. No optimizer, acceptance threshold, or FINUFFT setting is changed.
- Archived MSPB snapshots can now be analyzed offline through `mspb_pack.out
  pairstat`; see [MSPB_PACK.md](MSPB_PACK.md). This reuses the same adapter without
  running an optimization sweep or automatically reanalyzing historical data.

## Validation

`bash tests/run_pairstat_test.sh` compares every returned table against separate
direct calls to the existing library on independently assembled subsets: total plus
five occupied species, unequal radii, an empty species, 2D/3D rectangular and skewed
cells. It also checks incomplete/duplicate assignments, invalid parameters, and
report-stream errors. `python3 tests/test_cli.py` exercises CLI export and analysis.
`bash tests/run_export_test.sh` covers text/MSPB serialization independently.
