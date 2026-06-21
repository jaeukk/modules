# Module source provenance

This file records, per file copied into `modules/`, **which legacy copy was
chosen as canonical** and why. Policy: *prefer the subdirectory "module" copies*
(they are the leaner, decoupled versions), with explicit exceptions where the
subdir copy is a stub or is missing features that consumers need.

Legend: ✅ chosen source · ⚠️ exception to the default policy.

## core/
Source: **`cores/`** ✅ (newer than top-level; 2024 vs 2020-2022).

| File | From | Notes |
|------|------|-------|
| GeometryVector.{h,cpp} | cores/ | identical structure to top-level, newer date |
| PeriodicCellList.{h,cpp} | cores/ | cores/ slightly larger (1681 vs 1665 h-lines), newer |
| RandomGenerator.{h,cpp} | cores/ | newer |
| etc.{h,cpp} | cores/ | newer |

## pairstat/
Source: **`PairStat/`** ✅ (lean, decoupled from NonsphericalParticles).

| File | From | Notes |
|------|------|-------|
| PairCorrelation.{h,cpp} | PairStat/ | lean (808 vs 1068 cpp-lines); top-level extras dropped |
| StructureFactor.{h,cpp} | PairStat/ | lean (260 vs 1287 cpp-lines); ⚠️ drops nonspherical/spectral-density used by NonsphericalParticles/, FirstPassageTime/ — those consumers stay on legacy headers or get features ported behind USE_SHAPES |
| NumberVariance.{h,cpp} | PairStat/ | lean (296 vs 1524 cpp-lines) |
| Computation.{h,cpp} | PairStat/ | only exists in PairStat/ |
| ContactStatistics.{h,cpp} | PairStat/ | only exists in PairStat/ |

## potential/
Source: **`Potential/`** ✅ (newer, includes Move.h & MonteCarlo.h).

| File | From | Notes |
|------|------|-------|
| Potential.{h,cpp} | Potential/ | nearly identical, newer |
| MD_System.{h,cpp} | Potential/ | nearly identical, newer |
| StructureOptimization.{h,cpp} | Potential/ | nearly identical, newer |
| MonteCarlo.h, Move.h, CollectiveCoordinatePotential.h, AnalyzeData.h | Potential/ | present only in Potential/ |
| LJPotential.{h,cpp}, MC_System.{h,cpp} | top-level | no subdir copy |

## models/
Default policy would pick `Models/`, but:

| File | From | Notes |
|------|------|-------|
| GenerateConfigs.{h,cpp} | **top-level** | ⚠️ EXCEPTION: Models/ copy is a 129-line stub vs 479-line real implementation |
| RandomSequentialAddition.{h,cpp} | top-level | no clean subdir module |
| HardSphereFluids.{h,cpp} | top-level | no clean subdir module |
| RandomOrganization.{h,cpp} | top-level | no clean subdir module |
