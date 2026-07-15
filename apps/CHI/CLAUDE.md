# apps/CHI — known issues and review notes (auto-generated audit)

Date: 2026-07-15. Depends on: core + pairstat. `SpectralDensity.{h,cpp}` and the
empty `VolumeFractionVariance.{h,cpp}` were removed from apps/CHI; CHI now
consumes `pairstat/SpectralDensity` (+ `pairstat/StructureFactor` for `Mtilde`),
linking `pairstat/obj/SpectralDensity.o` and `pairstat/obj/StructureFactor.o`.
As of 2026-07-15 `LinealPathFunction.{h,cpp}` was also moved into pairstat (it is
a two-phase-media microstructure descriptor, a sibling of SpectralDensity); CHI
now includes it via `<LinealPathFunction.h>` and links
`pairstat/obj/LinealPathFunction.o`. The only remaining CHI own source is
`main.cpp` (plus the header `ComputationPacking.h`).

Normal `make` and strict bar (`-Wall -Wextra -Wpedantic -std=c++17 -O2
-fopenmp`) are both at 0 warnings for the apps/CHI TU (`main.cpp`); the only
strict-bar warnings come from core/ headers (`SpherePacking.h`, `ClusterMerger.h`
-Wunused-parameter) — out of scope.

## Unresolved review findings (not fixed — surgical-change rule)

- main.cpp:139,147,359 — `vpComputations` holds raw `new` pointers, deleted only
  on the `Exit` command; the EOF path (`return 2`) leaks them. Suggested:
  `std::vector<std::unique_ptr<ComputationPacking>>`.
- main.cpp:45 — user-supplied filename is the printf format; the argument
  `start_idx + i` is `size_t`, so a `%d` in the input format is UB. Suggested:
  pass `static_cast<int>(this->start_idx + i)` (or document `%zu`).
- main.cpp:67 — `ReadConfigPack::rng` member is never used. Suggested: remove.
- main.cpp:222-223 — `c.p.NumConfig()-StartIdx` underflows (unsigned) when
  `StartIdx > NumConfig()`. Suggested: check `StartIdx >= c.p.NumConfig()` first.
- (Relocated to `pairstat/LinealPathFunction.cpp` on 2026-07-15 — findings kept
  here for CHI's context; they now live in pairstat.)
- pairstat/LinealPathFunction.cpp:24-29,225,247 — boundary inconsistency: guard
  `*iter < size()` but index `*iter-1`, so the length equal to `size()`
  (valid index `size()-1`) is dropped. Suggested: `if (*iter <= v.size())`.
- pairstat/LinealPathFunction.cpp:141-146 — empty `if/else` blocks (dead code).
- pairstat/LinealPathFunction.cpp:152-178 — `if(false)` sanity-check block (dead code).
- Reuse: `ChordLengthDistribution`/`LinealPathFunction` (now in pairstat) share
  setup+summary blocks; `LinealPath`/`ChordLength` classes in ComputationPacking.h
  duplicate their ctor/Write. Suggested: factor a common driver when next touched.

## Out-of-directory warnings / pre-existing errors

Strict compiles surface -Wunused-parameter warnings from core/ headers
(`SpherePacking.h:127`, `ClusterMerger.h:65`); out of scope. The normal build
links cleanly against prebuilt core and pairstat objects.
