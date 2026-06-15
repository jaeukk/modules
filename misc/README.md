# misc

Small header-only support pieces that several apps share but that don't belong
to any single computational module. No `main`, no compiled objects of its own —
the headers are pulled in (via `-I$(misc)`) by consumers; the corresponding
`.cpp` implementations (where they exist) still live in the legacy tree.

| File | Purpose |
|------|---------|
| `Discrepancy.h` | `MergeSort` / `Merge_Sort` ordering helpers, `Crop`, van der Corput sequences, star-discrepancy templates. |
| `ParameterSet.h` | Abstract `ParameterSet` base + `CombinedParameterization` for pair-potential optimisation. |
| `Interfaces.h` | The legacy `cli()` command dispatcher and the `*CLI()` entry-point declarations (kept as reference; individual apps extract the command they need). |

## Dependencies
- **Lower modules:** `core` (GeometryVector, PeriodicCellList, etc), `potential`
  (`ParameterSet.h` → Potential.h), `visualization` (`Interfaces.h` → Plots.h,
  DisplaySpheres.h).
- **External:** GSL (via core).

## Notes
- Header-only: there is no makefile here. Consumers add `-I$(misc)` (and, because
  `Interfaces.h` includes the visualization headers, usually `-I$(visualization)`
  too).
