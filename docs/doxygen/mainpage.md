# Modules Physics Library

API documentation generated directly from the current C++ and Python source files
in this repository. Existing documentation comments are included, and undocumented
declarations are exposed for navigation. Descriptions have not been invented for
undocumented behavior.

## Explore the code

- **Classes:** class hierarchy, members, inheritance, and collaboration diagrams.
- **Files:** declarations, include dependencies, and linked source listings.
- **Search:** look up symbols such as `GetKs`, `PeriodicCellList`,
  `ShiftedCCPotential`, and `MultiStealthyPotential`.
- **Function pages:** caller and callee diagrams where Doxygen resolves references.

## Repository organization

| Area | Purpose |
| --- | --- |
| core | Geometry, periodic configurations, packing, configuration I/O, and utilities |
| pairstat | Pair statistics, structure factors, and spectral density |
| potential | Potentials, dynamics, Monte Carlo, and optimization |
| models | Configuration and disordered-model generators |
| nonspherical | Shapes, form factors, and dispersions |
| phonons, lp, misc, visualization | Specialized algorithms and visualization support |
| apps | Application drivers and application-specific algorithms |
| examples | Example consumers |

## Scope and interpretation

The build includes tracked and nonignored untracked source files in these areas,
including working-tree edits. Test directories and standalone test/smoke source
files are excluded, along with Graft caches, pilot tooling, external dependencies,
and ignored scratch code. Examples and application drivers remain included.

Doxygen parses the original source independently of Graft. The generated diagrams
are navigation aids, not a complete compiler-resolved call graph. Conditional
compilation uses Doxygen's default preprocessing configuration; optional build
variants and virtual dispatch may have incomplete coverage. Large diagrams are
limited to 60 nodes and depth 3.

The build's `provenance.json` records source hashes, Git state, generator versions,
and timestamp. `warnings.log` retains parser and documentation warnings. Rebuild
after source changes with `python docs/doxygen/build.py`.
