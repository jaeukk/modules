# Generate the API documentation

Open `generated/html/index.html` in a browser. This is actual Doxygen output from
the current repository sources, including uncommitted changes. The generated
HTML, XML, diagrams, logs, and provenance are local build artifacts excluded from
Git by this directory's `.gitignore`.

From the repository root:

```bash
python docs/doxygen/build.py
```

The script requires Doxygen and Graphviz (`dot`). It uses `DOXYGEN_BIN`, Doxygen
on PATH, or the isolated `/tmp/modules-doxygen-tool/usr/bin/doxygen`, in that
order. The temporary fallback also sets its own runtime-library search path.
Temporary tools may disappear after cleanup or reboot; a normal Doxygen install
or a custom `DOXYGEN_BIN` keeps this configuration reusable.

The initial environment uses Ubuntu's Doxygen 1.9.8 and Graphviz 2.43.0. Doxygen
and missing libraries were downloaded as Ubuntu packages and extracted under
`/tmp/modules-doxygen-tool` without installing system packages. No package
maintainer scripts were executed. The downloaded packages are in
`/tmp/modules-doxygen-packages`.

The script selects tracked and nonignored untracked C/C++ and Python sources
from the repository's library, application, and example directories. It excludes
directories named `test` or `tests`, source basenames `test`, `tests`, or `smoke`,
prefixes `test_` or `smoke_`, and suffixes `_test` or `_tests` (case-insensitive).
This excludes the FINUFFT test helpers and `core/smoke.cpp`. Examples and
application drivers remain included. Graft-generated cards and pilot scripts
are also excluded. Documentation comes from original code comments;
Doxygen resolves its own symbols and diagrams independently of Graft.

Search, source browsing, class hierarchies, include dependencies, and caller/callee
diagrams are enabled. `.h` files are explicitly mapped to C++. Diagram rendering
uses two workers, with graphs capped at 60 nodes and depth 3. Conditional build
variants are not enumerated, and static call diagrams may miss virtual dispatch.
The settings follow the [Doxygen configuration reference](https://www.doxygen.nl/manual/config.html).

Inspect `generated/warnings.log` for existing malformed or ambiguous source
comments and unresolved references. A successful build does not mean every
comment or call edge is correct. `generated/provenance.json` records source
hashes, Git state, tool versions, duration, exit status, and whether any source
changed during generation. No source comments are rewritten by this workflow.

The HTML is a multi-file site: keep its entire directory together when copying
it. Rebuilding clears the generated HTML and XML directories before generation,
so removed inputs cannot leave stale pages or search entries. The source files
are not deleted. `provenance.json` lists `excluded_test_sources` for each build.

## Historical initial build (before test exclusion)

The October 1, 2026 build completed in 75.21 seconds with exit status 0 and no
source changes. It generated 1,011 HTML pages, 5,102 SVG diagrams, and XML indexes
covering 138 classes, 57 structs, and 17 namespaces. Both `Force` and `AllForce`
are present for `ShiftedCCPotential` and `MultiStealthyPotential`. Local target
files linked from the landing page, file/class indexes, and those two class pages
were checked successfully.

The generator reported 909 warnings from existing documentation and graph-size
limits, plus 14 graph hash conflicts for repeated diagram output names. These
conflicts affect diagrams for symbols including `HSF`, `PrefixNumberSuffix`,
`ReadConfigPack`, and `RSA_Gen`; those diagrams should not be assumed reliable.
The HTML/API documentation remains available. Source-comment cleanup and
disambiguating duplicated symbols are separate work; original code was preserved.
Current build results are recorded in `generated/validation.json` and
`generated/warnings.log`; their counts supersede the historical figures above.

The initial effective configuration also contained the obsolete, harmless
`HTML_TIMESTAMP` setting; it has been removed from the reusable Doxyfile.
