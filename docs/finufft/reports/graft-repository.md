# Repository, Graft, and baseline discovery

Completed September 22, 2026. This bounded discovery succeeds: full clone, source-state provenance, genuine installed Graft operations, and a reproducible direct baseline are available. No FINUFFT backend or full optimization experiments were implemented or run. Later implementation must use Graft again during editing; this discovery index alone does not satisfy that later requirement.

## Provenance and boundaries

- Source repository: `/home/jaeukk/30_Codes/modules`; target application: `apps/CCOptimization`.
- Full clone: `../source/modules`, created with `git clone --no-hardlinks` from the source. Both have HEAD `ccdd9285b5f5cac37d37766bfa37f98cfd84ee3c`, all 47 reachable commits, and the clone is not shallow. Full structure is retained, not a reduced dependency extraction. Clone origin is the local source; original upstream is `https://github.com/jaeukk/modules.git`.
- Source has 17 tracked modified files and 16 untracked files. `provenance.json`, `source-status.txt`, and `source-dirty.patch` record them; untracked bytes are hashed in provenance. The clone uses committed HEAD, without applying unrelated dirty edits. The only dirty direct dependency header, `core/RandomGenerator.h`, differs in two comment spellings; the source remains untouched.
- No applicable on-disk AGENTS.md was found in the source ancestor chain, source tree, or workspace search. The supplied global instruction requires American English. Read BRIEF.md and control/COORDINATOR.md completely, and control/STATE.json. Read CCOptimization/CLAUDE.md for historical audit context; no historical build was counted as current evidence.
- Only `repository/` and generated clone build/index files were created. No tracked application or dependency source changed. Current clone status is only untracked `graft/` (build artifacts are already ignored). No pushes, publication, paid model calls, or global package changes.

## Build and executable baseline

Run from any directory:

```bash
bash /home/jaeukk/10_Projects/10_Active/multi-stealthy/_multiagent/tasks/2026-09-22-graft-finufft/repository/run-baseline.sh
```

`build.log` records a successful clean CCOptimization build (`make -j2`, exit 0). `baseline.cpp` exercises actual ShiftedCCPotential and MultiStealthyPotential with the same single-species configuration. `baseline.log` records successful compilation and execution, exit 0. Environment: GCC 13.3.0, GSL 2.8, NLopt 2.10.0, C++17/OpenMP; existing numerical libraries are read-only at `/home/jaeukk/install`. See `environment.json`. The baseline script currently names that installation explicitly, so another machine must substitute its dependency prefix. Node dependencies are isolated in the task's Graft checkout.

The baseline uses N=16, dimension 2, box side 4, RNG seed 20260922, four independent wavevectors `(1,0), (0,1), (1,1), (1,-1)` scaled by `2π/4`, unit weights, zero contact strength, and two OpenMP threads. Results:

```text
energy=7.7645917434367346
force_l2=7.9615103628771919
multi_energy=7.7645917434367346
force_max_difference=0
rho[0]=-6.9103317573827798,4.104958119187943
rho[1]=-0.24416008449275836,-5.3249737796220353
rho[2]=2.7392938155755857,2.5969308711680679
rho[3]=1.6974975103165226,-3.7531138381188067
```

This establishes execution and interface consistency, not finite-difference correctness, FINUFFT agreement, convergence, or either scientific gate. No long-running CLI defaults were invoked. Commands have explicit timeouts; build parallelism and numerical thread counts were capped at two. The reusable baseline additionally sets OMP_THREAD_LIMIT=2.

## Genuine Graft use and limitations

Inspected official [Graft documentation](https://github.com/trailhq/Graft) and actual source at commit `597d82f8c2f6b1dfa2daf547073a8424b5535eb2`, version 0.19.0. Full local checkout is `Graft/`; `npm ci --prefix repository/Graft` installed 74 packages and built the CLI, recorded in `graft-install.log`. Node 24.11.1 and npm 11.6.2 were used. The upstream package names itself `@nanonets/graft`; this was confirmed from the requested repository rather than inferred from the name.

The structural build ran successfully against the full clone with scope filters `apps/CCOptimization`, `core`, `pairstat`, and `potential`. It parsed 51/51 files and generated 867 nodes, 575 edges, and 31 cards. `graft-build.log` is the output. The graph lives at `../source/modules/graft/.graph/wiring.json`. Build did not change .gitignore, .ignore, agent configurations, or global settings. DO_NOT_TRACK=1 and CI=1 were set; no --deep, provider key, or paid endpoint was used. `graft check` reports synchronized structure with the meaning tier unbuilt, as expected.

Concrete retrieval artifacts:

| Operation | Evidence and contribution |
|---|---|
| `skeleton potential/CollectiveCoordinatePotential.h` | `graft-skeleton.txt`: quickly locates cached rho, energy, configuration setters, and derived variants. |
| `callers GetRho --json` | `graft-callers-rho.json`: identifies several definitions and explicitly reports ambiguous-name undercount. |
| `map --json` | `graft-map.txt`: scoped directory/symbol overview. |
| `grep 'GetRho\|AllForce\|CCPotential' --in apps/CCOptimization --json` | `graft-grep-rho.txt`: actual source hits expose wrappers and species dispatch. |
| `callers AllForce --json` | `graft-callers-allforce.txt`: caller exploration, subject to ambiguity limits. |
| `skeleton apps/CCOptimization/MultiStealthyPotential.cpp` | `graft-skeleton-multi.txt`: source method locations. |
| `grep 'GetKs\|SetConfiguration\|AllForce' --in potential/StructureOptimization.cpp --json` | `graft-optimizer-calls.json`: locates optimizer evaluation and gradient calls. |
| `check` | `graft-check.txt`: freshness evidence. |

`graft-operations.json` records argv and exit codes for the main retrieval batch. Graft is useful for orientation and symbol-grouped source retrieval here, but is not a complete C++ call graph:

1. `Graft/src/graph/generic.ts:46` maps `.h` to the C grammar; C++ uses `.cpp/.cc/.cxx/.hpp/.hh`. Our header skeleton omits base Force and other methods, labels many remaining methods as functions, and misidentifies constructors. A successful parse count does not prove semantic completeness.
2. Ambiguous names such as GetRho and AllForce cause cross-file edges to be dropped. Empty caller results do not mean no callers. Verify with Graft grep and actual source/makefiles.
3. No clangd was found on PATH, so compiler-backed `--lsp` was not attempted. No parser patches or source renaming were made.
4. Only the four scoped directories are indexed; other repository consumers require a wider graph or repository search before shared-header changes.
5. `graft-audit.json` records npm's one high-severity transitive js-yaml finding; no automatic dependency rewrite was performed. Runtime is used only for this local, bounded structural analysis.

Reusable CLI: `bash repository/graft.sh <command> ...`. It caps execution at 120 seconds and disables telemetry and automatic ignore-file edits. During implementation run scoped Graft grep/skeleton/callers before edits, and rerun build/check plus affected-caller queries afterward, saving new output under a distinct implementation evidence directory. It refreshes against working-tree changes by default. Keep one owner for graph rebuilds; coordinate rather than concurrently mutating the same cache. This worker can be redispatched to support that phase.

## Dependency map and numerical interfaces

| Layer | Files / responsibilities |
|---|---|
| CLI and setup | `apps/CCOptimization/main.cpp`: single/multi driver selection, initialization, modes, weights, optimizer choices and output. `makefile` compiles main.cpp and MultiStealthyPotential.cpp. DebugSimple.cpp is excluded scratch. |
| Species wrapper | `MultiStealthyPotential.{h,cpp}`: maps global particles to species configurations, wraps ShiftedCCPotential, adds soft contacts. `RepulsiveCCPotential.h` wraps the one-component contact case. |
| Numerical CC kernel | `potential/CollectiveCoordinatePotential.h`: ShiftedCCPotential plus weighted, logarithmic, S2, multiconfiguration, and Monte Carlo variants. |
| Objective/gradient consumers | `potential/StructureOptimization.cpp`: calls SetConfiguration, Energy, AllForce; converts Cartesian forces to fractional-coordinate derivatives using basis dot products and a minus sign. `potential/MD_System.cpp` consumes Potential energy/force for dynamics. |
| Geometry and modes | `core/GeometryVector`, `PeriodicCellList`, `SpherePacking`, `ConfigurationIO`, `RandomGenerator`, `etc`; GetKs is in `core/PeriodicCellList.cpp:509`. |
| Analysis | `pairstat/StructureFactor` and `PairCorrelation`; output analysis should retain direct summation as an independent check. |
| External/build | `config.mk` detects GSL/NLopt; core/pairstat/potential makefiles build objects. CCOptimization links specified objects, although dependency `objs` targets compile additional files. |

Base semantics verified in source: `rho(k)=sum_j exp(+i k·r_j)` with absolute Cartesian wavevectors; energy is `sum_k V_k |rho(k)|² / cell_volume + Shift`. Force is negative Cartesian gradient, with factor `2 V_k / cell_volume`. Public RhoReal/RhoImag align exactly with public constraints order. SetConfiguration invalidates rho. Direct AddConstraint does not invalidate an already-computed rho, so integration must handle/cache constraint mutation deliberately. Configuration lifetime and mode ordering are part of the interface.

GetKs excludes zero and chooses one representative of each opposite pair using the first nonzero reciprocal integer's sign. Preserve the actual returned constraints rather than rebuilding an approximate radial count. Multi CLI accepts cutoff chi normalized against total density, then reports realized species chi as `M/[d(N_s-1)]`; equal input cutoffs and species counts need explicit checking. Component CC energies and forces receive `u0 * N_s/N` weights. Soft repulsion uses radius sums, an overlap penalty, and e0; finite optimization energy does not impose hard non-overlap. These facts inform later model validation; no feasibility claim is made.

## Recommended integration ownership

Assign one implementation owner `potential/CollectiveCoordinatePotential.h`, any new FINUFFT kernel files, `apps/CCOptimization/MultiStealthyPotential.{h,cpp}`, `RepulsiveCCPotential.h`, main.cpp, and necessary makefile/config changes. Keep validation scripts, independent references, and benchmark manifests owned by validation. Shared-header derived variants must retain correct direct behavior unless explicitly supported.

The fastest safe boundary is selectable rho and all-particle-force evaluation behind the existing potential interface, keeping direct as default/reference. Preserve units, positive Fourier sign, volume normalization, per-species weights, mode ordering, cache invalidation, and negative-gradient convention. Public Force must remain consistent with AllForce. Critically, MultiStealthyPotential::AllForce currently loops its own Force, which calls each species Force individually: merely accelerating the species AllForce method leaves the multi path on particle-by-particle evaluation. Aggregate/cache species forces or explicitly route through each species AllForce in the integration phase, with thread-safe invalidation and bounded threads. Existing nested OpenMP and first-use caches merit review before parallel accelerated calls.

No reconstruction of particle positions appears in this active target path: reciprocal-to-particle force evaluation is the needed operation, not an inverse position reconstruction. Backend implementation, numerical Gate A, paired optimization Gate B, and five-component studies remain assigned future work. There is no blocker to beginning that implementation once the coordinator's literature/manual gate is met; Graft's header/caller incompleteness is a documented tooling limitation requiring source verification.
