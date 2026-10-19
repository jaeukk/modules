# Options A and B: species scheduling

These are temporary comparison names. Both use the existing C++ GetMultiCCO /
NLopt path and the same potential formulas and file formats.

- **A:** the original concurrent species scheduler. Internal FINUFFT threads stay
  at one. Species count plus one is a starting allocation, not a prescribed
  setting. Three species use four outer threads in this comparison.
- **B:** species are evaluated sequentially, each completing density and force
  work before the next. FINUFFT requests all `ParallelNumThread` threads, capped
  by `OMP_THREAD_LIMIT`. An additional asynchronous worker computes contact terms
  concurrently. As requested, contact time-shares the allocated CPUs; no FINUFFT
  thread is reserved for it. Four FINUFFT threads plus contact can make five
  runnable threads. **Use process affinity or scheduler CPU binding to limit CPU
  capacity to four.** The benchmark runner does this for both options.

Select A with `finufft`, or B with `finufft-sequential`, in a multicomponent
`ground` or `MD` block. In C++, select FINUFFT as usual, then call
`SetSequentialFinufftSpecies(true)` for B or `false` for A. The default remains A.
Returning to A restores one internal FINUFFT thread. The single-component GetCCO
CLI is unchanged.

B keeps the existing deterministic Type 1 particle partition, rather than adding
an independent block policy in core. Its blocks execute sequentially with internal
threading, and all its species are sequential. The bounded 2D comparison cases
are expected to be small enough for a single density block; the partition remains
input-dependent. B uses one batched Type 2 block for eligible base species: all
spatial components enter one execution with `ntransf = dimension`. It does not use
the force gate fitted to concurrent one-thread tasks. Unsupported inputs,
derived force implementations, and failed transforms retain direct force sums.
Compatible plans persist, and point state is refreshed with `setpts`.
Base species execute outside an enclosing OpenMP team. Derived density callbacks
retain their own virtual implementation and are responsible for any nested
parallelism they introduce.

With no contact work, no asynchronous worker starts. A one-thread effective
budget runs contact and species serially. Species exceptions restore the per-species
outer thread setting; destruction of the asynchronous future joins the contact
worker before its captured buffers are destroyed. Contact exceptions propagate
through `get()`. Species and contact results have separate buffers. Final force
assembly runs only after contact finishes; energy slots are summed in index order.

Internal multithreading can change floating-point summation order. B preserves
the scientific formulas but does not promise bitwise repetition, an identical
optimizer trajectory, or identical acceptance outcomes. The bounded comparison
uses direct-reference checks; it is not a new ground-state stopping protocol.

## Validation

```sh
make -C apps/CCOptimization USE_FINUFFT=1 test-sequential-species
make -C apps/CCOptimization USE_FINUFFT=1 soft_core_stealthy.out
python3 apps/CCOptimization/tests/threading/test_cli.py
```

The focused suite checks 1D–3D, skew cells, unequal and empty species, no contact,
cold/warm caches, particle motion, separate and fused energy/force orders, toggling
A/B, Type 1/2 failures, plan reuse/destruction, batched components, the internal
thread setting, sequential transform execution, finite-difference gradients, and
an exception during concurrent contact work. It also runs with
`OMP_THREAD_LIMIT=1`. Approximate comparisons use scaled tolerance 2e-10.
The existing direct-only CLI/export/archive/PairStat tests remain applicable.

## Performance comparison, after the implementation commit

```sh
make -C apps/CCOptimization USE_FINUFFT=1 bench_threading.out
python3 apps/CCOptimization/tests/threading/compare.py \
  --output /tmp/cco-options-comparison-NEW
```

The runner selects four distinct cores when the machine exposes enough, binds
itself and all children to those CPUs, and runs policies sequentially. It records
the commit, dirty state, executable hash, hardware, affinity, environment, raw
measurements, actual backends, internal thread settings, and commands.

Cases use three species, 300/600/1,000 particles per species, 2D, unit number
density, nominal per-species chi 0.30, packing fraction 0.70, contact strength 100,
and FINUFFT tolerance 1e-12. All species have contiguous particle IDs. Seeds 17–19
are shared between policies, with policy and size order rotated.

The evaluator benchmark performs 400 warm `SetConfiguration + EnergyAndAllForce`
calls per batch, three repetitions per seed. The optimization benchmark calls
`RelaxStructure_NLOPT`, as GetMultiCCO does, with default history 100 and a cap of
1,000 evaluations per coordinate-scaling phase. It includes that optimizer's
configuration handling and progress I/O, but excludes fixture setup and the
post-run independent direct check. Each process gets a separate output directory.
Both initial and final energy/force results are checked against direct sums.

CPU utilization is process CPU time / wall time / four. CPU accounting includes
contact and FINUFFT workers, as well as synchronization. Research datasets are not read or changed. Warm evaluation timing isolates
scheduling; complete optimizer trajectories can differ between policies, so a
bounded optimization is not a measurement of time to an accepted ground state.
