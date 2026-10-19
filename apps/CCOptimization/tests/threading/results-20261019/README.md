# October 19 local comparison evidence

These files preserve measurements from implementation commit 856d6d1, before the
production history default changed. JSON contains the original metadata, raw
measurements, endpoint hashes, and summaries. Absolute paths in metadata describe
the original environment, not a required installation path. No executable or
packing dataset is committed. See [the final report](../../../THREADING_RESULTS.md).

- `ab-small`: three species, 300/600/1,000 particles per species; 72 runs.
- `ab-4000`: the same comparison at 4,000 particles per species; 24 runs.
- `history-contact100`: Option A, 4,000 particles per species, phi 0.70,
  contact strength 100, histories 40/100, 18 bounded unsuccessful runs.
- `ground-contact1`: Option A, phi 0.10, **contact strength 1**, histories
  40/100, one seed with two repetitions each; four accepted ground states.
  Its Python validator independently checks saved coordinates, not optimization.
  The separate checks at strength 100 do not make this a strength-100 timing test.
- `endpoint-repeatability`: bounded phi-0.70 A/B endpoint replays. Neither
  schedule produced ground states in these replays.

## Reproduction

Build the committed comparison tool and run its documented runner with
`--history 100` for the original A/B settings. The historical 4,000-particle run
used the same binary with size arguments changed to 4,000. The direct executable
accepts `eval|opt A|B N_PER_SPECIES THREADS SEED EVALUATIONS REPETITION [HISTORY]`.
Bind the process to four physical cores and use the environment in metadata.
Always use a new output directory. `EVALUATIONS` limits each of three optimization
phases, not the whole run. The runner now defaults to history 40.

For either archived history harness, from `apps/CCOptimization`:

```sh
make USE_FINUFFT=1 -f makefile \
  -f tests/threading/results-20261019/history-contact100/history.mk \
  HISTORY_OUT=/tmp/cco-history100-reproduction \
  /tmp/cco-history100-reproduction/history.out
mkdir /tmp/cco-history100-run
cd /tmp/cco-history100-run
# Set affinity and environment as recorded in metadata before invoking:
/tmp/cco-history100-reproduction/history.out opt A 4000 4 17 1000 0 40
```

Use the other `history.mk` and a different `HISTORY_OUT` for `ground-contact1`;
its original cap was 10,000 per phase. Fixtures next to each harness retain the
measured contact strength and phi. The makefiles have been made portable; the
measurement metadata retains the original source and binary hashes. Rebuilt
executables and different hardware may produce different times and endpoints.
NLopt optimization remains entirely in C++.

Full MSPB snapshots and per-run logs are retained only in the original local
`/tmp/cco-*20261019` directories. This archive preserves their hashes and measured
validation results. Without the snapshots, coordinate validation requires a new
run; it cannot be reconstructed from hashes alone.
