# MSPB snapshot collections

`MSPB_Pack` manages an ordered ensemble in `<prefix>.MSPBPack`. It imports existing
`.mspb` files, restores snapshots by zero-based index, and exports individual files
without changing their MSPB v1 bytes. Every snapshot retains its cell, particle IDs,
species labels, declared species count (including empty species), relative
coordinates, and per-particle radii. Different snapshots may have different cells,
particle counts, dimensions, and species counts.

All code is local to `apps/CCOptimization`. `SpherePackingPack` supplies a useful API
pattern but is not generalized here: its native binary records do not store particle
IDs or species labels, and changing it would affect core consumers. An MSPB snapshot
can provide a `SpherePacking` view for existing algorithms, while retaining metadata
in the original snapshot.

## Build and use

From `apps/CCOptimization`:

```sh
make                         # Builds optimizer and mspb_pack.out.
make test-mspb-pack           # Archive API and end-to-end CLI tests.

# Append files in the specified order. Repeating the command appends duplicates.
./mspb_pack.out pack ensemble run_Success_0.mspb run_Success_1.mspb
./mspb_pack.out inspect ensemble

# Restore records 0 and 1 as restored/snapshot_0000.mspb and _0001.mspb.
mkdir -p restored
./mspb_pack.out unpack ensemble 0 2 restored/snapshot
```

Use a prefix, without `.MSPBPack`, for the archive argument; import arguments are
complete `.mspb` filenames. Parent directories must exist. Ranges must fit entirely
inside the archive. Export and analysis refuse existing output files. A multi-file
import commits each valid file in turn; if a later file fails, earlier appends remain.
Use explicitly ordered filenames or zero-padded names when shell globs determine order.

## C++ API

```cpp
#include "MSPB_Pack.h"

MSPB_Pack pack("ensemble");       // Opens existing archive; missing means empty.
pack.ImportFile("snapshot.mspb");
MSPBSnapshot snapshot = pack.GetConfig(0);
SpherePacking spheres = snapshot.ToSpherePacking();
// snapshot.particles[i].id and .species retain metadata absent from SpherePacking.
pack.AddConfig(snapshot);
pack.ExportFile(1, "copy.mspb");
const auto count = pack.NumConfig();
// pack.Clear();                  // Explicitly discards every stored record.
```

`Open(prefix)` attaches/reopens an archive. `NumConfig()` returns the handle's cached
count. `MSPBSnapshot::Read(stream)` reads one exact MSPB v1 snapshot from a seekable
binary stream; `Write(stream)` writes that format. The reader rejects unsupported
versions/flags, dimensions outside the core's supported range, byte-count mismatches,
truncation/trailing data, duplicate IDs, out-of-range species labels, singular cells,
nonfinite numbers, and nonpositive radii. IDs need not equal array indices. A valid
cell with zero particles can be stored, including zero declared species.

`ToSpherePacking()` constructs a core view with the original basis and radii. Core
normalizes coordinates into the periodic unit cell; the snapshot itself preserves
all coordinate bits, including finite coordinates outside that cell. As with other
core algorithms, numerical box scales must be representable by the core's double
arithmetic. Archive storage does not depend on building a core cell list.

## Direct optimizer output

In the existing `soft_core_stealthy.out multi` random-initialization `ground` workflow,
enter this before `run`:

```text
combined
export pack
run
```

Each accepted configuration is appended to `<prefix>_Success.MSPBPack`. This selects
archive output instead of separate `.mspb` files; `export binary` still writes those
individual files. `export txt` and `export off` retain their behavior. Legacy success
and restart packs remain enabled. MD and the incomplete filename-input workflow are
not extended.

When resuming, the MSPB archive count must match the legacy success count before
optimization starts. To enable archive output on a previous binary-export run, first
import **all** existing successes in success-index order into `<prefix>_Success`:

```sh
./mspb_pack.out pack run_Success run_Success_0.mspb run_Success_1.mspb
```

The two output stores are not one transaction. A failed export or an interruption
between legacy storage and the archive append can leave unequal counts; the next run
refuses to resume with `export pack` until the missing snapshots are supplied or the
run is restarted with a fresh prefix. Count matching checks alignment, not content
identity; do not attach an unrelated archive merely because its count matches.

## Offline PairStat pipeline

```sh
mkdir -p analysis
./mspb_pack.out pairstat ensemble 0 2 analysis/stats 3.0 0.2 0.4 0.2 1000 17
# Outputs analysis/stats_0000.txt and analysis/stats_0001.txt.
```

Arguments after the report prefix are `Kmax dK Rmax dR samples seed`. Each snapshot
feeds the existing `MultiComponentPairStat` adapter and `pairstat/` estimators, giving
combined and per-species structure factors, number variances, and volume-fraction
variances. The CLI uses one OpenMP thread. These are per-snapshot reports, not ensemble
averages. Declared empty species appear in report metadata with zero particles and
no tables. Empty total packings and dimensions other than 2D/3D cannot be analyzed.
See [PAIRSTAT_PIPELINE.md](PAIRSTAT_PIPELINE.md) for estimator assumptions and units.

This provides the requested analysis path within CCOptimization. It does not launch
or modify `apps/PairStatCLI`, whose current input interfaces do not preserve MSPB
species/radii metadata.

## Archive layout and persistence

All integers are fixed-width, unsigned, little-endian. There is no native struct
padding, compression, or checksum.

| Field | Encoding |
| --- | --- |
| Archive signature | 8 bytes: `MSPP\r\n\x1a\n` |
| Version | uint32, value 1 |
| Flags | uint32, value 0 |
| Each record | uint64 payload length, followed by one complete MSPB v1 snapshot |

MSPB payload layout is specified in [EXPORT_FORMAT.md](EXPORT_FORMAT.md). An archive
of N snapshots has 16 bytes of header plus 8 bytes of framing per snapshot, in
addition to the original file bytes. It reduces file count, not payload size.

There is no companion index file: `Open` scans framing to construct an in-memory
offset index, and `GetConfig` validates the requested payload. `inspect` validates
all payloads. Parsing and serialization hold one snapshot in memory; they do not
load the whole ensemble. Existing malformed or truncated archives are rejected,
not treated as empty or silently repaired.

Use one writer. A changed file size invalidates another handle's reads/appends until
it is reopened; same-size replacements and concurrent writers are not detected or
synchronized. Successful appends close/check the stream before updating the handle.
Detected write failures attempt to truncate back to the prior size. This is not an
fsync-backed durable transaction: interruption or power loss can leave a partial
last record, which `Open` rejects. Keep source snapshots/backups for recovery.

## Validation

- Independent binary fixtures and the existing optimizer exporter round-trip
  byte-for-byte through import, random access, and export.
- Tests cover 2D/3D skew cells, arbitrary IDs, unequal radii, empty species/packings,
  append/reopen/clear, index errors, stale handles, malformed headers/records,
  and rejection of invalid appends without modifying existing records.
- CLI tests cover direct optimization output and resume, legacy-to-archive migration,
  output collisions, restoration, and offline statistics against reports computed
  directly on accepted in-memory configurations.
