# Configuration export formats

`MultiStealthyPotential::ExportConfiguration(stream, format)` exports the latest
configuration supplied to the potential. Refresh it with `SetConfiguration(final)`
after optimization. The format is explicit, not selected by the filename.

- `ExportFormat::Text` (default), extension `.txt`: comment-prefixed cell metadata
  and a tab-separated particle table. Coordinates are relative; radii are Cartesian.
- `ExportFormat::Binary`, extension `.mspb`: **Multi-Stealthy Packing Binary**,
  version 1. The format preserves the cell, particle IDs, species labels, relative
  coordinates, and individual radii. It is a configuration snapshot, not an
  optimizer checkpoint: potential constraints and optimizer state are not stored.

```cpp
potential.SetConfiguration(final_configuration);
std::ofstream out("ground_state.mspb", std::ios::binary);
potential.ExportConfiguration(out, MultiStealthyPotential::ExportFormat::Binary);
out.close();
if (!out) throw std::runtime_error("Could not finish writing ground_state.mspb");
```

## MSPB version 1 layout

All integers are unsigned, fixed-width, little-endian. All floating-point values
are IEEE-754 binary64, also little-endian. Fields are contiguous with no padding,
compression, or checksum. The byte order is independent of the writer's machine.

| Field | Encoding |
|---|---|
| Signature | 8 bytes: `4d 53 50 42 0d 0a 1a 0a` (`MSPB\r\n\x1a\n`) |
| Version | uint32, value 1 |
| Dimension D | uint32 |
| Particle count N | uint64 |
| Declared species count | uint64, including empty species |
| Flags | uint32, value 3: bit 0 = relative coordinates; bit 1 = Cartesian radii |
| Basis vectors | D × D binary64 values: vector 0's Cartesian components, then vector 1, etc. |
| Particle records | N records, in global particle-index order, as specified below |

Each particle record contains uint64 ID, uint32 species label, D binary64 relative
coordinates, and one binary64 Cartesian radius. Labels are zero-based. Cartesian
positions are the linear combination of basis vectors weighted by the relative
coordinates. No square or orthogonal-cell assumption is made.

Total size: `36 + 8*D*D + N*(20 + 8*D)` bytes. Readers should reject unknown
versions/flags, invalid dimensions/counts/labels, truncated or trailing data, and
nonfinite values before constructing a configuration. A production reader and collection API are provided by `MSPBSnapshot::Read` and
`MSPB_Pack`; see [MSPB_PACK.md](MSPB_PACK.md).

The existing text default and old archived sweep files are unchanged by adding
binary support. Binary output streams must be opened with `std::ios::binary`.

Run regression tests with `bash tests/run_export_test.sh` from this directory.

## Multicomponent CLI

In `soft_core_stealthy.out multi`, choose `ground`, then enter `export txt` or
`export binary` before `run`. `export off` is the default, preserving legacy runs.
Each newly accepted configuration is exported as `<prefix>_Success_<index>.txt`
or `.mspb`, with a zero-based index in the success pack. These files contain all
components regardless of the legacy `combined` setting. Export runs directly on
the accepted in-memory configuration; it does not reread the configuration pack.
Existing configuration-pack/restart output is retained. MD export is not added.

For collection output, use `export pack` to append accepted configurations to
`<prefix>_Success.MSPBPack`. The archive supports indexed retrieval, import/export
of individual snapshots, and offline PairStat analysis; see [MSPB_PACK.md](MSPB_PACK.md)
for API, CLI examples, layout, and resume requirements.
