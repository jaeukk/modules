# CPmanager — ConfigurationPack and SpherePackingPack manager

A command-line utility to inspect and convert ensembles stored by `core`:

| Container | Data file | Index file | Particle data |
| --- | --- | --- | --- |
| `ConfigurationPack` | `<prefix>.ConfigPack` | `<prefix>.ConfigPackIndex` | Point configurations |
| `SpherePackingPack` | `<prefix>.SpherePack` | `<prefix>.SpherePackIndex` | Sphere positions and individual radii |

Both use the existing `Open`, `NumConfig`, `AddConfig`, `GetConfig`, and
`Clear` API. CPmanager shares its inspect/import/export implementation between
the two types, selecting the matching text I/O overload. SpherePack is binary
ensemble storage, without additional compression. Its native binary layout is
intended for compatible builds/platforms, not as a portable interchange format.

Depends on **`core` only**.

## Build

```sh
make            # -> CPmanager.out
make test       # binary restoration and CLI round-trip tests (requires Python 3)
make clean
```

## Usage

```sh
./CPmanager.out <option> [parameters...]
```

`<option>` selects the mode; the remaining arguments are mode-specific:

### `0` — CheckCP : print statistics of a pack

```sh
./CPmanager.out 0 <ConfigPack prefix>
```

Prints the number of configurations, the spatial dimension, and the
ensemble-averaged particle count and number density.

### `1` — CP2txt : unpack to text files

```sh
./CPmanager.out 1 <ConfigPack prefix> <start index> <num configs>
```

Writes `<num configs>` configurations (starting at `<start index>`) out of the
pack to text files named `<ConfigPack prefix>__%04d.txt` (e.g. `myPack__0000.txt`,
`myPack__0001.txt`, …). The count is clamped to the configurations available.

### `2` — txt2CP : pack text files into a ConfigPack

```sh
./CPmanager.out 2 <txt prefix> <start index> <num configs> <save prefix>
```

Reads text configurations named by `<txt prefix>` (a `printf`-style pattern
containing `%d`, filled with `start index … start index + num configs − 1`) and
**appends** them to the ConfigPack `<save prefix>`, creating it if needed.
The pattern is a filename **prefix**: `.txt` is appended by the reader.
Reading stops early at the first missing file. A valid box with no particles
is stored as an empty configuration.

### `3`, `4`, `5` — SpherePack inspection, export, and import

```sh
./CPmanager.out 3 <SpherePack prefix>
./CPmanager.out 4 <SpherePack prefix> <start index> <num packings>
./CPmanager.out 5 <txt prefix pattern> <start index> <num packings> <save prefix>
```

These mirror modes `0`, `1`, and `2`, using `SpherePackingPack`. Mode `3`
reports count, the first packing's dimension, average particle count, and average
number density. Mode `4` writes `<prefix>__%04d.txt`; mode `5` appends to the
`.SpherePack` / `.SpherePackIndex` pair. Export counts are clamped to available
records. Indices and counts must be non-negative integers.

Sphere text uses the existing `WriteConfiguration(SpherePacking, ...)` format:
a dimension line, one basis-vector row per dimension, then particle rows with
**four Cartesian coordinate columns followed by radius**. Unused coordinate
columns are zero. Text export uses 17 significant digits. Binary restoration
preserves stored doubles exactly; text round trips can incur floating-point
rounding when converting between Cartesian and relative coordinates.

## Examples

Run these commands from `apps/CPmanager` after building. Input patterns omit
`.txt`, which the reader appends automatically. Create the output directory
before importing, and use fresh pack prefixes when starting a new ensemble.

### Point configurations: import, inspect, and restore

Given `snap_0.txt` through `snap_49.txt` in the point-configuration text format:

```sh
mkdir -p data
./CPmanager.out 2 'snap_%d' 0 50 data/snaps
# Creates data/snaps.ConfigPack and data/snaps.ConfigPackIndex.
./CPmanager.out 0 data/snaps
./CPmanager.out 1 data/snaps 0 10
# Writes data/snaps__0000.txt ... data/snaps__0009.txt.
```

### Sphere packings: import, inspect, and restore

Given `sphere_0000.txt` through `sphere_0009.txt` with positions and radii:

```sh
mkdir -p data
./CPmanager.out 5 'sphere_%04d' 0 10 data/spheres
# Creates data/spheres.SpherePack and data/spheres.SpherePackIndex.
./CPmanager.out 3 data/spheres
./CPmanager.out 4 data/spheres 0 10
# Writes data/spheres__0000.txt ... data/spheres__0009.txt.
```

For example, a valid `sphere_0000.txt` containing two spheres in a 2D unit box is:

```text
2
1 0 0 0
0 1 0 0
0.25 0.25 0 0 0.1
0.75 0.75 0 0 0.15
```

To append `sphere_0010.txt` through `sphere_0014.txt`, then export those records:

```sh
./CPmanager.out 5 'sphere_%04d' 10 5 data/spheres
./CPmanager.out 4 data/spheres 10 5
# Writes data/spheres__0010.txt ... data/spheres__0014.txt.
```

To repack the first ten restored files into a separate ensemble:

```sh
./CPmanager.out 5 'data/spheres__%04d' 0 10 data/spheres_restored
./CPmanager.out 3 data/spheres_restored
```

## Notes

- The text I/O uses `core`'s `WriteConfiguration` / `ReadConfiguration`; the pack
  I/O uses `GetConfig`, `AddConfig`, and `NumConfig` on the selected container.
- Both `txt2CP` and `txt2SP` **append** — point it at a fresh `<save prefix>` to avoid mixing
  ensembles.

## Validation

`make test` checks exact binary restoration of box vectors, relative coordinates,
and radii after reopening, appending, and clearing. The CLI test independently
checks text → SpherePack → text → SpherePack → text with polydisperse 2D/3D skew
boxes and an empty packing (numeric tolerance `2e-14`), export range clamping,
invalid arguments, and the original ConfigPack modes. All test data is temporary.
