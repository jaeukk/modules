# CPmanager — ConfigurationPack file manager

A small command-line utility to **manage ConfigPack files** — the binary
container (`core`'s `ConfigurationPack`) that stores an *ensemble* of point
configurations / sphere packings in a single pair of files
(`<prefix>.ConfigPack` index + data). CPmanager lets you inspect a pack and
convert between the packed binary form and plain per-configuration text files.

Depends on **`core` only**.

## Build

```sh
make            # -> CPmanager.out
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
pack to text files named `<ConfigPack prefix>__%04d` (e.g. `myPack__0000`,
`myPack__0001`, …). The count is clamped to the configurations available.

### `2` — txt2CP : pack text files into a ConfigPack

```sh
./CPmanager.out 2 <txt prefix> <start index> <num configs> <save prefix>
```

Reads text configurations named by `<txt prefix>` (a `printf`-style pattern
containing `%d`, filled with `start index … start index + num configs − 1`) and
**appends** them to the ConfigPack `<save prefix>`, creating it if needed.
Reading stops early at the first missing/empty file.

## Examples

```sh
# inspect a pack
./CPmanager.out 0 data/2D_RSA_0.5

# export the first 10 configurations to text
./CPmanager.out 1 data/2D_RSA_0.5 0 10        # -> data/2D_RSA_0.5__0000 ...

# build a pack from text files snap_0.txt ... snap_49.txt
./CPmanager.out 2 snap_%d.txt 0 50 data/snaps
```

## Notes

- The text I/O uses `core`'s `WriteConfiguration` / `ReadConfiguration`; the pack
  I/O uses `ConfigurationPack::{GetConfig, AddConfig, NumConfig}`.
- `txt2CP` **appends** — point it at a fresh `<save prefix>` to avoid mixing
  ensembles.
