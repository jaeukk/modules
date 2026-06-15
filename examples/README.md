# examples

Reference **consumer** of the modular library — the pattern to copy when
migrating a real subproject.

| File | Demonstrates |
|------|--------------|
| `g2_sk_demo.cpp` | Build a point pattern with `core`, compute g(r) and S(k) with `pairstat`, end-to-end. |

## Build & run

```sh
make && ./g2_sk_demo.out
```

For a Poisson pattern the demo prints g₂ and S(k) tables (≈ 1 away from k→0).

## The migration recipe (see this makefile)

1. `include ../config.mk` at the top of your makefile.
2. Use `-I$(core) -I$(pairstat) …` include paths (never `../File.h` or absolute
   `/home/jaeukk/...` paths).
3. In sources, include module headers with angle brackets: `#include <StructureFactor.h>`.
4. Link against the module objects, building them first via
   `make -C $(core) objs` etc.
5. Do **not** define library globals like `Verbosity` in your `main` — the
   `core` module owns them now.
