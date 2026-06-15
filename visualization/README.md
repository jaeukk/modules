# visualization

Rendering / plotting **declarations** for configurations and pair statistics.
Header-only here: the implementations rely on heavy optional libraries (MathGL
for 2D plots, OpenGL/GLU/GLUT for 3D sphere display) and still live in the
legacy tree, so this module only exposes the headers that consumers compile
against (via `-I$(visualization)`). Linking a binary that actually draws needs
those legacy `.cpp` files plus MathGL/OpenGL.

| File | Purpose |
|------|---------|
| `Plots.h` | 2D plots: `Plot`, directional structure factor, packings, function plots (MathGL / Grace backends). |
| `DisplaySpheres.h` | Interactive 3D OpenGL display of spheres / configurations / structure factor, and movies. |

## Dependencies
- **Lower modules:** `core` (PeriodicCellList), `misc` (`Plots.h` → ParameterSet.h).
- **External (link only):** MathGL, OpenGL/GLU/GLUT — guarded behind the heavy
  includes; not required to *compile* a consumer that merely declares these.

## Notes
- Header-only: no makefile. This is the `viz` tier from the plan, renamed
  `visualization`.
