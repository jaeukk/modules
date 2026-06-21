#==============================================================================
# modules/visualization/makefile   (plotting / 3-D display, leaf module)
#
# Implementations migrated from the top-level legacy tree to join the headers
# that were already relocated here (Plots.h, DisplaySpheres.h):
#   math3d.{h,cpp}      standalone 3-D vector/matrix math (no external deps)
#   glframe.h           OpenGL camera frame (needs GL/freeglut, USE_OPENGL)
#   Plots.cpp           MathGL plotting of g2 / S(k) / configs (USE_MATHGL)
#   DisplaySpheres.cpp  freeglut interactive 3-D viewer    (USE_OPENGL)
#
# Depends on: core (etc), pairstat (StructureFactor), potential (Potential),
# misc (ParameterSet, via Plots.h). External (optional): MathGL, OpenGL+freeglut.
#
# HEAVY OPTIONAL LIBS: the actual drawing code is #ifdef-guarded behind
# USE_MATHGL / USE_OPENGL; neither MathGL nor freeglut is installed on this box,
# so by default those branches compile out (math3d.o + small Plots.o /
# DisplaySpheres.o stubs).  Enable with e.g.  make USE_MATHGL=1 USE_OPENGL=1
# (add -I/-L for mgl2 and freeglut as needed).
#
#   make / make objs   # compile the module objects
#   make clean
#==============================================================================
include ../config.mk

ODIR := obj/

SRCS := math3d.cpp Plots.cpp DisplaySpheres.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I$(potential) -I$(misc) -I.

# Legacy graphics toggles (separate from config.mk's USE_CGAL/VORO/TJ).
USE_MATHGL ?= 0
USE_OPENGL ?= 0
PREPRO += $(if $(filter 1,$(USE_MATHGL)),-DUSE_MATHGL) \
          $(if $(filter 1,$(USE_OPENGL)),-DUSE_OPENGL)

.PHONY: all objs clean
all: objs
objs: $(OBJS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
