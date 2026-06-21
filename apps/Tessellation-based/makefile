#==============================================================================
# modules/apps/Tessellation-based/makefile   (tiling / space-filling-packing)
#
# Migrated from the legacy top-level tiling sources.  The legacy monolithic
# TilingProblem.cpp was split into:
#   SpacefillingPacking.{h,cpp}  space-filling multiscale sphere packing
#   Computation.{h,cpp}          the Computation analysis-task hierarchy
#   main.cpp                     TilingProblemCLI + the "Tiling" entry point
# plus Voronoi.{h,cpp} (relocated as-is) and PolyPacking.{h,cpp}
# (the MakeMono2Poly* helpers extracted from the retired etc1.{h,cpp}).
#
# Depends on: core, pairstat (StructureFactor/NumberVariance/PairCorrelation),
# models (GenerateConfigs/RSA/HardSphereFluids), misc (Discrepancy, only under
# USE_VORO) and visualization (DisplaySpheres). External: GSL.
#
# HEAVY OPTIONAL LIBS: the Voronoi analyses need CGAL (USE_CGAL), voro++
# (USE_VORO) and MathGL (USE_MATHGL) — none are installed on this box, so by
# default those branches are #ifdef-guarded out.  This makefile compiles the
# objects (the compile check); linking a fully-functional binary needs those
# libraries.  Enable e.g. with:  make USE_CGAL=1 USE_VORO=1
#
#   make / make objs   # compile the module objects
#   make clean
#==============================================================================
include ../../config.mk

ODIR := obj/

SRCS := SpacefillingPacking.cpp Computation.cpp PolyPacking.cpp Voronoi.cpp main.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I$(models) -I$(misc) -I$(visualization) -I.

# MathGL is gated by USE_MATHGL in the legacy sources (not the config.mk flags);
# surface it here so `make USE_MATHGL=1 ...` works once MathGL is installed.
USE_MATHGL ?= 0
PREPRO += $(if $(filter 1,$(USE_MATHGL)),-DUSE_MATHGL)

.PHONY: all objs clean
all: objs
objs: $(OBJS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
