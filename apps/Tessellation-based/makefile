#==============================================================================
# modules/apps/Tessellation-based/makefile   (tiling / space-filling-packing)
#
# Migrated from the legacy top-level tiling sources.  The legacy monolithic
# TilingProblem.cpp was split into:
#   SpacefillingPacking.{h,cpp}  space-filling multiscale sphere packing
#   Computation.{h,cpp}          the Computation analysis-task hierarchy
#   main.cpp                     TilingProblemCLI + the "Tiling" entry point
#   TessellationAnalysisCLI.cpp  main.cpp decoupled from the models generators
#                                (reads configs from files; no RSA_Gen/HSF/HSF_Gen)
# plus Voronoi.{h,cpp} (relocated as-is) and PolyPacking.{h,cpp}
# (the MakeMono2Poly* helpers extracted from the retired etc1.{h,cpp}).
#
# HEAVY OPTIONAL LIBS: the Voronoi analyses need CGAL (USE_CGAL), voro++
# (USE_VORO) and MathGL (USE_MATHGL).  By default those branches are
# #ifdef-guarded out and `make objs` only compiles the objects.  With voro++
# installed (voro++-dev; auto-detected via VORO_PREFIX in config.mk) build a
# runnable decoupled CLI with:
#
#   make USE_VORO=1                 # -> TessellationCLI.out (config-reading CLI)
#   make objs                       # compile the library objects only
#   make clean
#
# main.cpp (the generator-coupled CLI) additionally needs models/HSF_Gen.h, which
# has not been migrated yet, so it is not part of the default build.
#==============================================================================
include ../../config.mk

ODIR   := obj/
# The default (and only) built binary is the ANALYSIS-ONLY CLI: it is compiled
# from TessellationAnalysisCLI.cpp (config-reading; decoupled from the models
# generators). The generator-coupled driver main.cpp is intentionally NOT built.
TARGET := TessellationCLI.out

LIB_SRCS := SpacefillingPacking.cpp Computation.cpp PolyPacking.cpp Voronoi.cpp
LIB_OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(LIB_SRCS))
# Analysis-only CLI object (initial configs loaded from files, no internal
# RSA_Gen / HardSphereFluids / HSF_Gen generation).
CLI_OBJ  := $(ODIR)TessellationAnalysisCLI.o

INC := -I$(core) -I$(pairstat) -I$(models) -I$(visualization) -I. $(VORO_INC)

# MathGL is gated by USE_MATHGL in the legacy sources (not the config.mk flags).
USE_MATHGL ?= 0
PREPRO += $(if $(filter 1,$(USE_MATHGL)),-DUSE_MATHGL)

# Objects of the lower modules the linked CLI needs.
CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o \
             $(core)/obj/SpherePacking.o $(core)/obj/ConfigurationIO.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
PAIR_OBJS := $(pairstat)/obj/StructureFactor.o $(pairstat)/obj/SpectralDensity.o \
             $(pairstat)/obj/NumberVariance.o $(pairstat)/obj/VolumeFractionVariance.o \
             $(pairstat)/obj/PairCorrelation.o
MODEL_OBJS := $(models)/obj/RandomSequentialAddition.o $(models)/obj/GenerateConfigs.o

.PHONY: all objs clean deps
all: $(TARGET)

objs: $(LIB_OBJS)

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(pairstat) objs
	$(MAKE) -C $(models) objs

$(TARGET): $(CLI_OBJ) $(LIB_OBJS) deps
	$(CC) -o $@ $(CLI_OBJ) $(LIB_OBJS) $(CORE_OBJS) $(PAIR_OBJS) $(MODEL_OBJS) $(VORO_LIB) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
