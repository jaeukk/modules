#==============================================================================
# modules/apps/CCOptimization/makefile   (migrated consumer)
#
# Collective-coordinate optimization (soft-core stealthy, etc.): generate ground
# states of CC potentials via optimization + MD.
# Depends on: core, pairstat (StructureFactor + Mtilde, PairCorrelation),
#             potential (CollectiveCoordinatePotential, MD_System,
#             StructureOptimization, Potential).  External: GSL, NLopt.
#
# NOTE: only main.cpp is compiled (DebugSimple.cpp is standalone scratch, and the
# Voronoi/TilingProblem/CCInversionProblem/DisplaySpheres includes in main.cpp
# are inside dead #if 0 blocks -- matching the legacy makefile's object list).
#
#   make && ./soft_core_stealthy.out
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := soft_core_stealthy.out

SRCS := main.cpp MultiStealthyPotential.cpp MultiComponentPairStat.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I$(potential) -I. $(NLOPT_INC)

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o $(core)/obj/ConfigurationIO.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
PAIR_OBJS := $(pairstat)/obj/StructureFactor.o $(pairstat)/obj/PairCorrelation.o \
             $(pairstat)/obj/NumberVariance.o $(pairstat)/obj/VolumeFractionVariance.o
POT_OBJS  := $(potential)/obj/Potential.o $(potential)/obj/StructureOptimization.o \
             $(potential)/obj/MD_System.o

.PHONY: all clean deps
all: $(TARGET)

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(pairstat) objs
	$(MAKE) -C $(potential) objs

$(TARGET): $(OBJS) deps
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(PAIR_OBJS) $(POT_OBJS) $(LFLAGS) $(NLOPT_LIB)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
