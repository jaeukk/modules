#==============================================================================
# modules/apps/inverse/makefile   (CC-inversion application)
#
# Collective-coordinate inversion problem: reconstruct point configurations from
# their collective coordinates.  This is an APPLICATION (uses the library
# modules), placed under apps/ per the project layout.
#
# Contents:
#   CCInversion.h          CC-inversion potentials (FrankInvertCCPotential*),
#                          SimpleHarmonic*, HS_pot.
#   CCInversionProblem.*   the inversion-problem drivers (FrankInvCCUniqueTest*,
#                          CCInversionProblem()).
#   LatticeUtil.h          GetLattice_config + Find_LeastCommonDivider, relocated
#                          here from the retired etc1.h (only consumer).
#   main.cpp               CC-inversion entry point (the "CCInversion" command
#                          extracted from the legacy Interfaces.h cli()).
#
# Depends on: core, pairstat (StructureFactor), potential (Potential,
# StructureOptimization), and the support module misc (Discrepancy -> MergeSort,
# ParameterSet, Interfaces).  External: GSL, NLopt.
# The "CCInversion" command (CCInversionProblem1, used by main) does NOT need the
# visualization implementations, so a runnable executable links without MathGL/
# OpenGL.
#
#   make / make all    # build CCInversion.out
#   make objs          # compile this app's objects only (no link)
#   make clean
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := CCInversion.out

SRCS := main.cpp CCInversionProblem.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

# Discrepancy.cpp (MergeSort etc.) lives in the misc module; compile it here.
MISC_OBJS := $(ODIR)Discrepancy.o

INC := -I$(core) -I$(pairstat) -I$(potential) -I$(misc) -I$(visualization) -I. $(NLOPT_INC)

# Object files of the lower modules this app links against.
CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
PAIRSTAT_OBJS := $(pairstat)/obj/StructureFactor.o
POTENTIAL_OBJS := $(potential)/obj/Potential.o $(potential)/obj/StructureOptimization.o

.PHONY: all objs clean core-objs pairstat-objs potential-objs
all: $(TARGET)

objs: $(OBJS) $(MISC_OBJS)

core-objs:
	$(MAKE) -C $(core) objs
pairstat-objs:
	$(MAKE) -C $(pairstat) objs
potential-objs:
	$(MAKE) -C $(potential) objs

$(TARGET): $(OBJS) $(MISC_OBJS) core-objs pairstat-objs potential-objs
	$(CC) -o $@ $(OBJS) $(MISC_OBJS) $(CORE_OBJS) $(PAIRSTAT_OBJS) $(POTENTIAL_OBJS) $(LFLAGS) $(NLOPT_LIB)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

# misc-module source compiled into this app's obj dir.
$(ODIR)Discrepancy.o: $(misc)/Discrepancy.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
