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
# StructureOptimization), and the header-only support modules misc (Discrepancy,
# ParameterSet, Interfaces) and visualization (Plots, DisplaySpheres).
# External: GSL, NLopt.
# The CCInversionProblem() driver references Plots/Interfaces declarations
# (compile-safe headers); building a runnable executable would additionally need
# the visualization implementations (MathGL/OpenGL).  This makefile compiles the
# objects (the compile check); it does not link an executable yet.
#
#   make / make objs   # compile the module objects
#   make clean
#==============================================================================
include ../../config.mk

ODIR := obj/

SRCS := main.cpp CCInversionProblem.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I$(potential) -I$(misc) -I$(visualization) -I. $(NLOPT_INC)

.PHONY: all objs clean
all: objs
objs: $(OBJS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
