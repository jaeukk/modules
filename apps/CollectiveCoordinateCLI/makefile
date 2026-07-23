#==============================================================================
# modules/apps/CollectiveCoordinateCLI/makefile   (migrated consumer)
#
# Ge Zhang's stealthy-hyperuniform ground-state command-line interface.
# Generates and relaxes collective-coordinate ("stealthy") ground states via
# gradient optimization, molecular dynamics, MD annealing and Wang-Landau MC.
# Adapted from the legacy c++/CollectiveCoordinateCLI.cpp (the 'Debug' command
# and its Solvers/Degeneracy tail were dropped in migration).
#
# Depends on: core, pairstat (StructureFactor, PairCorrelation),
#             potential (Potential, StructureOptimization, MD_System, MC_System,
#             CollectiveCoordinatePotential), models (RandomSequentialAddition),
#             misc (Interfaces/ParameterSet headers) and the visualization
#             headers.  The real MathGL/OpenGL plotters are unavailable here, so
#             viz_stubs.cpp provides app-local no-op stand-ins (see README.md).
# External: GSL, NLopt.
#
#   make               # build CollectiveCoordinateCLI.out
#   make objs          # compile this app's objects only (no link)
#   make clean
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := CollectiveCoordinateCLI.out

SRCS := main.cpp CollectiveCoordinateCLI.cpp viz_stubs.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I$(potential) $(models_INC) -I$(misc) \
       -I$(visualization) -I. $(NLOPT_INC)

# Object files of the lower modules this app links against.
CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o \
             $(core)/obj/ConfigurationIO.o $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
PAIR_OBJS := $(pairstat)/obj/StructureFactor.o $(pairstat)/obj/PairCorrelation.o
POT_OBJS  := $(potential)/obj/Potential.o $(potential)/obj/StructureOptimization.o \
             $(potential)/obj/MD_System.o $(potential)/obj/MC_System.o
MODEL_OBJS := $(models)/obj/RandomSequentialAddition.o

.PHONY: all objs clean deps
all: $(TARGET)

objs: $(OBJS)

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(pairstat) objs
	$(MAKE) -C $(potential) objs
	$(MAKE) -C $(models) objs

$(TARGET): $(OBJS) deps
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(PAIR_OBJS) $(POT_OBJS) $(MODEL_OBJS) $(LFLAGS) $(NLOPT_LIB)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
