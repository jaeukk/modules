#==============================================================================
# modules/phonons/makefile   (phonon-frequency module, leaf)
#
# Migrated from the top-level legacy tree:
#   PhononFrequency.{h,cpp}    dynamical matrix + phonon dispersion of a crystal
#                              (diagonalised with GSL complex eigensolvers)
#   PlotPhononFrequencies.cpp  band-structure plotting along k-paths (MathGL)
#
# Depends on: core (GeometryVector/PeriodicCellList), potential (PairPotential),
# GSL. External (optional): MathGL for the plotting (USE_MATHGL).
#
# No main() — this is a library module (objects only).  The plotting code is
# #ifdef-guarded behind USE_MATHGL; MathGL is not installed here, so by default
# PlotPhononFrequencies.o compiles with the plot bodies stubbed out.
#
#   make / make objs   # compile the module objects
#   make clean
#==============================================================================
include ../config.mk

ODIR := obj/

SRCS := PhononFrequency.cpp PlotPhononFrequencies.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(potential) -I.

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
