#==============================================================================
# modules/potential/makefile
#
# Potentials, dynamics & optimization:
#   Potential, LJPotential, MD_System, MC_System, MonteCarlo (h), Move (h),
#   StructureOptimization, CollectiveCoordinatePotential (h), AnalyzeData (h).
# Depends on: core.  External: GSL, NLopt, OpenMP.
#
# This is a library module (no main); `make` compiles all objects, which is the
# self-consistency / against-core compile check.
#
#   make / make objs   # compile module objects
#   make clean
#==============================================================================
include ../config.mk

ODIR := obj/

SRCS := Potential.cpp LJPotential.cpp MD_System.cpp MC_System.cpp \
        StructureOptimization.cpp LuksanPLIS.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(potential) $(NLOPT_INC)

.PHONY: all objs clean
all: objs
objs: $(OBJS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

# LuksanPLIS.cpp promises NLopt's arithmetic bit for bit: never contract a*b+c into a fused multiply-add,
# whatever architecture flags a build adds.
$(ODIR)LuksanPLIS.o: CFLAGS += -ffp-contract=off

# No per-object header tracking here, so every object depends on every header of this module: a changed
# class layout (for example a new virtual function in Potential.h) must rebuild all of them.
$(OBJS): $(wildcard *.h)

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
