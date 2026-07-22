#==============================================================================
# modules/models/makefile
#
# Point-pattern generation models, organized into per-model subdirs:
#   ./            GenerateConfigs (ConfigGen base, PoissonPtGen, predefined
#                 lattices, numK_points) + Lattices/Lattice.cpp (LatticeGen impl,
#                 #included into GenerateConfigs.o)
#   RSA/          RandomSequentialAddition (voxel-split RSA) + RSA_Gen.h
#   HSF/          HardSphereFluids (MC hard spheres) + HSF_Gen.h
#   Pert/         GenPerturb + GenPtsDefect (header-only)
#   Lattices/     Lattice.h/.cpp (LatticeGen)
#   RandOrg/      RandomOrganization (Berthier active-passive model)
# Depends on: core.  Optional: pairstat (S(k)-driven generation).
#
# The compiled objects keep flat basenames in obj/ (GenerateConfigs.o,
# RandomSequentialAddition.o, HardSphereFluids.o, RandomOrganization.o) so that
# consumers linking $(models)/obj/*.o are unaffected by the subdir layout.
# Flat header includes (<RSA_Gen.h>, <Lattice.h>, ...) come from $(models_INC),
# defined in config.mk.
#
#   make / make objs                # compile module objects (core only)
#   make WITH_SK=1                  # also enable StructureFactor-based features
#   make clean
#==============================================================================
include ../config.mk

ODIR := obj/

SRCS := GenerateConfigs.cpp RSA/RandomSequentialAddition.cpp \
        HSF/HardSphereFluids.cpp RandOrg/RandomOrganization.cpp
OBJS := $(addprefix $(ODIR),$(notdir $(SRCS:.cpp=.o)))

# Resolve the flat object basenames back to their subdir sources.
vpath %.cpp . RSA HSF RandOrg

INC := $(models_INC) -I$(core)

# Opt-in: structure-factor-driven configuration generation (pulls in pairstat).
WITH_SK ?= 0
ifeq ($(WITH_SK),1)
INC    += -I$(pairstat)
EXTRA  := -DSTRUCTUREFACTOR_INCLUDED
endif

.PHONY: all objs clean
all: objs
objs: $(OBJS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(EXTRA) $(INC) -o $@ $<

# GenerateConfigs.cpp #includes the Lattices/Lattice.cpp translation-unit fragment
# (LatticeGen impl) and, via GenerateConfigs.h, the split subdir headers.  Declare
# those as explicit prerequisites so the object is rebuilt when they change.
$(ODIR)GenerateConfigs.o: GenerateConfigs.h Lattices/Lattice.cpp Lattices/Lattice.h Pert/GenPerturb.h

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
