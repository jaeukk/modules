#==============================================================================
# modules/models/makefile
#
# Point-pattern generation models:
#   GenerateConfigs (lattices, perturbed lattices, vacancies, paracrystals),
#   RandomSequentialAddition (RSA), HardSphereFluids, RandomOrganization.
# Depends on: core.  Optional: pairstat (S(k)-driven generation).
#
# Sources here are the TOP-LEVEL copies (the Models/ subdir GenerateConfigs is a
# stub -- see ../MODULE_SOURCES.md).
#
#   make / make objs                # compile module objects (core only)
#   make WITH_SK=1                  # also enable StructureFactor-based features
#   make clean
#==============================================================================
include ../config.mk

ODIR := obj/

SRCS := GenerateConfigs.cpp RandomSequentialAddition.cpp HardSphereFluids.cpp \
        RandomOrganization.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(models)

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

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
