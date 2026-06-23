#==============================================================================
# modules/apps/weighted/makefile   (migrated consumer)
#
# A simplified PairStatisticsCLI for WEIGHTED / MARKED point patterns:
# weighted structure factor, weighted g2, and local variance of a marked
# configuration.  Self-contained — ships its own Computation / PairCorrelation /
# LocalVariance / WeightedStructureFactor (distinct from the pairstat module's
# same-named files).  Depends on: core only (+ header-only Boost, OpenMP).
# (Legacy weighted/ hardcoded core = /home/jaeukk/30_Codes/new_c++/cores/ )
#
#   make && ./weighted.out
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := weighted.out

SRCS := main.cpp Computation.cpp LocalVariance.cpp PairCorrelation.cpp \
        WeightedStructureFactor.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I.

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o

.PHONY: all clean core-objs
all: $(TARGET)

core-objs:
	$(MAKE) -C $(core) objs

$(TARGET): $(OBJS) core-objs
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
