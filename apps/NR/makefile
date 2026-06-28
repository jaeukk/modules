#==============================================================================
# modules/apps/NR/makefile   (migrated consumer)
#
# Number-restricted / nearest-neighbour distributions of point patterns.
# Depends on: core, pairstat (NumberVariance), models (GenerateConfigs + RSA +
# HardSphereFluids).   Legacy NR/ reached the library via ../*.cpp.
#
#   make && ./NR.out
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := NR.out

# main.cpp + NR_distribution.cpp + examples.cpp (example routines).
# distance-metrics.cpp stays excluded (standalone scratch, not wired into NR).
SRCS := main.cpp NR_distribution.cpp examples.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I$(models) -I.

CORE_OBJS  := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o $(core)/obj/ConfigurationIO.o \
              $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
PAIR_OBJS  := $(pairstat)/obj/NumberVariance.o
MODEL_OBJS := $(models)/obj/GenerateConfigs.o $(models)/obj/HardSphereFluids.o \
              $(models)/obj/RandomSequentialAddition.o

.PHONY: all clean deps
all: $(TARGET)

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(pairstat) objs
	$(MAKE) -C $(models) objs

$(TARGET): $(OBJS) deps
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(PAIR_OBJS) $(MODEL_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
