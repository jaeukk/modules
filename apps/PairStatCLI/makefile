#==============================================================================
# modules/apps/PairStatCLI/makefile   (migrated consumer)
#
# PairStatisticsCLI: ensemble pair-statistics driver for point configurations
# (g(r), S(k), number variance, contact statistics).  Self-contained app split
# out of the pairstat module: ships the Computation driver class + cli.cpp main.
# Depends on: core + pairstat (links their object files).
#
#   make && ./PairStatCLI.out
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := PairStatCLI.out

SRCS := cli.cpp Computation.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I.

# Objects of the lower modules this app links against.
CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
PAIRSTAT_OBJS := $(pairstat)/obj/PairCorrelation.o $(pairstat)/obj/StructureFactor.o \
                 $(pairstat)/obj/NumberVariance.o $(pairstat)/obj/ContactStatistics.o

.PHONY: all clean core-objs pairstat-objs
all: $(TARGET)

core-objs:
	$(MAKE) -C $(core) objs

pairstat-objs:
	$(MAKE) -C $(pairstat) objs

$(TARGET): $(OBJS) core-objs pairstat-objs
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(PAIRSTAT_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
