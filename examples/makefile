#==============================================================================
# modules/examples/makefile
#
# Proof consumer that builds purely against the modular library (core +
# pairstat).  This is the reference pattern for migrating real subprojects:
# include ../config.mk, use -I$(module) include paths, and link the module
# objects built by `make -C $(module) objs`.
#
#   make && ./g2_sk_demo.out
#==============================================================================
include ../config.mk

ODIR   := obj/
TARGET := g2_sk_demo.out

INC := -I$(core) -I$(pairstat)

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
PAIRSTAT_OBJS := $(pairstat)/obj/PairCorrelation.o \
                 $(pairstat)/obj/StructureFactor.o \
                 $(pairstat)/obj/NumberVariance.o

.PHONY: all clean deps
all: $(TARGET)

# Build the module objects we link against.
deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(pairstat) objs

$(TARGET): g2_sk_demo.cpp deps | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $(ODIR)g2_sk_demo.o g2_sk_demo.cpp
	$(CC) -o $@ $(ODIR)g2_sk_demo.o $(CORE_OBJS) $(PAIRSTAT_OBJS) $(LFLAGS)

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
