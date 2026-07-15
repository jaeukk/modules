#==============================================================================
# modules/apps/CHI/makefile   (migrated consumer)
#
# Two-phase-media descriptors of packings: lineal-path function L(z), spectral
# density chi(k), volume-fraction variance.  Depends on: core only.
# (Legacy CHI/ used a BROKEN core path: /home/jaeukk/30_Codes/new_c++/cores/ )
#
#   make && ./CHI.out
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := CHI.out

SRCS := main.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I.

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o $(core)/obj/ConfigurationIO.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o

PAIRSTAT_OBJS := $(pairstat)/obj/SpectralDensity.o $(pairstat)/obj/StructureFactor.o $(pairstat)/obj/LinealPathFunction.o

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
