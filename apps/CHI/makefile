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

SRCS := main.cpp LinealPathFunction.cpp SpectralDensity.cpp \
        VolumeFractionVariance.cpp
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
