#==============================================================================
# modules/apps/FirstPassageTime/makefile   (migrated consumer)
#
# First-passage-time / mean-survival-time computations on dispersions of
# nonspherical inclusions.  Depends on: core, nonspherical.
# Legacy reached the library via ../*.h and ../NonsphericalParticles/*.h.
#
#   make && ./FirstPassageTime.out
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := FirstPassageTime.out

SRCS := main.cpp FirstPassageTime.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(nonspherical) -I.

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o $(core)/obj/ConfigurationIO.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
NS_OBJS   := $(nonspherical)/obj/ParticleShapes.o $(nonspherical)/obj/FormFactors.o \
             $(nonspherical)/obj/Dispersions.o

.PHONY: all clean deps
all: $(TARGET)

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(nonspherical) objs

$(TARGET): $(OBJS) deps
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(NS_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
