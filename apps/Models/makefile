#==============================================================================
# modules/apps/Models/makefile   (migrated consumer)
#
# CLI that generates ensembles of selected disordered models (Poisson, RSA, HSF,
# Perturbations) and writes them as a ConfigurationPack or SpherePacking text
# files.  Adapted from the legacy c++/Models/CLI.cpp.
# Depends on: core, models (GenerateConfigs + RSA + HardSphereFluids).
#
#   make && ./models.out -h
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := models.out

SRCS := main.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

# $(models_INC) expands to the flat per-model include paths (config.mk).
INC := -I$(core) $(models_INC) -I.

CORE_OBJS  := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o \
              $(core)/obj/ConfigurationIO.o $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
MODEL_OBJS := $(models)/obj/GenerateConfigs.o $(models)/obj/HardSphereFluids.o \
              $(models)/obj/RandomSequentialAddition.o

.PHONY: all clean deps
all: $(TARGET)

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(models) objs

$(TARGET): $(OBJS) deps
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(MODEL_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
