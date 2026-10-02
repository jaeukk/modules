#==============================================================================
# modules/apps/Statistics/tests/makefile
#
# Validation tools for gn.out (not part of the app):
#   gen.out          generator of 1D perturbed lattices / Poisson patterns (test inputs)
#   check_url_gn.py  compares gn.out output with the exact g_n of the 1D URL (no build)
# Depends on: core, models.
#
#   make && ./gen.out PL < input
#==============================================================================
include ../../../config.mk

ODIR   := obj/
TARGET := gen.out

INC := -I$(core) $(models_INC)

CORE_OBJS  := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o \
              $(core)/obj/ConfigurationIO.o $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
MODEL_OBJS := $(models)/obj/GenerateConfigs.o

.PHONY: all clean deps
all: $(TARGET)

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(models) objs

$(TARGET): $(ODIR)gen.o deps
	$(CC) -o $@ $(ODIR)gen.o $(CORE_OBJS) $(MODEL_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
