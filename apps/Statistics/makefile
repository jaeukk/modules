#==============================================================================
# modules/apps/Statistics/makefile   (migrated consumer)
#
# n-body correlation functions g_n (n = 2, 3, 4) and the tau_n order metric of
# 1D point patterns, plus hole / overlap statistics (gn.out), and a small
# generator of 1D perturbed lattices / Poisson patterns (gen.out).
# Migrated from the legacy new_c++/Statistics/ consumer.
# Depends on: core (gn.out); core + models (gen.out).
#
#   make && ./gn.out gn_1D < input
#==============================================================================
include ../../config.mk

ODIR    := obj/
TARGET  := gn.out
TARGET2 := gen.out

SRCS  := main.cpp gn.cpp tau_n.cpp
OBJS  := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))
OBJS2 := $(ODIR)gen.o

INC := -I$(core) $(models_INC) -I.

CORE_OBJS  := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o \
              $(core)/obj/ConfigurationIO.o $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
MODEL_OBJS := $(models)/obj/GenerateConfigs.o

.PHONY: all clean deps
all: $(TARGET) $(TARGET2)

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(models) objs

$(TARGET): $(OBJS) deps
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(LFLAGS)

$(TARGET2): $(OBJS2) deps
	$(CC) -o $@ $(OBJS2) $(CORE_OBJS) $(MODEL_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR)main.o $(ODIR)gn.o: gn.h
$(ODIR)main.o $(ODIR)tau_n.o: tau_n.h

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET) $(TARGET2)
