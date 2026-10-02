#==============================================================================
# modules/apps/Statistics/makefile   (migrated consumer)
#
# n-body correlation functions g_n (n = 2, 3, 4) and the tau_n order metric of
# 1D point patterns, plus hole / overlap statistics (gn.out).
# Migrated from the legacy new_c++/Statistics/ consumer.
# Depends on: core.  Validation tools live in tests/ (make -C tests).
#
#   make && ./gn.out gn_1D < input
#==============================================================================
include ../../config.mk

ODIR    := obj/
TARGET  := gn.out

SRCS  := main.cpp gn.cpp tau_n.cpp
OBJS  := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I.

CORE_OBJS  := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o \
              $(core)/obj/ConfigurationIO.o $(core)/obj/RandomGenerator.o $(core)/obj/etc.o

.PHONY: all clean deps
all: $(TARGET)

deps:
	$(MAKE) -C $(core) objs

$(TARGET): $(OBJS) deps
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR)main.o $(ODIR)gn.o: gn.h
$(ODIR)main.o $(ODIR)tau_n.o: tau_n.h

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
