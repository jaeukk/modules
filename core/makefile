#==============================================================================
# modules/core/makefile
#
# Foundation module: GeometryVector, PeriodicCellList, RandomGenerator, etc.
# Builds the object files (for consumers to link) and a small smoke-test
# executable from smoke.cpp.
#
#   make          # build objects + smoke.out
#   make objs     # build object files only (what other modules link against)
#   make clean
#==============================================================================
include ../config.mk

ODIR   := obj/
TARGET := smoke.out

# Public translation units of this module.
SRCS := GeometryVector.cpp PeriodicCellList.cpp SpherePacking.cpp \
        ConfigurationIO.cpp RandomGenerator.cpp etc.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

# This module only needs its own headers on the include path.
INC  := -I$(core)

.PHONY: all objs clean
all: $(TARGET)

objs: $(OBJS)

$(TARGET): $(OBJS) $(ODIR)smoke.o
	$(CC) -o $@ $^ $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
