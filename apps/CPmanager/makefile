#==============================================================================
# modules/apps/CPmanager/makefile   (migrated consumer; legacy "EXC_formats")
#
# ConfigPack file-format manager: convert ConfigPack <-> txt-format point
# configurations.  Depends on: core only.
#
#   make && ./CPmanager.out
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := CPmanager.out

SRCS := main.cpp
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
