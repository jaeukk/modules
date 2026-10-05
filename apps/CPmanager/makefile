#==============================================================================
# modules/apps/CPmanager/makefile   (migrated consumer; legacy "EXC_formats")
#
# ConfigPack/SpherePack manager: convert binary ensembles <-> text
# configurations and sphere packings.  Depends on: core only.
#
#   make && ./CPmanager.out
#==============================================================================
include ../../config.mk

ODIR   := obj/
TARGET := CPmanager.out

SRCS := main.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I.

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o $(core)/obj/ConfigurationIO.o \
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
	rm -f $(ODIR)*.o $(ODIR)binary_roundtrip.out $(TARGET)

.PHONY: test
test: $(TARGET) $(ODIR)binary_roundtrip.out
	@tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; ./$(ODIR)binary_roundtrip.out "$$tmp/pack"
	python3 tests/roundtrip.py ./$(TARGET)

$(ODIR)binary_roundtrip.o: tests/binary_roundtrip.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR)binary_roundtrip.out: $(ODIR)binary_roundtrip.o core-objs
	$(CC) -o $@ $< $(CORE_OBJS) $(LFLAGS)
