#==============================================================================
# modules/pairstat/makefile
#
# Pair-statistics module: PairCorrelation g(r), StructureFactor S(k),
# NumberVariance, Computation driver, ContactStatistics.
# Depends on: core.
#
#   make          # build module objects + pairstat.out (CLI from cli.cpp)
#   make objs     # module object files only
#   make clean
#==============================================================================
include ../config.mk

ODIR   := obj/
TARGET := pairstat.out

SRCS := PairCorrelation.cpp StructureFactor.cpp NumberVariance.cpp \
        Computation.cpp ContactStatistics.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

# Objects of the core module this module links against.
CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o

INC := -I$(core) -I$(pairstat)

.PHONY: all objs clean core-objs
all: $(TARGET)

objs: $(OBJS)

# Ensure the core module's objects exist before linking.
core-objs:
	$(MAKE) -C $(core) objs

$(TARGET): $(OBJS) $(ODIR)cli.o core-objs
	$(CC) -o $@ $(OBJS) $(ODIR)cli.o $(CORE_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
