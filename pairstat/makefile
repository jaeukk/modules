#==============================================================================
# modules/pairstat/makefile
#
# Pair-statistics module (library only): PairCorrelation g(r), StructureFactor
# S(k), SpectralDensity chi_V(k), NumberVariance, VolumeFractionVariance, and
# ContactStatistics.
# Depends on: core.
# The Computation driver + CLI now live in apps/PairStatCLI/.
#
#   make / make objs   # module object files
#   make clean
#==============================================================================
include ../config.mk

ODIR   := obj/

SRCS := PairCorrelation.cpp StructureFactor.cpp SpectralDensity.cpp \
        NumberVariance.cpp VolumeFractionVariance.cpp ContactStatistics.cpp \
        LinealPathFunction.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat)

.PHONY: all objs clean
all: objs

objs: $(OBJS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
