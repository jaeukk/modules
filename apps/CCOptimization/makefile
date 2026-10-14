#==============================================================================
# modules/apps/CCOptimization/makefile   (migrated consumer)
#
# Collective-coordinate optimization (soft-core stealthy, etc.): generate ground
# states of CC potentials via optimization + MD.
# Depends on: core, pairstat (StructureFactor + Mtilde, PairCorrelation),
#             potential (CollectiveCoordinatePotential, MD_System,
#             StructureOptimization, Potential).  External: GSL, NLopt.
#
# NOTE: only main.cpp is compiled (DebugSimple.cpp is standalone scratch, and the
# Voronoi/TilingProblem/CCInversionProblem/DisplaySpheres includes in main.cpp
# are inside dead #if 0 blocks -- matching the legacy makefile's object list).
#
#   make && ./soft_core_stealthy.out
#==============================================================================
include ../../config.mk

# USE_FINUFFT=1 compiles the optional FINUFFT backend (select it at run time with the
# `finufft` keyword). CCPOTENTIAL_USE_FINUFFT must match in every object of one program (ODR),
# so FINUFFT builds keep their objects in a separate directory.
USE_FINUFFT ?= 0
ifeq ($(USE_FINUFFT),1)
  ifeq ($(FINUFFT_PREFIX),)
    $(error USE_FINUFFT=1 but finufft.h was not found; set FINUFFT_PREFIX=/path)
  endif
  ODIR    := obj_finufft/
  FINUFFT := -DCCPOTENTIAL_USE_FINUFFT $(FINUFFT_INC)
  LIBS_FINUFFT := $(FINUFFT_LIB)
else
  ODIR    := obj/
endif
TARGET := soft_core_stealthy.out

SRCS := main.cpp MultiStealthyPotential.cpp MultiComponentPairStat.cpp MSPB_Pack.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I$(potential) -I. $(NLOPT_INC) $(FINUFFT)

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o $(core)/obj/ConfigurationIO.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o
PAIR_OBJS := $(pairstat)/obj/StructureFactor.o $(pairstat)/obj/PairCorrelation.o \
             $(pairstat)/obj/NumberVariance.o $(pairstat)/obj/VolumeFractionVariance.o
POT_OBJS  := $(potential)/obj/Potential.o $(potential)/obj/StructureOptimization.o \
             $(potential)/obj/MD_System.o

.PHONY: all clean deps
all: $(TARGET) mspb_pack.out

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(pairstat) objs
	$(MAKE) -C $(potential) objs

$(TARGET): $(OBJS) deps
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(PAIR_OBJS) $(POT_OBJS) $(LFLAGS) $(NLOPT_LIB) $(LIBS_FINUFFT)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -MMD -MP -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f obj/*.[od] obj_finufft/*.[od] $(TARGET) mspb_pack.out bench_finufft.out test_finufft_multi.out test_indexed_view.out

# Offline archive management and multicomponent analysis stay local to this app.
mspb_pack.out: $(ODIR)mspb_pack_main.o $(ODIR)MSPB_Pack.o $(ODIR)MultiStealthyPotential.o $(ODIR)MultiComponentPairStat.o deps
	$(CC) -o $@ $(filter %.o,$^) $(CORE_OBJS) $(PAIR_OBJS) $(POT_OBJS) $(LFLAGS) $(NLOPT_LIB) $(LIBS_FINUFFT)

.PHONY: test-mspb-pack
test-mspb-pack: all
	bash tests/run_mspb_pack_test.sh
	python3 tests/test_mspb_pack_cli.py

# Per-evaluation timing, backend agreement and optimization benchmark (see tests/bench_finufft.py).
bench_finufft.out: $(ODIR)bench_finufft.o $(ODIR)MultiStealthyPotential.o deps
	$(CC) -o $@ $(filter %.o,$^) $(CORE_OBJS) $(PAIR_OBJS) $(POT_OBJS) $(LFLAGS) $(NLOPT_LIB) $(LIBS_FINUFFT)

$(ODIR)bench_finufft.o: tests/bench_finufft.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -MMD -MP -o $@ $<

# Shared FINUFFT correctness oracle (Direct vs FINUFFT, repeatability, cache invalidation).
test_finufft_multi.out: $(ODIR)test_finufft_multi.o $(ODIR)MultiStealthyPotential.o deps
	$(CC) -o $@ $(filter %.o,$^) $(CORE_OBJS) $(PAIR_OBJS) $(POT_OBJS) $(LFLAGS) $(NLOPT_LIB) $(LIBS_FINUFFT)

$(ODIR)test_finufft_multi.o: tests/test_finufft_multi.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -MMD -MP -o $@ $<

.PHONY: test-finufft
test-finufft:
	$(MAKE) USE_FINUFFT=1 test_finufft_multi.out
	./test_finufft_multi.out

# Indexed species views versus the former per-species sub-Configurations (bitwise; Direct and FINUFFT).
test_indexed_view.out: $(ODIR)test_indexed_view.o $(ODIR)MultiStealthyPotential.o deps
	$(CC) -o $@ $(filter %.o,$^) $(CORE_OBJS) $(PAIR_OBJS) $(POT_OBJS) $(LFLAGS) $(NLOPT_LIB) $(LIBS_FINUFFT)

$(ODIR)test_indexed_view.o: tests/test_indexed_view.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -MMD -MP -o $@ $<

.PHONY: test-indexed-view
test-indexed-view:
	$(MAKE) USE_FINUFFT=1 test_indexed_view.out
	./test_indexed_view.out

# Header dependencies (CollectiveCoordinatePotential.h changes must rebuild every object).
-include $(wildcard $(ODIR)*.d)
