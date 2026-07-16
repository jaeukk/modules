#==============================================================================
# modules/apps/OrderMetric_PeriodicNetwork/makefile   (migrated consumer)
#
# Two independent drivers built from the legacy NonsphericalParticles project:
#
#   GetChivDispersions.out : orientationally-averaged spectral density chi_V(k)
#                            of nonspherical dispersions (main.cpp).
#                            Depends on: core, nonspherical.
#   NetworkBV.out          : local volume-fraction variance \overline{B}_V of
#                            periodic networks (bv_main.cpp + PeriodicNetworks).
#                            Depends on: core, nonspherical, pairstat, models.
#
# PeriodicNetworks.{h,cpp} guard members behind #if defined USE_SHAPES, so this
# app compiles with -DUSE_SHAPES (the nonspherical module is always shapes-on).
#
#   make               # builds both .out
#   make clean
#==============================================================================
include ../../config.mk

# PeriodicNetworks.cpp defines members its header only declares under USE_SHAPES.
CFLAGS += -DUSE_SHAPES

ODIR := obj/

INC := -I$(core) -I$(nonspherical) -I$(pairstat) -I$(models) -I.

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o $(core)/obj/SpherePacking.o \
             $(core)/obj/ConfigurationIO.o $(core)/obj/RandomGenerator.o $(core)/obj/ParameterInput.o $(core)/obj/etc.o
NS_OBJS   := $(nonspherical)/obj/ParticleShapes.o $(nonspherical)/obj/FormFactors.o \
             $(nonspherical)/obj/Dispersions.o $(nonspherical)/obj/SpectralDensity_Dispersion.o
PS_OBJS   := $(pairstat)/obj/StructureFactor.o $(pairstat)/obj/SpectralDensity.o $(pairstat)/obj/NumberVariance.o \
             $(pairstat)/obj/VolumeFractionVariance.o
MOD_OBJS  := $(models)/obj/GenerateConfigs.o

.PHONY: all clean deps
all: GetChivDispersions.out NetworkBV.out

deps:
	$(MAKE) -C $(core) objs
	$(MAKE) -C $(nonspherical) objs
	$(MAKE) -C $(pairstat) objs
	$(MAKE) -C $(models) objs

# Target 1: chi_V(k) of dispersions -- core + nonspherical only.
GetChivDispersions.out: $(ODIR)main.o deps
	$(CC) -o $@ $(ODIR)main.o $(CORE_OBJS) $(NS_OBJS) $(LFLAGS)

# Target 2: \overline{B}_V of periodic networks -- core + nonspherical + pairstat + models.
NetworkBV.out: $(ODIR)bv_main.o $(ODIR)PeriodicNetworks.o deps
	$(CC) -o $@ $(ODIR)bv_main.o $(ODIR)PeriodicNetworks.o $(CORE_OBJS) $(NS_OBJS) $(PS_OBJS) $(MOD_OBJS) $(LFLAGS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o GetChivDispersions.out NetworkBV.out
