#==============================================================================
# modules/apps/FFH/makefile   (migrated consumer)
#
# Fast-Fourier-transform homogenization of two-phase media.
# Depends on: core, pairstat (Mtilde; GetKs comes from core), FFTW3.
#
# FFTW3: this machine has no system fftw3.h, so we point FFTW_PREFIX at the conda
# 'pmp' env which ships fftw3.h + libfftw3{,_omp}.so. Override on another machine:
#   make FFTW_PREFIX=/path/with/include/fftw3.h
# Set BEFORE including config.mk so config.mk's `?=` keeps this value.
#==============================================================================
FFTW_PREFIX ?= /home/jaeukk/miniconda3/envs/pmp

include ../../config.mk

ODIR   := obj/
TARGET := GetEffEps.out

SRCS := main.cpp Tensors.cpp FFH.cpp TwoPhaseDispersion.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(pairstat) -I. $(FFTW_INC)

CORE_OBJS := $(core)/obj/GeometryVector.o $(core)/obj/PeriodicCellList.o \
             $(core)/obj/RandomGenerator.o $(core)/obj/etc.o

ifeq ($(FFTW_PREFIX),)
$(warning FFH: fftw3.h not found -- install libfftw3-dev or set FFTW_PREFIX)
endif

.PHONY: all clean core-objs
all: $(TARGET)

core-objs:
	$(MAKE) -C $(core) objs

$(TARGET): $(OBJS) core-objs
	$(CC) -o $@ $(OBJS) $(CORE_OBJS) $(LFLAGS) $(FFTW_LIB)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o $(TARGET)
