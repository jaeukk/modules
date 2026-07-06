#==============================================================================
# config.mk  --  shared build configuration for the modules/ tree
#
# Every module and every consumer makefile should `include` this file.  It
# provides:
#   * ROOT / MODS                 absolute paths to the repo and modules/ tree
#   * per-module include vars     $(core) $(pairstat) $(potential) ...
#   * auto-detected library paths $(GSL_PREFIX) $(NLOPT_PREFIX)
#   * canonical compiler flags    $(CFLAGS) $(LFLAGS)
#   * optional-feature toggles    USE_OPENMP / USE_CGAL / USE_VORO / USE_TJ
#
# Goal: a subproject no longer hardcodes /home/jaeukk/... paths.  The library
# prefix is discovered from a candidate list, so the SAME makefile works on the
# HOME laptop (~/lib/gsl), this WSL box (~/install) and the DELLA cluster
# (/usr/local) without edits.  Override any variable on the command line, e.g.
#   make GSL_PREFIX=/opt/gsl
#
# Author: Jaeuk Kim  --  modularization pass
#==============================================================================

#----- Repo layout ------------------------------------------------------------
# Resolve ROOT as the directory that contains this config.mk, then go up one
# (config.mk lives in modules/).  Works regardless of where make is invoked.
MODS := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
ROOT := $(abspath $(MODS)/..)

#----- Per-module include directories -----------------------------------------
core         := $(MODS)/core
pairstat     := $(MODS)/pairstat
models       := $(MODS)/models
potential    := $(MODS)/potential
nonspherical := $(MODS)/nonspherical
inverse      := $(MODS)/inverse
tiling       := $(MODS)/tiling
viz          := $(MODS)/viz
# Header-only support modules split out of the apps:
#   misc          : Discrepancy, ParameterSet, Interfaces (the legacy cli() dispatcher)
#   visualization : Plots, DisplaySpheres (declarations; impls need MathGL/OpenGL)
misc          := $(MODS)/misc
visualization := $(MODS)/visualization
# Leaf library modules migrated from the top-level legacy tree:
#   phonons : PhononFrequency (+ PlotPhononFrequencies); core + potential + GSL
#   lp      : TorquatoJiao linear-programming jamming; needs glpk + Gurobi (USE_TJ)
phonons       := $(MODS)/phonons
lp            := $(MODS)/lp

#----- Library prefix auto-detection ------------------------------------------
# First candidate whose include/gsl/gsl_matrix.h exists wins.
GSL_CANDIDATES   := $(HOME)/install $(HOME)/lib/gsl /usr/local /usr
GSL_PREFIX       ?= $(firstword $(foreach d,$(GSL_CANDIDATES),\
                      $(if $(wildcard $(d)/include/gsl/gsl_matrix.h),$(d))))

NLOPT_CANDIDATES := $(HOME)/install /usr/local /usr
NLOPT_PREFIX     ?= $(firstword $(foreach d,$(NLOPT_CANDIDATES),\
                      $(if $(wildcard $(d)/include/nlopt.h),$(d))))

# FFTW3 (needed by FFH). Requires the dev header fftw3.h, not just libfftw3.so.3.
FFTW_CANDIDATES  := $(HOME)/install /usr/local /usr
FFTW_PREFIX      ?= $(firstword $(foreach d,$(FFTW_CANDIDATES),\
                      $(if $(wildcard $(d)/include/fftw3.h),$(d))))

# voro++ (optional; used by apps/Tessellation-based with USE_VORO=1). Requires the
# dev headers voro++/voro++.hh (the voro++-dev package), not just libvoro++.so.
VORO_CANDIDATES  := $(HOME)/install /usr/local /usr
VORO_PREFIX      ?= $(firstword $(foreach d,$(VORO_CANDIDATES),\
                      $(if $(wildcard $(d)/include/voro++/voro++.hh),$(d))))

ifeq ($(GSL_PREFIX),)
  $(warning config.mk: could not auto-detect GSL; set GSL_PREFIX=/path manually)
endif

#----- External include/link flags --------------------------------------------
GSL_INC   := -I$(GSL_PREFIX)/include
GSL_LIB   := -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
NLOPT_INC := $(if $(NLOPT_PREFIX),-I$(NLOPT_PREFIX)/include)
NLOPT_LIB := $(if $(NLOPT_PREFIX),-L$(NLOPT_PREFIX)/lib -lnlopt)
FFTW_INC  := $(if $(FFTW_PREFIX),-I$(FFTW_PREFIX)/include)
FFTW_LIB  := $(if $(FFTW_PREFIX),-L$(FFTW_PREFIX)/lib -lfftw3 -lfftw3_omp)
# voro++ (header included as <voro++/voro++.hh>). Like NLOPT/FFTW, these are set
# when the prefix is found; a consumer references them only under USE_VORO=1.
VORO_INC  := $(if $(VORO_PREFIX),-I$(VORO_PREFIX)/include)
VORO_LIB  := $(if $(VORO_PREFIX),-L$(VORO_PREFIX)/lib -lvoro++)
# rpath so the produced binary finds the shared libs without LD_LIBRARY_PATH.
# Built with ifneq (not $(if ...)) because the commas in -Wl,-rpath,DIR would
# otherwise be parsed as extra arguments to the make $(if) function.
RPATH :=
ifneq ($(GSL_PREFIX),)
RPATH += -Wl,-rpath,$(GSL_PREFIX)/lib
endif
ifneq ($(NLOPT_PREFIX),)
RPATH += -Wl,-rpath,$(NLOPT_PREFIX)/lib
endif

#----- Optional features (override on the command line to enable) --------------
USE_OPENMP ?= 1
USE_CGAL   ?= 0
USE_VORO   ?= 0
USE_TJ     ?= 0

OMP_FLAG := $(if $(filter 1,$(USE_OPENMP)),-fopenmp)
PREPRO   := $(if $(filter 1,$(USE_CGAL)),-DUSE_CGAL) \
            $(if $(filter 1,$(USE_VORO)),-DUSE_VORO) \
            $(if $(filter 1,$(USE_TJ)),-DUSE_TJ)

#----- Compiler ---------------------------------------------------------------
# make pre-defines CC=cc, so use := (override on the command line still wins).
CC    := g++
STD   ?= -std=c++17
DEBUG ?= -g
OPT   ?= -O3

# Base flags shared by every module.  A module adds its own -I$(core) etc.
COMMON_INC := $(GSL_INC)
BASE_LIBS  := -lm $(GSL_LIB)

CFLAGS ?= -Wall $(STD) $(DEBUG) $(OPT) $(OMP_FLAG) $(PREPRO) $(COMMON_INC) -c
LFLAGS ?= -Wall $(STD) $(DEBUG) $(OPT) $(OMP_FLAG) $(BASE_LIBS) $(RPATH)
