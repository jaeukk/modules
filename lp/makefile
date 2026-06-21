#==============================================================================
# modules/lp/makefile   (linear-programming jamming, leaf module)
#
# Migrated from the top-level legacy tree:
#   TorquatoJiao.{cpp,h}   the Torquato-Jiao adaptive-shrinking-cell LP jamming
#                          algorithm for hard-sphere packings (TJJammer)
#   gurobi_c.h, gurobi_c++.h   vendored Gurobi C/C++ SDK headers (third-party,
#                          copied verbatim from the legacy tree; needed only
#                          under USE_TJ)
#
# Depends on: core (etc/PeriodicCellList), GSL.  External (USE_TJ only): glpk +
# Gurobi + Boost.Filesystem.
#
# The ENTIRE TorquatoJiao.cpp is wrapped in `#ifdef USE_TJ`, and TorquatoJiao.h
# provides a stub TJJammer when USE_TJ is off — so by default (USE_TJ=0, which is
# the config.mk default) this compiles to an essentially empty object.  glpk and
# Gurobi are not installed on this box, so a USE_TJ=1 build is currently blocked.
# Enable (once those libs exist) with:  make USE_TJ=1
#
#   make / make objs   # compile the module object
#   make clean
#==============================================================================
include ../config.mk

ODIR := obj/

SRCS := TorquatoJiao.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

# -I. picks up the vendored gurobi_c*.h siblings under USE_TJ.
INC := -I$(core) -I.

.PHONY: all objs clean
all: objs
objs: $(OBJS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
