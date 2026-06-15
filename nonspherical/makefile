#==============================================================================
# modules/nonspherical/makefile
#
# Nonspherical-particle shapes and dispersions (packings of nonspherical
# inclusions). Library module (no main).
#   ParticleShapes : shape descriptors (spheres, ellipses, polygons, ...).
#   Dispersions    : SpherePacking-like container of nonspherical particles + I/O.
# Depends on: core.
#
#   make / make objs   # compile module objects
#   make clean
#==============================================================================
include ../config.mk

ODIR := obj/

SRCS := ParticleShapes.cpp FormFactors.cpp Dispersions.cpp
OBJS := $(patsubst %.cpp,$(ODIR)%.o,$(SRCS))

INC := -I$(core) -I$(nonspherical)

.PHONY: all objs clean
all: objs
objs: $(OBJS)

$(ODIR)%.o: %.cpp | $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

clean:
	rm -f $(ODIR)*.o
