# Include after apps/CCOptimization/makefile; keep generated artifacts outside the repo.
HISTORY_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
HISTORY_OUT ?= /tmp/cco-history-contact100-reproduction
$(HISTORY_OUT):
	mkdir -p $@
$(HISTORY_OUT)/history.o: $(HISTORY_DIR)history.cpp $(HISTORY_DIR)fixtures.h | $(HISTORY_OUT)
	$(CC) $(CFLAGS) $(INC) -I$(HISTORY_DIR) -o $@ $<
$(HISTORY_OUT)/history.out: $(HISTORY_OUT)/history.o $(ODIR)MultiStealthyPotential.o deps
	$(CC) -o $@ $(filter %.o,$^) $(CORE_OBJS) $(PAIR_OBJS) $(POT_OBJS) $(LFLAGS) $(NLOPT_LIB) $(LIBS_FINUFFT) -Wl,--wrap=nlopt_optimize
