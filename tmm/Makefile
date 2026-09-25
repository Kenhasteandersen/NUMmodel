# Build the NUM model coupled to the TMM (PETSc) driver.
#
#   source <venv>/bin/activate          # venv holding petsc + tmmlib
#   make -C tmm
#
# PETSC_DIR and TMM_DIR are discovered from the active venv if not set.

PETSC_DIR ?= $(shell python -c "from petsc import get_petsc_dir; print(get_petsc_dir())")
TMM_DIR   ?= $(shell python -c "from tmmlib import get_tmm_dir; print(get_tmm_dir())")

# NUM's f_* entry points live in the matlab wrapper library, not the core one.
NUM_BUILD ?= $(CURDIR)/../build/Fortran
NUM_LIB   := NUMmodel_matlab

GFORTRAN_DIR := $(dir $(shell gfortran -print-file-name=libgfortran.dylib))

CC      := mpicc
FC      := mpif90
INCS    := -I$(PETSC_DIR)/include -I$(TMM_DIR)/include -I$(TMM_DIR)/include/tmm
CFLAGS  := -O2 $(INCS)

# The driver, plus the stubs we do not override. external_forcing_num.c takes
# the place of the tmm_external_forcing.c stub.
TMM_SRCS  := tmm_main.c tmm_external_bc.c tmm_monitor.c tmm_misfit.c
TMM_OBJS  := $(TMM_SRCS:%.c=%.o)
OBJS      := $(TMM_OBJS) external_forcing_num.o insolation.o

LDFLAGS := -L$(TMM_DIR)/lib -ltmm -Wl,-rpath,$(TMM_DIR)/lib \
           -L$(PETSC_DIR)/lib -lpetsc -Wl,-rpath,$(PETSC_DIR)/lib \
           -L$(NUM_BUILD) -l$(NUM_LIB) -Wl,-rpath,$(NUM_BUILD) \
           -L$(GFORTRAN_DIR) -lgfortran -Wl,-rpath,$(GFORTRAN_DIR)

numtmm: $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

%.o: $(TMM_DIR)/src/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

external_forcing_num.o: external_forcing_num.c
	$(CC) -c $< -o $@ $(CFLAGS)

# Shared with the MOPS example; TMM does not ship it in libtmm.
insolation.o: $(TMM_SRC_REPO)/models/mops3.0/code/src/insolation.F
	$(FC) -c $< -o $@

clean:
	rm -f *.o numtmm

.PHONY: clean
