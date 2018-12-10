#
# Makefile for DSM_f90 1.0
#
# Compiler: gfortran
# 
# Author: Zhichao Li @ AAQRL, Washington University in St. Louis
#

# main program name
MPROG = DSM

# object files *.o
ODIR = ./obj

# compiler flags required
FLAGS = -Wall -g

# project dependent
_OBJS = preci.mod \
        globals.mod \
        calbeta_common.mod \
        aux_data_funcs.mod \
        funcs.mod \
        quad_glegq.mod \
        dtwodq_mod.mod \
        calbeta_mod.mod \
        ot_mod.mod \
        xdoteq_mod.mod \
        rk45_init0.mod \
        discsecpb.mod

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = gfortran

# compiling
$(ODIR)/%.mod : %.f90
	$(CC) $(FLAGS) -c -o $@ $<

# linking
$(MPROG) : $(OBJS)
	$(CC) -o $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.mod ./*.mod *~ core $(MPROG)
