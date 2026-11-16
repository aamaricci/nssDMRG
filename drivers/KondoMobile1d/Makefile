FC=mpif90

DMRG_INC:=$(shell pkg-config --cflags dmrg)
DMRG_LIB:=$(shell pkg-config --libs dmrg)

FLAG=-cpp -D_MPI -D_DBLE -O2 -ffree-line-length-none
TARGET=$(HOME)/.bin/dmrg_kondo_mobile_1d

.PHONY: all clean

all: $(TARGET)

$(TARGET): dmrg_kondo_mobile_1d.f90
	$(FC) $(FLAG) $< -o $@ $(DMRG_INC) $(DMRG_LIB)

clean:
	rm -f *.mod *.o
