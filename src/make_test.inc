FC=mpif90

TYPE=DBLE
FLAG=-cpp -D_$(TYPE) -O0 -g -fcheck=all -fbacktrace -fbounds-check \
     -ffree-line-length-none

.PHONY: test geometry clean

test: geometry

geometry:
	$(FC) $(FLAG) -D_TEST DMRG_GEOMETRY.f90 -o $(HOME)/.bin/testDMRG_GEOMETRY

clean:
	rm -f *.mod *.o *~
