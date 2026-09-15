FC=mpif90
PLAT=gnu

LIB_DMRG=dmrg

GLOB_INC+=$(shell pkg-config --cflags ${LIB_DMRG})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_DMRG})

GLOB_INC+=$(shell pkg-config --cflags scifor)# libtcmalloc libprofiler)
GLOB_LIB+=$(shell pkg-config --libs   scifor)# libtcmalloc libprofiler)


ifeq ($(PLAT),intel)
FFLAG=-O1 -ftz
DFLAG=-p -O0 -g
AFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
FPPSERIAL =-fpp -D_
FPPMPI =-fpp -D_	
endif

ifeq ($(PLAT),gnu)
#FFLAG = -O1 -ffast-math -march=native -funroll-loops -ffree-line-length-none -g
FFLAG = -O3 -funroll-loops -ffree-line-length-none 
DFLAG = -w -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
AFLAG = -w -O0 -p -g  -fbacktrace -fwhole-file -fcheck=all -fbounds-check -fsanitize=address -fdebug-aux-vars -Wall -Waliasing -Wsurprising -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wdo-subscript -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wno-align-commons -Wno-overwrite-recursive -Wno-tabs -Wreal-q-constant -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime
FPPSERIAL =-cpp -D_
FPPMPI =-cpp -D_MPI
endif


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

define colorecho	
	@tput setaf $2
	@tput bold
	@echo $1
	@tput sgr0
endef

all: spin hubbard


spin: FLAG:=${FFLAG} ${FPPMPI}
spin: 
	@echo ""
	$(call colorecho,"compiling dmrg_spin_1d.f90 ", 1)
	$(FC) -o Spin1d_SU2/dmrg_spin_1d Spin1d_SU2/dmrg_spin_1d.f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created MPI dmrg_spin_1d in Spin1d_SU2/", 2)
	@echo ""

hubbard: FLAG:=${FFLAG} ${FPPMPI}
hubbard: 
	@echo ""	
	$(call colorecho,"compiling dmrg_hm_1d.f90 ", 1)
	$(FC) -o Hubbard1d/dmrg_hm_1d Hubbard1d/dmrg_hm_1d.f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created MPI dmrg_hm_1d in  Hubbard1d/", 2)
	@echo ""	

mpi: FLAG:=${FFLAG} ${FPPMPI}
mpi: spin hubbard

serial: FLAG:=${FFLAG} ${FPPSERIAL}
serial: spin hubbard


debug: FLAG:=${DFLAG} ${FPPMPI} -D_DEBUG
debug: spin hubbard

serial_debug: FLAG:=${DFLAG} ${FPPSERIAL} -D_DEBUG
serial_debug: spin hubbard



clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~



#########################################################################
