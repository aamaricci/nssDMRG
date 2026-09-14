FC=mpif90

SF_INC:=$(shell pkg-config --cflags scifor)
SF_LIB:=$(shell pkg-config --libs scifor)

TYPE=DBLE

FLAG= -cpp -D_DEBUG -D_$(TYPE) -O0 -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace -fbounds-check  -ffree-line-length-none -fPIC -w -fallow-argument-mismatch -fopenmp


define VER
MODULE VERSION
  implicit none
  character(len=41),parameter :: git_code_version = "$(shell git rev-parse HEAD)"	
END MODULE VERSION
endef
export VER


SRC=..
DOTS=../BLOCKS
FOCK=../FOCK
LIST=../LIST
MATRIX=../MATRIX


test: sites blocks


sites: MYOBJ=VERSION.o INPUT_VARS.o AUX_FUNCS.o MATRIX_SPARSE.o QN_TUPLE_BASIS.o LIST_OPERATORS.o LIST_SECTORS.o HLOCAL.o 
sites: VERSION.o INPUT_VARS.o AUX_FUNCS.o MATRIX_SPARSE.o QN_TUPLE_BASIS.o LIST_OPERATORS.o LIST_SECTORS.o HLOCAL.o 
	${FC} $(FLAG) -D_TEST ${MYOBJ} SITES.f90 -o $(HOME)/.bin/testSITES ${SF_INC} ${SF_LIB}

blocks: MYOBJ=VERSION.o INPUT_VARS.o AUX_FUNCS.o HLOCAL.o MATRIX_SPARSE.o MATRIX_BLOCKS.o QN_TUPLE_BASIS.o LIST_SECTORS.o LIST_OPERATORS.o LIST_OMATRICES.o SITES.o
blocks: VERSION.o INPUT_VARS.o AUX_FUNCS.o HLOCAL.o MATRIX_SPARSE.o MATRIX_BLOCKS.o QN_TUPLE_BASIS.o LIST_SECTORS.o LIST_OPERATORS.o LIST_OMATRICES.o SITES.o
	${FC} $(FLAG) -D_TEST ${MYOBJ} BLOCKS.f90 -o $(HOME)/.bin/testBLOCKS ${SF_INC} ${SF_LIB}

define DMRGCONF
 IMEASURE=F                                    !Bool to perform measurements. T for post-processing.
 TS=-1.000000000,-1.000000000                  !Hopping amplitudes
 MH=0.d0                                       !Crystal field splittings
 LAMBDA=0.d0                                   !off-diagonal amplitude
 DMRGTYPE=i                                    !Set the DMRG algorithm: 'i,I'=Infinite, 'f,F'=Finite. Default i
 PBCDMRG=F                                     !Set the Periodic(T) or Open(F) boundary conditions. Default F (OBC)
 LDMRG=20                                      !iDMRG steps to take=max length of the SB.
 MDMRG=3000                                     !Number of states for truncation. If 0 use Edmrg as threshold.
 EDMRG=0.d0                                    !Threshold energy for truncation. If 0d0 use fixed Mdmrg.
 NSWEEP=1                                      !Number of DMRG sweep to take for finite DMRG algorithm (min 1).
 MSWEEP=1                                      !!list of states for each sweep in a finite DMRG algorithm.
 ESWEEP=0.d0                                   !!list of error threshold for each sweep in a finite DMRG algorithm.
 DMRG_MODE=normal                              !Flag to set DMRG type: normal=normal, superc=superconductive, nonsu2=broken SU(2)
 QNDIM=2                                       !Total  conserved abelian quantum numbers to consider.
 QNTYPE=local                                  !Type of conserved Quantum Numbers, local: q in [0,1], global: Q in [-Ldmrg,Ldmrg]
 DMRG_QN=5.000000000E-01,5.000000000E-01       !Target Sector QN in units specified by QNtype
 NORB=1                                        !Number of impurity orbitals.
 ULOC=2.000000000                              !Values of the local interaction per orbital
 UST=0.d0                                      !Value of the inter-orbital interaction term
 JH=0.d0                                       !Hunds coupling
 JX=0.d0                                       !S-E coupling, Jxy Heisenberg
 JP=0.d0                                       !P-H coupling, Jz  Heisenberg
 HFMODE=F                                      !Flag to set the Hartree form of the interaction (n-1/2).
 DEG_EVALS_THRESHOLD=1.000000000E-01           !threshold for degenerate eigenvalues of rho
 VERBOSE=6                                     !Verbosity level: 0=almost nothing --> 5:all. Really: all
 SPARSE_H=F                                    !Sparse H*v: True = allocate sparse; False=direct product using QN decomposition
 LANC_METHOD=arpack                            !select the lanczos method: ARPACK (default), LANCZOS (T=0 only)
 LANC_NEIGEN=1                                 !Number of states per SB sector to be determined.
 LANC_NCV_FACTOR=5                             !Set the size of the Arpack block (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
 LANC_NCV_ADD=0                                !Set the size of the Arpack block (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
 LANC_NITER=512                                !Number of Lanczos iteration in spectrum determination.
 LANC_TOLERANCE=1.000000000E-12                !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.
 LANC_DIM_THRESHOLD=1                          !Dimension threshold for Lapack use.
 SAVE_UMAT=T                                   !Logical flag to save Rotation matrices, default=T.
 SAVE_BLOCK=T                                  !Logical flag to trigger block storage. DEBUG enforces T
 SAVE_ALL_BLOCKS=F                             !Logical flag to save all blocks (T) or just the last (F, default). DEBUG enforces T
 STORE_UMAT=b                                  !Chr flag to save (and retrieve) Rotation matrices: b=Block, f=File
 UMAT_FILE=umat                                !Name suffix of the stored rotation matrices, used to measure in DMRG.
 BLOCK_FILE=block                              !Name prefix of the stored block file, used to restart DMRG.
 LOGFILE=6                                     !LOG unit.
endef

$(file > DMRG.conf,$(DMRGCONF))

VERSION.o:
	@echo "$$VER" > $(SRC)/VERSION.f90
	$(FC) $(FLAG) $(SF_INC) -c $(SRC)/VERSION.f90

INPUT_VARS.o: 
	$(FC) $(FLAG) $(SF_INC) -c $(SRC)/INPUT_VARS.f90

AUX_FUNCS.o: 
	$(FC) $(FLAG) $(SF_INC) -c $(SRC)/AUX_FUNCS.f90

HLOCAL.o: 
	$(FC) $(FLAG) $(SF_INC) -c $(FOCK)/HLOCAL.f90

MATRIX_SPARSE.o: 
	$(FC) $(FLAG) $(SF_INC) -c $(MATRIX)/MATRIX_SPARSE.f90

MATRIX_BLOCKS.o: 
	$(FC) $(FLAG) $(SF_INC) -c $(MATRIX)/MATRIX_BLOCKS.f90

QN_TUPLE_BASIS.o: 
	$(FC) $(FLAG) $(SF_INC) -c $(LIST)/QN_TUPLE_BASIS.f90

LIST_SECTORS.o: 
	$(FC) $(FLAG) $(SF_INC) -c $(LIST)/LIST_SECTORS.f90

LIST_OPERATORS.o:
	$(FC) $(FLAG) $(SF_INC) -c $(LIST)/LIST_OPERATORS.f90

LIST_OMATRICES.o:
	$(FC) $(FLAG) $(SF_INC) -c $(LIST)/LIST_OMATRICES.f90

SITES.o:
	$(FC) $(FLAG) $(SF_INC) -c $(DOTS)/SITES.f90

clean: 
	@echo "Cleaning:"
	@rm -fv *.mod *.o *~ 
	@echo ""


.f90.o:	
	$(FC) $(FLAG) $(SF_INC) -c $< 
