MODULE INPUT_VARS
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str,set_store_size,to_lower
  USE VERSION
  implicit none

  !input variables
  !=========================================================
  character(len=1)             :: DMRGtype
  !Set the DMRG algorithm: 'i,I'=Infinite, 'f,F'=Finite. Default i
  logical                      :: PBCdmrg
  !Set the Periodic(T) or Open(F) boundary conditions. Default F (OBC)
  integer                      :: Ldmrg,iNlat,fNlat
  !# of iDMRG steps to take, Ldmrg=max length of the SB.
  integer                      :: Mdmrg
  !# of states to retain at truncation. If 0 use Edmrg as threshold.  
  real(8)                      :: Edmrg
  !Threshold energy used to evaluate the number of states to keep.
  !If 0d0 use fixed Mdmrg.
  real(8)                      :: QNtruncation_error
  !Threshold error for the Quantum Number truncation
  integer                      :: QNtruncation_dim
  !Threshold dimension for the Quantum Number truncation   
  integer                      :: QNdim
  !Number of conserved Quantum Numbers to consider:
  character(len=12)            :: QNtype
  !Type of conserved Quantum Numbers, local: q in [0,1], global: Q in [-Ldmrg,Ldmrg] 
  real(8),allocatable          :: Dmrg_QN(:)
  !Desired Target Quantum Numbers: size(DMRG_QN)=QNdim
  character(len=7)             :: Dmrg_mode            !
  !Flag to set the DMRG mode: normal[,superc,nonsu2]
  integer                      :: Nsweep
  !# of DMRG sweep to take for finite DMRG algorithm.
  real(8),allocatable          :: Esweep(:)
  !list of threshold energies at each sweep in the finite DMRG algorithm.
  integer,allocatable          :: Msweep(:)
  !list of states to optimize at each sweep in the finite DMRG algorithm.
  integer                      :: Norb
  !# of orbitals
  integer                      :: Nspin=2
  !Nspin=# spins
  ! integer                    :: SUn
  ! !Spin representation 
  real(8),allocatable          :: Uloc(:)
  !local interactions
  real(8)                      :: Ust
  !intra-orbitals interactions
  real(8)                      :: Jh
  !J_Hund: Hunds' coupling constant 
  real(8)                      :: Jx
  !J_X: coupling constant for the spin-eXchange interaction term
  real(8)                      :: Jp
  !J_P: coupling constant for the Pair-hopping interaction term
  real(8)                      :: xmu=0d0
  !chemical potential
  real(8)                      :: temp
  !temperature
  real(8)                      :: deg_evals_threshold
  !threshold for degenerate eigenvalues of rho
  real(8)                      :: eps
  !broadening
  real(8)                      :: wini,wfin
  !frequency range
  logical                      :: HFmode
  !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)                      :: cutoff
  !cutoff for spectral summation
  real(8)                      :: gs_threshold
  !Energy threshold for ground state degeneracy loop up
  logical,allocatable          :: gf_flag(:)
  !evaluate Green's functions for Norb+1
  logical,allocatable          :: chispin_flag(:)
  !evaluate spin susceptibility for Norb+1
  integer                      :: filling
  !Total number of allowed electrons
  logical                      :: sparse_H
  !flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE,
  !or direct on-the-fly H*v product (mem++, cpu--
  integer                      :: verbose
  !Flag to control verbosity of the library
  ! character(len=12)            :: lanc_method
  ! !select the lanczos method to be used in the determination of the spectrum.
  ! !ARPACK (default), LANCZOS (T=0 only) 
  real(8)                      :: lanc_tolerance
  !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos. 
  integer                      :: lanc_niter
  !Max number of Lanczos iterations
  integer                      :: lanc_ngfiter
  !Max number of iteration in resolvant tri-diagonalization
  integer                      :: lanc_ncv_factor
  !Set the size of the block used by Lanczos-Arpack as:
  !(Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer                      :: lanc_ncv_add
  !Adds up to the size of the block to prevent it to become too small
  !(Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer                      :: lanc_neigen
  !Number of required eigenvalues per sector in the SuperBlock
  integer                      :: lanc_dim_threshold
  !Min dimension threshold to use Lanczos determination of the
  !spectrum rather than Lapack based exact diagonalization.
  logical                      :: lanc_v0_dble
  !Set the initial vector used in Arpack to be identical to 1d0/N. This is an experimental feature
  !used to enforce: i) the same starting point, ii) Hermiticity breaking in CMPLX code
  logical                      :: save_umat
  !flag to save Rotation matrices, default=T
  logical                      :: save_block
  !flag to trigger block storage to file
  logical                      :: save_all_blocks
  !flag to save all blocks or just the last one (true if Build_type=DEBUG)
  character(len=1)             :: store_umat
  !index flag to save (and retrieve) Rotation matrices
  character(len=:),allocatable :: umat_file
  !Name prefix of the stored Umatrix file, used to restart DMRG.
  character(len=:),allocatable :: block_file
  !Name prefix of the stored block file, used to restart DMRG.

  !
  !Some parameters for function dimension:
  !=========================================================
  integer                      :: Lmats
  integer                      :: Lreal
  integer                      :: Ltau

  !LOG
  !=========================================================
  integer,save         :: LOGfile


  !THIS IS JUST A RELOCATED GLOBAL VARIABLE
  character(len=200)   :: input_file=""


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine read_input(INPUTunit)
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    character(len=256) :: block_file_
    character(len=256) :: umat_file_
    logical          :: master=.true.
    integer          :: i,rank=0,add,dim
#ifdef _MPI
    if(check_MPI())then
       master=get_Master_MPI(MPI_COMM_WORLD)
       rank  =get_Rank_MPI(MPI_COMM_WORLD)
    endif
#endif
    !
    !Store the name of the input file:
    input_file=str(INPUTunit)
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(DMRGtype,"DMRGtype",INPUTunit,&
         default='i',&
         comment="Set the DMRG algorithm: 'i,I'=Infinite, 'f,F'=Finite. Default i")
    select case(to_lower(DMRGtype))
    case('i','f');continue
    case default;stop "ERROR DMRG: unsupported DMRGtype: DMRGtype !=['i','f']"
    end select

    call parse_input_variable(PBCdmrg,"PBCdmrg",INPUTunit,&
         default=.false.,&
         comment="Set the Periodic(T) or Open(F) boundary conditions. Default F (OBC)")

    call parse_input_variable(Ldmrg,"Ldmrg",INPUTunit,&
         default=5,&
         comment="iDMRG steps to take=max length of the SB.")
    if(PBCdmrg)then
       if(mod(Ldmrg,2)==0)then
          Ldmrg=Ldmrg+1 !Ldmrg = 2*k+1 for some k=1,...
          write(LOGfile,*)"PBC + Ldmrg%2==0 => Shift Ldmrg to Odd number:",Ldmrg
       endif
    endif
    iNlat=2*Ldmrg+2
    fNlat=2*Ldmrg


    call parse_input_variable(Mdmrg,"Mdmrg",INPUTunit,&
         default=0,&
         comment="Number of states for truncation. If 0 use Edmrg as threshold.")

    call parse_input_variable(Edmrg,"Edmrg",INPUTunit,&
         default=1d-8,&
         comment="Threshold energy for truncation. If 0d0 use fixed Mdmrg.")

    call parse_input_variable(Nsweep,"Nsweep",INPUTunit,default=1,&
         comment="Number of DMRG sweep to take for finite DMRG algorithm (min 1).")

    allocate(Msweep(Nsweep),Esweep(Nsweep))
    call parse_input_variable(Msweep,"Msweep",INPUTunit,&
         default=(/(Mdmrg,i=1,Nsweep )/),&
         comment="!list of states for each sweep in a finite DMRG algorithm.")
    call parse_input_variable(Esweep,"Esweep",INPUTunit,&
         default=(/(Edmrg,i=1,Nsweep )/),&
         comment="!list of error threshold for each sweep in a finite DMRG algorithm.")

    call parse_input_variable(dmrg_mode,"DMRG_MODE",INPUTunit,&
         default='normal',&
         comment="Flag to set DMRG type: normal=normal, superc=superconductive, nonsu2=broken SU(2)")

     call parse_input_variable(QNtruncation_error,"QNtruncation_error",INPUTunit,&
         default=0d0,&
         comment="Threshold error on the the symmetry sector RDM truncation: W_k=sum_{i\in QN(k) w_i < QNtruncation_error are excluded")

     call parse_input_variable(QNtruncation_dim,"QNtruncation_dim",INPUTunit,&
         default=1,&
         comment="Threshold dimension on the symmetry sector RDM truncation: dim(QN(k)<QNtruncation_dim are excluded.")

    call parse_input_variable(QNdim,"QNdim",INPUTunit,&
         default=1,&
         comment="Total  conserved abelian quantum numbers to consider.")
    allocate(Dmrg_QN(QNdim))    
    call parse_input_variable(QNtype,"QNtype",INPUTunit,&
         default="local",&
         comment="Type of conserved Quantum Numbers, local: q in [0,1], global: Q in [-Ldmrg,Ldmrg]")
    call parse_input_variable(DMRG_QN,"DMRG_QN",INPUTunit,&
         default=(/(0d0,i=1,QNdim )/),&
         comment="Target Sector QN in units specified by QNtype")

    call parse_input_variable(Norb,"NORB",INPUTunit,&
         default=1,&
         comment="Number of impurity orbitals.")

    !>Interaction parameters:
    allocate(Uloc(Norb))
    call parse_input_variable(uloc,"ULOC",INPUTunit,&
         default=(/( 0d0,i=1,size(Uloc) )/),&
         comment="Values of the local interaction per orbital")

    call parse_input_variable(ust,"UST",INPUTunit,&
         default=0d0,&
         comment="Value of the inter-orbital interaction term")

    call parse_input_variable(Jh,"JH",INPUTunit,&
         default=0d0,&
         comment="Hunds coupling")

    call parse_input_variable(Jx,"JX",INPUTunit,&
         default=0d0,&
         comment="S-E coupling, Jxy Heisenberg")

    call parse_input_variable(Jp,"JP",INPUTunit,&
         default=0d0,&
         comment="P-H coupling, Jz  Heisenberg")

    call parse_input_variable(hfmode,"HFMODE",INPUTunit,&
         default=.true.,&
         comment="Flag to set the Hartree form of the interaction (n-1/2).")

    call parse_input_variable(deg_evals_threshold,"DEG_EVALS_THRESHOLD",INPUTunit,&
         default=1.d-1,&
         comment="threshold for degenerate eigenvalues of rho")
    !
    call parse_input_variable(verbose,"VERBOSE",INPUTunit,&
         default=6,&
         comment="Verbosity level: 0=almost nothing --> 5:all. Really: all")
    !
    call parse_input_variable(sparse_H,"SPARSE_H",INPUTunit,&
         default=.false.,&
         comment="Sparse H*v: True = allocate sparse; False=direct product using QN decomposition")
    !
    !> Lanczos parameters:
    ! call parse_input_variable(lanc_method,"LANC_METHOD",INPUTunit,default="arpack",&
    !      comment="select the lanczos method: ARPACK (default), LANCZOS (T=0 only)")
    call parse_input_variable(lanc_neigen,"LANC_NEIGEN",INPUTunit,default=2,&
         comment="Number of states per SB sector to be determined.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=10,&
         comment="Set the size of the Arpack block (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=0,&
         comment="Set the size of the Arpack block (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,&
         comment="Number of Lanczos iteration in spectrum determination.")
    ! call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,&
    !      comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=1d-12,&
         comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,&
         default=1024,comment="Dimension threshold for Lapack use.")
    call parse_input_variable(lanc_v0_dble,"LANC_V0_dble",INPUTunit,default=.false.,&
         comment="Set the initial vector in Arpack to 1d0/N. Experimental feature ensuring i) same starting point, ii) fix Hermiticity breaking in CMPLX code")
    !
    !>Save Blocks tuning:
    call parse_input_variable(save_umat,"SAVE_UMAT",INPUTunit,default=.true.,&
         comment="Logical flag to save Rotation matrices, default=T.")
    call parse_input_variable(save_block,"SAVE_BLOCK",INPUTunit,default=.true.,&
         comment="Logical flag to trigger block storage. DEBUG enforces T ")
    !
    call parse_input_variable(save_all_blocks,"SAVE_ALL_BLOCKS",INPUTunit,default=.false.,&
         comment="Logical flag to save all blocks (T) or just the last (F, default). DEBUG enforces T")
    !
    call parse_input_variable(store_umat,"STORE_UMAT",INPUTunit,default='b',&
         comment="Chr flag to save (and retrieve) Rotation matrices: b=Block, f=File")
    !
    call parse_input_variable(umat_file_,"UMAT_FILE",INPUTunit,default='umat',&
         comment="Name suffix of the stored rotation matrices, used to measure in DMRG.")
    umat_file=str(umat_file_)
    call parse_input_variable(block_file_,"BLOCK_FILE",INPUTunit,default='block',&
         comment="Name prefix of the stored block file, used to restart DMRG.")
    block_file=str(block_file_)
    !
    !
    !    
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")
    !
#ifdef _DEBUG
    save_block      =.true.
    save_all_blocks =.true.
#endif
    !
    store_umat=to_lower(store_umat)
    select case(store_umat)
    case('b','f');continue
    case default;stop "Input error: store_umat!=[b,f]"
    end select
    !
    call set_store_size(1)
    !
    if(master)then
       call print_input()
       call save_input(INPUTunit)
       call scifor_version()
       call code_version(git_code_version)
    endif
    !
  end subroutine read_input



  ! call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=2,comment="Number of spin degeneracy")
  ! call parse_input_variable(filling,"FILLING",INPUTunit,default=0,comment="Total number of allowed electrons")

  !> Local parameters:
  ! call parse_input_variable(temp,"TEMP",INPUTunit,default=0d0,comment="temperature")
  ! call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,&
  !      comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")


  ! call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
  ! call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
  ! call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")

  !> Aux parameters (GF)
  ! call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
  ! call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
  ! call parse_input_variable(Lmats,"LMATS",INPUTunit,default=4096,comment="Number of Matsubara frequencies.")
  ! call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
  ! call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1024,comment="Number of imaginary time points.")
  !
  ! allocate(gf_flag(Norb))
  ! allocate(chispin_flag(Norb))
  ! call parse_input_variable(gf_flag,"GF_FLAG",INPUTunit,&default=(/( .false.,i=1,size(gf_flag) )/),comment="Flag to activate Greens functions calculation")
  ! call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=(/( .false.,i=1,size(chispin_flag) )/),comment="Flag to activate spin susceptibility calculation.")
  ! !

  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE INPUT_VARS
