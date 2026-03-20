MODULE DMRG_GLOBAL
  USE SCIFOR
  USE INPUT_VARS
  USE AUX_FUNCS
  USE MATRIX_GRAPH
  USE MATRIX_SPARSE, id=>sp_eye
  USE MATRIX_BLOCKS
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE SITES
  USE BLOCKS
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  ! abstract interface
  !    function UserHconnect(left,right) 
  !      USE BLOCKS
  !      implicit none
  !      type(block)                          :: left
  !      type(block)                          :: right
  !      real(8),dimension(:,:),allocatable :: UserHconnect ![Nso,Nso]
  !    end function UserHconnect
  ! end interface
  ! procedure(UserHconnect),pointer,public  :: Hmodel=>null()
  ! !

#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable :: HopH
#else
  real(8),dimension(:,:),allocatable    :: HopH
#endif




  abstract interface
     subroutine sparse_HxV(Nloc,v,Hv)
       integer                          :: Nloc
#ifdef _CMPLX
       complex(8),dimension(Nloc)       :: v
       complex(8),dimension(Nloc)       :: Hv
#else
       real(8),dimension(Nloc)          :: v
       real(8),dimension(Nloc)          :: Hv
#endif
     end subroutine sparse_HxV
  end interface
  procedure(sparse_HxV),pointer,public  :: spHtimesV_p=>null()



  type(sparse_matrix)                            :: spHsb
  !
  real(8)                                        :: truncation_error_left,truncation_error_right
  character(len=:),allocatable                   :: suffix
  real(8),dimension(:),allocatable               :: target_Qn,current_target_QN
  integer                                        :: current_L
  type(block)                                    :: init_left,init_right
  logical                                        :: init_called=.false.
	logical                                       :: dims_set=.false.
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable          :: gs_vector
#else
  real(8),dimension(:,:),allocatable             :: gs_vector
#endif
  real(8),dimension(:),allocatable               :: gs_energy
  real(8),dimension(:),allocatable               :: rho_left_evals
  real(8),dimension(:),allocatable               :: rho_right_evals
  type(blocks_matrix)                            :: psi_left
  type(blocks_matrix)                            :: psi_right
  type(blocks_matrix)                            :: rho_left
  type(blocks_matrix)                            :: rho_right
  !GLOBAL LEFT & RIGHT & DOT 
  type(block)                                    :: left,right
  type(site),dimension(:),allocatable            :: dot
  !
  !SUPERBLOCK SHARED THINGS
  integer,dimension(:),allocatable               :: sb_states
  type(sectors_list)                             :: sb_sector
  !
  integer                                        :: Mstates
  real(8)                                        :: Estates
  !
  !Block to Global Map M and its inverse
  !b2g_Map: maps BLOCK indices/labels set during iDMRG construction TO GLOBAL indices 1...N
  !g2b_Map: maps GLOBAL indices TO BLOCK indices/labels as obtained during iDMRG construction
  ! e.g.
  ! PBC:  [7,5,3,1,2,4,6]
  !       [1,2,3,4,5,6,7]
  !
  ! OBC:  [1,2,3,4,5,6,7]
  !       [1,2,3,4,5,6,7]
  !
  ! maps:
  !
  !block index:           1 -> global index           4
  !block index:           2 -> global index           5
  !block index:           3 -> global index           3
  !block index:           4 -> global index           6
  !block index:           5 -> global index           2
  !block index:           6 -> global index           7
  !block index:           7 -> global index           1
  !
  !global index:           1 -> block index           7
  !global index:           2 -> block index           5
  !global index:           3 -> block index           3
  !global index:           4 -> block index           1
  !global index:           5 -> block index           2
  !global index:           6 -> block index           4
  !global index:           7 -> block index           6
  !
  integer,dimension(:),allocatable               :: b2gMap,g2bMap
  !
  !Memory pool for HxV direct product 
  type(sparse_matrix),allocatable,dimension(:)   :: Hleft,Hright
  type(sparse_matrix),allocatable,dimension(:,:) :: A,B
  integer,dimension(:),allocatable               :: Dls,Drs,Offset
  integer,dimension(:,:),allocatable             :: RowOffset,ColOffset,isb2jsb
  integer,dimension(:,:),allocatable             :: IsHconjg


  !Profiling stuff
  character(len=100)                             :: prof_file="dmrg_profile.out"
  real(8),allocatable,dimension(:)               :: prof_times
  character(len=32),allocatable,dimension(:)     :: prof_names
  !
  real(8)                                        :: t0,t1

  real(8)                                        :: t_dmrg_step      !_MAIN
  real(8)                                        :: t_enlarge_blocks !_CONNECT
  real(8)                                        :: t_connect_blocks !_CONNECT
  real(8)                                        :: t_sb_get_states  !_SUPERBLOCK
  real(8)                                        :: t_sb_diag        !_SUPERBLOCK
  real(8)                                        :: t_setup_sb_sparse!_SUPERBLOC_SETUP
  real(8)                                        :: t_setup_sb_direct!_SUPERBLOC_SETUP
  real(8)                                        :: t_hxv_sparse     !_SUPERBLOC_SETUP
  real(8)                                        :: t_hxv_direct     !_SUPERBLOC_SETUP
  real(8)                                        :: t_hxv_1LxHR      !_SUPERBLOC_SETUP
  real(8)                                        :: t_hxv_HLx1R      !_SUPERBLOC_SETUP
  real(8)                                        :: t_hxv_AxB        !_SUPERBLOC_SETUP
  real(8)                                        :: t_hxv_B          !_SUPERBLOC_SETUP
  real(8)                                        :: t_hxv_A          !_SUPERBLOC_SETUP  
  real(8)                                        :: t_rdm_get        !_RDM
  real(8)                                        :: t_rdm_diag       !_RDM
  real(8)                                        :: t_rdm_renorm     !_RDM
  !
  real(8)                                        :: rdcd_sb_dim
  real(8)                                        :: full_sb_dim
  real(8)                                        :: t_a2av           !_GLOBAL
  real(8)                                        :: t_sctr           !_GLOBAL
  real(8)                                        :: t_gthr           !_GLOBAL
  real(8)                                        :: t_agthr          !_GLOBAL
  real(8)                                        :: kb_vecDim
  real(8)                                        :: kb_sent_a2av
  real(8)                                        :: kb_recv_a2av
  real(8)                                        :: kb_sent_sctr
  real(8)                                        :: kb_recv_sctr
  real(8)                                        :: kb_sent_gthr
  real(8)                                        :: kb_recv_gthr
  real(8)                                        :: kb_sent_agthr
  real(8)                                        :: kb_recv_agthr
  real(8)                                        :: kb_agthr_sb_states
  real(8)                                        :: kb_sb_setup_bcast
  !
  real(8)                                        :: D_kb_size
  real(8)                                        :: C_kb_size
  real(8)                                        :: I_kb_size  
  real(8)                                        :: DATA_kb_size
  !
  integer                                        :: NumOp





  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                        :: MpiComm_Global=MPI_COMM_NULL
  integer                                        :: MpiComm=MPI_COMM_NULL
#else
  integer                                        :: MpiComm_Global=0
  integer                                        :: MpiComm=0

#endif
  logical                                        :: MpiStatus=.false.
  logical                                        :: MpiMaster=.true.
  integer                                        :: MpiRank=0
  integer                                        :: MpiSize=1
  integer,allocatable,dimension(:)               :: mpiDls
  integer,allocatable,dimension(:)               :: mpiDrs
  integer,allocatable,dimension(:)               :: mpiDl,mpiDr
  integer,allocatable,dimension(:)               :: mpiOffset
  integer,allocatable,dimension(:)               :: mpiSBComm ! sub-communicators for each sb sector
  integer,allocatable,dimension(:)               :: mpiNactive ! number of active ranks in each sb sector








#ifdef _MPI
  interface scatter_vector_MPI
     module procedure  :: scatter_vector_MPI_v1
     module procedure  :: scatter_vector_MPI_v2
  end interface scatter_vector_MPI

  interface gather_vector_MPI
     module procedure  :: gather_vector_MPI_v1
     module procedure  :: gather_vector_MPI_v2
  end interface gather_vector_MPI

  public :: scatter_vector_MPI
  public :: gather_vector_MPI
  public :: allgather_vector_MPI
#endif


contains


  



  
  



  !##################################################################
  !##################################################################
  !                   MPI AUX FUNCTIONS
  !##################################################################
  !##################################################################
  subroutine dmrg_set_MpiComm()
#ifdef _MPI
    integer :: ierr
    integer :: typesize
    MpiComm_Global = MPI_COMM_WORLD
    MpiComm        = MPI_COMM_WORLD
    MpiStatus      = .true.
    MpiSize        = get_Size_MPI(MpiComm_Global)
    MpiRank        = get_Rank_MPI(MpiComm_Global)
    MpiMaster      = get_Master_MPI(MpiComm_Global)
    !
    !Returns the byte dimension for a single element of MPI_DOUBLE_PRECISION/COMPLEX type
    call MPI_Type_size(MPI_DOUBLE_COMPLEX, typesize, ierr);C_kb_size =dble(typesize)/1024d0
    call MPI_Type_size(MPI_DOUBLE_PRECISION,typesize,ierr);D_kb_size =dble(typesize)/1024d0
    call MPI_Type_size(MPI_INTEGER, typesize, ierr)       ;I_kb_size =dble(typesize)/1024d0
#ifdef _CMPLX
    DATA_kb_size=C_kb_size
#else
    DATA_kb_size=D_kb_size
#endif
#ifdef _DEBUG
    if(MpiMaster)write(Logfile,*)"DEBUG: dmrg_set_MpiComm - setting MPI comm"
#endif
#endif
  end subroutine dmrg_set_MpiComm

  subroutine dmrg_del_MpiComm()
#ifdef _MPI    
    MpiComm_Global = MPI_UNDEFINED
    MpiComm        = MPI_UNDEFINED
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
#ifdef _DEBUG
    if(MpiMaster)write(Logfile,*)"DEBUG: dmrg_del_MpiComm - deleting MPI comm"
#endif
    MpiMaster      = .true.
#endif
  end subroutine dmrg_del_MpiComm



  !##################################################################
  !##################################################################
  !              MPI AUXILIARY FUNCTIONS
  !##################################################################
  !##################################################################
  !
  !####################################################################
  !               ALL-2-ALL-V VECTOR MPI TRANSPOSITION 
  !####################################################################
#ifdef _MPI
  subroutine vector_transpose_MPI(nrow,qcol,a,ncol,qrow,b,sub_comm)
    !
    integer                             :: nrow !Global number of rows 
    integer                             :: ncol !Global number of columns
    integer                             :: qrow !Local number of rows on each thread
    integer                             :: qcol !Local number of columns on each thread
    integer                             :: sub_comm !Sub-communicator for active ranks
#ifdef _CMPLX    
    complex(8)                          :: a(nrow,qcol) ! Input vector to be transposed
    complex(8)                          :: b(ncol,qrow) ! Output vector :math:`b = v^T`
    complex(8),dimension(:),allocatable :: Vtmp
#else
    real(8)                             :: a(nrow,qcol) ! Input vector to be transposed
    real(8)                             :: b(ncol,qrow) ! Output vector :math:`b = v^T`
    real(8),dimension(:),allocatable    :: Vtmp
#endif
    integer,allocatable,dimension(:,:)  :: send_counts,send_offset
    integer,allocatable,dimension(:,:)  :: recv_counts,recv_offset
    integer                             :: counts,Ntot
    integer                             :: i,j,irank,ierr
    !
    if (sub_comm == MPI_COMM_NULL) return
	  !
		!
		!Round robin distribution of columns: hardcoded for now, can be easily changed to block distribution if needed	
    counts = Nrow/MpiSize
    Ntot   = Ncol/MpiSize
    if(mod(Ncol,MpiSize)/=0)Ntot=Ntot+1
    !
    allocate(send_counts(0:MpiSize-1,Ntot));send_counts=0
    allocate(send_offset(0:MpiSize-1,Ntot));send_offset=0
    allocate(recv_counts(0:MpiSize-1,Ntot));recv_counts=0
    allocate(recv_offset(0:MpiSize-1,Ntot));recv_offset=0
    !
	 !Round robin distribution of rows: hardcoded for now, can be easily changed to block distribution if needed

    do i=1,qcol
       do irank=0,MpiSize-1
          if(irank < mod(Nrow,MpiSize))then
             send_counts(irank,i) = counts+1
          else
             send_counts(irank,i) = counts
          endif
       enddo
    enddo
    !
    do i=1,Ntot
       call MPI_AllToAll(&
            send_counts(0:,i),1,MPI_INTEGER,&
            recv_counts(0:,i),1,MPI_INTEGER,&
            MpiComm,ierr) !not sure if this is the right communicator to use here, but it should be fine since only active ranks have non-zero counts
    enddo
    !
    do i=1,Ntot
       do irank=1,MpiSize-1
          send_offset(irank,i) = send_counts(irank-1,i) + send_offset(irank-1,i)
       enddo
    enddo
    !
    !Get the irank=0 elements, i.e. first entries:
    recv_offset(0,1) = 0
    do i=2,Ntot
       recv_offset(0,i) = sum(recv_counts(0,:i-1))
    enddo
    !the rest of the entries:
    do i=1,Ntot
       do irank=1,MpiSize-1
          recv_offset(irank,i) = recv_offset(irank-1,i) + sum(recv_counts(irank-1,:))
       enddo
    enddo
    !
    !
    t0=t_start()
    do j=1,Ntot
       !Fix issue with empty columns arising from having few MPI nodes
       if(j<=size(A,2))then
          Vtmp = A(:,j)            !automatic allocation
       else
          allocate(Vtmp(0))
       endif

#ifdef _CMPLX
       call MPI_AllToAllV(& ! A(:,j),send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
            Vtmp,send_counts(:,j),send_offset(:,j),MPI_DOUBLE_COMPLEX,&
            B(:,:),recv_counts(:,j),recv_offset(:,j),MPI_DOUBLE_COMPLEX,&
            sub_comm,ierr)
#else
       call MPI_AllToAllV(& ! A(:,j),send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
            Vtmp,send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
            B(:,:),recv_counts(:,j),recv_offset(:,j),MPI_DOUBLE_PRECISION,&
            sub_comm,ierr)
#endif
       deallocate(Vtmp)
    enddo
    t_a2av=t_a2av+t_stop()
    !
    call local_transpose(b,ncol,qrow)
    !
    !
    !Accumulate onto global variable (initialized or reset to 0 somewhere else)
    kb_sent_a2av = kb_sent_a2av + sum(send_counts)*DATA_kb_size
    kb_recv_a2av = kb_recv_a2av + sum(recv_counts)*DATA_kb_size
    !
    return
  end subroutine vector_transpose_MPI
  !
  subroutine local_transpose(mat,nrow,ncol)
    integer                         :: nrow,ncol
#ifdef _CMPLX
    complex(8),dimension(Nrow,Ncol) :: mat
#else
    real(8),dimension(Nrow,Ncol)    :: mat
#endif
    mat = transpose(reshape(mat,[Ncol,Nrow]))
  end subroutine local_transpose
  !=========================================================








  !! Scatter V into the arrays Vloc on each thread: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine scatter_vector_MPI_v1(MpiComm,v,vloc)
    integer                          :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:)          :: v    !size[N]
    complex(8),dimension(:)          :: vloc !size[Nloc]
#else
    real(8),dimension(:)             :: v    !size[N]
    real(8),dimension(:)             :: vloc !size[Nloc]
#endif
    integer                          :: i,k,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: pCounts,pOffset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>4)write(Logfile,*)"DEBUG: d_scatter_vector_MPI: scatter v into vloc"
#endif
    !
    if( MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null )return
    ! stop "scatter_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "scatter_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(pCounts(0:MpiSize-1)) ; pCounts=0
    allocate(pOffset(0:MpiSize-1)) ; pOffset=0
    !
    Vloc=0d0
    do k=1,size(sb_sector)
       !
       !Get Counts;
       pCounts=0
       pOffset=0
       call MPI_AllGather(mpiDl(k),1,MPI_INTEGER,pCounts,1,MPI_INTEGER,MpiComm,MpiIerr)
       kb_recv_sctr = kb_recv_sctr + sum(pCounts)*DATA_kb_size
       !
       !Get Offset:
       pOffset(0)=0
       do i=1,MpiSize-1
          pOffset(i) = pOffset(i-1) + pCounts(i-1)
       enddo
       !
       if(MpiMaster)then
          v_start = 1 + Offset(k)
          v_end = Drs(k)*Dls(k)+OffSet(k)
       else
          v_start = 1
          v_end = 1
       endif
       !
       vloc_start = 1 + mpiOffset(k)
       vloc_end   = mpiDl(k)+mpiOffSet(k)
       !
       t0=t_start()
#ifdef _CMPLX
       call MPI_Scatterv(V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_COMPLEX,&
            Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
#else
       call MPI_Scatterv(V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_PRECISION,&
            Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_PRECISION,0,MpiComm,MpiIerr)
#endif
       t_sctr=t_sctr+t_stop()
    enddo
    !
    ! call Sum_MPI(MpiComm,kb_recv_sctr)
    ! call Max_MPI(MpiCOmm,t_sctr)
    kb_sent_sctr=kb_recv_sctr
    !
    return
  end subroutine scatter_vector_MPI_v1

  subroutine scatter_vector_MPI_v2(MpiComm,v,vloc)
    integer                   :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:,:) :: v    !size[N,Neigen]
    complex(8),dimension(:,:) :: vloc !size[Nloc,Neigen]
#else
    real(8),dimension(:,:)    :: v    !size[N,Neigen]
    real(8),dimension(:,:)    :: vloc !size[Nloc,Neigen]
#endif
    integer                   :: N,Nloc,Neigen,i
    !
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>4)write(Logfile,*)"DEBUG scatter_vector_MPI: scatter many v"
#endif
    N      = size(v,1)
    Nloc   = size(vloc,1)
    Neigen = size(vloc,2)
    if( size(v,2) < Neigen ) stop "scatter_vector_MPI error: size(v,2) < Neigen"
    !
    do i=1,Neigen
       call scatter_vector_MPI_v1(MpiComm,v(:,i),vloc(:,i))
    end do
    !
    return
  end subroutine scatter_vector_MPI_v2



  !! gather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine gather_vector_MPI_v1(MpiComm,vloc,v)
    integer                          :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:)          :: vloc !size[Nloc]
    complex(8),dimension(:)          :: v    !size[N]
#else
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
#endif
    integer                          :: i,k,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: pCounts,pOffset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>4)write(Logfile,*)"DEBUG d_gather_basis_MPI: gather  v"
#endif
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    !stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(pCounts(0:MpiSize-1)) ; pCounts=0
    allocate(pOffset(0:MpiSize-1)) ; pOffset=0
    !
    V = 0d0
    do k=1,size(sb_sector)
       !
       !Get Counts;
       pCounts=0
       pOffset=0
       call MPI_AllGather(mpiDl(k),1,MPI_INTEGER,pCounts,1,MPI_INTEGER,MpiComm,MpiIerr)
       kb_sent_gthr = kb_sent_gthr + sum(pCounts)*DATA_kb_size
       !
       !Get Offset:
       pOffset(0)=0
       do i=1,MpiSize-1
          pOffset(i) = pOffset(i-1) + pCounts(i-1)
       enddo
       !
       vloc_start = 1 + mpiOffset(k)
       vloc_end   = mpiDl(k)+mpiOffSet(k)
       !
       if(MpiMaster)then
          v_start = 1 + Offset(k)
          v_end = Drs(k)*Dls(k)+OffSet(k)
       else
          v_start = 1
          v_end = 1
       endif
       !
       t0=t_start()
#ifdef _CMPLX
       call MPI_Gatherv(Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_COMPLEX,&
            V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
#else
       call MPI_Gatherv(Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_PRECISION,&
            V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_PRECISION,0,MpiComm,MpiIerr)
#endif
       t_gthr=t_gthr+t_stop()
    enddo
    !
    ! call Sum_MPI(MpiComm,kb_sent_gthr)
    ! call Max_MPI(MpiComm,t_gthr)
    kb_recv_gthr=kb_sent_gthr
    return
  end subroutine gather_vector_MPI_v1

  subroutine gather_vector_MPI_v2(MpiComm,vloc,v)
    integer                   :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:,:) :: vloc !size[Nloc,Neigen]
    complex(8),dimension(:,:) :: v    !size[N,Neigen]
#else
    real(8),dimension(:,:)    :: vloc !size[Nloc,Neigen]
    real(8),dimension(:,:)    :: v    !size[N,Neigen]
#endif
    integer                   :: N,Nloc,Neigen,i
    !
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>4)write(Logfile,*)"DEBUG gather_vector_MPI: gather many v"
#endif
    N      = size(v,1)
    Nloc   = size(vloc,1)
    Neigen = size(vloc,2)
    if( size(v,2) < Neigen ) stop "gather_vector_MPI error: size(v,2) < Neigen"
    !
    do i=1,Neigen
       call gather_vector_MPI_v1(MpiComm,vloc(:,i),v(:,i))
    end do
    !
    return
  end subroutine gather_vector_MPI_v2




  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine allgather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
#ifdef _CMPLX
    complex(8),dimension(:)          :: vloc !size[Nloc]
    complex(8),dimension(:)          :: v    !size[N]
#else
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
#endif
    integer                          :: i,k,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: pCounts,pOffset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(MpiMaster.AND.verbose>4)write(Logfile,*)"DEBUG d_allgather_basis_MPI: allgather v"
#endif
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    ! stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N    = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "allgather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(pCounts(0:MpiSize-1)) ; pCounts=0
    allocate(pOffset(0:MpiSize-1)) ; pOffset=0
    !
    V = 0d0
    do k=1,size(sb_sector)
       !
       !Get Counts;
       pCounts=0
       pOffset=0
       call MPI_AllGather(mpiDl(k),1,MPI_INTEGER,pCounts,1,MPI_INTEGER,MpiComm,MpiIerr)
       kb_recv_agthr = kb_recv_agthr + sum(pCounts)*DATA_kb_size
       !
       !Get Offset:
       pOffset(0)=0
       do i=1,MpiSize-1
          pOffset(i) = pOffset(i-1) + pCounts(i-1)
       enddo
       !
       vloc_start = 1 + mpiOffset(k)
       vloc_end   = mpiDl(k)+mpiOffSet(k)
       !
       v_start = 1 + Offset(k)
       v_end = Drs(k)*Dls(k)+OffSet(k)
       !
       t0 = t_start()
#ifdef _CMPLX
       call MPI_AllGatherv(Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_COMPLEX,&
            V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_COMPLEX,MpiComm,MpiIerr)
#else
       call MPI_AllGatherv(Vloc(vloc_start:vloc_end),mpiDl(k),MPI_DOUBLE_PRECISION,&
            V(v_start:v_end),pCounts,pOffset,MPI_DOUBLE_PRECISION,MpiComm,MpiIerr)
#endif
       t_agthr=t_agthr+t_stop()
    enddo
    !
    return
  end subroutine allgather_vector_MPI

#endif









  !##################################################################
  !##################################################################
  !                   PROFILE AUX FUNCTIONS
  !##################################################################
  !##################################################################
  subroutine reset_profile()
    !
    t0                = 0d0
    t1                = 0d0
    t_a2av            = 0d0 !_GLOBAL
    t_sctr            = 0d0 !_GLOBAL
    t_gthr            = 0d0 !_GLOBAL
    t_agthr           = 0d0 !_GLOBAL
    t_dmrg_step       = 0d0 !_MAIN
    t_enlarge_blocks  = 0d0 !_CONNECT
    t_connect_blocks  = 0d0 !_CONNECT
    t_sb_get_states   = 0d0 !_SUPERBLOCK
    t_sb_diag         = 0d0 !_SUPERBLOCK
    t_setup_sb_sparse = 0d0 !_SUPERBLOC_SETUP
    t_setup_sb_direct = 0d0 !_SUPERBLOC_SETUP
    t_hxv_direct      = 0d0 !_SUPERBLOC_SETUP
    t_hxv_1LxHR       = 0d0 !_SUPERBLOC_SETUP
    t_hxv_HLx1R       = 0d0 !_SUPERBLOC_SETUP
    t_hxv_AxB         = 0d0 !_SUPERBLOC_SETUP
    t_hxv_B           = 0d0 !_SUPERBLOC_SETUP
    t_hxv_A           = 0d0 !_SUPERBLOC_SETUP
    t_rdm_get         = 0d0 !_RDM
    t_rdm_diag        = 0d0 !_RDM
    t_rdm_renorm      = 0d0 !_RDM
    kb_vecDim         = 0d0
    kb_sent_a2av      = 0d0
    kb_recv_a2av      = 0d0
    kb_sent_sctr      = 0d0
    kb_recv_sctr      = 0d0
    kb_sent_gthr      = 0d0
    kb_recv_gthr      = 0d0    
    kb_sent_agthr     = 0d0
    kb_recv_agthr     = 0d0
    kb_agthr_sb_states= 0d0
    kb_sb_setup_bcast = 0d0
    !
    if(allocated(prof_times))deallocate(prof_times)
    if(allocated(prof_names))deallocate(prof_names)
  end subroutine reset_profile





  subroutine set_profile
    integer :: unit,i,Ierr,q
#ifdef _PROFILE
#ifdef _MPI
    if(MpiStatus)then
       call Sum_MPI(MpiComm,kb_vecDim)
       kb_vecDim=kb_vecDim/MpiSize
       !
       call Max_MPI(MpiComm,t_enlarge_blocks)
       !
       call Max_MPI(MpiComm,t_sb_get_states)
       call Sum_MPI(MpiComm,kb_agthr_sb_states)
       !
       call Max_MPI(MpiComm,t_sb_diag)
       call Max_MPI(MpiComm,t_sctr)
       call Sum_MPI(MpiComm,kb_sent_sctr)
       call Sum_MPI(MpiComm,kb_recv_sctr)
       !
       call Max_MPI(MpiComm,t_setup_sb_sparse)
       call Max_MPI(MpiComm,t_hxv_sparse)
       !
       call Max_MPI(MpiComm,t_setup_sb_direct)
       call Sum_MPI(MpiComm,kb_sb_setup_bcast)
       !
       call Max_MPI(MpiComm,t_hxv_direct)
       call Max_MPI(MpiComm,t_hxv_1LxHR)
       call Max_MPI(MpiComm,t_hxv_HLx1R)
       call Max_MPI(MpiComm,t_hxv_AxB)
       call Max_MPI(MpiComm,t_hxv_B)
       call Max_MPI(MpiComm,t_hxv_A)
       !

       call Max_MPI(MpiComm,t_a2av)
       call Sum_MPI(MpiComm,kb_sent_a2av)
       call Sum_MPI(MpiComm,kb_recv_a2av)
       !
       call Max_MPI(MpiComm,t_rdm_get)
       call Max_MPI(MpiComm,t_gthr)
       call Sum_MPI(MpiComm,kb_sent_gthr)
       call Sum_MPI(MpiComm,kb_recv_gthr)
       !
       call Max_MPI(MpiComm,t_rdm_diag)
       call Max_MPI(MpiComm,t_rdm_renorm)
       !
       call Max_MPI(MpiComm,t_dmrg_step)
    endif
#endif
    !
    call append_profile("total_length_sb",dble(left%length + right%length))
    call append_profile("reduced_sb_dim",rdcd_sb_dim)
    call append_profile("kb_vecDim_av",kb_vecDim)
    call append_profile("full_sb_dim",full_sb_dim)
    call append_profile("time_enlarge_blocks",t_enlarge_blocks) !>= t_connect_blocks
    !
    call append_profile("time_sb_get_states",t_sb_get_states)
    call append_profile("kb_agthr_sb_states",kb_agthr_sb_states)
    !
    call append_profile("time_sb_diag",t_sb_diag)
    call append_profile("kb_sent_sctr",kb_sent_sctr)
    call append_profile("kb_recv_sctr",kb_recv_sctr)
    !
    call append_profile("time_setup_sb_sparse",t_setup_sb_sparse) !>= t_connect_blocks
    call append_profile("time_hxv_sparse",t_hxv_sparse)
    !
    call append_profile("time_setup_sb_direct",t_setup_sb_direct)
    call append_profile("kb_sb_setup_bcast",kb_sb_setup_bcast)    
    !
    call append_profile("time_hxv_direct",t_hxv_direct)
    call append_profile("time_hxv_1LxHR",t_hxv_1LxHR)
    call append_profile("time_hxv_HLx1R",t_hxv_HLx1R)
    call append_profile("time_hxv_AxB",t_hxv_AxB)
    call append_profile("time_hxv_B",t_hxv_B)
    call append_profile("time_hxv_A",t_hxv_A)
    !
    call append_profile("time_a2av",t_a2av)
    call append_profile("kb_sent_a2av",kb_sent_a2av)
    call append_profile("kb_recv_a2av",kb_recv_a2av)
    !
    call append_profile("HxV_Nop",dble(NumOp))
    !
    call append_profile("time_rdm_get",t_rdm_get)
    call append_profile("kb_sent_gthr",kb_sent_gthr)
    call append_profile("kb_recv_gthr",kb_recv_gthr)
    call append_profile("time_rdm_diag",t_rdm_diag)
    call append_profile("time_rdm_renorm",t_rdm_renorm)
    !
    call append_profile("time_dmrg_step",t_dmrg_step)
    !
    if(MpiMaster)then
       unit=fopen("full_profile_"//to_lower(DMRGtype)//"DMRG.out",append=.true.)
       write(unit,*)"# STEP:",left%length
       do i=1,size(prof_times)
          write(unit,*)prof_names(i),prof_times(i)
       enddo
       write(unit,*)""
       close(unit)
    endif
#endif
    call reset_profile()
    !
  end subroutine set_profile






  subroutine append_profile(string,time)
    character(len=*) :: string
    real(8)          :: time
    if(MpiMaster)call append(prof_names,string)
    if(MpiMaster)call append(prof_times,time)
  end subroutine append_profile






  !##################################################################
  !##################################################################
  !                   AUX FUNCTIONS
  !##################################################################
  !##################################################################
  function sector_qn_dim(self,qn) result(dim)
    type(sectors_list)   :: self
    real(8),dimension(:) :: qn
    integer              :: dim
    dim = 0
    if(.not.self%has_qn(qn))return
    dim =  size(self%map(qn=qn))
  end function sector_qn_dim





  subroutine dmrg_graphic(label,link)
    integer                   :: label
    character(len=*),optional :: link
    character(len=1)          :: ilink
    integer                   :: i,N,Mleft,Mright,LMleft,LMright,index,Ltot
    character(:),allocatable  :: Ldot,Rdot
    real(8)                   :: eps=1d-6
    integer                   :: M=100
    character(len=3)         :: br,bl,tr,tl
    !
    ! turn each integer (byte) into a 1-byte character and concatenate
    bl = char(int(Z'E2')) // char(int(Z'94')) // char(int(Z'94')) !└
    br = char(int(Z'E2')) // char(int(Z'94')) // char(int(Z'98')) !┘
    tl = char(int(Z'E2')) // char(int(Z'94')) // char(int(Z'8C')) !┌
    tr = char(int(Z'E2')) // char(int(Z'94')) // char(int(Z'90')) !┐
    !
    !
    ilink='n';if(present(link))ilink=to_lower(link(1:1))
    !
    Ltot = Ldmrg
    Ldot = bold_green('=')
    Rdot = bold_red('=')
    !
    N = int(Ltot/(M+eps))+1
    !
    ! call execute_command_line("clear")
    do i=1,3
       write(LOGfile,*)""
    enddo
    Mleft  = int(left%length/(N+eps))+1
    Mright = int(right%length/(N+eps))+1
    index  = nint(mod(dble(left%length),N+eps))
    LMleft = Ltot/N-Mleft
    LMright= Ltot/N-Mright
    !
    write(LOGfile,"(A,2I4,2x)",advance="no")"left; right=",left%length+1,right%length+1
    if(PBCdmrg)then
       !
       write(LOGfile,"(*(A1))",advance="no")(" ",i=1,LMleft-1)
       select case(iLink)
       case default
          write(LOGfile,"(A3,A1)",advance="no")tl,"|"
          write(LOGfile,"(*(A))",advance="no")(trim(Ldot),i=1,Mleft)
          write(LOGfile,"(A)",advance="no")bold_green("*")//bold("|")//bold_red("*")
          write(LOGfile,"(*(A))",advance="no")(trim(Rdot),i=1,Mright)
          write(LOGfile,"(A1,A3)",advance='no')"|",tr
       case('p')
          write(LOGfile,"(A,A1,A)",advance="no")tl,"|",bold_green("*")
          write(LOGfile,"(*(A))",advance="no")(trim(Ldot),i=1,Mleft)
          write(LOGfile,"(A)",advance="no")bold("|")
          write(LOGfile,"(*(A))",advance="no")(trim(Rdot),i=1,Mright)
          write(LOGfile,"(A,A1,A)",advance="no")bold_red("*"),"|",tr
       end select
       ! write(LOGfile,"(*(A))",advance="no")(" ",i=1,LMright-1)
       if(Ltot<=M)then
          write(LOGfile,"(A1)",advance='yes')""
       else
          write(LOGfile,"(2x,I3,2x,A,1x,A)",advance='yes')index,Ldot//":"//str(N)//&
               "- ;",Rdot//":"//str(N)//"-"
       endif
       write(LOGfile,"(A,10x)",advance="no")"            "
       write(LOGfile,"(*(A1))",advance="no")(" ",i=1,LMleft-1)
       write(LOGfile,"(A)",advance="no")bl
       write(LOGfile,"(*(A))",advance="no")(".",i=1,Mleft+Mright+5)
       write(LOGfile,"(A)",advance="no")br
       ! write(LOGfile,"(*(A))",advance="no")(" ",i=1,LMright-1)
       write(LOGfile,"(A1)",advance='yes')""
       !
    else
       !
       select case(label)
       case default
          write(LOGfile,"(A1)",advance="no")"|"
          write(LOGfile,"(*(A))",advance="no")(" ",i=1,LMleft)
          write(LOGfile,"(*(A))",advance="no")(trim(Ldot),i=1,Mleft)
          write(LOGfile,"(A)",advance="no")bold_green("*")//bold("|")//bold_red("*")
          write(LOGfile,"(*(A))",advance="no")(trim(Rdot),i=1,Mright)
          write(LOGfile,"(*(A))",advance="no")(" ",i=1,LMright)
          write(LOGfile,"(A1)",advance='no')"|"
       case(1)
          write(LOGfile,"(A1)",advance="no")"|"
          write(LOGfile,"(*(A))",advance="no")(trim(Ldot),i=1,Mleft)
          write(LOGfile,"(A)",advance="no")bg_green(">")//"|"//bold_red("*")
          write(LOGfile,"(*(A))",advance="no")(trim(Rdot),i=1,Mright)
          write(LOGfile,"(A1)",advance="no")"|"
       case(2)
          write(LOGfile,"(A1)",advance="no")"|"
          write(LOGfile,"(*(A))",advance="no")(trim(Ldot),i=1,Mleft)
          write(LOGfile,"(A)",advance="no")bold_green("*")//"|"//bg_red("<")
          write(LOGfile,"(*(A))",advance="no")(trim(Rdot),i=1,Mright)
          write(LOGfile,"(A1)",advance="no")"|"
       end select
       if(Ltot<=M)then
          write(LOGfile,"(A1)",advance='yes')""
       else
          write(LOGfile,"(2x,I3,2x,A,1x,A)",advance='yes')&
               index,Ldot//"->"//str(N)//"= ;",Rdot//"->"//str(N)//"-"
       endif

    endif
    !
    call wait(250)
    !
  end subroutine dmrg_graphic







END MODULE DMRG_GLOBAL





