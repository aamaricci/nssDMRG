MODULE DMRG_SUPERBLOCK
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_SETUP
  USE MPI
  implicit none
  private


  public :: sb_get_states
  public :: sb_get_shares
  public :: sb_diag

  public :: sb_build_Hv
  public :: sb_vecDim_Hv
  public :: sb_delete_Hv

  integer                          :: ispin
  integer                          :: iorb,jorb
  integer                          :: io,jo



contains




  !-----------------------------------------------------------------!
  ! purpose: build the list of states compatible with the specified
  ! quantum numbers
  !SUPERBLOCK SHARED THINGS:
  ! integer,dimension(:),allocatable   :: sb_states
  ! type(sectors_list)                 :: sb_sector
  !-----------------------------------------------------------------!
  subroutine sb_get_states()
    integer                          :: ql,ipr,rDim,total_states
    integer                          :: ir,il,istate,unit,k,Nsl,i0,i1
    integer                          :: Istart,Istep,Ierr
    real(8),dimension(:),allocatable :: left_qn,right_qn
    integer,dimension(:),allocatable :: left_map,right_map
    integer,dimension(:),allocatable :: Astates,Istates
    integer,dimension(:),allocatable :: Nl,Nr,Offset,Nk
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock get states"
#endif
    !
    if(MpiMaster)call start_timer("Build SB states")
    t0=t_start()
    !
    !INIT SB STATES OBJECTS:
    call sb_del_states()
    !
    Nsl = size(left%sectors(1))
    if(Nsl==0)then
       if(MpiMaster)write(LOGfile,*)"sb_get_states error: Nsl==0, blocks seem empty."
       stop
    endif
    !ipr = max(2,Nsl/10)!; if(ipr==0)ipr = 2
    allocate(Nl(Nsl),Nr(Nsl),Offset(Nsl),Nk(Nsl))
    Nl=0
    Nr=0
    Nk=0
    Offset=0
    !
    if(MpiMaster)rDim=right%Dim
#ifdef _MPI
    if(MpiStatus)call Bcast_MPI(MpiComm,rDim)
#endif
    !
    !Pre-computation loop:
    total_states = 0
    do ql = 1, Nsl
       left_qn  = left%sectors(1)%qn(index=ql)
       right_qn = current_target_qn - left_qn
       if (.not. right%sectors(1)%has_qn(right_qn)) cycle
       left_map  = left%sectors(1)%map(qn=left_qn)
       right_map = right%sectors(1)%map(qn=right_qn)
       Nl(ql) = size(left_map)
       Nr(ql) = size(right_map)
       Nk(ql) = Nl(ql) * Nr(ql)
       total_states = total_states + Nk(ql)
    end do
    !Get index Offset among sectors
    Offset(1) = 0
    do ql = 2, Nsl
       Offset(ql) = Offset(ql-1) + Nk(ql-1)
    end do
    !Check sanity
    if(total_states==0)then
       if(MpiMaster) call stop_timer("Build SB states")
       t_sb_get_states=t_stop()
       stop "sb_get_states ERROR: total_states=0. There are no SB states."
    else
       if(MpiMaster)write(LOGfile,*)"Total States:",total_states
    endif
    !
    !Main work here:
    Istart=1
    Istep =1
#ifdef _MPI
    if(MpiStatus)then
       Istart=1+MpiRank
       Istep =MpiSize
    endif
#endif
    if(allocated(sb_states))deallocate(sb_states)
    allocate(sb_states(total_states));sb_states = 0
    allocate(Astates(total_states))  ;Astates   = 0
    do ql=1,Nsl
       left_qn   = left%sectors(1)%qn(index=ql)
       right_qn  = current_target_qn - left_qn
       if(.not.right%sectors(1)%has_qn(right_qn))cycle
       left_map  = left%sectors(1)%map(qn=left_qn)
       right_map = right%sectors(1)%map(qn=right_qn)
       !
       !There is no real advantage in doing this MPI
       do il=istart,Nl(ql),istep
          i0 = 1+(il-1)*Nr(ql) + Offset(ql)
          i1 = il*Nr(ql)       + Offset(ql)
          sb_states( i0 : i1 ) = right_map(:) + (left_map(il)-1)*rDim
          Astates( i0 : i1 ) = (/( k, k=1, Nr(ql) )/) + i0-1 !== Offset(ql) + (il-1)*Nr(ql)
       enddo
       if(MpiMaster.AND.Nl(ql)>10)&           !.AND.mod(ql,ipr)==0)&
            write(LOGfile,*)"ql:"//str(ql)//"/"//str(Nsl)//" N(ql):"//str(Nl(ql))
    enddo
    !
    !Bulk MPI reduction
#ifdef _MPI
    if (MpiStatus) then
       call MPI_ALLREDUCE(MPI_IN_PLACE, sb_states, total_states, &
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE, Astates, total_states,   &
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
    end if
#endif
    !
    ! Finalize sector states calculation
    do ql = 1, Nsl
       if (Nk(ql) == 0) cycle
       left_qn  = left%sectors(1)%qn(index=ql)
       i0 = 1     + Offset(ql)
       i1 = Nk(ql)+ Offset(ql)
       call sb_sector%appends(qn=left_qn, istates=Astates(i0 : i1) )
    enddo
    deallocate(Astates)
    !
    if(MpiMaster)call stop_timer("Build SB states")
    t_sb_get_states=t_stop()
    !
  end subroutine sb_get_states






  subroutine sb_del_states()
    !    
    if(allocated(sb_states))deallocate(sb_states)
    call sb_sector%free()
    !
  end subroutine sb_del_states





  subroutine sb_get_shares()
    integer :: unit,q
    call sb_build_dims()
    if(MpiMaster)then
       unit=fopen("sb_shares_"//to_lower(DMRGtype)//"DMRG.out",append=.true.)
       write(unit,*)"# STEP:",left%length
       do q=1,size(sb_sector)
#ifdef _MPI
          if(MpiStatus)then
             write(unit,*)q,mpiDls(q)*Drs(q),mpiDls(q),Drs(q)
          else
             write(unit,*)q,Dls(q)*Drs(q),Dls(q),Drs(q)
          endif
#else
          write(unit,*)q,Dls(q)*Drs(q),Dls(q),Drs(q)
#endif
       enddo
       write(unit,*)""
       flush(unit)
       close(unit)
    endif
    call Barrier_MPI(MpiComm)
    call sb_delete_dims()
  end subroutine sb_get_shares


  !-----------------------------------------------------------------!
  ! Purpose: Diagonalize the SuperBlock problem.
  !-----------------------------------------------------------------!
  subroutine sb_diag()
    integer                               :: m_sb
    integer                               :: Nitermax,Neigen,Nblock
    real(8),dimension(:),allocatable      :: evals
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hsb,gs_tmp
#else
    real(8),dimension(:,:),allocatable    :: Hsb,gs_tmp
#endif
    integer                               :: vecDim,Nloc,m_tmp
    logical                               :: exist,lanc_solve,fMpi
    real(8) :: EH
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock diagonalization"
#endif
    !
    m_sb = size(sb_states)
    if(m_sb==0)stop "sb_diag ERROR: size(sb_states)==0"
    !
    !Set Lanczos params
    Neigen   = min(m_sb,lanc_neigen)
    Nitermax = min(m_sb,lanc_niter)
    Nblock   = min(m_sb,lanc_ncv_factor*Lanc_Neigen+lanc_ncv_add)
    !
    !Decide how to operate on H_sb
    lanc_solve = .true.
    if(Neigen==m_sb)lanc_solve=.false.
    if(m_sb  <= lanc_dim_threshold)lanc_solve=.false.
    !
    !Allocate EigPairs:
    if(allocated(gs_energy))deallocate(gs_energy)
    allocate(gs_energy(Neigen));gs_energy=zero
    !
    !
    if(lanc_solve)then          !Use (P)-Arpack
       !
       call sb_build_Hv()
       !
       if(allocated(gs_vector))deallocate(gs_vector)
       vecDim = sb_vecDim_Hv()
       !
       allocate(gs_vector(vecDim,Neigen));gs_vector=zero
       !
       if(MpiMaster)call start_timer("Diag H_sb")
       t0=t_start()
#ifdef _MPI
       if(MpiStatus)then
          !This condition protect against problems too small compared to MpiSize:
          !Solve serial and scatter the result
          call sb_check_Hv(fMPI)
          if(fMpi)then
             allocate(gs_tmp(m_sb,Neigen))
             if(MpiMaster)call sp_eigh(spHtimesV_p,gs_energy,gs_tmp,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,&
                  iverbose=(verbose>4),NumOp=NumOp)
             call Bcast_MPI(MpiComm,gs_energy)
             call sb_build_dims()
             call scatter_vector_MPI(MpiComm,gs_tmp,gs_vector)
             call sb_delete_dims()
             deallocate(gs_tmp)
          else
             call sp_eigh(MpiComm,spHtimesV_p,gs_energy,gs_vector,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,&
                  iverbose=(verbose>4),NumOp=NumOp)
          endif
       else
          call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
               Nblock,&
               Nitermax,&
               tol=lanc_tolerance,&
               iverbose=(verbose>4),NumOp=NumOp)
       endif
#else
       call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
            Nblock,&
            Nitermax,&
            tol=lanc_tolerance,&
            iverbose=(verbose>4),NumOp=NumOp)
#endif
       if(MpiMaster)call stop_timer("Diag H_sb")
       t_sb_diag=t_stop()
       !        allocate(gs_tmp,mold=gs_vector)
       !        call spHtimesV_p(vecDim,gs_vector(:,1),gs_tmp(:,1))
       !        EH = dot_product(gs_vector(:,1),  gs_tmp(:,1))
       ! #ifdef _MPI
       !        if(MpiStatus)call Bcast_MPI(MpiComm,EH)
       ! #endif
       !        write(200,*)current_L,EH/current_L
       !
    else !use LAPACK
       !
       call sb_build_Hv(Hsb)
       !
       if(allocated(gs_vector))deallocate(gs_vector)
       vecDim = sb_vecDim_Hv()
       allocate(gs_vector(vecDim,Neigen));gs_vector=zero
       allocate(evals(m_sb))
       !
       if(MpiMaster)call start_timer("Diag H_sb")
       t0=t_start()
       if(MpiMaster)call eigh(Hsb,evals)
       if(MpiMaster)call stop_timer("Diag H_sb")
       t_sb_diag=t_stop()
       !
#ifdef _MPI
       if(MpiStatus)then
          call Bcast_MPI(MpiComm,evals)
          call sb_build_dims()
          call scatter_vector_MPI(MpiComm,Hsb(:,1:Neigen),gs_vector)
          gs_energy(1:Neigen)   = evals(1:Neigen)
          call sb_delete_dims()
       else
          gs_vector(:,1:Neigen) = Hsb(:,1:Neigen)
          gs_energy(1:Neigen)   = evals(1:Neigen)
       endif
#else
       gs_vector(:,1:Neigen) = Hsb(:,1:Neigen)
       gs_energy(1:Neigen)   = evals(1:Neigen)          
#endif
       !
       deallocate(Hsb,evals)
       !
    endif
    !
    !Free Memory
    call sb_delete_Hv()
    !
  end subroutine sb_diag









  !##################################################################
  !         SETUP THE SUPERBLOCK HAMILTONIAN PROBLEM
  ! . if Hmat: returb H^SB as dense matrix there for Lapack use
  ! . if sparse_H = T: build H^SB as sparse matrix
  ! . if sparse_H = F: setup H^SB terms and blocks for H*v procedure
  !##################################################################
  subroutine sb_build_Hv(Hmat)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable,optional :: Hmat
#else
    real(8),dimension(:,:),allocatable,optional    :: Hmat
#endif
    real(8),dimension(:),allocatable               :: qn
    integer                                        :: q,m_sb,Nsb,irank,vecDim
    logical :: fMPI
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: SuperBlock build H*v"
#endif
    !
    if(.not.allocated(sb_states))stop "build_Hv_superblock ERROR: sb_states not allocated"
    m_sb = size(sb_states)
    Nsb  = size(sb_sector)
    !
#ifdef _MPI
    if(allocated(Dls))deallocate(Dls)
    if(allocated(Drs))deallocate(Drs)
    if(allocated(mpiDls))deallocate(mpiDls)
    allocate(Dls(Nsb),Drs(Nsb),mpiDls(Nsb))
    do q=1,Nsb
       qn        = sb_sector%qn(index=q)
       Dls(q)    = sector_qn_dim(left%sectors(1),qn)
       Drs(q)    = sector_qn_dim(right%sectors(1),current_target_qn - qn)
       mpiDls(q) = Dls(q)/MpiSize
       if(MpiRank< mod(Dls(q),MpiSize))mpiDls(q) = mpiDls(q)+1
    enddo
#ifdef _DEBUG
    if(MpiStatus.AND.verbose>4.AND.(MpiComm/=Mpi_Comm_Null).AND.MpiSize>=1)then
       if(MpiMaster)write(*,*)"         mpiRank,          mpiDls        -        mpiDl        -      mpiL    -   mpiOffset"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)write(*,*)MpiRank,mpiDls,"-",&
               Drs(:)*mpiDls(:),"-",sum(Drs(:)*mpiDls(:)),"-",&
               (sum(Drs(1:q-1)*mpiDls(1:q-1)),q=1,Nsb)
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
    deallocate(Dls,Drs,mpiDls)
#endif
    !
    !IF PRESENT HMAT: get SB_H sparse > dump it to dense Hmat > return
    if(present(Hmat))then
       if(allocated(Hmat))deallocate(Hmat)
       allocate(Hmat(m_sb,m_sb));Hmat=zero
       !Nullify HxV function pointer:
       spHtimesV_p => null()
       !
       !>Build Sparse Hsb: (Root only)
       if(MpiMaster)then
          call Setup_SuperBlock_Sparse() !<- no MPI here
          !
          !Dump Hsb to dense matrix as required:
          write(LogFile,"(A)")"LAPACK: spHsb.dump(Hmat)"
          call spHsb%dump(Hmat)
       endif
    else
       !
       !Build SuperBLock HxV operation: stored or direct
       select case(sparse_H)
       case(.true.)
          call Setup_SuperBlock_Sparse() !<- no MPI here yet
          spHtimesV_p => spMatVec_sparse_main
          !
       case(.false.)
          call Setup_SuperBlock_Direct() !<- SETUP MPI here
          spHtimesV_p => spMatVec_direct_main
#ifdef _MPI
          if(MpiStatus)then
             call sb_check_Hv(fMPI)
             if(fMPI)then  !this is true for all nodes at once see sb_vecDim_Hv
                write(LogFile,"(A)")"ARPACK: using SERIAL over MPI as MpiSize > N"
                spHtimesV_p => spMatVec_direct_main
             else
                spHtimesV_p => spMatVec_MPI_direct_main
             endif
          endif
#endif
          !
       end select
    endif
    !
  end subroutine sb_build_Hv





  subroutine sb_delete_Hv()
    integer :: i,j
    !
    spHtimesV_p => null()
    !
    call spHsb%free()
    if(allocated(Hleft))then
       do concurrent(i=1:size(Hleft))
          call Hleft(i)%free()
       enddo
       deallocate(Hleft)
    endif
    if(allocated(Hright))then
       do concurrent(i=1:size(Hright))
          call Hright(i)%free()
       enddo
       deallocate(Hright)
    endif
    if(allocated(A))then
       do concurrent(i=1:size(A,1),j=1:size(A,2))
          call A(i,j)%free()
       enddo
       deallocate(A)
    endif
    if(allocated(B))then
       do concurrent(i=1:size(B,1),j=1:size(B,2))
          call B(i,j)%free()
       enddo
       deallocate(B)
    endif
    !
    call sb_delete_dims()
    if(allocated(isb2jsb))deallocate(isb2jsb)
    if(allocated(IsHconjg))deallocate(IsHconjg)
    if(allocated(RowOffset))deallocate(RowOffset)
    if(allocated(ColOffset))deallocate(ColOffset)
    !
  end subroutine sb_delete_Hv






  function sb_vecDim_Hv() result(vecDim)
    integer                          :: vecDim           !vector or vector chunck dimension
    real(8),dimension(:),allocatable :: qn
    integer                          :: q,Nsb,i,unit
    !
    Nsb  = size(sb_sector)
    if(allocated(Dls))deallocate(Dls)
    if(allocated(Drs))deallocate(Drs)
    if(allocated(mpiDls))deallocate(mpiDls)
    if(allocated(mpiDl))deallocate(mpiDl)
    allocate(Dls(Nsb),Drs(Nsb),mpiDls(Nsb),mpiDl(Nsb))
    do q=1,Nsb
       qn     = sb_sector%qn(index=q)
       Dls(q) = sector_qn_dim(left%sectors(1),qn)
       Drs(q) = sector_qn_dim(right%sectors(1),current_target_qn - qn)
    enddo
#ifdef _MPI
    if(MpiStatus)then
       do q=1,Nsb
          mpiDls(q) = Dls(q)/MpiSize
          if(MpiRank < mod(Dls(q),MpiSize))mpiDls(q) = mpiDls(q)+1
          mpiDl(q)  = Drs(q)*mpiDls(q)
       enddo
    else
       do q=1,Nsb
          mpiDls(q)= Dls(q)
          mpiDl(q)  = Drs(q)*mpiDls(q)
       enddo
       if(sum(mpiDl)/=size(sb_states))stop "sb_vecDim_Hv error: no MPI but vecDim != m_sb"
    endif
#else
    do q=1,Nsb
       mpiDls(q)= Dls(q)
       mpiDl(q)  = Drs(q)*mpiDls(q)
    enddo
    if(sum(mpiDl)/=size(sb_states))stop "sb_vecDim_Hv error: no MPI but vecDim != m_sb"
#endif
    !
    vecDim = sum(mpiDl)
    kb_vecDim = vecDim*DATA_kb_size
  end function sb_vecDim_Hv





  subroutine sb_check_Hv(anyZero)
    integer :: vecDim           !vector or vector chunck dimension
    integer :: ierr
    logical :: hasZero,anyZero
    !
    anyZero=.false.
    !
#ifdef _MPI
    if(MpiStatus)then
       vecDim  = sb_vecDim_Hv()
       hasZero = (vecDim == 0)  !T if a rank has vecDim==0
       call MPI_ALLREDUCE(hasZero, anyZero, 1, MPI_LOGICAL, MPI_LOR, MpiComm, ierr)
       !anyZero=T if at least 1 nodes has vecDim==0
    end if
#endif
  end subroutine sb_check_Hv

END MODULE DMRG_SUPERBLOCK













