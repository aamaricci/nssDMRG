MODULE DMRG_SUPERBLOCK
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_SETUP
  USE MPI
  implicit none
  private


  public :: sb_build_dims
  public :: sb_delete_dims
  public :: sb_get_states
  public :: sb_get_shares
  public :: sb_diag
  public :: sb_set_current_qn
  public :: sb_build_Hv
  public :: sb_vecDim_Hv
  public :: sb_delete_Hv

  integer :: ierr
  integer :: ispin
  integer :: iorb,jorb
  integer :: io,jo



contains



  !##############################################################
  !##############################################################
  !       BUILD THE SUPERBLOCK DIMENSIONS (MPI SHARES)
  !##############################################################
  !##############################################################
  subroutine sb_build_dims(quiet)
    logical,optional                 :: quiet  
    integer                          :: Nsb
    integer                          :: isb
    real(8),dimension(:),allocatable :: qn
    integer                          :: color,n_active,unit
    logical                          :: quiet_
    !
    quiet_=.true.;if(present(quiet))quiet_=quiet
    !
    Nsb  = size(sb_sector)
    if(Nsb==0)stop "sb_build_dims error: size(sb_sector)==0"
    !
    call sb_delete_dims()
    !
    allocate(Dls(Nsb))
    allocate(Drs(Nsb))
    allocate(Offset(Nsb))
#ifdef _MPI  
    allocate(mpiDls(Nsb))
    allocate(mpiDrs(Nsb))
    allocate(mpiDl(Nsb))
    allocate(mpiDr(Nsb))
    allocate(mpiOffset(Nsb))
    allocate(mpiSBCOMM(Nsb))
    allocate(mpiNactive(Nsb)) 
    mpiOffset=0
    mpiSBCOMM=Mpi_Comm_Null
    mpiNactive=MpiSize
#endif
    !
    Offset=0
    do isb=1,Nsb
      qn      = sb_sector%qn(index=isb)
      Dls(isb)= sector_qn_dim(left%sectors(1),qn)
      Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
      Offset(isb)=sum(Dls(1:isb-1)*Drs(1:isb-1))
      !
#ifdef _MPI
      ! MPI VARS SETUP 
      mpiDls(isb) = Dls(isb)/MpiSize
      mpiDrs(isb) = Drs(isb)/MpiSize
      if(MpiRank < mod(Dls(isb),MpiSize))mpiDls(isb) = mpiDls(isb)+1
      if(MpiRank < mod(Drs(isb),MpiSize))mpiDrs(isb) = mpiDrs(isb)+1
      mpiDl(isb) = Drs(isb)*mpiDls(isb)
      mpiDr(isb) = mpiDrs(isb)*Dls(isb)
      !
      mpiOffset(isb)=sum(Drs(1:isb-1)*mpiDls(1:isb-1))
      !
      !Get how many active ranks for this sector: 
      n_active = max( min(Dls(isb),MpiSize) , min(Drs(isb),MpiSize) )
      !Here we SPLIT the communicators because we do not want the idle ranks to 
      !be completely shut off the computation but we want them to wait for the active ranks to 
      !finish and then move on together to the next sector.
      if(n_active >= MpiSize )then
         mpiNactive(isb) = MpiSize
         mpiSBCOMM(isb)  = MpiComm   !All ranks are active in this sector
      else
         mpiNactive(isb) = n_active
         color           = MPI_UNDEFINED
         if(MpiRank < mpiNactive(isb))color = 1
         call MPI_Comm_split(MpiComm_Global, color, MpiRank, mpiSBCOMM(isb), ierr)
         if(.not.quiet_)then
            if(MpiMaster)then
               if(verbose>4)write(LOGfile,"(A,I6,A,I6,A,I6,A,I6,A,I6)")"Reducing N_cpu to Nactive:",n_active,&
               " for sector ",isb,"/",Nsb," with Dls=",Dls(isb)," and Drs=",Drs(isb)
               unit=fopen("sb_active_"//to_lower(DMRGtype)//"DMRG.dmrg",append=.true.)
               write(unit,*)left%length,isb,Dls(isb),Drs(isb),MpiSize,mpiNactive(isb)
               flush(unit)
               close(unit)
            endif
#ifdef _DEBUG         
            if(mpisbcomm(isb) == mpi_comm_null) then
               write(logfile,*) "rank ", mpirank, " is correctly idle (comm is null)"
            else
               write(logfile,*) "rank ", mpirank, " is active in comm ", mpisbcomm(isb)
            endif
#endif                     
         endif
      endif		
#endif
    enddo
    !
    dims_set=.true.
    !	
  end subroutine sb_build_dims

  
  subroutine sb_delete_dims
    integer :: q
    if(allocated(Dls))       deallocate(Dls)
    if(allocated(Drs))       deallocate(Drs)
    if(allocated(Offset))    deallocate(Offset)
#ifdef _MPI
    if(allocated(mpiDls))    deallocate(mpiDls)
    if(allocated(mpiDrs))    deallocate(mpiDrs)
    if(allocated(mpiDl))     deallocate(mpiDl)
    if(allocated(mpiDr))     deallocate(mpiDr)
    if(allocated(mpiOffset)) deallocate(mpiOffset)
    if(allocated(mpiNactive))deallocate(mpiNactive)
    if(allocated(mpiSBComm))then
      do q=1, size(mpiSBComm)
        if(mpiSBComm(q) /= MPI_COMM_NULL .AND. mpiSBComm(q) /= MpiComm_Global) &
        call MPI_Comm_free(mpiSBComm(q), ierr)
      end do
      deallocate(mpiSBComm)
    endif
    MpiComm_Global=MpiComm
#endif
    !
    dims_set=.false.
    !
  end subroutine sb_delete_dims








    
  !##############################################################
  !##############################################################
  !            BUILD THE SUPERBLOCK STATES LIST
  !##############################################################
  !##############################################################
  !shared objects:
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
    integer :: unit,q,irank
    call sb_build_dims(quiet=.true.)
    if(MpiMaster)then
      unit=fopen("sb_shares_"//to_lower(DMRGtype)//"DMRG.dmrg",append=.true.)
      write(unit,*)"# STEP:",left%length
      do q=1,size(sb_sector)
#ifdef _MPI
         if(MpiStatus)then
            write(unit,*)q,Drs(q),Dls(q),mpiDls(q),MpiSize
         else
            write(unit,*)q,Drs(q),Dls(q)
         endif
#else
         write(unit,*)q,Drs(q),Dls(q)
#endif
      enddo
      write(unit,*)""
      flush(unit)
      close(unit)
    endif
    call sb_delete_dims()
  end subroutine sb_get_shares






  !##############################################################
  !##############################################################
  !            Diagonalize the SuperBlock problem               !
  !##############################################################
  !##############################################################
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
#ifdef _CMPLX
    complex(8) :: uno=(1.d0,0.d0)
#else
    real(8)    :: uno=1.d0
#endif    

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
                  vrandom=(.not.lanc_v0_dble),&
                  iverbose=(verbose>4),NumOp=NumOp)
             call Bcast_MPI(MpiComm,gs_energy)
             call scatter_vector_MPI(MpiComm,gs_tmp,gs_vector)
             deallocate(gs_tmp)
          else
             call sp_eigh(MpiComm,spHtimesV_p,gs_energy,gs_vector,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,&
                  vrandom=(.not.lanc_v0_dble),&
                  iverbose=(verbose>4),NumOp=NumOp)
          endif
       else
          call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
               Nblock,&
               Nitermax,&
               tol=lanc_tolerance,&
               vrandom=(.not.lanc_v0_dble),&
               iverbose=(verbose>4),NumOp=NumOp)
       endif
#else
       call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
            Nblock,&
            Nitermax,&
            tol=lanc_tolerance,&
            vrandom=(.not.lanc_v0_dble),&
            iverbose=(verbose>4),NumOp=NumOp)
#endif
       if(MpiMaster)call stop_timer("Diag H_sb")
       t_sb_diag=t_stop()
       !
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
          call scatter_vector_MPI(MpiComm,Hsb(:,1:Neigen),gs_vector)
          gs_energy(1:Neigen)   = evals(1:Neigen)
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
  !##################################################################
  !         SETUP THE SUPERBLOCK HAMILTONIAN PROBLEM
  !##################################################################
  ! . if Hmat: return H^SB as dense matrix there for Lapack use
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
    !
    call sb_build_dims(quiet=.false.)
    !    
#ifdef _MPI
#ifdef _DEBUG
    if(MpiStatus.AND.verbose>4.AND.(MpiComm/=Mpi_Comm_Null).AND.MpiSize>=1)then
        if(MpiMaster)then 
          write(LOGfile,"(A6,A3,A9)",advance="no")"Irank"," | ","mpiDls(q)" 
          write(LOGfile,"(*(A9))",advance="no")(str(q),q=2,Nsb)
          write(LOGfile,"(A3,A6,A3,A12)",advance="no")" | ", "mpiDl"," | ","mpiOffset(q)"
          write(LOGfile,"(*(A9))",advance="yes")(str(q),q=2,Nsb)
        endif 
        do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)then
            write(LOGfile,"(I6,A3)",advance="no")irank," | "
            write(LOGfile,"(*(I9))",advance="no")(mpiDls(q),q=1,Nsb)
            write(LOGfile,"(A3,I6,A3)",advance="no")" | ",sum(Drs(:)*mpiDls(:))," | "
            write(LOGfile,"(*(I9))",advance="yes")(sum(Drs(1:q-1)*mpiDls(1:q-1)),q=1,Nsb)
          endif
        enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
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
    !
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
    if(.not.dims_set)stop "sb_vecDim_Hv error: dimensions not set. Call sb_build_dims() first."
    !
    Nsb  = size(sb_sector)
#ifdef _MPI
    if(MpiStatus)then
      vecDim = sum(mpiDl)
    else
      vecDim = sum(Drs*Dls)
      if(vecDim/=size(sb_states))stop "sb_vecDim_Hv error: no MPI but vecDim != m_sb"
    endif
#else
    vecDim = sum(Drs*Dls)
    if(vecDim/=size(sb_states))stop "sb_vecDim_Hv error: no MPI but vecDim != m_sb"
#endif
    kb_vecDim = vecDim*DATA_kb_size
  end function sb_vecDim_Hv

!   function sb_vecDim_Hv() result(vecDim)
!     integer                          :: vecDim           !vector or vector chunck dimension
!     real(8),dimension(:),allocatable :: qn
!     integer                          :: q,Nsb,i,unit
!     !
!     Nsb  = size(sb_sector)
!     if(allocated(Dls))deallocate(Dls)
!     if(allocated(Drs))deallocate(Drs)
!     if(allocated(mpiDls))deallocate(mpiDls)
!     if(allocated(mpiDl))deallocate(mpiDl)
!     allocate(Dls(Nsb),Drs(Nsb),mpiDls(Nsb),mpiDl(Nsb))
!     do q=1,Nsb
!        qn     = sb_sector%qn(index=q)
!        Dls(q) = sector_qn_dim(left%sectors(1),qn)
!        Drs(q) = sector_qn_dim(right%sectors(1),current_target_qn - qn)
!     enddo
! #ifdef _MPI
!     if(MpiStatus)then
!        do q=1,Nsb
!           mpiDls(q) = Dls(q)/MpiSize
!           if(MpiRank < mod(Dls(q),MpiSize))mpiDls(q) = mpiDls(q)+1
!           mpiDl(q)  = Drs(q)*mpiDls(q)
!        enddo
!     else
!        do q=1,Nsb
!           mpiDls(q)= Dls(q)
!           mpiDl(q)  = Drs(q)*mpiDls(q)
!        enddo
!        if(sum(mpiDl)/=size(sb_states))stop "sb_vecDim_Hv error: no MPI but vecDim != m_sb"
!     endif
! #else
!     do q=1,Nsb
!        mpiDls(q)= Dls(q)
!        mpiDl(q)  = Drs(q)*mpiDls(q)
!     enddo
!     if(sum(mpiDl)/=size(sb_states))stop "sb_vecDim_Hv error: no MPI but vecDim != m_sb"
! #endif
!     !
!     vecDim = sum(mpiDl)
!     kb_vecDim = vecDim*DATA_kb_size
!   end function sb_vecDim_Hv





  !##################################################################
  !##################################################################
  !                   AUX FUNCTIONS
  !##################################################################
  !##################################################################
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


  subroutine sb_set_current_qn()
    current_L         = left%length + right%length
    select case(str(to_lower(QNtype(1:1))))
    case default;stop "DMRG_MAIN error: QNtype != [local,global]"
    case("l")
       current_target_QN = int(target_qn*current_L*Norb)
    case("g")
       current_target_QN = min(current_L,int(target_qn*Norb)) !to check
    end select
  end subroutine sb_set_current_qn



END MODULE DMRG_SUPERBLOCK













