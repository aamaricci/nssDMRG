MODULE DMRG_SUPERBLOCK
  USE ISO_C_BINDING, only: c_double,c_int,c_int64_t
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK_SETUP
  USE DMRG_SUPERBLOCK_HXV
  USE SF_SP_LINALG, only: primme_matmat_d,primme_matmat_z,PRIMME_DEFAULT_MIN_MATVECS
#ifdef _MPI
   USE MPI
   USE SF_MPI
#endif
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
  public :: sb_save_measure_state
  public :: sb_load_measure_state
  public :: sb_save_measure_blocks
  public :: sb_load_measure_blocks

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
       if(MpiMaster)then
          write(LOGfile,*)"No retained superblock states for target:",current_target_QN
          write(LOGfile,*)"Length:",current_L," density:",target_density," offset:",target_offset
          write(LOGfile,*)"Check charge integrality, physical bounds and retained block sectors."
       endif
       error stop "sb_get_states ERROR: requested QN sector is empty"
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
    complex(8),allocatable                :: wf_seed_block(:,:)
#else
    real(8),dimension(:,:),allocatable    :: Hsb,gs_tmp
    real(8),allocatable                   :: wf_seed_block(:,:)
#endif
    integer                               :: vecDim,Nloc,m_tmp
    logical                               :: exist,lanc_solve,fMPI
    logical                               :: primme_failed
#ifdef _MPI
    logical                               :: primme_failed_global
    integer                               :: primme_mpierr
#endif
    integer(c_int64_t)                    :: primme_nblock,primme_nitermax
    integer(c_int64_t)                    :: primme_init_size,primme_numop,primme_niter,primme_nconv
    integer(c_int)                        :: primme_iexit
#ifdef _CMPLX
    procedure(primme_matmat_z),pointer    :: primme_matmat_p
#else
    procedure(primme_matmat_d),pointer    :: primme_matmat_p
#endif
#ifdef _CMPLX
    complex(8),allocatable                :: wf_seed(:)
#else
    real(8),allocatable                   :: wf_seed(:)
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
#ifdef _MPI
    ! PARPACK needs at least one local basis element on every process. When
    ! the full superblock is smaller than the MPI task count, use the dense
    ! LAPACK path directly instead of falling back to a serial ARPACK solve.
    if(MpiStatus.and.m_sb<MpiSize)then
       lanc_solve=.false.
       if(MpiMaster)write(LogFile,"(A,I0,A,I0,A)") &
            "LAPACK: superblock dimension ",m_sb," is smaller than MPI size ",MpiSize, &
            "; skipping PARPACK."
    endif
#endif
    !
    !Allocate EigPairs:
    if(allocated(gs_energy))deallocate(gs_energy)
    allocate(gs_energy(Neigen));gs_energy=zero
    !
    !
    if(lanc_solve)then
       !
       call sb_build_Hv()
       call sb_prepare_wf_guess(wf_seed,wf_seed_block)
       !
       if(allocated(gs_vector))deallocate(gs_vector)
       vecDim = sb_vecDim_Hv()
       !
       allocate(gs_vector(vecDim,Neigen));gs_vector=zero
       !
       if(MpiMaster)call start_timer("Diag H_sb "//to_upper(str(lanc_method)))
       t0=t_start()
       select case(lanc_method)
       case("arpack")
#ifdef _MPI
          if(MpiStatus)then
             call sb_check_Hv(fMPI)
             if(fMPI)then
                allocate(gs_tmp(m_sb,Neigen))
                if(MpiMaster)call sp_eigh(spHtimesV_p,gs_energy,gs_tmp,&
                     Nblock,&
                     Nitermax,&
                     tol=lanc_tolerance,&
                     iverbose=(verbose>4),NumOp=NumOp)
                call Bcast_MPI(MpiComm,gs_energy)
                call scatter_vector_MPI(MpiComm,gs_tmp,gs_vector)
                deallocate(gs_tmp)
             else
                call sp_eigh(MpiComm,spHtimesV_p,gs_energy,gs_vector,&
                     Nblock,&
                     Nitermax,&
                     tol=lanc_tolerance,v0=wf_seed,&
                     iverbose=(verbose>4),NumOp=NumOp)
             endif
          else
             call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,v0=wf_seed,&
                  iverbose=(verbose>4),NumOp=NumOp)
          endif
#else
          call sp_eigh(spHtimesV_p,gs_energy,gs_vector,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,v0=wf_seed,&
                  iverbose=(verbose>4),NumOp=NumOp)
#endif
       case("primme")
          primme_nblock=int(Nblock,c_int64_t)
          primme_nitermax=int(Nitermax,c_int64_t)
          primme_init_size=0_c_int64_t
          primme_numop=0_c_int64_t
          primme_niter=0_c_int64_t
          primme_nconv=0_c_int64_t
          primme_iexit=0_c_int
          primme_failed=.false.
#ifdef _MPI
          if(MpiStatus)then
             call sb_check_Hv(fMPI)
             if(fMPI)then
                allocate(gs_tmp(m_sb,Neigen))
                gs_tmp=zero
                if(MpiMaster)then
                   primme_matmat_p=>sb_primme_matmat_serial
                   if(allocated(wf_seed_block))then
                      primme_init_size=int(min(size(wf_seed_block,2),Neigen),c_int64_t)
                      gs_tmp(:,1:primme_init_size)=wf_seed_block(:,1:primme_init_size)
                   endif
                   call sp_eigh_primme(gs_energy,gs_tmp,primme_matmat_p,&
                        method=PRIMME_DEFAULT_MIN_MATVECS,tol=lanc_tolerance,&
                        Nitermax=primme_nitermax,Nblock=primme_nblock,&
                        initSize=primme_init_size,iexit=primme_iexit,&
                        NumOp=primme_numop,Niter=primme_niter,Nconv=primme_nconv)
                endif
                call Bcast_MPI(MpiComm,gs_energy)
                call scatter_vector_MPI(MpiComm,gs_tmp,gs_vector)
                deallocate(gs_tmp)
             else
                primme_matmat_p=>sb_primme_matmat_mpi
                if(allocated(wf_seed_block))then
                   primme_init_size=int(min(size(wf_seed_block,2),Neigen),c_int64_t)
                   gs_vector(:,1:primme_init_size)=wf_seed_block(:,1:primme_init_size)
                endif
                call sp_eigh_primme(MpiComm,gs_energy,gs_vector,primme_matmat_p,&
                     method=PRIMME_DEFAULT_MIN_MATVECS,tol=lanc_tolerance,&
                     Nitermax=primme_nitermax,Nblock=primme_nblock,&
                     initSize=primme_init_size,iexit=primme_iexit,&
                     NumOp=primme_numop,Niter=primme_niter,Nconv=primme_nconv)
             endif
          else
#endif
             primme_matmat_p=>sb_primme_matmat_serial
             if(allocated(wf_seed_block))then
                primme_init_size=int(min(size(wf_seed_block,2),Neigen),c_int64_t)
                gs_vector(:,1:primme_init_size)=wf_seed_block(:,1:primme_init_size)
             endif
             call sp_eigh_primme(gs_energy,gs_vector,primme_matmat_p,&
                  method=PRIMME_DEFAULT_MIN_MATVECS,tol=lanc_tolerance,&
                  Nitermax=primme_nitermax,Nblock=primme_nblock,&
                  initSize=primme_init_size,iexit=primme_iexit,&
                  NumOp=primme_numop,Niter=primme_niter,Nconv=primme_nconv)
#ifdef _MPI
          endif
#endif
          NumOp=int(primme_numop)
          primme_failed=(primme_iexit/=0)
#ifdef _MPI
          if(MpiStatus)then
             call MPI_Allreduce(primme_failed,primme_failed_global,1,MPI_LOGICAL,MPI_LOR, &
                  MpiComm,primme_mpierr)
             if(primme_mpierr/=MPI_SUCCESS) &
                  error stop "sb_diag: MPI_Allreduce of PRIMME status failed"
             primme_failed=primme_failed_global
          endif
#endif
          if(primme_failed)then
             if(MpiMaster)write(LOGfile,*) &
                  "PRIMME failed with status:",primme_iexit," converged states:",primme_nconv
             error stop "sb_diag ERROR: PRIMME did not complete successfully"
          endif
       case default
          error stop "sb_diag ERROR: unsupported LANC_METHOD"
       end select
       if(MpiMaster)call stop_timer("Diag H_sb")
       t_sb_diag=t_stop()
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
    if(allocated(wf_trial))deallocate(wf_trial)
    !
    if(lanc_solve.and.lanc_check_residual)call sb_check_ritz_residuals()
    !Free Memory
    call sb_delete_Hv()
    !
  end subroutine sb_diag


  ! Adapt the scalar superblock HxV callback to the PRIMME serial MatMat API.
  ! This keeps the serial and empty-rank fallback on the existing HxV path.
  subroutine sb_primme_matmat_serial(n,nvec,vin,vout)
    integer,intent(in) :: n,nvec
#ifdef _CMPLX
    complex(c_double),intent(in)  :: vin(n,nvec)
    complex(c_double),intent(out) :: vout(n,nvec)
#else
    real(c_double),intent(in)     :: vin(n,nvec)
    real(c_double),intent(out)    :: vout(n,nvec)
#endif
    integer :: iv
    !
    if(nvec<1)error stop "sb_primme_matmat_serial: empty vector block"
    do iv=1,nvec
       call spHtimesV_p(n,vin(:,iv),vout(:,iv))
    enddo
  end subroutine sb_primme_matmat_serial


#ifdef _MPI
  ! Adapt the MPI superblock callback to the PRIMME distributed MatMat API.
  ! Lazy HxV filters are shared across columns; cached and sparse paths reuse
  ! their existing scalar callback for each column.
  subroutine sb_primme_matmat_mpi(n,nvec,vin,vout)
    integer,intent(in) :: n,nvec
#ifdef _CMPLX
    complex(c_double),intent(in)  :: vin(n,nvec)
    complex(c_double),intent(out) :: vout(n,nvec)
#else
    real(c_double),intent(in)     :: vin(n,nvec)
    real(c_double),intent(out)    :: vout(n,nvec)
#endif
    integer :: iv
    !
    if(nvec<1)error stop "sb_primme_matmat_mpi: empty vector block"
    if(.not.MpiStatus)error stop "sb_primme_matmat_mpi: MPI is not active"
    if(direct_H_lazy.and..not.sparse_H)then
       call spMatMat_MPI_lazy_main(n,vin,vout)
    else
       do iv=1,nvec
          call spHtimesV_p(n,vin(:,iv),vout(:,iv))
       enddo
    endif
  end subroutine sb_primme_matmat_mpi
#endif


  ! Prepare the distributed target-state seed. Rank zero owns wf_trial; the
  ! eigensolver receives a local vector for ARPACK or a local block for PRIMME.
  subroutine sb_prepare_wf_guess(wf_seed,wf_seed_block)
#ifdef _CMPLX
    complex(8),allocatable,intent(out)      :: wf_seed(:),wf_seed_block(:,:)
    complex(8),allocatable                 :: wf_loc(:,:)
#else
    real(8),allocatable,intent(out)        :: wf_seed(:),wf_seed_block(:,:)
    real(8),allocatable                    :: wf_loc(:,:)
#endif
    logical                                :: wf_available,fMPI
    integer                                :: mpierr,nvalid
    ! Collectively decide whether a root-owned trial block is available.
    wf_available=.false.
    fMPI=.false.
    nvalid=0
    if(MpiMaster)then
       wf_available=allocated(wf_trial)
       nvalid=wf_trial_nvalid
    endif
#ifdef _MPI
    if(MpiStatus)then
       call MPI_Bcast(wf_available,1,MPI_LOGICAL,0,MpiComm,mpierr)
       if(mpierr/=MPI_SUCCESS)error stop "sb_prepare_wf_guess: MPI_Bcast availability failed"
       call MPI_Bcast(wf_trial_correlated,1,MPI_LOGICAL,0,MpiComm,mpierr)
       if(mpierr/=MPI_SUCCESS)error stop "sb_prepare_wf_guess: MPI_Bcast correlated flag failed"
       call MPI_Bcast(nvalid,1,MPI_INTEGER,0,MpiComm,mpierr)
       if(mpierr/=MPI_SUCCESS)error stop "sb_prepare_wf_guess: MPI_Bcast guess count failed"
    endif
#endif
    if(.not.wf_available.or.nvalid==0)return
#ifdef _MPI
    ! Empty-row ranks use a serial fallback. PRIMME can receive the complete
    ! root-owned block there; ARPACK retains its usual unseeded fallback.
    if(MpiStatus)then
       call sb_check_Hv(fMPI)
       if(fMPI)then
          if(lanc_method=="primme".and.lanc_target_prediction.and.wf_trial_correlated)then
             if(MpiMaster)allocate(wf_seed_block,source=wf_trial(:,1:nvalid))
          endif
          return
       endif
    endif
#endif
    ! Rank zero supplies the complete block; other ranks allocate a send-buffer
    ! placeholder, then every column follows the normal local SB vector layout.
    if(.not.allocated(wf_trial))then
       allocate(wf_trial(1,nvalid))
       wf_trial=zero
    endif
    allocate(wf_loc(sb_vecDim_Hv(),nvalid))
#ifdef _MPI
    if(MpiStatus)then
       call scatter_vector_MPI(MpiComm,wf_trial(:,1:nvalid),wf_loc)
    else
       wf_loc=wf_trial(:,1:nvalid)
    endif
#else
    wf_loc=wf_trial(:,1:nvalid)
#endif
    if(lanc_target_prediction.and.wf_trial_correlated)then
       select case(lanc_method)
       case("arpack")
          ! Input validation guarantees one requested state for ARPACK.
          allocate(wf_seed,source=wf_loc(:,1))
       case("primme")
          allocate(wf_seed_block,source=wf_loc)
       end select
    endif
  end subroutine sb_prepare_wf_guess


  subroutine sb_check_ritz_residuals()
    integer :: k,nloc
    logical :: empty_rank
    real(8) :: local_norm,total_norm
#ifdef _CMPLX
    complex(8),allocatable :: hv(:)
#else
    real(8),allocatable :: hv(:)
#endif
    call sb_check_Hv(empty_rank)
    if(empty_rank)then
       if(MpiMaster)write(LOGfile,*)"Ritz residual check skipped: serial fallback"
       return
    endif
    nloc=size(gs_vector,1)
    allocate(hv(nloc))
    do k=1,size(gs_vector,2)
       call spHtimesV_p(nloc,gs_vector(:,k),hv)
       local_norm=sum(abs(hv-gs_energy(k)*gs_vector(:,k))**2)
       total_norm=local_norm
#ifdef _MPI
       if(MpiStatus)call AllReduce_MPI(MpiComm,local_norm,total_norm)
#endif
       if(MpiMaster)write(LOGfile,'(A,I4,2ES24.15)') &
            "Ritz check: state, energy, residual ",k,gs_energy(k),sqrt(total_norm)
    enddo
  end subroutine sb_check_ritz_residuals









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
    logical                                        :: fMPI
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
          if(direct_H_lazy)then
             spHtimesV_p => spMatVec_lazy_main
          else
             spHtimesV_p => spMatVec_cache_main
          endif
#ifdef _MPI
          if(MpiStatus)then
             call sb_check_Hv(fMPI)
             if(fMPI)then
                write(LogFile,"(A)")"ARPACK: using SERIAL over MPI as MpiSize > N"
                if(direct_H_lazy)then
                   spHtimesV_p => spMatVec_lazy_main
                else
                   spHtimesV_p => spMatVec_cache_main
                endif
             else
                if(direct_H_lazy)then
                   spHtimesV_p => spMatVec_MPI_lazy_main
                else
                   spHtimesV_p => spMatVec_MPI_cache_main
                endif
             endif
          endif
#endif
          !
       end select
    endif
    !
  end subroutine sb_build_Hv





  subroutine sb_delete_Hv()
    spHtimesV_p => null()
    call spHsb%free()
    call free_superblock_setup_state()
    call sb_delete_dims()
  end subroutine sb_delete_Hv










  !##################################################################
  !##################################################################
  !           SAVE/LOAD MEASURE STATE
  !includes:
  ! * gs_energy
  ! * gs_vector
  ! * sb_states
  ! * sb_sector
  !
  !##################################################################
  !##################################################################  
  subroutine sb_save_measure_state()
   integer                                :: unit,Nstates,Nenergy,Nrow,Ncol,Nqn
   character(len=:),allocatable           :: file
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: gs_full
#else
    real(8),dimension(:,:),allocatable    :: gs_full
#endif
    !
    if(.not.allocated(sb_states))return
    if(.not.allocated(gs_vector))return
    if(.not.allocated(gs_energy))return
    !
    !Build parallel dims (which have been deleted by sb_delete_Hv)
    call sb_build_dims(quiet=.true.)
    Nstates=size(sb_states)
    Nenergy=size(gs_energy)
    Ncol=size(gs_vector,2)
    if(MpiStatus)then
       if(MpiMaster)then
          allocate(gs_full(Nstates,Ncol))
       else
          allocate(gs_full(1,Ncol))
       endif
       gs_full=zero
#ifdef _MPI
       call gather_vector_MPI(MpiComm,gs_vector,gs_full)
#endif
    else
       allocate(gs_full(Nstates,Ncol))
       gs_full = gs_vector
    endif
    !
    if(MpiMaster)then
       file = str(measure_file)//".restart"
       open(free_unit(unit),file=str(file))
       write(unit,*)"#NSSDMRG_MEASURE_STATE"
       write(unit,*)MpiSize
       write(unit,*)left%length,right%length,current_L
       Nqn=size(current_target_qn)
       write(unit,*)Nqn
       write(unit,*)current_target_qn
       write(unit,*)Nstates
       write(unit,*)sb_states
       call sb_sector%write(unit=unit)
       write(unit,*)Nenergy
       write(unit,*)gs_energy
       Nrow=size(gs_full,1)
       Ncol=size(gs_full,2)
       write(unit,*)Nrow,Ncol
       close(unit)
       !
       file = str(measure_file)//".gs.restart"
       open(free_unit(unit),file=str(file))
       write(unit,*)gs_full
       close(unit)
    endif
    call sb_delete_dims()
    deallocate(gs_full)
  end subroutine sb_save_measure_state


  subroutine sb_load_measure_state(found)
    logical,optional                      :: found
    logical                               :: file_exists
    integer                               :: unit,Ncores
    integer                               :: left_length,right_length,stored_current_L
    integer                               :: Nstates,Nenergy,Nrow,Ncol,Nqn
    character(len=64)                     :: magic
    character(len=:),allocatable          :: file
    character(len=:),allocatable          :: vector_file
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: gs_full
#else
    real(8),dimension(:,:),allocatable    :: gs_full
#endif
    !
    file        = str(measure_file)//".restart"
    vector_file = str(measure_file)//".gs.restart"
    inquire(file=str(file),exist=file_exists)
    if(.not.file_exists)then
       file        = str(measure_restart_file)//".restart"
       vector_file = str(measure_restart_file)//".gs.restart"
       inquire(file=str(file),exist=file_exists)
    endif
    if(present(found))found=file_exists
    if(.not.file_exists)then !it will stops in the calling routine
       if(MpiMaster)then
          write(LOGfile,*)"sb_load_measure_state: no metadata file at ",str(measure_file)//".restart"
          write(LOGfile,*)"sb_load_measure_state: no metadata file at ",str(measure_restart_file)//".restart"
       endif
       return
    endif
    !Metadata and vector must belong to the same checkpoint prefix.
    inquire(file=str(vector_file),exist=file_exists)
    if(.not.file_exists)then
       if(MpiMaster)write(LOGfile,*)"sb_load_measure_state: missing ground-state vector ",str(vector_file)
       error stop 1
    endif
    !
    !Read text:
    open(free_unit(unit),file=str(file))
    read(unit,*)magic
    if(str(magic)/="#NSSDMRG_MEASURE_STATE")stop "sb_load_measure_state error: unsupported file format"
    !
    !Read MpiSize_Ncores
    read(unit,*)Ncores
    if(Ncores/=MpiSize)then   
         if(MpiMaster)write(LOGfile,*) "sb_load_measure_state warning: MPI size mismatch"
    endif
    !
    !Read block dimensions and current_L
    read(unit,*)left_length,right_length,stored_current_L
    if(left_length/=left%length)stop "sb_load_measure_state error: left block length mismatch"
    if(right_length/=right%length)stop "sb_load_measure_state error: right block length mismatch"
    current_L = stored_current_L
    !
    !Read number of quantum numbers (Left) + QN target
    read(unit,*)Nqn
    if(allocated(current_target_qn))deallocate(current_target_qn)
    allocate(current_target_qn(Nqn))
    read(unit,*)current_target_qn
    !
    !Read SB_States
    read(unit,*)Nstates
    if(allocated(sb_states))deallocate(sb_states)
    allocate(sb_states(Nstates))
    read(unit,*)sb_states
    !
    !Read SB_Sector
    call sb_sector%read(unit=unit)
    read(unit,*)Nenergy
    !
    !Read GS_energy and GS_vector
    if(allocated(gs_energy))deallocate(gs_energy)
    allocate(gs_energy(Nenergy))
    read(unit,*)gs_energy
    read(unit,*)Nrow,Ncol
    close(unit)
    call sb_build_dims(quiet=.true.)
    if(allocated(gs_vector))deallocate(gs_vector)
    allocate(gs_vector(sb_vecDim_Hv(),Ncol))
    if(MpiStatus)then
       if(MpiMaster)then
          allocate(gs_full(Nrow,Ncol))
          open(free_unit(unit),file=str(vector_file))
          read(unit,*)gs_full
          close(unit)
       else
          allocate(gs_full(1,Ncol));gs_full=zero
       endif
#ifdef _MPI
       call scatter_vector_MPI(MpiComm,gs_full,gs_vector)
#endif
    else
       allocate(gs_full(Nrow,Ncol))
       open(free_unit(unit),file=str(vector_file))
       read(unit,*)gs_full
       close(unit)
       gs_vector = gs_full
    endif
    call sb_delete_dims()
    deallocate(gs_full)
  end subroutine sb_load_measure_state


  subroutine sb_save_measure_blocks()
    call left%save(suffix_dmrg('left',type='i')//".restart",gzip=.false.,&
         include_omatrices=.false.,file_prefix=str(measure_block_file))
    call right%save(suffix_dmrg('right',type='i')//".restart",gzip=.false.,&
         include_omatrices=.false.,file_prefix=str(measure_block_file))
  end subroutine sb_save_measure_blocks



  subroutine sb_load_measure_blocks()
   logical                      :: left_file,right_file
   character(len=:),allocatable :: left_suffix,right_suffix,prefix
   !
   left_suffix  = suffix_dmrg('left',type='i')//".restart"
   right_suffix = suffix_dmrg('right',type='i')//".restart"
   !
   !Inquire for left and right measure_block_file files: it normal or restart file
   prefix=str(measure_block_file)
   inquire(file=str(prefix)//str(left_suffix),exist=left_file)
   inquire(file=str(prefix)//str(right_suffix),exist=right_file)
   if(.not.(left_file.and.right_file))then
      prefix=str(measure_block_restart_file)
      inquire(file=str(prefix)//str(left_suffix),exist=left_file)
      inquire(file=str(prefix)//str(right_suffix),exist=right_file)
   endif
   !
   !IF files have been found: load, ELSE raise a WARNING
   if(left_file.and.right_file)then
      call left%load(str(left_suffix),file_prefix=str(prefix))
      call right%load(str(right_suffix),file_prefix=str(prefix))
   elseif(left_file.neqv.right_file)then
      if(MpiMaster)write(LOGfile,*)"Init_Measure_DMRG WARNING: incomplete measurement block restart pair found."
   endif
 end subroutine sb_load_measure_blocks










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
    real(8),allocatable :: requested(:),lower(:),upper(:)
    current_L         = left%length + right%length
    !L is the total number of sites in the two enlarged blocks.
    !Calculate each target as INT(L*rho)+offset.
    !INT discards the fractional part: INT(2.7)=2 and INT(-2.7)=-2.
    !For example, 0.33333333*6 is slightly below 2 and gives 1.
    !Keep the offset as supplied: spin QNs can be half-integer.
    requested = dble(int(target_density*current_L)) + target_offset
    current_target_QN=requested
    !In the first growth steps, an offset may ask for more particles or spin
    !than the small chain can contain. Temporarily limit the target to what fits.
    !Example: rho_up=0.5, offset_up=3.
    !  L=4: request 5 up electrons, but only 4 fit; use 4.
    !  L=6: request 6 up electrons; all 6 fit, so use the requested target.
    !This correction requires the same local QN limits on every site.
    !At L >= 2*Ldmrg, use the requested target without correcting it.
    !Its physical validity was already checked in init_dmrg for standard bases.
    if(target_uniform_bounds.and.current_L<2*Ldmrg)then
       lower=current_L*target_site_min
       upper=current_L*target_site_max
       current_target_QN=max(lower,min(upper,requested))
       if(any(current_target_QN/=requested).and.MpiMaster)then
          write(LOGfile,*)"Warm-up QN clipping at length:",current_L
          write(LOGfile,*)"Requested target:",requested
          write(LOGfile,*)"Effective target:",current_target_QN
       endif
    endif
    !A target inside these limits may still be unavailable in the saved blocks.
    !For example, DMRG truncation may have removed the states needed to build it.
    !sb_get_states checks that states exist for the target we actually use.
    !Do not change the target to match the QNs left after truncation:
    !that would change the physical sector we are calculating.
    !With different local limits, leave the target unchanged and let the SB
    !construction check whether it can be built.
  end subroutine sb_set_current_qn



END MODULE DMRG_SUPERBLOCK











