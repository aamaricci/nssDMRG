MODULE DMRG_GF
  USE DMRG_GLOBAL
  USE DMRG_SUPERBLOCK
  USE DMRG_SUPERBLOCK_SETUP
  USE DMRG_MEASURE
  implicit none
  private

  public :: GF_DMRG
  public :: Green_Op_DMRG
  public :: Green_Obc_Mode_DMRG
  public :: Lanczos_Green_DMRG
  public :: Build_GF_Targets_DMRG
  public :: Clear_GF_Targets_DMRG

  interface GF_DMRG
     module procedure :: Green_Op_DMRG
  end interface GF_DMRG

contains

  subroutine Build_GF_Targets_DMRG(step_type)
    character(len=*),intent(in) :: step_type
    type(sparse_matrix)         :: Op,OpDgr
    real(8),dimension(:),allocatable :: dq
    integer                     :: imode
    character(len=64)           :: label
    !
    call Clear_GF_Targets_DMRG()
    if(.not.gf_target)return
    if(to_lower(str(step_type(1:1)))/='f')return
    if(PBCdmrg)stop "Build_GF_Targets_DMRG error: GF targeting v1 supports OBC fDMRG only"
    if(gf_target_nmodes<=0)stop "Build_GF_Targets_DMRG error: GF_TARGET_NMODES must be > 0"
#ifdef _MPI
    if(MpiStatus)stop "Build_GF_Targets_DMRG error: GF targeting is not MPI-enabled yet"
#endif
    !
    call setup_obc_b2g_map()
    call ensure_gf_target_umat()
    Op = get_gf_target_operator()
    OpDgr = hconjg(Op)
    !
    do imode=1,gf_target_nmodes
       if(gf_target_modes(imode)<1.OR.gf_target_modes(imode)>current_L)&
            stop "Build_GF_Targets_DMRG error: GF_TARGET_MODES index outside [1,L]"
       if(gf_target_weight_hole>0d0)then
          write(label,'(A,I0)')"hole_n",gf_target_modes(imode)
          dq = infer_operator_qn_shift(Op)
          call append_obc_mode_target(Op,gf_target_modes(imode),gf_target_weight_hole,str(label),dq)
       endif
       if(gf_target_weight_particle>0d0)then
          write(label,'(A,I0)')"particle_n",gf_target_modes(imode)
          dq = infer_operator_qn_shift(OpDgr)
          call append_obc_mode_target(OpDgr,gf_target_modes(imode),gf_target_weight_particle,str(label),dq)
       endif
    enddo
    call Op%free()
    call OpDgr%free()
  end subroutine Build_GF_Targets_DMRG


  subroutine Clear_GF_Targets_DMRG()
    integer :: i
    if(.not.allocated(dmrg_targets))return
    do i=1,size(dmrg_targets)
       if(allocated(dmrg_targets(i)%qn))deallocate(dmrg_targets(i)%qn)
       if(allocated(dmrg_targets(i)%sb_states))deallocate(dmrg_targets(i)%sb_states)
       if(allocated(dmrg_targets(i)%vector))deallocate(dmrg_targets(i)%vector)
       call dmrg_targets(i)%sb_sector%free()
    enddo
    deallocate(dmrg_targets)
  end subroutine Clear_GF_Targets_DMRG


  subroutine append_obc_mode_target(Op,mode,weight,label,dq)
    type(sparse_matrix),intent(in)       :: Op
    integer,intent(in)                   :: mode
    real(8),intent(in)                   :: weight
    character(len=*),intent(in)          :: label
    real(8),dimension(:),intent(in)      :: dq
    integer,dimension(:),allocatable     :: old_sb_states,old_Dls,old_Drs,old_Offset
    type(sectors_list)                   :: old_sb_sector
    real(8),dimension(:),allocatable     :: old_target_qn,new_target_qn
#ifdef _CMPLX
    complex(8),dimension(:),allocatable  :: phi,tmp
#else
    real(8),dimension(:),allocatable     :: phi,tmp
#endif
    real(8)                              :: coeff,kn,nrm,pi_
    integer                              :: pos
    !
    old_target_qn = current_target_qn
    old_sb_states = sb_states
    old_sb_sector = sb_sector
    call sb_build_dims(quiet=.true.)
    old_Dls       = Dls
    old_Drs       = Drs
    old_Offset    = Offset
    call sb_delete_dims()
    !
    new_target_qn = old_target_qn + dq
    current_target_qn = new_target_qn
    call sb_get_states()
    call sb_build_dims(quiet=.true.)
    allocate(phi(sum(Dls*Drs)));phi=zero
    pi_ = acos(-1d0)
    kn  = pi_*dble(mode)/dble(current_L+1)
    do pos=1,current_L
       coeff = sqrt(2d0/dble(current_L+1))*sin(kn*dble(pos))
       if(abs(coeff)<epsilon(1d0))cycle
       tmp = apply_op_between_qn_sectors(Op,pos,gs_vector(:,1),&
            old_sb_sector,old_target_qn,old_Dls,old_Drs,old_Offset,&
            sb_sector,current_target_qn,Dls,Drs,Offset)
       phi = phi + coeff*tmp
       if(allocated(tmp))deallocate(tmp)
    enddo
    nrm = sqrt(max(0d0,real_scalar_product(phi,phi)))
    if(nrm>max(lanc_tolerance,epsilon(1d0)))then
       phi = phi/nrm
       call append_target_state(label,weight,current_target_qn,sb_states,sb_sector,phi)
    elseif(MpiMaster)then
       write(LOGfile,*)"Build_GF_Targets_DMRG warning: skipping zero-norm target "//str(label)
    endif
    !
    call sb_delete_dims()
    if(allocated(sb_states))deallocate(sb_states)
    sb_states = old_sb_states
    sb_sector = old_sb_sector
    current_target_qn = old_target_qn
    call old_sb_sector%free()
    if(allocated(phi))deallocate(phi)
  end subroutine append_obc_mode_target


  subroutine append_target_state(label,weight,qn,states,sectors,vec)
    character(len=*),intent(in)       :: label
    real(8),intent(in)                :: weight
    real(8),dimension(:),intent(in)   :: qn
    integer,dimension(:),intent(in)   :: states
    type(sectors_list),intent(in)     :: sectors
#ifdef _CMPLX
    complex(8),dimension(:),intent(in) :: vec
#else
    real(8),dimension(:),intent(in)    :: vec
#endif
    type(dmrg_target_state),dimension(:),allocatable :: tmp
    integer :: n
    !
    n=0
    if(allocated(dmrg_targets))n=size(dmrg_targets)
    allocate(tmp(n+1))
    if(n>0)tmp(1:n)=dmrg_targets
    tmp(n+1)%label  = label
    tmp(n+1)%weight = weight
    tmp(n+1)%qn = qn
    tmp(n+1)%sb_states = states
    tmp(n+1)%sb_sector = sectors
    tmp(n+1)%vector = vec
    call move_alloc(tmp,dmrg_targets)
  end subroutine append_target_state


  function get_gf_target_operator() result(Op)
    type(sparse_matrix) :: Op
    character(len=:),allocatable :: key,keyn
    key  = trim(gf_target_op)
    keyn = trim(gf_target_op)//"_n"
    if(dot(1)%operators%has_key(key))then
       Op = dot(1)%operators%op(key)
    elseif(dot(1)%operators%has_key(keyn))then
       Op = dot(1)%operators%op(keyn)
    else
       stop "get_gf_target_operator error: GF_TARGET_OP not found in local site operator list"
    endif
  end function get_gf_target_operator


  subroutine setup_obc_b2g_map()
    integer :: i
    if(allocated(b2gMap))deallocate(b2gMap)
    allocate(b2gMap(current_L))
    b2gMap = (/(i,i=1,current_L)/)
  end subroutine setup_obc_b2g_map


  subroutine ensure_gf_target_umat()
    character(len=:),allocatable :: default_suffix,current_suffix
    if(.not.MpiMaster)return
    if(size(left%omatrices)<=1)then
       default_suffix = suffix_dmrg('left',type='i')//".restart"
       current_suffix = suffix_dmrg('left')//".restart"
       call left%load_umat(str(default_suffix),left%length)
       if(size(left%omatrices)<=1.and.str(current_suffix)/=str(default_suffix))&
            call left%load_umat(str(current_suffix),left%length)
    endif
    if(size(right%omatrices)<=1)then
       default_suffix = suffix_dmrg('right',type='i')//".restart"
       current_suffix = suffix_dmrg('right')//".restart"
       call right%load_umat(str(default_suffix),right%length)
       if(size(right%omatrices)<=1.and.str(current_suffix)/=str(default_suffix))&
            call right%load_umat(str(current_suffix),right%length)
    endif
  end subroutine ensure_gf_target_umat

  !##################################################################
  ! Compute G_i(z)=<psi0|O_i^dag (z+E0-H)^(-1) O_i|psi0>.
  ! The caller supplies z, e.g. z=omega+i*eta for the retarded GF.
  !##################################################################
  function Green_Op_DMRG(Op,pos,z,niter,qn_shift) result(Gz)
    type(sparse_matrix),intent(in)        :: Op
    integer,intent(in)                    :: pos
    complex(8),dimension(:),intent(in)    :: z
    integer,optional,intent(in)           :: niter
    real(8),dimension(:),optional,intent(in) :: qn_shift
    complex(8),dimension(size(z))         :: Gz
#ifdef _CMPLX
    complex(8),dimension(:),allocatable   :: phi
#else
    real(8),dimension(:),allocatable      :: phi
#endif
    integer,dimension(:),allocatable      :: old_sb_states
    integer,dimension(:),allocatable      :: old_Dls,old_Drs,old_Offset
    type(sectors_list)                    :: old_sb_sector
    real(8),dimension(:),allocatable      :: old_target_qn,target_qn,op_qn_shift
    integer                               :: niter_
    !
    call Init_Measure_DMRG(" GF")
    if(.not.allocated(gs_vector))stop "Green_Op_DMRG error: gs_vector is not allocated"
    if(.not.allocated(gs_energy))stop "Green_Op_DMRG error: gs_energy is not allocated"
    !
    old_target_qn = current_target_qn
    if(present(qn_shift))then
       op_qn_shift = qn_shift
    else
       op_qn_shift = infer_operator_qn_shift(Op)
    endif
    if(size(op_qn_shift)/=size(old_target_qn))&
         stop "Green_Op_DMRG error: operator QN shift has wrong dimension"
    target_qn = old_target_qn + op_qn_shift
    !
    if(qn_is_zero(op_qn_shift))then
       phi = Apply_Op_DMRG(Op,pos,gs_vector(:,1))
    else
#ifdef _MPI
       if(MpiStatus)stop "Green_Op_DMRG error: sector-changing GF seed is not MPI-enabled yet"
#endif
       old_sb_states = sb_states
       old_sb_sector = sb_sector
       old_Dls       = Dls
       old_Drs       = Drs
       old_Offset    = Offset
       !
       current_target_qn = target_qn
       call sb_get_states()
       call sb_build_dims(quiet=.true.)
       phi = apply_op_between_qn_sectors(Op,pos,gs_vector(:,1),&
            old_sb_sector,old_target_qn,old_Dls,old_Drs,old_Offset,&
            sb_sector,current_target_qn,Dls,Drs,Offset)
    endif
    !
    niter_ = lanc_ngfiter
    if(present(niter))niter_=niter
    !
    call sb_build_Hv()
    Gz = Lanczos_Green_DMRG(phi,z,gs_energy(1),niter_)
    call sb_delete_Hv()
    !
    if(.not.qn_is_zero(op_qn_shift))then
       if(allocated(sb_states))deallocate(sb_states)
       sb_states = old_sb_states
       sb_sector = old_sb_sector
       current_target_qn = old_target_qn
       call old_sb_sector%free()
    endif
    call End_Measure_DMRG()
    !
    if(allocated(phi))deallocate(phi)
    !
  end function Green_Op_DMRG


  function Green_Obc_Mode_DMRG(Op,mode,z,niter,dagger) result(Gz)
    type(sparse_matrix),intent(in)        :: Op
    integer,intent(in)                    :: mode
    complex(8),dimension(:),intent(in)    :: z
    integer,optional,intent(in)           :: niter
    logical,optional,intent(in)           :: dagger
    complex(8),dimension(size(z))         :: Gz
    type(sparse_matrix)                   :: Oseed
    integer,dimension(:),allocatable      :: old_sb_states
    type(sectors_list)                    :: old_sb_sector
    real(8),dimension(:),allocatable      :: old_target_qn,dq
    integer                               :: niter_
    logical                               :: dagger_
    !
    dagger_=.false.;if(present(dagger))dagger_=dagger
    call Init_Measure_DMRG(" GF OBC mode")
    if(.not.allocated(gs_vector))stop "Green_Obc_Mode_DMRG error: gs_vector is not allocated"
    if(.not.allocated(gs_energy))stop "Green_Obc_Mode_DMRG error: gs_energy is not allocated"
#ifdef _MPI
    if(MpiStatus)stop "Green_Obc_Mode_DMRG error: OBC mode GF is not MPI-enabled yet"
#endif
    if(PBCdmrg)stop "Green_Obc_Mode_DMRG error: OBC mode GF requires PBCdmrg=F"
    !
    call setup_obc_b2g_map()
    call ensure_gf_target_umat()
    if(dagger_)then
       Oseed = hconjg(Op)
    else
       Oseed = Op
    endif
    dq = infer_operator_qn_shift(Oseed)
    old_target_qn = current_target_qn
    old_sb_states = sb_states
    old_sb_sector = sb_sector
    call Clear_GF_Targets_DMRG()
    call append_obc_mode_target(Oseed,mode,1d0,"gf_obc_mode",dq)
    if(.not.allocated(dmrg_targets))then
       Gz=cmplx(0d0,0d0,kind=8)
    else
       current_target_qn = dmrg_targets(1)%qn
       if(allocated(sb_states))deallocate(sb_states)
       sb_states = dmrg_targets(1)%sb_states
       sb_sector = dmrg_targets(1)%sb_sector
       niter_ = lanc_ngfiter
       if(present(niter))niter_=niter
       call sb_build_Hv()
       Gz = Lanczos_Green_DMRG(dmrg_targets(1)%vector,z,gs_energy(1),niter_)
       call sb_delete_Hv()
    endif
    !
    if(allocated(sb_states))deallocate(sb_states)
    sb_states = old_sb_states
    sb_sector = old_sb_sector
    current_target_qn = old_target_qn
    call old_sb_sector%free()
    call Clear_GF_Targets_DMRG()
    call Oseed%free()
    call End_Measure_DMRG()
  end function Green_Obc_Mode_DMRG


  function Lanczos_Green_DMRG(seed,z,E0,niter) result(Gz)
#ifdef _CMPLX
    complex(8),dimension(:),intent(in)    :: seed
    complex(8),dimension(:),allocatable   :: q,qold,w
#else
    real(8),dimension(:),intent(in)       :: seed
    real(8),dimension(:),allocatable      :: q,qold,w
#endif
    complex(8),dimension(:),intent(in)    :: z
    real(8),intent(in)                    :: E0
    integer,intent(in)                    :: niter
    complex(8),dimension(size(z))         :: Gz
    real(8),dimension(:),allocatable      :: alpha,beta
    real(8)                               :: norm2,bnorm,tol
    integer                               :: Nloc,nmax,nlan,j
    !
    if(.not.associated(spHtimesV_p))stop "Lanczos_Green_DMRG error: H*v is not initialized"
    Nloc = size(seed)
    nmax = min(max(1,niter),Nloc)
    tol  = max(lanc_tolerance,epsilon(1d0))
    !
    allocate(alpha(nmax),beta(nmax))
    allocate(q(Nloc),qold(Nloc),w(Nloc))
    alpha=0d0
    beta =0d0
    qold =zero
    !
    norm2 = real_scalar_product(seed,seed)
    if(norm2<=tol)then
       Gz=cmplx(0d0,0d0,kind=8)
       call cleanup()
       return
    endif
    q = seed/sqrt(norm2)
    !
    nlan = 0
    do j=1,nmax
       call spHtimesV_p(Nloc,q,w)
       if(j>1)w = w - beta(j-1)*qold
       alpha(j) = real_scalar_product(q,w)
       w = w - alpha(j)*q
       bnorm = sqrt(max(0d0,real_scalar_product(w,w)))
       nlan = j
       if(bnorm<=tol)exit
       beta(j) = bnorm
       qold = q
       q = w/beta(j)
    enddo
    !
    Gz = norm2*continued_fraction(z,E0,alpha,beta,nlan)
    call cleanup()
    !
  contains
    subroutine cleanup()
      if(allocated(alpha))deallocate(alpha)
      if(allocated(beta))deallocate(beta)
      if(allocated(q))deallocate(q)
      if(allocated(qold))deallocate(qold)
      if(allocated(w))deallocate(w)
    end subroutine cleanup
  end function Lanczos_Green_DMRG


  function continued_fraction(z,E0,alpha,beta,nlan) result(cf)
    complex(8),dimension(:),intent(in) :: z
    real(8),intent(in)                 :: E0
    real(8),dimension(:),intent(in)    :: alpha,beta
    integer,intent(in)                 :: nlan
    complex(8),dimension(size(z))      :: cf
    complex(8)                         :: den
    integer                            :: iw,j
    !
    cf=cmplx(0d0,0d0,kind=8)
    if(nlan<1)return
    do iw=1,size(z)
       den = z(iw) + cmplx(E0-alpha(nlan),0d0,kind=8)
       do j=nlan-1,1,-1
          den = z(iw) + cmplx(E0-alpha(j),0d0,kind=8) - beta(j)*beta(j)/den
       enddo
       cf(iw) = 1d0/den
    enddo
  end function continued_fraction


  function apply_op_between_qn_sectors(Op,pos,v,from_sector,from_target,from_Dls,from_Drs,from_Offset,&
       to_sector,to_target,to_Dls,to_Drs,to_Offset) result(Ov)
    type(sparse_matrix),intent(in)       :: Op
    integer,intent(in)                   :: pos
    type(sectors_list),intent(in)        :: from_sector,to_sector
    real(8),dimension(:),intent(in)      :: from_target,to_target
    integer,dimension(:),intent(in)      :: from_Dls,from_Drs,from_Offset
    integer,dimension(:),intent(in)      :: to_Dls,to_Drs,to_Offset
#ifdef _CMPLX
    complex(8),dimension(:),intent(in)   :: v
    complex(8),dimension(:),allocatable  :: Ov
    complex(8)                           :: val
#else
    real(8),dimension(:),intent(in)      :: v
    real(8),dimension(:),allocatable     :: Ov
    real(8)                              :: val
#endif
    type(sparse_matrix)                  :: Oi,Ok
    real(8),dimension(:),allocatable     :: qn_to,qn_from,qr_to,qr_from
    integer,dimension(:),allocatable     :: rows,cols
    character(len=1)                     :: label
    integer                              :: L,R,N,kt,kf,ir,il,jcol,jl
    integer                              :: i,j
    !
    L = left%length
    R = right%length
    N = L+R
    if(pos<1.OR.pos>N)stop "apply_op_between_qn_sectors error: Pos not in [1,Ldmrg]"
    label='l'; if(pos>L)label='r'
    !
    Oi = Build_Op_DMRG(Op,pos)
    Oi = Advance_Op_DMRG(Oi,pos)
    allocate(Ov(sum(to_Dls*to_Drs)))
    Ov=zero
    !
    do kt=1,size(to_sector)
       qn_to = to_sector%qn(index=kt)
       qr_to = to_target - qn_to
       do kf=1,size(from_sector)
          qn_from = from_sector%qn(index=kf)
          qr_from = from_target - qn_from
          select case(label)
          case("l")
             if(.not.qn_equal(qr_to,qr_from))cycle
             rows = left%sectors(1)%map(qn=qn_to)
             cols = left%sectors(1)%map(qn=qn_from)
             Ok = sp_filter(Oi,rows,cols)
             do ir=1,to_Drs(kt)
                do il=1,to_Dls(kt)
                   i = ir + (il-1)*to_Drs(kt) + to_Offset(kt)
                   do jcol=1,Ok%row(il)%Size
                      val = Ok%row(il)%vals(jcol)
                      jl  = Ok%row(il)%cols(jcol)
                      j   = ir + (jl-1)*from_Drs(kf) + from_Offset(kf)
                      Ov(i) = Ov(i) + val*v(j)
                   enddo
                enddo
             enddo
          case("r")
             if(.not.qn_equal(qn_to,qn_from))cycle
             rows = right%sectors(1)%map(qn=qr_to)
             cols = right%sectors(1)%map(qn=qr_from)
             Ok = sp_filter(Oi,rows,cols)
             do il=1,to_Dls(kt)
                do ir=1,to_Drs(kt)
                   i = ir + (il-1)*to_Drs(kt) + to_Offset(kt)
                   do jcol=1,Ok%row(ir)%Size
                      val = Ok%row(ir)%vals(jcol)
                      jl  = Ok%row(ir)%cols(jcol)
                      j   = jl + (il-1)*from_Drs(kf) + from_Offset(kf)
                      Ov(i) = Ov(i) + val*v(j)
                   enddo
                enddo
             enddo
          end select
          call Ok%free()
       enddo
    enddo
    !
    call Oi%free()
  end function apply_op_between_qn_sectors


  function infer_operator_qn_shift(Op) result(dq)
    type(sparse_matrix),intent(in)    :: Op
    real(8),dimension(:),allocatable  :: dq
    real(8),dimension(:),allocatable  :: qi,qj,dqij
    integer                           :: i,j,jj,iq
    logical                           :: found
    !
    if(size(dot)<1)stop "infer_operator_qn_shift error: dot is not allocated"
    allocate(dq(size(current_target_qn)))
    dq=0d0
    found=.false.
    do i=1,Op%Nrow
       qi = local_state_qn(i)
       do jj=1,Op%row(i)%Size
          if(Op%row(i)%vals(jj)==zero)cycle
          j = Op%row(i)%cols(jj)
          qj = local_state_qn(j)
          dqij = qi-qj
          if(.not.found)then
             dq=dqij
             found=.true.
          elseif(.not.qn_equal(dq,dqij))then
             stop "infer_operator_qn_shift error: operator mixes QN sectors; pass qn_shift or decompose Op"
          endif
       enddo
    enddo
    if(.not.found)dq=0d0
  contains
    function local_state_qn(state) result(qn)
      integer,intent(in)                 :: state
      real(8),dimension(:),allocatable   :: qn
      integer,dimension(:),allocatable   :: states
      integer                            :: isec
      do isec=1,size(dot(1)%sectors(1))
         states = dot(1)%sectors(1)%map(index=isec)
         if(any(states==state))then
            qn = dot(1)%sectors(1)%qn(index=isec)
            return
         endif
      enddo
      stop "infer_operator_qn_shift error: local state not found in site sectors"
    end function local_state_qn
  end function infer_operator_qn_shift


  function qn_equal(a,b) result(equal)
    real(8),dimension(:),intent(in) :: a,b
    logical                         :: equal
    equal = size(a)==size(b)
    if(equal)equal = all(abs(a-b)<1d-10)
  end function qn_equal


  function qn_is_zero(a) result(is_zero)
    real(8),dimension(:),intent(in) :: a
    logical                         :: is_zero
    is_zero = all(abs(a)<1d-10)
  end function qn_is_zero


#ifdef _CMPLX
  function real_scalar_product(a,b) result(dot)
    complex(8),dimension(:),intent(in) :: a,b
    complex(8)                         :: local_dot,global_dot
    real(8)                            :: dot
    local_dot = dot_product(a,b)
    global_dot = local_dot
#ifdef _MPI
    if(MpiStatus)call AllReduce_MPI(MpiComm,local_dot,global_dot)
#endif
    dot = real(global_dot,8)
  end function real_scalar_product
#else
  function real_scalar_product(a,b) result(dot)
    real(8),dimension(:),intent(in) :: a,b
    real(8)                         :: local_dot,dot
    local_dot = dot_product(a,b)
    dot = local_dot
#ifdef _MPI
    if(MpiStatus)call AllReduce_MPI(MpiComm,local_dot,dot)
#endif
  end function real_scalar_product
#endif

END MODULE DMRG_GF
