MODULE DMRG_RDM  
  USE DMRG_GLOBAL
  USE DMRG_SUPERBLOCK 
  implicit none
  private


  public :: sb_get_rdm
  public :: renormalize_block
  public :: wf_build_trial, wf_prepare_center, wf_reset

#ifdef _CMPLX
  complex(8),allocatable :: wf_old_psi(:),wf_ul(:,:),wf_ur(:,:),wf_center(:,:),wf_correlated(:,:,:,:)
#else
  real(8),allocatable :: wf_old_psi(:),wf_ul(:,:),wf_ur(:,:),wf_center(:,:),wf_correlated(:,:,:,:)
#endif
  integer,allocatable :: wf_old_states(:)
  real(8),allocatable :: wf_old_target(:)
  real(8),allocatable :: wf_qn_l(:,:),wf_qn_r(:,:)
  integer :: wf_old_nr=0

contains

  subroutine wf_reset()
    if(allocated(wf_old_psi))deallocate(wf_old_psi)
    if(allocated(wf_old_states))deallocate(wf_old_states)
    if(allocated(wf_old_target))deallocate(wf_old_target)
    if(allocated(wf_ul))deallocate(wf_ul)
    if(allocated(wf_ur))deallocate(wf_ur)
    if(allocated(wf_center))deallocate(wf_center)
    if(allocated(wf_correlated))deallocate(wf_correlated)
    if(allocated(wf_qn_l))deallocate(wf_qn_l)
    if(allocated(wf_qn_r))deallocate(wf_qn_r)
    if(allocated(wf_trial))deallocate(wf_trial)
    wf_trial_correlated=.false.
  end subroutine wf_reset



  !##################################################################
  !       GET REDUCED DENSITY MATRIX 
  !##################################################################
  subroutine sb_get_rdm()
    integer                               :: isb,igs,Ngs,wf_est_m
    integer                               :: Nleft,Nright,m_sb,Neig
    real(8)                               :: gs_weight
    real(8),dimension(:),allocatable      :: sb_qn,qn
    integer,dimension(:),allocatable      :: sb_map
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: v_state
    complex(8),dimension(:,:),allocatable :: rho
#else
    real(8),dimension(:,:),allocatable    :: v_state
    real(8),dimension(:,:),allocatable    :: rho
#endif
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: get RDM"
#endif
    !
    call sb_build_dims(quiet=.true.)
    !
    m_sb = size(sb_states)
    Neig = size(gs_vector,2)
    allocate(v_state(m_sb,Neig));v_state=zero
#ifdef _MPI
    if(MpiStatus)then
       call gather_vector_MPI(MpiComm,gs_vector,v_state)
    else
       v_state = gs_vector   
    endif
#else
    v_state = gs_vector    
#endif
    Ngs=get_gs_degeneracy()
    if((wf_predict_diag.or.lanc_gs_prediction).and.MpiMaster)then
       if(allocated(wf_old_psi))deallocate(wf_old_psi)
       if(allocated(wf_old_states))deallocate(wf_old_states)
       if(allocated(wf_old_target))deallocate(wf_old_target)
       wf_est_m=Mstates
       if(wf_est_m<=0)wf_est_m=max(left%dim,right%dim)
       if(int(m_sb,8)*int(wf_est_m,8)**2<=50000000_8)then
          allocate(wf_old_psi,source=v_state(:,1))
          allocate(wf_old_states,source=sb_states)
          allocate(wf_old_target,source=current_target_qn)
          wf_old_nr=right%dim
       else
          if(allocated(wf_correlated))deallocate(wf_correlated)
          ! The diagnostic uses a direct contraction. Keep it out of
          ! production-sized runs until a sectorwise contraction exists.
          if(allocated(wf_center))deallocate(wf_center)
          if(allocated(wf_qn_l))deallocate(wf_qn_l)
          if(allocated(wf_qn_r))deallocate(wf_qn_r)
          write(LOGfile,*)"WF trial skipped: projected-center work limit"
       endif
    endif
    gs_weight=1d0/dble(Ngs)
    !
    !
    t0=t_start()
    call rho_left%free()
    call rho_right%free()
    !
    if(MpiMaster)then
       call start_timer("Get \rho")
       do isb=1,size(sb_sector)
          sb_qn   = sb_sector%qn(index=isb)
          sb_map  = sb_sector%map(index=isb)
          Nleft   = size(left%sectors(1)%map(qn=sb_qn))
          Nright  = size(right%sectors(1)%map(qn=(current_target_qn - sb_qn)))
          if(Nleft*Nright==0)cycle
          !
          !LEFT:
          qn  = sb_qn
          rho = build_density_matrix(Nleft,Nright,v_state(:,1),sb_map,'left')
          rho = gs_weight*rho
          do igs=2,Ngs
             rho = rho + gs_weight*build_density_matrix(&
                  Nleft,Nright,v_state(:,igs),sb_map,'left')
          enddo
          call rho_left%append(rho,qn=qn,map=left%sectors(1)%map(qn=qn))
          !
          !RIGHT:
          qn  = current_target_qn-sb_qn
          rho = build_density_matrix(Nleft,Nright,v_state(:,1),sb_map,'right') 
          rho = gs_weight*rho
          do igs=2,Ngs
             rho = rho + gs_weight*build_density_matrix(&
                  Nleft,Nright,v_state(:,igs),sb_map,'right')
          enddo
          call rho_right%append(rho,qn=qn,map=right%sectors(1)%map(qn=qn))
          !
       enddo
       call stop_timer("Get \rho")
    endif
    t_rdm_get=t_stop()
    !

    if(MpiMaster)then
      call start_timer("Diag rho LEFT")
      t0=t_start()
      call rho_left%eigh(sort=.true.,reverse=.true.)
      t_rdm_diag=t_rdm_diag+t_stop()
      call stop_timer("Diag rho LEFT")
      if(allocated(rho_left_evals))deallocate(rho_left_evals)
      allocate(rho_left_evals, source=rho_left%evals())
      where(rho_left_evals<0d0)rho_left_evals=0d0
      !
      call start_timer("Diag rho RIGHT")
      t0=t_start()
      call rho_right%eigh(sort=.true.,reverse=.true.)
      t_rdm_diag=t_rdm_diag+t_stop()
      call stop_timer("Diag rho RIGHT")
      if(allocated(rho_right_evals))deallocate(rho_right_evals)
      allocate(rho_right_evals, source=rho_right%evals())
      where(rho_right_evals<0d0)rho_right_evals=0d0
    endif  
    !
    call sb_delete_dims()
    if(allocated(rho))deallocate(rho)
    if(allocated(sb_map))deallocate(sb_map)
    if(allocated(sb_qn))deallocate(sb_qn)
    if(allocated(qn))deallocate(qn)
  end subroutine sb_get_rdm









  !##################################################################
  !              BUILD REDUCED DENSITY MATRIX 
  !##################################################################
  function build_density_matrix(Nleft,Nright,psi,map,direction) result(rho)
    integer                               :: Nleft,Nright
    integer,dimension(nleft*nright)       :: map
    character(len=*)                      :: direction
    integer                               :: il,ir,i,j
#ifdef _CMPLX
    complex(8),dimension(:)               :: psi
    complex(8),dimension(:,:),allocatable :: rho
    complex(8),dimension(nleft,nright)    :: psi_tmp
#else
    real(8),dimension(:)                  :: psi
    real(8),dimension(:,:),allocatable    :: rho
    real(8),dimension(nleft,nright)       :: psi_tmp
#endif
    !
    !
    if(allocated(rho))deallocate(rho)
    !
    !These two give the same results
    !psi_tmp = transpose(reshape(psi(map), [nright,nleft]))
    do concurrent(il=1:Nleft,ir=1:Nright)
       i = map(ir + (il-1)*Nright)
       psi_tmp(il,ir) = psi(i)
    enddo
    !
#ifdef _CMPLX
    select case(to_lower(str(direction)))
    case ('left','l')
       allocate(rho(nleft,nleft));rho=zero
       rho  = matmul( psi_tmp,  conjg(transpose(psi_tmp)) )
    case ('right','r')
       allocate(rho(nright,nright));rho=zero
       rho  = transpose(matmul( conjg(transpose(psi_tmp)), psi_tmp  ))
    end select
    if(any(abs(rho-conjg(transpose(rho)))>1d-12))&
         stop "build_density_matrix error: rho not Hermitian"
#else
    select case(to_lower(str(direction)))
    case ('left','l')
       allocate(rho(nleft,nleft));rho=zero
       rho  = matmul(psi_tmp,  transpose(psi_tmp) )
    case ('right','r')
       allocate(rho(nright,nright));rho=zero
       rho  = matmul(transpose(psi_tmp), psi_tmp  )
    end select
    if(any(abs(rho-transpose(rho))>1d-12))&
         stop "build_density_matrix error: rho not Symmetric"
#endif
  end function build_density_matrix






  !##################################################################
  !   RENORMALIZE THE L/R BLOCKS USING RDM
  !##################################################################
  subroutine renormalize_block(label,mtr)
    character(len=*),intent(in)      :: label
    integer                          :: mtr
    integer                          :: m_left,m_right
    integer                          :: m_s,m_e
    integer                          :: m_max
    integer                          :: j_
    integer                          :: m_err(1)
    real(8)                          :: e_,err
    real(8),dimension(:),allocatable :: qn
    integer                          :: i,j,r,l,im,unit
    type(tbasis)                     :: left_basis,right_basis
    type(sparse_matrix)              :: trRho_left,trRho_right
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: renormalize block "//str(label)
#endif
    !
    m_left  = left%dim
    m_right = right%dim
    !

    if(allocated(qn))deallocate(qn)
    allocate(qn, mold=current_target_qn)
    !
    select case(to_lower(str(label(1:1))))
    case ("l","s")!!!!!      LEFT
       mtr = m_left
       if(left%length+right%length==2)return
       !
       if(MpiMaster)then
          !Build Truncated Density Matrices:
          call start_timer("Renormalize "//to_lower(str(label)));t0=t_start()
          m_max = min(m_left,size(rho_left_evals))
          if(Estates>0d0)then
             m_err = minloc(abs(1d0-cumulate(rho_left_evals)-Estates))
             m_s   = min(m_err(1),m_max)
             if(Mstates>0.and.m_s>Mstates)then
                if(MpiMaster)write(LOGfile,*)"renormalize_block WARNING: Mdmrg caps Edmrg left states from",m_s,"to",Mstates
                m_s = Mstates
             endif
          elseif(Mstates>0)then
             m_s   = min(Mstates,m_max)
          else
             stop "Mdmrg=Edmrg=0 can not fix a threshold for the RDM"
          endif
          !
          e_=rho_left_evals(m_s)
          j_=m_s
          do i=j_+1,size(rho_left_evals)
             err = abs(e_-rho_left_evals(i))/e_
             if(err<=deg_evals_threshold.and.m_s<m_max)m_s=m_s+1
          enddo
          !>truncation-rotation matrices:
          truncation_error_left  = 1d0 - sum(rho_left_evals(1:m_s))
          trRho_left             = rho_left%sparse(m_left,m_s)
          if((wf_predict_diag.or.lanc_gs_prediction).and.allocated(wf_old_psi))then
             if(allocated(wf_ul))deallocate(wf_ul)
             wf_ul=trRho_left%as_matrix()
          endif
          !>Store all the rotation/truncation matrices:
          if(block_umat_cache)call left%put_omat(str(left%length),trRho_left)
          if(save_umat)call left%save_omat(str(left%length),trRho_left,&
               suffix_dmrg('left')//".restart",append=(left%length/=2),gzip=.false.)
       endif
#ifdef _MPI
       if(MpiStatus)call Bcast_MPI(MpiComm,m_s)
       if(MpiStatus)call trRho_left%bcast()
#endif       
       !>Renormalize Blocks:
       call operators_renormalization(left,trRho_left)
       ! call left%renormalize(trRho_left)
       if(MpiMaster)then
          call stop_timer("Renormalize "//to_lower(str(label)))
          t_rdm_renorm=t_rdm_renorm+t_stop()
       endif
       !
       !
       !>Prepare output and update basis state
       do im=1,m_s
          if(MpiMaster)qn = rho_left%qn(m=im)
#ifdef _MPI
          if(MpiStatus)call Bcast_Mpi(MpiComm,qn)
#endif
          call left_basis%append( qn=qn )
       enddo
       call left%set_basis(basis=left_basis)
       !
       Mtr = m_s
       !
       !
#ifdef _DEBUG
       if(MpiMaster.AND.verbose>5)then
          unit=fopen("lambdas_L_"//str(left%length)//".dat")       
          do i=1,size(rho_left_evals)
             err = abs(e_-rho_left_evals(i))/e_
             write(unit,*)i,rho_left_evals(i),err,1d0-sum(rho_left_evals(1:i))
             if(i==m_s)write(unit,*)" "
          enddo
          close(unit)
       endif
#endif
       call left_basis%free()
       call trRho_left%free()
       !
    case ("r","e")!!!!!      RIGHT
       mtr  = m_right
       if(left%length+right%length==2)return
       !
       if(MpiMaster)then
          !Build Truncated Density Matrices:
          call start_timer("Renormalize "//to_lower(str(label)));t0=t_start()
          m_max = min(m_right,size(rho_right_evals))
          if(Estates>0d0)then
             m_err = minloc(abs(1d0-cumulate(rho_right_evals)-Estates))
             m_e   = min(m_err(1),m_max)
             if(Mstates>0.and.m_e>Mstates)then
                if(MpiMaster)write(LOGfile,*)"renormalize_block WARNING: Mdmrg caps Edmrg right states from",m_e,"to",Mstates
                m_e = Mstates
             endif
          elseif(Mstates>0)then
             m_e   = min(Mstates,m_max)
          else
             stop "Mdmrg=Edmrg=0 can not fix a threshold for the RDM"
          endif
          !
          e_=rho_right_evals(m_e)
          j_=m_e
          do i=j_+1,size(rho_right_evals)
             err = abs(e_-rho_right_evals(i))/e_
             if(err<=deg_evals_threshold.and.m_e<m_max)m_e=m_e+1
          enddo
          !>truncation-rotation matrices:
          truncation_error_right = 1d0 - sum(rho_right_evals(1:m_e))
          trRho_right            = rho_right%sparse(m_right,m_e)
          if((wf_predict_diag.or.lanc_gs_prediction).and.allocated(wf_old_psi))then
             if(allocated(wf_ur))deallocate(wf_ur)
             wf_ur=trRho_right%as_matrix()
          endif
          !>Store all the rotation/truncation matrices:
          if(block_umat_cache)call right%put_omat(str(right%length),trRho_right)
          if(save_umat)call right%save_omat(str(right%length),trRho_right,&
               suffix_dmrg('right')//".restart",append=(right%length/=2),gzip=.false.)
       endif
#ifdef _MPI
       if(MpiStatus)call Bcast_MPI(MpiComm,m_e)
       if(MpiStatus)call trRho_right%bcast()
#endif
       !>Renormalize Blocks:
       call operators_renormalization(right,trRho_right)
       ! call right%renormalize(trRho_right)
       if(MpiMaster)then
          call stop_timer("Renormalize "//to_lower(str(label)));
          t_rdm_renorm=t_rdm_renorm+t_stop()
       endif
       !
       !
       !>Prepare output and update basis state
       do im=1,m_e
          if(MpiMaster)qn = rho_right%qn(m=im)
#ifdef _MPI
          if(MpiStatus)call Bcast_Mpi(MpiComm,qn)
#endif
          call right_basis%append( qn=qn )
       enddo
       call right%set_basis(basis=right_basis)
       !
       Mtr = m_e
       !
#ifdef _DEBUG
       if(MpiMaster.AND.verbose>5)then
          unit = fopen("lambdas_R_"//str(left%length)//".dat")       
          do i=1,size(rho_right_evals)
             err = abs(e_-rho_right_evals(i))/e_
             write(unit,*)i,rho_right_evals(i),err,1d0-sum(rho_right_evals(1:i))
             if(i==m_e)write(unit,*)" "
          enddo
          close(unit)
       endif
#endif
       call right_basis%free()
       call trRho_right%free()
    case default;stop "renormalize block: wrong label, not in [left-sys|right-env]"
    end select
    !
    return
  end subroutine renormalize_block


  subroutine wf_prepare_center()
    type(tbasis) :: basis_l,basis_r
    integer :: k,i,j,a,b
    real(8) :: wf_tbegin,wf_tend
#ifdef _CMPLX
    complex(8) :: amplitude
#else
    real(8) :: amplitude
#endif
    if((.not.wf_predict_diag.and..not.lanc_gs_prediction).or..not.MpiMaster)return
    if(.not.allocated(wf_old_psi))return
    if(.not.allocated(wf_ul).or..not.allocated(wf_ur))return
    call cpu_time(wf_tbegin)
    call wf_predict_obc()
    call cpu_time(wf_tend)
    if(wf_predict_diag)write(LOGfile,*)"WF correlated construction CPU seconds:",wf_tend-wf_tbegin
    if(allocated(wf_center))deallocate(wf_center)
    allocate(wf_center(size(wf_ul,2),size(wf_ur,2)))
    wf_center=zero
    ! C = U_L^dagger Psi U_R^*. The right conjugation follows the
    ! right-density-matrix convention in build_density_matrix.
    do k=1,size(wf_old_states)
       i=(wf_old_states(k)-1)/wf_old_nr+1
       j=mod(wf_old_states(k)-1,wf_old_nr)+1
       do a=1,size(wf_center,1)
          do b=1,size(wf_center,2)
#ifdef _CMPLX
             amplitude=conjg(wf_ul(i,a))*conjg(wf_ur(j,b))
#else
             amplitude=wf_ul(i,a)*wf_ur(j,b)
#endif
             wf_center(a,b)=wf_center(a,b)+amplitude*wf_old_psi(k)
          enddo
       enddo
    enddo
    deallocate(wf_old_psi,wf_old_states,wf_ul,wf_ur)
    call left%get_basis(basis_l)
    call right%get_basis(basis_r)
    if(allocated(wf_qn_l))deallocate(wf_qn_l)
    if(allocated(wf_qn_r))deallocate(wf_qn_r)
    wf_qn_l=basis_l%dump()
    wf_qn_r=basis_r%dump()
    call basis_l%free();call basis_r%free()
  end subroutine wf_prepare_center


  ! Same-edge McCulloch contraction. All old virtual bases are those
  ! actually used to construct the current enlarged blocks.
  subroutine wf_predict_obc()
    integer :: ml,mr,dl,dr,nl,nr,k,i,j,a,b,ss,tt
    real(8),allocatable :: singular(:)
#ifdef _CMPLX
    complex(8),allocatable :: psi(:,:),f(:,:),g(:,:),u(:,:),vh(:,:),ci(:,:),tmp(:,:),fv(:,:),gv(:,:)
#else
    real(8),allocatable :: psi(:,:),f(:,:),g(:,:),u(:,:),vh(:,:),ci(:,:),tmp(:,:),fv(:,:),gv(:,:)
#endif
    if(allocated(wf_correlated))deallocate(wf_correlated)
    if(PBCdmrg.or..not.allocated(wf_center))return
    ml=size(wf_center,1);mr=size(wf_center,2)
    dl=size(wf_ul,1)/ml;dr=size(wf_ur,1)/mr
    if(dl/=dr) return
    nl=size(wf_ul,2);nr=size(wf_ur,2)
    allocate(psi(size(wf_ul,1),size(wf_ur,1)));psi=zero
    do k=1,size(wf_old_states)
       i=(wf_old_states(k)-1)/wf_old_nr+1
       j=mod(wf_old_states(k)-1,wf_old_nr)+1
       psi(i,j)=wf_old_psi(k)
    enddo
#ifdef _CMPLX
    f=matmul(conjg(transpose(wf_ul)),psi)
    g=matmul(psi,conjg(wf_ur))
#else
    f=matmul(transpose(wf_ul),psi)
    g=matmul(psi,wf_ur)
#endif
    allocate(u(ml,ml),vh(mr,mr),singular(min(ml,mr)),ci(mr,ml))
    call svd(wf_center,singular,u,vh)
    ci=zero
    do k=1,size(singular)
       if(singular(k)<=1d-10*maxval(singular))cycle
       do i=1,mr
          do j=1,ml
#ifdef _CMPLX
             ci(i,j)=ci(i,j)+conjg(vh(k,i))*conjg(u(j,k))/singular(k)
#else
             ci(i,j)=ci(i,j)+vh(k,i)*u(j,k)/singular(k)
#endif
          enddo
       enddo
    enddo
    allocate(wf_correlated(nl,dl,dr,nr),fv(nl,mr),gv(ml,nr))
    do ss=1,dl
       fv=f(:,1+(ss-1)*mr:ss*mr)
       tmp=matmul(fv,ci)
       do tt=1,dr
          do a=1,ml
             gv(a,:)=g(tt+(a-1)*dl,:)
          enddo
          wf_correlated(:,ss,tt,:)=matmul(tmp,gv)
       enddo
    enddo
  end subroutine wf_predict_obc

  subroutine wf_build_trial(link)
    character(len=1),intent(in) :: link
    type(tbasis) :: basis_l,basis_r,site_l,site_r
    integer :: dL,dR,mL,mR,k,il,ir,sl,sr,xl,xr,pL,pR,pSL,pSR,phase,npair
    real(8),allocatable :: ql(:),qr(:),qsl(:),qsr(:),delta(:)
    real(8) :: norm
    wf_trial_correlated=.false.
    if(allocated(wf_trial))deallocate(wf_trial)
    if((.not.wf_predict_diag.and..not.lanc_gs_prediction).or..not.MpiMaster)return
    if(.not.allocated(wf_center))return
    if(.not.allocated(wf_qn_l).or..not.allocated(wf_qn_r))return
    mL=size(wf_center,1); mR=size(wf_center,2)
    dL=dot(left%length)%dim; dR=dot(left%length+2)%dim
    if(dL/=dR.or.left%dim/=mL*dL.or.right%dim/=mR*dR)return
    call left%get_basis(basis_l)
    call right%get_basis(basis_r)
    call dot(left%length)%get_basis(site_l)
    call dot(left%length+2)%get_basis(site_r)
    allocate(delta,source=current_target_qn-wf_old_target)
    allocate(wf_trial(size(sb_states)));wf_trial=zero
    npair=0
    do sl=1,dL
       do sr=1,dR
          qsl=site_l%qn(sl); qsr=site_r%qn(sr)
          if(all(abs(qsl+qsr-delta)<1d-8))npair=npair+1
       enddo
    enddo
    if(npair==0)then
       call basis_l%free();call basis_r%free()
       call site_l%free();call site_r%free()
       deallocate(wf_trial)
       return
    endif
    do k=1,size(sb_states)
       xl=(sb_states(k)-1)/right%dim+1
       xr=mod(sb_states(k)-1,right%dim)+1
       select case(link)
       case('n')
          il=(xl-1)/dL+1; sl=mod(xl-1,dL)+1
          sr=(xr-1)/mR+1; ir=mod(xr-1,mR)+1
          ql=basis_l%qn(xl); qr=basis_r%qn(xr)
          qsl=site_l%qn(sl); qsr=site_r%qn(sr)
       case('p')
          sl=(xl-1)/mL+1; il=mod(xl-1,mL)+1
          ir=(xr-1)/dR+1; sr=mod(xr-1,dR)+1
          ql=basis_l%qn(xl); qr=basis_r%qn(xr)
          qsl=site_l%qn(sl); qsr=site_r%qn(sr)
       end select
       if(allocated(wf_correlated).and..not.PBCdmrg)then
          wf_trial(k)=wf_correlated(il,sl,sr,ir)
          wf_trial_correlated=.true.
          cycle
       endif
       if(.not.all(abs(qsl+qsr-delta)<1d-8))cycle
       if(any(abs(ql-qsl-wf_qn_l(il,:))>1d-8))error stop "WF trial: left tensor order/QN mismatch"
       if(any(abs(qr-qsr-wf_qn_r(ir,:))>1d-8))error stop "WF trial: right tensor order/QN mismatch"
       ! Fermion parity of a tensor product is additive modulo two.
       ! The old block parities follow by subtracting the site QNs.
       pL=mod(nint(sum(ql-qsl)),2); pR=mod(nint(sum(qr-qsr)),2)
       pSL=mod(nint(sum(qsl)),2); pSR=mod(nint(sum(qsr)),2)
       select case(link)
       case('n') ! L,R,sL,sR -> L,sL,sR,R
          phase=1-2*mod(pR*(pSL+pSR),2)
       case('p') ! L,R,sL,sR -> sL,L,R,sR
          phase=1-2*mod(pSL*(pL+pR),2)
       end select
       wf_trial(k)=phase*wf_center(il,ir)/sqrt(dble(npair))
    enddo
    norm=sqrt(sum(abs(wf_trial)**2))
    if(norm>1d-14)then
       if(wf_predict_diag)write(LOGfile,"(A,A1,A,I0,A,ES14.6)") &
            "WF trial link=",link," local pairs=",npair," raw norm=",norm
       wf_trial=wf_trial/norm
    else
       deallocate(wf_trial)
    endif
    call basis_l%free();call basis_r%free()
    call site_l%free();call site_r%free()
  end subroutine wf_build_trial






  !+------------------------------------------------------------------+
  !PURPOSE:  Rernmalize operators in self.operators using Urho:
  ! O' = Urho^+ . O. Urho
  ! Note: dq(O') = dq(O) because renormalization acts within the same QN
  ! sector subspaces. Updating the existing key preserves its dq metadata.
  !+------------------------------------------------------------------+
  subroutine operators_renormalization(self,Urho)
    type(block)                  :: self
    type(sparse_matrix)          :: Urho
    integer                      :: i,N,M  !N=self%dim,M=truncated dimension
    type(sparse_matrix)          :: Op,rOp
    character(len=:),allocatable :: key
    !
    !0. Master checks dimensions:
    N = Urho%Nrow
    M = Urho%Ncol
    if(MpiMaster.AND.(N/=self%dim)) stop "self.renormalize error: size(Umat,1) != self.dim"
    !
    do i=1,size(self%operators)
       !1. Master retrieve operator and check dimensions
       if(MpiMaster)Op  = self%operators%op(index=i)
#ifdef _MPI
       if(MpiStatus)call Op%bcast()
#endif
       if( any( [Op%Nrow,Op%Ncol] /= [N,N] ) ) &
            stop "self.renormalize error: shape(Op) != [N,N] N=size(Rho,1)"
       !2. All nodes rotate&truncate
       rOp = rotate_and_truncate(Op)
       !3. Master store renormalized operator
       if(MpiMaster)then
          key  = self%operators%key(index=i)
          call self%update_op(str(key),rOp)
       endif
    enddo
    self%dim = M
    !
    call Op%free()
    call rOp%free()
    !
  contains
    !
    !Udgr.O.U: [M,N].[N,N].[N,M]=[M,M]
    function rotate_and_truncate(Op) result(RotOp)
      type(sparse_matrix) :: Op
      type(sparse_matrix) :: RotOp
#ifdef _MPI
      if(MpiStatus)then
         RotOp = (Urho%dgr().pm.Op).pm.Urho
      else
         RotOp = matmul( matmul(Urho%dgr(),Op), Urho)
      endif
#else
      RotOp = matmul( matmul(Urho%dgr(),Op), Urho)
#endif
    end function rotate_and_truncate
    !
  end subroutine operators_renormalization


END MODULE DMRG_RDM
