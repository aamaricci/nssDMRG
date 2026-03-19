MODULE DMRG_RDM  
  USE DMRG_GLOBAL
  USE DMRG_SUPERBLOCK 
  implicit none
  private


  public :: sb_get_rdm
  public :: renormalize_block

contains



  !##################################################################
  !       GET REDUCED DENSITY MATRIX 
  !##################################################################
  subroutine sb_get_rdm()
    integer                               :: isb
    integer                               :: Nleft,Nright,m_sb,Neig
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
    call sb_build_dims()
    !
#ifdef _MPI
    if(MpiStatus)then
       m_sb = size(sb_states)
       Neig = size(gs_vector,2)
       allocate(v_state(m_sb,Neig));v_state=zero
       call gather_vector_MPI(MpiComm,gs_vector,v_state)
    endif
#endif
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
          qn  = sb_qn
#ifdef _MPI
          if(MpiStatus)then
             rho = build_density_matrix(Nleft,Nright,v_state(:,1),sb_map,'left')
          else
             rho = build_density_matrix(Nleft,Nright,gs_vector(:,1),sb_map,'left')
          endif
#else
          rho = build_density_matrix(Nleft,Nright,gs_vector(:,1),sb_map,'left')
#endif
          call rho_left%append(rho,qn=qn,map=left%sectors(1)%map(qn=qn))
          !
          !
          qn  = current_target_qn-sb_qn
#ifdef _MPI
          if(MpiStatus)then
             rho = build_density_matrix(Nleft,Nright,v_state(:,1),sb_map,'right') 
          else
             rho = build_density_matrix(Nleft,Nright,gs_vector(:,1),sb_map,'right') 
          endif
#else
          rho = build_density_matrix(Nleft,Nright,gs_vector(:,1),sb_map,'right') 
#endif
          call rho_right%append(rho,qn=qn,map=right%sectors(1)%map(qn=qn))
          !
       enddo
       call stop_timer("Get \rho")
    endif
    t_rdm_get=t_stop()
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
    if(any(abs(rho-conjg(transpose(rho)))/=0d0))&
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
          call start_timer("Diag rho "//to_lower(str(label)));
          t0=t_start()
          call rho_left%eigh(sort=.true.,reverse=.true.)
          call stop_timer("Diag rho "//to_lower(str(label)));t_rdm_diag=t_rdm_diag+t_stop()
          !
          if(allocated(rho_left_evals))deallocate(rho_left_evals)
          allocate(rho_left_evals, source=rho_left%evals())
          !
          !Build Truncated Density Matrices:
          call start_timer("Renormalize "//to_lower(str(label)));t0=t_start()
          if(Mstates/=0)then
             m_s   = min(Mstates,m_left,size(rho_left_evals))
          elseif(Estates/=0d0)then
             m_err = minloc(abs(1d0-cumulate(rho_left_evals)-Estates))
             m_s   = m_err(1)
          else
             stop "Mdmrg=Edmrg=0 can not fix a threshold for the RDM"
          endif
          !
          e_=rho_left_evals(m_s)
          j_=m_s
          do i=j_+1,size(rho_left_evals)
             err = abs(e_-rho_left_evals(i))/e_
             if(err<=deg_evals_threshold)m_s=m_s+1
          enddo
          !>truncation-rotation matrices:
          truncation_error_left  = 1d0 - sum(rho_left_evals(1:m_s))
          trRho_left             = rho_left%sparse(m_left,m_s)
          !>Store all the rotation/truncation matrices:
          if(save_umat)then
             select case(store_umat)
             case('b');call left%put_omat(str(left%length),trRho_left,'')
             case('f');call left%write_omat(str(left%length),trRho_left,'',suffix_dmrg('left')//".restart")
             end select
          endif
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
       call rho_left%free()
       !
    case ("r","e")!!!!!      RIGHT
       mtr  = m_right
       if(left%length+right%length==2)return
       !
       if(MpiMaster)then
          call start_timer("Diag rho "//to_lower(str(label)));t0=t_start()
          call rho_right%eigh(sort=.true.,reverse=.true.)
          call stop_timer("Diag rho "//to_lower(str(label)));t_rdm_diag=t_rdm_diag+t_stop()
          !
          if(allocated(rho_right_evals))deallocate(rho_right_evals)
          allocate(rho_right_evals, source=rho_right%evals())
          !
          !Build Truncated Density Matrices:
          call start_timer("Renormalize "//to_lower(str(label)));t0=t_start()
          if(Mstates/=0)then
             m_e   = min(Mstates,m_right,size(rho_right_evals))
          elseif(Estates/=0d0)then
             m_err = minloc(abs(1d0-cumulate(rho_right_evals)-Estates))
             m_e   = m_err(1)
          else
             stop "Mdmrg=Edmrg=0 can not fix a threshold for the RDM"
          endif
          !
          e_=rho_right_evals(m_e)
          j_=m_e
          do i=j_+1,size(rho_right_evals)
             err = abs(e_-rho_right_evals(i))/e_
             if(err<=deg_evals_threshold)m_e=m_e+1
          enddo
          !>truncation-rotation matrices:
          truncation_error_right = 1d0 - sum(rho_right_evals(1:m_e))
          trRho_right            = rho_right%sparse(m_right,m_e)
          !>Store all the rotation/truncation matrices:
          if(save_umat)then
             select case(store_umat)
             case('b');call right%put_omat(str(right%length),trRho_right,'')
             case('f');call right%write_omat(str(right%length),trRho_right,'',suffix_dmrg('right')//".restart")
             end select
          endif
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
       call rho_right%free()
    case default;stop "renormalize block: wrong label, not in [left-sys|right-env]"
    end select
    !
    return
  end subroutine renormalize_block






  !+------------------------------------------------------------------+
  !PURPOSE:  
  !+------------------------------------------------------------------+
  subroutine operators_renormalization(self,Urho)
    type(block)                  :: self
    type(sparse_matrix)          :: Urho
    integer                      :: i,N,M  !N=self%dim,M=truncated dimension
    type(sparse_matrix)          :: Op,rOp
    character(len=:),allocatable :: key,type
    !
    !0. Master checks dimensions:
    N = Urho%Nrow
    M = Urho%Ncol
    if(MpiMaster.AND.(N/=self%dim)) stop "self.renormalize error: size(Umat,1) != self.dim"
    !
    do i=1,size(self%operators)
       !1. Master retrieve operator and check dimensions
       if(MpiMaster)Op  = self%operators%op(index=i)
       if(MpiStatus)call Op%bcast()
       if( any( [Op%Nrow,Op%Ncol] /= [N,N] ) ) &
            stop "self.renormalize error: shape(Op) != [N,N] N=size(Rho,1)"
       !2. All nodes rotate&truncate
       rOp = rotate_and_truncate(Op)
       !4. Master store renormalized operator
       if(MpiMaster)then
          key  = self%operators%key(index=i)
          type = self%operators%type(index=i)
          call self%put_op(str(key),rOp, type)
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










!     call start_timer("Diag Rho")
!     call rho_left%eigh(sort=.true.,reverse=.true.)
!     call rho_right%eigh(sort=.true.,reverse=.true.)
!     rho_left_evals  = rho_left%evals()
!     rho_right_evals = rho_right%evals()
!     call stop_timer()
!     !
!     call start_timer("Renormalize Blocks + Setup Basis")
!     !Build Truncated Density Matrices:
!     if(Mstates/=0)then
!        m_s = min(Mstates,m_left,size(rho_left_evals))
!        m_e = min(Mstates,m_right,size(rho_right_evals))       
!     elseif(Estates/=0d0)then
!        m_err = minloc(abs(1d0-cumulate(rho_left_evals)-Estates))
!        m_s   = m_err(1)
!        m_err = minloc(abs(1d0-cumulate(rho_right_evals)-Estates))
!        m_e   = m_err(1)
!     else
!        stop "Mdmrg=Edmrg=0 can not fix a threshold for the RDM"
!     endif
!     !
!     e_=rho_left_evals(m_s)
!     j_=m_s
!     do i=j_+1,size(rho_left_evals)
!        err=abs(rho_left_evals(i)-e_)/e_
!        if(err<=1d-1)m_s=m_s+1
!     enddo
!     e_=rho_right_evals(m_e)
!     j_=m_e
!     do i=j_+1,size(rho_right_evals)          
!        err=abs(rho_right_evals(i)-e_)/e_
!        if(err<=1d-1)m_e=m_e+1
!     enddo
!     !
!     truncation_error_left  = 1d0 - sum(rho_left_evals(1:m_s))
!     truncation_error_right = 1d0 - sum(rho_right_evals(1:m_e))
!     !
!     !>truncation-rotation matrices:
!     trRho_left  = rho_left%sparse(m_left,m_s)
!     trRho_right = rho_right%sparse(m_right,m_e)
!     !
!     !
!     !>Store all the rotation/truncation matrices:
!     call left%put_omat(str(left%length),trRho_left)
!     call right%put_omat(str(right%length),trRho_right)
!     !
!     !>Renormalize Blocks:
!     call left%renormalize(as_matrix(trRho_left))
!     call right%renormalize(as_matrix(trRho_right))
!     !
!     !>Prepare output and update basis state
!     do im=1,m_s
!        call left_basis%append( qn=rho_left%qn(m=im) )
!     enddo
!     do im=1,m_e
!        call right_basis%append( qn=rho_right%qn(m=im) )
!     enddo
!     call left%set_basis(basis=left_basis)
!     call right%set_basis(basis=right_basis)
!     !
! #ifdef _DEBUG
!     unit     = fopen("lambdas_L_"//str(left%length)//".dat")       
!     do i=1,m_s
!        write(unit,*)i,rho_left_evals(i),floor(log10(abs(rho_left_evals(i))))
!     enddo
!     write(unit,*)" "
!     do i=m_s+1,size(rho_left_evals)
!        write(unit,*)i,rho_left_evals(i),floor(log10(abs(rho_left_evals(i))))
!     enddo
!     close(unit)
!     !
!     unit     = fopen("lambdas_R_"//str(left%length)//".dat")       
!     do i=1,m_e
!        write(unit,*)i,rho_right_evals(i),floor(log10(abs(rho_right_evals(i))))
!     enddo
!     write(unit,*)" "
!     do i=m_s+1,size(rho_right_evals)
!        write(unit,*)i,rho_right_evals(i),floor(log10(abs(rho_right_evals(i))))
!     enddo
!     close(unit)
!     !
!     call trRho_left%show(file="TrRho_L_"//str(left%length)//".dat")
!     call trRho_right%show(file="TrRho_R_"//str(right%length)//".dat")
!     !
!     call rho_left%show(file="NewBasis_L_"//str(left%length)//".dat")
!     call rho_right%show(file="NewBasis_R_"//str(right%length)//".dat")
! #endif
!     !
!     call stop_timer()
!     write(LOGfile,"(A,2ES24.15)")"Truncation Errors                    :",&
!          truncation_error_left,truncation_error_right
!     call left_basis%free()
!     call right_basis%free()
!     call trRho_left%free()
!     call trRho_right%free()
!     call rho_left%free()
!     call rho_right%free()
