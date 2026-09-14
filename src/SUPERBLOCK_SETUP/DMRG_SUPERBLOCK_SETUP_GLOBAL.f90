MODULE DMRG_SUPERBLOCK_SETUP_GLOBAL
  USE DMRG_GLOBAL
  implicit none
  private

  integer,public                         :: tNso=0
  type(sparse_matrix),allocatable,public :: Hleft(:),Hright(:)
  type(sparse_matrix),allocatable,public :: A(:,:),B(:,:)
  type(tstates),allocatable,public       :: SBleft_states(:),SBright_states(:)
  type(tstates),allocatable,public       :: SBleft_maps(:),SBright_maps(:)
  type(sparse_matrix),public             :: Lazy_Hl,Lazy_Hr
  type(sparse_matrix),allocatable,public :: Lazy_Sl_n(:),Lazy_Sr_n(:)
  type(sparse_matrix),allocatable,public :: Lazy_Sl_p(:),Lazy_Sr_p(:)
  type(sparse_matrix),allocatable,public :: Lazy_Cr_n(:),Lazy_Cr_p(:)
  type(sparse_matrix),allocatable,public :: Lazy_CdgP_n(:),Lazy_CdgP_p(:)
  real(8),allocatable,public             :: Lazy_dql_n(:,:),Lazy_dql_p(:,:)
  real(8),allocatable,public             :: Lazy_dqr_n(:,:),Lazy_dqr_p(:,:)
  integer,allocatable,public             :: RowOffset(:,:),ColOffset(:,:),isb2jsb(:,:)
  integer,allocatable,public             :: IsHconjg(:,:)


  public :: setup_sector_filter_maps
  !
  public :: free_superblock_setup_state
  !
  public :: filter_left_operator
  public :: filter_right_operator
  public :: free_lazy_operators



  
contains





  !##################################################################
  !              FREE SUPERBLOCK SETUP STATE
  !    free the memory allocated during setup of the SB operators
  !##################################################################
  subroutine free_superblock_setup_state()
    call free_sparse_vector(Hleft)
    call free_sparse_vector(Hright)
    call free_sparse_matrix_array(A)
    call free_sparse_matrix_array(B)
    if(allocated(SBleft_states))deallocate(SBleft_states)
    if(allocated(SBright_states))deallocate(SBright_states)
    if(allocated(SBleft_maps))deallocate(SBleft_maps)
    if(allocated(SBright_maps))deallocate(SBright_maps)
    call free_lazy_operators()
    if(allocated(RowOffset))deallocate(RowOffset)
    if(allocated(ColOffset))deallocate(ColOffset)
    if(allocated(isb2jsb))deallocate(isb2jsb)
    if(allocated(IsHconjg))deallocate(IsHconjg)
    tNso=0
  end subroutine free_superblock_setup_state

  subroutine free_sparse_vector(vector)
    type(sparse_matrix),allocatable,intent(inout) :: vector(:)
    integer :: i
    if(.not.allocated(vector))return
    do i=1,size(vector)
       call vector(i)%free()
    enddo
    deallocate(vector)
  end subroutine free_sparse_vector



  subroutine free_sparse_matrix_array(matrix)
    type(sparse_matrix),allocatable,intent(inout) :: matrix(:,:)
    integer :: i,j
    if(.not.allocated(matrix))return
    do j=1,size(matrix,2)
       do i=1,size(matrix,1)
          call matrix(i,j)%free()
       enddo
    enddo
    deallocate(matrix)
  end subroutine free_sparse_matrix_array



  subroutine free_lazy_operators()
    call Lazy_Hl%free()
    call Lazy_Hr%free()
    call free_sparse_vector(Lazy_Sl_n)
    call free_sparse_vector(Lazy_Sr_n)
    call free_sparse_vector(Lazy_Sl_p)
    call free_sparse_vector(Lazy_Sr_p)
    call free_sparse_vector(Lazy_CdgP_n)
    call free_sparse_vector(Lazy_Cr_n)
    call free_sparse_vector(Lazy_CdgP_p)
    call free_sparse_vector(Lazy_Cr_p)
    if(allocated(Lazy_dql_n))deallocate(Lazy_dql_n)
    if(allocated(Lazy_dqr_n))deallocate(Lazy_dqr_n)
    if(allocated(Lazy_dql_p))deallocate(Lazy_dql_p)
    if(allocated(Lazy_dqr_p))deallocate(Lazy_dqr_p)
  end subroutine free_lazy_operators








  !##################################################################
  !               SETUP SECTOR FILTER MAPS
  !  you get the Left/Right states for each sector from the SB list. 
  !##################################################################
  subroutine setup_sector_filter_maps()
    integer :: isb,istate
    if(allocated(SBleft_maps))deallocate(SBleft_maps)
    if(allocated(SBright_maps))deallocate(SBright_maps)
    allocate(SBleft_maps(size(sb_sector)),SBright_maps(size(sb_sector)))
    do isb=1,size(sb_sector)
       allocate(SBleft_maps(isb)%states(left%Dim))
       allocate(SBright_maps(isb)%states(right%Dim))
       SBleft_maps(isb)%states=0
       SBright_maps(isb)%states=0
       do istate=1,size(SBleft_states(isb)%states)
          SBleft_maps(isb)%states(SBleft_states(isb)%states(istate))=istate
       enddo
       do istate=1,size(SBright_states(isb)%states)
          SBright_maps(isb)%states(SBright_states(isb)%states(istate))=istate
       enddo
    enddo
  end subroutine setup_sector_filter_maps

  function filter_left_operator(Op,irow,icol) result(Op_q)
    type(sparse_matrix),intent(in) :: Op
    integer,intent(in) :: irow,icol
    type(sparse_matrix) :: Op_q
    Op_q=sp_filter(Op,SBleft_states(irow)%states,SBleft_maps(icol)%states,&
         size(SBleft_states(icol)%states))
  end function filter_left_operator

  function filter_right_operator(Op,irow,icol) result(Op_q)
    type(sparse_matrix),intent(in) :: Op
    integer,intent(in) :: irow,icol
    type(sparse_matrix) :: Op_q
    Op_q=sp_filter(Op,SBright_states(irow)%states,SBright_maps(icol)%states,&
         size(SBright_states(icol)%states))
  end function filter_right_operator



END MODULE DMRG_SUPERBLOCK_SETUP_GLOBAL
