MODULE DMRG_SUPERBLOCK_SETUP_SPARSE
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private

  public :: Setup_SuperBlock_Sparse
  !
contains





  !##################################################################
  !              SETUP THE SUPERBLOCK HAMILTONIAN
  !                      SPARSE MODE
  !    H^SB = H^L x 1^R  + 1^L x H^R + H^LR
  !    H^LR = sum_p A_p x B_p
  ! 
  ! * sparse: get the sparse global SB Hamiltonian spHsb
  !##################################################################
  !POSSIBLY INCLUDE MPI HERE... this is probably a less efficient version
  !note that only the SB Hamiltonian needs to be constructed in parallel form
  !the other operators (small) are stored by each cpu.
  !In principle one could store any sparse matrix in parallel and build H^SB as MPI too.
  subroutine Setup_SuperBlock_Sparse()
    integer                      :: m_left,m_right
    character(len=:),allocatable :: type
    type(sparse_matrix)          :: H2
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Setup SB Sparse"
#endif
    !
    if(MpiMaster)call start_timer("Setup SB Sparse")
    t0=t_start()
    !
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Sparse ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Sparse ERROR: Missing right.H operator in the list"
    !
    type=str(left%type())
    if(type/=str(right%type()))&
         stop "Setup_SuperBlock_Sparse ERROR: left.Type != right.Type"
    !
    m_left = left%dim
    m_right= right%dim
    !
    select case(to_lower(type(1:1)))
    case default;stop "Setup_SuperBlock_Sparse ERROR: wrong left/right.Type"
    case ("s")
       H2 = connect_spin_blocks(left,right,sb_states,link="n")
       if(PBCdmrg)H2 = H2 + connect_spin_blocks(left,right,sb_states,link="p")
    case ("f","e")
       H2 = connect_fermion_blocks(left,right,sb_states,link="n")
       if(PBCdmrg)H2 = H2 + connect_fermion_blocks(left,right,sb_states,link="p")
    end select
    !
    spHsb= H2 & 
         + sp_kron(left%operators%op("H"),id(m_right),sb_states) &
         + sp_kron(id(m_left),right%operators%op("H"),sb_states) 
    !
    if(MpiMaster)call stop_timer("Setup SB Sparse")
    t_setup_sb_sparse=t_stop()
    t_connect_blocks=0d0
    !
  end subroutine Setup_SuperBlock_Sparse


  

END MODULE DMRG_SUPERBLOCK_SETUP_SPARSE









