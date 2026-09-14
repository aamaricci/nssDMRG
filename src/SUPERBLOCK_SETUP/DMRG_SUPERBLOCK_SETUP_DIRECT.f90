MODULE DMRG_SUPERBLOCK_SETUP_DIRECT
  USE DMRG_GLOBAL
  USE DMRG_SUPERBLOCK_SETUP_DIRECT_SPIN
  USE DMRG_SUPERBLOCK_SETUP_DIRECT_FERMION
  implicit none
  private


  public :: Setup_SuperBlock_Direct


contains


  !##################################################################
  !         SETUP THE SUPERBLOCK HAMILTONIAN using DIRECT method
  !    i.e. apply H to wavefunction on the fly.
  !    memory is used wither in Lazy mode (block operators are retrieved at 
  !    each step) or in Cache mode (block operators are stored in memory).
  !##################################################################
  subroutine Setup_SuperBlock_Direct()
    character(len=:),allocatable :: type
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Setup SB Direct"
#endif
    !Some checks:
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
    type=str(left%type())
    if(type/=str(right%type()))&
         stop "Setup_SuperBlock_Direct ERROR: left.Type != right.Type"

    !         
    t0=t_start()
    select case(to_lower(type(1:1)))
    case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.type"//str(type(1:1))
    case ("s")    ;call Setup_SuperBlock_Spin_Direct()
    case ("f","e");call Setup_SuperBlock_Fermion_Direct()
    end select
    t_setup_sb_direct=t_stop()
    !
  end subroutine Setup_SuperBlock_Direct




END MODULE DMRG_SUPERBLOCK_SETUP_DIRECT
