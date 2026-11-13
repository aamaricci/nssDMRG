MODULE DMRG_GEOMETRY
  implicit none
  private

  ! Return the map from block positions to iDMRG growth labels.
  ! For OBC the order is the identity. For PBC, block_length must be odd
  ! and growth starts at the central position, alternating right and left:
  ! for block_length=7 the returned map is [7,5,3,1,2,4,6].
  ! Invalid lengths terminate with error stop. The routine is serial and
  ! has no MPI side effects.
  public :: Build_Block_Growth_Map

  ! Return the map from iDMRG growth labels to physical positions for a
  ! chain containing an even number of sites. Consecutive entries form the
  ! left/right site pair consumed at each symmetric growth step. In PBC,
  ! Nsites/2 must be odd. Invalid sizes terminate with error stop.
  public :: Build_DMRG_Growth_Map

contains


  function Build_Block_Growth_Map(block_length,pbc) result(block_to_growth)
    integer,intent(in)               :: block_length
    logical,intent(in)               :: pbc
    integer,dimension(:),allocatable :: block_to_growth
    integer                          :: center,i
    !
    ! Validate the block geometry.
    if(block_length<=0)&
         error stop 'Build_Block_Growth_Map: block_length must be positive'
    if(pbc.and.mod(block_length,2)==0)&
         error stop 'Build_Block_Growth_Map: PBC block_length must be odd'
    !
    ! Construct the block-position to growth-label map.
    allocate(block_to_growth(block_length))
    if(.not.pbc)then
       block_to_growth=[(i,i=1,block_length)]
       return
    endif
    center=(block_length+1)/2
    block_to_growth(center)=1
    do i=1,(block_length-1)/2
       block_to_growth(center+i)=2*i
       block_to_growth(center-i)=2*i+1
    enddo
    !
  end function Build_Block_Growth_Map


  function Build_DMRG_Growth_Map(Nsites,pbc) result(growth_to_physical)
    integer,intent(in)               :: Nsites
    logical,intent(in)               :: pbc
    integer,dimension(:),allocatable :: growth_to_physical
    integer,dimension(:),allocatable :: block_to_growth
    integer                          :: block_length,growth,position
    !
    ! Validate the complete-chain geometry.
    if(Nsites<=0.or.mod(Nsites,2)/=0)&
         error stop 'Build_DMRG_Growth_Map: Nsites must be positive and even'
    block_length=Nsites/2
    if(pbc.and.mod(block_length,2)==0)&
         error stop 'Build_DMRG_Growth_Map: PBC Nsites/2 must be odd'
    !
    ! Invert the block map and pair every left position with its reflection
    ! in the right block.
    block_to_growth=Build_Block_Growth_Map(block_length,pbc)
    allocate(growth_to_physical(Nsites))
    do position=1,block_length
       growth=block_to_growth(position)
       growth_to_physical(2*growth-1)=position
       growth_to_physical(2*growth)=Nsites+1-position
    enddo
    !
  end function Build_DMRG_Growth_Map


END MODULE DMRG_GEOMETRY




!##################################################################
!##################################################################
!                          TEST
!##################################################################
!##################################################################
#ifdef _TEST
program testDMRG_GEOMETRY
  USE DMRG_GEOMETRY
  implicit none
  integer,dimension(:),allocatable :: map
  !
  map=Build_Block_Growth_Map(7,.true.)
  call assert_map(map,[7,5,3,1,2,4,6],'PBC block map')
  map=Build_DMRG_Growth_Map(10,.true.)
  call assert_map(map,[3,8,4,7,2,9,5,6,1,10],'PBC physical map')
  map=Build_DMRG_Growth_Map(8,.false.)
  call assert_map(map,[1,8,2,7,3,6,4,5],'OBC physical map')
  print *,'DMRG_GEOMETRY TEST: PASS'

contains


  subroutine assert_map(actual,expected,message)
    integer,dimension(:),intent(in) :: actual,expected
    character(len=*),intent(in)     :: message
    !
    if(size(actual)/=size(expected))&
         error stop 'testDMRG_GEOMETRY: '//trim(message)//' size'
    if(any(actual/=expected))&
         error stop 'testDMRG_GEOMETRY: '//trim(message)//' value'
    !
  end subroutine assert_map


end program testDMRG_GEOMETRY
#endif
