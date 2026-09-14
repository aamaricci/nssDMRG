MODULE LIST_OMATRICES
  USE SCIFOR, only: str,free_unit
  USE MATRIX_SPARSE
  implicit none
  private

  character(len=*),parameter :: omatrices_list_header = "NSSDMRG_OMATRICES_LIST"
  integer,parameter          :: omatrices_list_format = 2

  type omattype
     integer                      :: index=0
     character(len=:),allocatable :: ckey
     type(sparse_matrix)          :: mat
     type(omattype),pointer       :: next=>null()
  end type omattype



  !U MATRICES LIST
  type omatrices_list
     integer               :: size=0
     type(omattype),pointer :: root=>null()
   contains
     procedure,pass :: free     => free_omatrices_list
     procedure,pass :: put      => put_omatrices_list
     procedure,pass :: append   => append_omatrices_list
     procedure,pass :: load     => load_omatrices_list
     procedure,pass :: dump     => dump_op_omatrices_list
     procedure,pass :: op       => get_op_omatrices_list
     procedure,pass :: key      => get_key_omatrices_list
     procedure,pass :: has_key  => has_key_omatrices_list
     procedure,pass :: is_valid => is_valid_omatrices_list
     procedure,pass :: show     => show_omatrices_list
     procedure,pass :: write    => write_omatrices_list
     procedure,pass :: read     => read_omatrices_list
  end type omatrices_list



  !GENERIC CONSTRUCTOR
  interface omatrices_list
     module procedure :: construct_from_single_operator
  end interface omatrices_list 



  !EQUALITY
  interface assignment(=)
     module procedure :: equality_omatrices_list
  end interface


  !INTRINSIC FUNCTION SIZE(OPERATORS_LIST)
  intrinsic :: size
  interface size
     module procedure :: size_omatrices_list
  end interface


  public :: omatrices_list
  public :: size
  public :: assignment(=)




contains




  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################     
  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor: given a key+operator
  !+------------------------------------------------------------------+
  function construct_from_single_operator(key,op) result(self)
    type(omatrices_list)                     :: self
    character(len=*),intent(in)              :: key
    type(sparse_matrix),intent(in)           :: op
    call self%free()
    allocate(self%root)
    call self%put(key,op)
  end function construct_from_single_operator


  !+------------------------------------------------------------------+
  !PURPOSE:  Free an operators_list (destructor) 
  !+------------------------------------------------------------------+
  recursive subroutine free_omatrices_list(self)
    class(omatrices_list),intent(inout) :: self
    type(omattype),pointer              :: p,c
    if(.not.associated(self%root))return
    do
       p=>self%root
       c=>p%next
       if(.not.associated(c))exit
       p%next=>c%next
       c%next=>null()
       call c%mat%free()
       if(allocated(c%ckey))deallocate(c%ckey)
       deallocate(c)
    enddo
    self%size=0
    deallocate(self%root)
    nullify(self%root,p,c)
  end subroutine free_omatrices_list







  !##################################################################
  !##################################################################
  !       PUT/LOAD/APPEND  - GET/DUMP OPERATORS IN/FROM A LIST
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Put a sparse matrix in the omatrices_list
  !+------------------------------------------------------------------+
  subroutine put_omatrices_list(self,key,op)
    class(omatrices_list),intent(inout)  :: self
    character(len=*),intent(in)          :: key
    type(sparse_matrix),intent(in)       :: op
    type(omattype),pointer               :: p,c
    logical                              :: found
    if(.not.associated(self%root))allocate(self%root)
    found=.false.
    p=>self%root
    c=>p%next
    do
       if(.not.associated(c))exit
       if(str(c%ckey)==str(key))then
          found=.true.
          exit
       endif
       p=>c
       c=>c%next
    enddo
    if(found)then
       c%ckey=str(key)
       c%mat=op
    else
       allocate(p%next)
       p%next%ckey=str(key)
       p%next%index=p%index+1
       p%next%mat=op
       p%next%next=>null()
       self%size=self%size+1
    endif
    nullify(p,c)
  end subroutine put_omatrices_list



  !+------------------------------------------------------------------+
  !PURPOSE:  Append == Put a sparse matrix as operator in the operators_list
  !+------------------------------------------------------------------+
  subroutine append_omatrices_list(self,key,op)
    class(omatrices_list),intent(inout) :: self
    character(len=*),intent(in)         :: key
    type(sparse_matrix),intent(in)      :: op
    call self%put(str(key),op)
  end subroutine append_omatrices_list




  !+------------------------------------------------------------------+
  !PURPOSE:  Load a dense matrix as operator in the omatrices_list
  !+------------------------------------------------------------------+
  subroutine load_omatrices_list(self,key,op)
    class(omatrices_list),intent(inout)      :: self
    character(len=*),intent(in)              :: key
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in)     :: op
#else
    real(8),dimension(:,:),intent(in)        :: op
#endif
    if(.not.associated(self%root))allocate(self%root)
    call self%put(key,as_sparse(op))
  end subroutine load_omatrices_list



  !+------------------------------------------------------------------+
  !PURPOSE: Dump operator of the omatrices_list as a dense matrix given a key 
  !+------------------------------------------------------------------+
  function dump_op_omatrices_list(self,key) result(matrix)
    class(omatrices_list),intent(inout)   :: self
    character(len=*),intent(in)           :: key
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: matrix
#else
    real(8),dimension(:,:),allocatable    :: matrix
#endif
    matrix = as_matrix( self%op(key=key) )  
  end function dump_op_omatrices_list





  !##################################################################
  !##################################################################
  !              RETRIEVE CONTENT: OP, KEY, DQ
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE: Return a sparse matrix of the list given:
  ! + key: the operator corresponding to the key value
  ! + indx: the operator corresponding to  the indx value
  !+------------------------------------------------------------------+
  function get_op_omatrices_list(self,key,index) result(op)
    class(omatrices_list),intent(in)     :: self
    character(len=*),intent(in),optional :: key
    integer,intent(in),optional          :: index
    type(sparse_matrix)                  :: op
    type(omattype),pointer               :: c
    logical                              :: found
    nullify(c)
    if(present(key).eqv.present(index))&
         stop "get_op_omatrices_list: specify exactly one of key or index"
    found=.false.
    if(associated(self%root))c=>self%root%next
    do while(associated(c))
       if(present(key))then
          found=str(c%ckey)==str(key)
       else
          found=c%index==index
       endif
       if(found)exit
       c=>c%next
    enddo
    if(.not.found)stop "get_op_omatrices_list: matrix not found"
    op=c%mat
    nullify(c)
  end function get_op_omatrices_list




  !+------------------------------------------------------------------+
  !PURPOSE: Return key of the operators_list  corresponding to:
  ! + indx: the given indx value
  !+------------------------------------------------------------------+  
  function get_key_omatrices_list(self,index) result(key)
    class(omatrices_list),intent(in) :: self
    integer,intent(in)               :: index
    character(len=:),allocatable     :: key
    type(omattype),pointer           :: c
    if(index<1.OR.index>self%size)stop "get_key_omatrices_list: index out of range"
    c=>self%root%next
    do while(c%index/=index)
       c=>c%next
    enddo
    key=str(c%ckey)
    nullify(c)
  end function get_key_omatrices_list

  
  !+------------------------------------------------------------------+
  !PURPOSE:  Returns True is key exists, False otherwise
  !+------------------------------------------------------------------+  
  function has_key_omatrices_list(self,key) result(found)
    class(omatrices_list),intent(in) :: self
    character(len=*),intent(in)      :: key
    logical                          :: found
    type(omattype),pointer           :: c
    found=.false.
    if(.not.associated(self%root))return
    c=>self%root%next
    do while(associated(c))
       if(str(c%ckey)==str(key))then
          found=.true.
          exit
       endif
       c=>c%next
    enddo
    nullify(c)
  end function has_key_omatrices_list



  !+------------------------------------------------------------------+
  !PURPOSE:  Check if operators_list is a valid one, ie the operators in the
  ! dictionary have the right dimensions:
  !+------------------------------------------------------------------+  
  function is_valid_omatrices_list(self) result(valid)
    class(omatrices_list),intent(in) :: self
    logical                          :: valid
    integer                          :: count
    type(omattype),pointer           :: c
    nullify(c)
    valid=.true.
    count=0
    if(.not.associated(self%root))then
       valid=self%size==0
       return
    endif
    c=>self%root%next
    do while(associated(c))
       count=count+1
       valid=valid.AND.c%index==count
       valid=valid.AND.allocated(c%ckey)
       valid=valid.AND.c%mat%status
       c=>c%next
    enddo
    valid=valid.AND.count==self%size
    nullify(c)
  end function is_valid_omatrices_list




  !##################################################################
  !##################################################################
  !              INTRISIC FUNCTIONS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Returns the size of given operators_list
  !+------------------------------------------------------------------+  
  function size_omatrices_list(self) result(n)
    type(omatrices_list),intent(in) :: self
    integer                         :: n
    n=self%size
  end function size_omatrices_list





  !##################################################################
  !##################################################################
  !               SHOW 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Pretty print an omatrices_list
  !+------------------------------------------------------------------+  
  subroutine show_omatrices_list(self,fmt,unit,file)
    class(omatrices_list),intent(in) :: self
    character(len=*),intent(in),optional :: fmt
    integer,intent(in),optional      :: unit
    character(len=*),intent(in),optional :: file
    integer                          :: unit_
    type(omattype),pointer           :: c
    nullify(c)
    unit_=6
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    write(unit_,"(A7,I12)")"Size :",self%size
    write(unit_,"(A18)")"------------------"
    if(associated(self%root))c=>self%root%next
    do while(associated(c))
       write(unit_,"(A7,I12)")"Index: ",c%index
       write(unit_,"(A7,A)")"Key  : ",str(c%ckey)
       call c%mat%display()
       write(unit_,*)""
       c=>c%next
    enddo
    if(present(file))close(unit_)
    nullify(c)
  end subroutine show_omatrices_list




  subroutine write_omatrices_list(self,file,unit)
    class(omatrices_list),intent(in) :: self
    character(len=*),intent(in),optional :: file
    integer,intent(in),optional      :: unit
    integer                          :: unit_
    type(omattype),pointer           :: c
    nullify(c)
    unit_=-1
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    if(unit_==-1)stop "write_omatrices_list: no file or unit"
    if(.not.self%is_valid())stop "write_omatrices_list: invalid list"
    write(unit_,"(A,1X,I0)")omatrices_list_header,omatrices_list_format
    write(unit_,*)self%size
    if(associated(self%root))c=>self%root%next
    do while(associated(c))
       write(unit_,*)str(c%ckey)
       call c%mat%write(unit=unit_)
       c=>c%next
    enddo
    if(present(file))close(unit_)
    nullify(c)
  end subroutine write_omatrices_list


  subroutine read_omatrices_list(self,file,unit)
    class(omatrices_list),intent(inout) :: self
    character(len=*),intent(in),optional :: file
    integer,intent(in),optional          :: unit
    integer                              :: unit_,ios,version,n,i
    character(len=512)                   :: line,header,key
    type(sparse_matrix)                  :: mat
    unit_=-1
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    if(unit_==-1)stop "read_omatrices_list: no file or unit"
    read(unit_,"(A)",iostat=ios)line
    if(ios/=0)stop "read_omatrices_list: unable to read header"
    read(line,*,iostat=ios)header,version
    if(ios/=0.OR.str(header)/=omatrices_list_header.OR.version/=omatrices_list_format)&
         stop "read_omatrices_list: incompatible format"
    call self%free()
    read(unit_,*,iostat=ios)n
    if(ios/=0.OR.n<0)stop "read_omatrices_list: invalid size"
    do i=1,n
       read(unit_,*,iostat=ios)key
       if(ios/=0)stop "read_omatrices_list: unable to read key"
       call mat%read(unit=unit_)
       call self%put(str(key),mat)
    enddo
    call mat%free()
    if(present(file))close(unit_)
  end subroutine read_omatrices_list




  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Equality between two operators_lists (full copy)
  !+------------------------------------------------------------------+  
  subroutine equality_omatrices_list(lhs,rhs)
    type(omatrices_list),intent(inout) :: lhs
    type(omatrices_list),intent(in)    :: rhs
    integer                            :: i
    call lhs%free()
    do i=1,size(rhs)
       call lhs%put(rhs%key(i),rhs%op(index=i))
    enddo
  end subroutine equality_omatrices_list

END MODULE LIST_OMATRICES
