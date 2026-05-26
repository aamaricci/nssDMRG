MODULE LIST_OPERATORS
  USE SCIFOR, only:str,free_unit
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  implicit none
  private

  type optype
     integer                      :: index=0
     character(len=:),allocatable :: ckey
     character(len=:),allocatable :: ctype
     type(sparse_matrix)          :: ope
     type(optype),pointer         :: next   =>null()
  end type optype


  !OPERATORS DICTIONARY 
  type operators_list
     integer              :: size=0
     type(optype),pointer :: root   =>null()
   contains
     procedure,pass :: free     => free_operators_list     !destructor
     procedure,pass :: put      => put_operators_list      !put sparse operator
     procedure,pass :: append   => append_operators_list   !put sparse operator
     procedure,pass :: show     => show_operators_list     !show operators_list to screen
     procedure,pass :: load     => load_operators_list     !load dense matrix operator
     procedure,pass :: dump     => dump_op_operators_list  !dump dense matrix operator
     procedure,pass :: get      => get_all_operators_list  !get {key,operator,type}
     procedure,pass :: op       => get_op_operators_list   !return operator given:key,indx
     procedure,pass :: key      => get_key_operators_list  !return key given: indx
     procedure,pass :: type     => get_type_operators_list  !return type given: key,indx
     procedure,pass :: keys     => keys_operators_list     !return all the keys
     procedure,pass :: types    => types_operators_list     !return all the types
     procedure,pass :: has_key  => has_key_operators_list  !True if key exists
     procedure,pass :: is_valid => is_valid_operators_list !True if operators_list is valid
     procedure,pass :: shape    => shape_operators_list
     procedure,pass :: write    => write_operators_list  !write operators list
     procedure,pass :: read     => read_operators_list  !read operators list
  end type operators_list


  !GENERIC CONSTRUCTOR
  interface operators_list
     module procedure :: construct_from_single_operator
     module procedure :: construct_from_array_operator
  end interface operators_list

  !EQUALITY 
  interface assignment(=)
     module procedure :: equality_operators_list
  end interface assignment(=)

  !INTRINSIC FUNCTION SIZE(OPERATORS_LIST)
  intrinsic :: size
  interface size
     module procedure :: size_operators_list
  end interface size

  !INTRINSIC FUNCTION SHAPE(OPERATORS_LIST)
  interface shape
     module procedure :: shape_operators_list
  end interface shape

  public :: operators_list
  public :: size
  public :: shape
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
  function construct_from_single_operator(key,op,type) result(self)
    type(operators_list)      :: self
    character(len=*)          :: key
    character(len=*),optional :: type
    type(sparse_matrix)       :: op
    character(len=16)         :: type_
    type_='';if(present(type))type_=str(type)
    call self%free()
    allocate(self%root)
    call self%put(key,op,type_)
  end function construct_from_single_operator


  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor: given a list of keys+operators
  !+------------------------------------------------------------------+
  function construct_from_array_operator(keys,ops,types) result(self)
    type(operators_list)                            :: self
    character(len=*),dimension(:)                   :: keys
    character(len=*),dimension(size(keys)),optional :: types
    type(sparse_matrix),dimension(size(keys))       :: ops
    character(len=16),dimension(size(keys))         :: types_
    integer                                         :: i
    types_=''
    if(present(types))then
       do i=1,size(keys)
          types_(i)=str(types(i))
       enddo
    endif
    call self%free()
    allocate(self%root)
    do i=1,size(keys)
       call self%put(keys(i),ops(i),types_(i))
    enddo
  end function construct_from_array_operator






  !+------------------------------------------------------------------+
  !PURPOSE:  Free an operators_list (destructor) 
  !+------------------------------------------------------------------+
  recursive subroutine free_operators_list(self)
    class(operators_list),intent(inout) :: self
    type(optype),pointer                :: p,c
    if(.not.associated(self%root))return
    do
       p => self%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next
       c%next => null()
       c%index=  0
       call c%ope%free()         !<- use sparse_matrix free procedure
       if(allocated(c%ckey))deallocate(c%ckey)
       if(allocated(c%ctype))deallocate(c%ctype)
       deallocate(c)
    enddo
    self%size=0
    self%root=>null()
    p=>null()
    c=>null()
  end subroutine free_operators_list





  !##################################################################
  !##################################################################
  !       PUT/LOAD/APPEND  - GET/DUMP OPERATORS IN/FROM A LIST
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Put a sparse matrix as operator in the operators_list
  !+------------------------------------------------------------------+
  subroutine put_operators_list(self,key,op,type)
    class(operators_list),intent(inout)  :: self
    character(len=*),intent(in)          :: key
    character(len=*),intent(in),optional :: type
    type(sparse_matrix),intent(in)       :: op
    type(optype),pointer                 :: p,c
    logical                              :: iadd
    character(len=16)                    :: type_
    !
    type_='';if(present(type))type_=str(type)
    !
    if(.not.associated(self%root))allocate(self%root)
    !
    iadd = .false.
    p => self%root
    c => p%next
    do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if ( (str(c%ckey)  == str(key)) ) then
          iadd = .true.
          exit
       endif
       p => c
       c => c%next
    end do
    !
    if(iadd)then                !KEY exists: update operator
       c%ckey         = str(key)
       c%ctype        = str(type)
       ! call c%ope%copy(op)
       c%ope          = op
    else                        !QN does not exist: create a new element
       allocate(p%next)
       p%next%ckey    = str(key)
       p%next%ctype   = str(type)
       p%next%index   = p%index+1
       p%next%ope     = op
       ! call p%next%ope%copy(op)
       if(.not.associated(c))then
          p%next%next => null()
       else
          p%next%next => c
       end if
       self%size      = self%size+1
    endif
    p=>null()
    c=>null()
  end subroutine put_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE:  Append == Put a sparse matrix as operator in the operators_list
  !+------------------------------------------------------------------+
  subroutine append_operators_list(self,key,op,type)
    class(operators_list),intent(inout)  :: self
    character(len=*),intent(in)          :: key
    character(len=*),intent(in),optional :: type
    type(sparse_matrix),intent(in)       :: op
    character(len=16)                    :: type_
    type_='';if(present(type))type_=str(type)
    call self%put(str(key),op,str(type_))
  end subroutine append_operators_list




  !+------------------------------------------------------------------+
  !PURPOSE:  Load a dense matrix as operator in the operators_list
  !+------------------------------------------------------------------+
  subroutine load_operators_list(self,key,op,type)
    class(operators_list),intent(inout)  :: self
    character(len=*),intent(in)          :: key
    character(len=*),intent(in),optional :: type
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in) :: op
#else
    real(8),dimension(:,:),intent(in)    :: op
#endif
    character(len=16)                    :: type_
    type_='';if(present(type))type_=str(type)    
    if(.not.associated(self%root))allocate(self%root)
    call self%put(key,as_sparse(op),type_)
  end subroutine load_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE: Dump operator of the operators_list as a dense matrix  given a key 
  !+------------------------------------------------------------------+
  function dump_op_operators_list(self,key) result(matrix)
    class(operators_list),intent(inout)   :: self
    character(len=*),intent(in)           :: key
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: matrix
#else
    real(8),dimension(:,:),allocatable    :: matrix
#endif
    matrix = as_matrix( self%op(key=key) )  
  end function dump_op_operators_list





  !##################################################################
  !##################################################################
  !              RETRIEVE CONTENT: OP, KEY
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Get {key,operator,type} of the list given an index
  !+------------------------------------------------------------------+
  subroutine get_all_operators_list(self,index,key,op,type)
    class(operators_list),intent(inout)   :: self
    integer,intent(in)                    :: index
    character(len=*),intent(out)          :: key
    character(len=*),intent(out),optional :: type
    type(sparse_matrix),intent(out)       :: op
    integer                               :: index_
    type(optype),pointer                  :: c
    logical                               :: ifound
    !
    index_=index
    if(index_>self%size.OR.index_<=0)stop "get_sectors_list: index !in [1,self.size]"    
    !
    ifound=.false.
    c => self%root%next
    do                            !traverse the list until KEY is found
       if(.not.associated(c))exit
       if(c%index == index_) then
          ifound=.true.
          exit          
       endif
       c => c%next
    end do
    if(.not.ifound)stop "get_op_operators_list error: not found"
    !
    op  = c%ope
    key = str(c%ckey)
    if(present(type))type= str(c%ctype)
    !
    c=>null()
  end subroutine get_all_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return operator of the list as sparse matrix given:
  ! + key: the operator corresponding to the key value
  ! + indx: the operator corresponding to  the indx value
  !+------------------------------------------------------------------+
  function get_op_operators_list(self,key,index) result(op)
    class(operators_list)     :: self
    character(len=*),optional :: key
    integer,optional          :: index
    type(sparse_matrix)       :: op
    integer                   :: index_
    type(optype),pointer      :: c
    logical                   :: ifound
    !
    index_=self%size;if(present(index))index_=index
    if(index_>self%size.OR.index_<=0)stop "get_op_operators_list: index !in [1,self.size]"
    if(.not.present(index).AND..not.present(key))stop "get_op_operators_list: no input given: use index=i OR key=str"
    ifound=.false.
    c => self%root%next
    loop:do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(present(key))then
          if (str(c%ckey) == str(key)) then
             ifound=.true.
             exit loop
          endif
       elseif(c%index == index_)then
          ifound=.true.
          exit
       endif
       c => c%next
    end do loop
    if(.not.ifound)then
       if(present(key))then
          print*,"get_op_operators_list error: key not found",key
       else
          print*,"get_op_operators_list error: index not found",index_
       endif
       stop
    endif
    !
    op = c%ope    
    !
    c=>null()
  end function get_op_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return key of the operators_list  corresponding to:
  ! + indx: the given indx value
  !+------------------------------------------------------------------+  
  function get_key_operators_list(self,index) result(key)
    class(operators_list)        :: self
    integer                      :: index
    character(len=:),allocatable :: key
    integer                      :: index_
    type(optype),pointer         :: c
    logical                      :: ifound
    !
    index_=index
    if(index_>self%size.OR.index_<=0)stop "get_key_operators_list: index !in [1,self.size]"
    !
    c => self%root%next
    do                            !traverse the list until index is found
       if(.not.associated(c))exit
       if(c%index == index_) then
          ifound=.true.
          exit
       endif
       c => c%next
    end do
    if(.not.ifound)stop "get_key error: not found"
    !
    key = str(c%ckey)
    !
    c=>null()
  end function get_key_operators_list


  !+------------------------------------------------------------------+
  !PURPOSE: Return operator type of the list given:
  ! + key  : the type corresponding to the key value
  ! + index: the type corresponding to  the indx value
  !+------------------------------------------------------------------+
  function get_type_operators_list(self,index,key) result(type)
    class(operators_list)        :: self
    integer                      :: index
    character(len=*),optional    :: key
    character(len=:),allocatable :: type
    integer                      :: index_
    type(optype),pointer         :: c
    logical                      :: ifound
    !
    index_=index
    if(index_>self%size.OR.index_<=0)stop "get_op_operators_list: index !in [1,self.size]"
    !
    ifound=.false.
    c => self%root%next
    loop:do                            !traverse the list until QN is found
       if(.not.associated(c))exit
       if(present(key))then
          if (str(c%ckey) == str(key)) then
             ifound=.true.
             exit loop
          endif
       elseif(c%index == index_)then
          ifound=.true.
          exit
       endif
       c => c%next
    end do loop
    if(.not.ifound)stop "get_type_operators_list error: not found"
    !
    type = str(c%ctype)
    !
    c=>null()
  end function get_type_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return all the keys in the operators_list
  !+------------------------------------------------------------------+  
  function keys_operators_list(self,len) result(keys)
    class(operators_list)                       :: self
    integer                                     :: len
    character(len=len),dimension(:),allocatable :: keys
    integer                                     :: i,Nsize
    Nsize=size(self)
    allocate(keys(Nsize))
    do i=1,Nsize
       keys(i) = str(self%key(i))
    enddo
  end function keys_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return all the keys in the operators_list
  !+------------------------------------------------------------------+  
  function types_operators_list(self,len) result(types)
    class(operators_list)                       :: self
    integer                                     :: len
    character(len=len),dimension(:),allocatable :: types
    integer                                     :: i,Nsize
    Nsize=size(self)
    allocate(types(Nsize))
    do i=1,Nsize
       types(i) = str(self%type(i))
    enddo
  end function types_operators_list




  !+------------------------------------------------------------------+
  !PURPOSE:  Check if operators_list is a valid one, ie the operators in the
  ! dictionary have the right dimensions:
  ! N = size(dim)
  ! N=1 => shape(op)=[dim(1),dim(1)]
  ! N>1 => shape(op)= [dim(1),dim(1)]&&...&&[dim(N),dim(N)].OR.[prod(dim),prod(dim)]
  !+------------------------------------------------------------------+
  function is_valid_operators_list(self,dim) result(bool)
    class(operators_list),intent(inout) :: self
    integer,optional                    :: dim
    integer                             :: dim_
    integer                             :: oshape(2),Ndim,i
    logical                             :: bool,bdim
    type(optype),pointer                :: c
    bool = .true.
    dim_ = 0;if(present(dim))dim_=dim
    c => self%root%next
    do 
       if(.not.associated(c))exit
       if(dim_==0)then
          oshape = [c%ope%Nrow,c%ope%Ncol]!shape(c%ope)
          dim_   = oshape(1)
       endif
       bool = bool.AND.(all([c%ope%Nrow,c%ope%Ncol] == [dim_,dim_]))
       c => c%next
    enddo
    c=>null()
  end function is_valid_operators_list




  !##################################################################
  !##################################################################
  !              ENUMERATOR & ITERATORS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Returns the size of given operators_list
  !+------------------------------------------------------------------+
  function size_operators_list(self) result(size)
    class(operators_list),intent(in) :: self
    integer                          :: size
    size = self%size
  end function size_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE:  Returns True is key exists, False otherwise
  !+------------------------------------------------------------------+
  function has_key_operators_list(self, key) result(bool)
    class(operators_list),intent(inout) :: self
    character(len=*),intent(in)         :: key
    logical                             :: bool
    type(optype),pointer                :: c
    !
    bool=.false.
    if(.not.associated(self%root))return
    c => self%root%next
    do                            !traverse the list until index is found
       if(.not.associated(c))exit
       if(str(c%ckey) == str(key)) then
          bool=.true.
          exit
       endif
       c => c%next
    end do
    c=>null()
  end function has_key_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE:  Returns the shape of the operators in the operators_list
  ! If valid list all operators have same shape so the first is fine. 
  !+------------------------------------------------------------------+
  function shape_operators_list(self) result(shape)
    class(operators_list),intent(inout) :: self
    integer,dimension(2)                :: shape
    type(optype),pointer                :: c
    logical                             :: bool
    bool = self%is_valid()
    if(.not.bool)stop "shape_operator_list: not a valid list"
    c => self%root%next
    shape = [c%ope%Nrow,c%ope%Ncol]
  end function shape_operators_list









  !##################################################################
  !##################################################################
  !               SHOW 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Pretty print an operators_list
  !+------------------------------------------------------------------+
  recursive subroutine show_operators_list(self,fmt,unit,file)
    class(operators_list),intent(inout) :: self
    character(len=*),optional           :: fmt
    integer,optional                    :: unit
    character(len=32)                   :: fmt_
    integer                             :: i,count=0
    type(optype),pointer                :: c
    character(len=*),optional           :: file
    integer                             :: unit_
    unit_=6
    if(present(unit))unit_=unit
    if(present(file))open(free_unit(unit_),file=str(file))
    !
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    !
    write(unit_,"(A6,I12)")"Size :",self%size
    write(unit_,"(A18)")"------------------"
    c => self%root%next
    do
       if(.not.associated(c))exit
       count=count+1
       write(unit_,"(A6,I12)")  "Index:",c%index
       write(unit_,"(A6,A)")"Key  :",str(c%ckey)
       write(unit_,"(A6,A)")"Type :",str(c%ctype)
       call c%ope%display()
       write(unit_,*)""
       c => c%next
    end do
    c=>null()
    if(present(file))close(unit_)
  end subroutine show_operators_list




  subroutine write_operators_list(self,file,unit)
    class(operators_list),intent(inout) :: self
    character(len=*),optional           :: file
    integer,optional                    :: unit
    integer                             :: i
    type(optype),pointer                :: c
    integer                             :: unit_
    !
    unit_=-1
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    if(unit_==-1)stop "write_operators_list error: no input +file or +unit given"
    !
    write(unit_,*)self%size
    if(.not.associated(self%root))then
       if(present(file))close(unit_)
       return
    endif
    c => self%root%next
    do
       if(.not.associated(c))exit
       write(unit_,*)str(c%ckey)
       if(str(c%ctype)=="")then
          write(unit_,*)"none"
       else
          write(unit_,*)str(c%ctype)
       endif
       call c%ope%write(unit=unit_)
       c => c%next
    end do
    c=>null()
    if(present(file))close(unit_)
  end subroutine write_operators_list





  subroutine read_operators_list(self,file,unit)
    class(operators_list),intent(inout) :: self
    character(len=*),optional           :: file
    integer,optional                    :: unit
    integer                             :: i,ListSize
    integer                             :: unit_
    character(len=512)                  :: key
    character(len=512)                  :: type
    type(sparse_matrix)                 :: ope
    !
    unit_=-1
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    if(unit_==-1)stop "read_operators_list error: no input +file or +unit given"
    !
    read(unit_,*)ListSize
    do i=1,ListSize
       read(unit_,*)key
       read(unit_,*)type
       if(str(type)=="none")type=""
       call ope%read(unit=unit_)
       call self%append(key=str(key),op=ope,type=str(type))
    end do
    if(present(file))close(unit_)
  end subroutine read_operators_list





  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Equality between two operators_lists
  !+------------------------------------------------------------------+
  subroutine equality_operators_list(A,B)
    type(operators_list),intent(inout) :: A
    type(operators_list),intent(in)    :: B
    integer                            :: i
    call A%free()
    do i=1,size(B)
       call A%put(B%key(index=i),B%op(index=i),B%type(index=i))
    enddo
  end subroutine equality_operators_list




END MODULE LIST_OPERATORS






!##################################################################
!##################################################################
!##################################################################
!##################################################################
!                          /_  __/ ____/ ___/_  __/
!                           / / / __/  \__ \ / /   
!                          / / / /___ ___/ // /    
!                         /_/ /_____//____//_/     
!##################################################################
!##################################################################
!##################################################################
!##################################################################
#ifdef _TEST
program testOPERATORS_TUPLE
  USE SCIFOR
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  implicit none
  type(operators_list)                  :: my_list,a_list
  type(operators_list)                  :: copy_list,clist(2)
  type(sparse_matrix)                   :: spSz,spSp,spH,spK,a,b,c
  integer                               :: i,j,n
  logical                               :: bool
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable :: mat
  complex(8),dimension(2,2),parameter   :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  complex(8),dimension(2,2),parameter   :: S0=pauli_0
  complex(8),dimension(2,2),parameter   :: Sz=pauli_z
  complex(8),dimension(2,2),parameter   :: Sx=pauli_x
  complex(8),dimension(2,2),parameter   :: Splus=reshape([zero,zero,one,zero],[2,2])
  complex(8),dimension(4,4)             :: Gamma13,Gamma03
#else
  real(8),dimension(:,:),allocatable    :: mat
  real(8),dimension(2,2),parameter      :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  real(8),dimension(2,2),parameter      :: S0=pauli_0
  real(8),dimension(2,2),parameter      :: Sz=pauli_z
  real(8),dimension(2,2),parameter      :: Sx=pauli_x
  real(8),dimension(2,2),parameter      :: Splus=reshape([zero,zero,one,zero],[2,2])
  real(8),dimension(4,4)                :: Gamma13,Gamma03
#endif
  character(len=10)                     :: key,type
  character(len=10),allocatable         :: keys(:)
  integer,parameter                     :: sec=500
  integer,parameter                     :: sec=1


  Gamma13=kron(Sx,Sz)
  Gamma03=kron(S0,Sz)


  print*,"TEST CONSTRUCTOR, PUT, SHOW, FREE"
  my_list = operators_list(&
       ['H0','Sz','Sp','P ','C '],&
       [sparse(Hzero),sparse(Sz),sparse(Splus),sparse(Sz),sparse(Sz)],&
       ['bose ','bose ','bose ','sign ','fermi'])
  call my_list%show()
  call my_list%free()
  call wait(sec)


  print*,"TEST LOAD matrices"
  call my_list%load("H0",Hzero,'b')
  call my_list%load("Sz",Sz,'s')
  call my_list%load("Sp",Splus,'fermi')
  print*,"TEST SHOW"
  call my_list%show()
  call my_list%free()
  call wait(sec)



  print*,"TEST (CONSTRUCT + )APPEND matrices"
  call my_list%append("H0",as_sparse(Hzero),'b')
  call my_list%append("Sz",as_sparse(Sz),'b')
  call my_list%append("Sp",as_sparse(Splus),'Bosonic')
  call my_list%show()
  print*,""
  call wait(sec)





  print*,"TEST RETRIEVE FUNCTIONALITIES"
  print*,"TEST .DUMP"
  print*,"Mat.allocated:",allocated(mat)
  print*,"dump Sp -> Mat"
  mat = my_list%dump("Sp")
  print*,"Mat.allocated:",allocated(mat)
  do i=1,size(mat,1)
     write(*,*)(mat(i,j),j=1,size(mat,2))
  enddo
  deallocate(mat)
  print*,""
  call wait(sec)


  print*,"TEST .GET"
  do i=1,size(my_list)
     call my_list%get(index=i,key=key,op=a,type=type)
     print*,i
     print*,key
     print*,type
     call a%show()
  enddo
  print*,""
  call wait(sec)


  print*,"TEST .KEY + .OP + ITERATION over index"
  do i=1,size(my_list)
     a = my_list%op(index=i)
     print*,i
     call a%show()
  enddo
  print*,""
  do i=1,size(my_list)
     key = my_list%key(index=i)
     a = my_list%op(key=key)
     print*,i,key
     call a%show
  enddo
  print*,""
  call wait(sec)


  print*,"TEST HAS_KEY"
  print*,"list has key Sz",my_list%has_key("Sz")
  print*,"list has key SZ",my_list%has_key("SZ")
  print*,""
  call wait(sec)



  print*,"TEST IS_VALID "
  print*,my_list%is_valid()
  print*,"is valid with dim=2"
  print*,my_list%is_valid(dim=2)
  print*,"is not valid with dim=3"
  print*,my_list%is_valid(dim=3)
  print*,"is not valid once appended s_0.x.s_3"
  call my_list%append("W",as_sparse(Gamma03),'b')
  print*,my_list%is_valid()
  call my_list%free
  print*,""
  call wait(sec)



  call my_list%append("H0",as_sparse(Hzero),'b')
  call my_list%append("Sz",as_sparse(Sz),'b')
  call my_list%load("Sp",Splus,'b')



  print*,"TEST DEEP COPY ="
  copy_list = my_list
  call copy_list%show()
  print*,copy_list%is_valid()
  print*,""
  call wait(sec)


  print*,"TEST my_list.o('key')"
  print*,"before a=empty"
  call a%free
  call a%show
  print*,"a = my_list%op('Sz')"
  a = my_list%op("Sz")
  print*,"a.print"
  call a%show
  print*,""
  call wait(sec)



  print*,"TEST ITERATION SIZE:"
  do i=1,size(my_list)
     a = my_list%op(index=i)
     print*,i,my_list%key(i),my_list%type(i)
     call a%show()
  enddo
  print*,""
  call wait(sec)

  print*,"TEST ITERATION KEYS:"
  keys = my_list%keys(len(keys))
  do i=1,size(keys)
     a = my_list%op(key=str(keys(i)))
     print*,i,str(keys(i))
     call a%show()
  enddo
  print*,""
  call wait(sec)



  print*,"TEST DEEP COPY '='"
  Gamma13=kron(Sx,Sz)
  Gamma03=kron(S0,Sz)
  call a_list%append("gamma13",as_sparse(Gamma13),'b')
  call a_list%append("gamma03",as_sparse(Gamma03),'b')
  call a_list%append("Gamma33",as_sparse(kron(Sz,Sz)),'b')

  clist(1) = my_list
  clist(2) = a_list

  call clist(1)%show()
  call clist(2)%show()
  print*,""
  call wait(sec)

  print*,"TEST WRITE/READ"
  call a_list%free()
  call a_list%append("gamma13",as_sparse(Gamma13),'b')
  call a_list%append("gamma03",as_sparse(Gamma03),'b')
  call a_list%append("Gamma33",as_sparse(kron(Sz,Sz)),'b')

  print*,"write:"
  call a_list%write(file="a_list_write.dat")

  print*,"show:"
  call a_list%show()

  print*,"free:"
  call a_list%free()

  print*,"read:"
  call a_list%read(file="a_list_write.dat")

  print*,"show again:"
  call a_list%show()

end program testOPERATORS_TUPLE
#endif
