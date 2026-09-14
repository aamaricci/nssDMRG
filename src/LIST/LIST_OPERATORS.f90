MODULE LIST_OPERATORS
  USE SCIFOR, only:str,free_unit
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  implicit none
  private

  character(len=*),parameter :: operators_list_header = "NSSDMRG_OPERATORS_LIST"
  integer,parameter          :: operators_list_format = 3

  type optype
     integer                      :: index=0
     character(len=:),allocatable :: ckey
     character(len=:),allocatable :: ctype
     type(sparse_matrix)          :: ope
     real(8),allocatable          :: dq(:)
     type(optype),pointer         :: next   =>null()
  end type optype


  !OPERATORS DICTIONARY 
  type operators_list
     integer              :: size=0
     type(optype),pointer :: root   =>null()
   contains
     procedure,pass :: free         => free_operators_list     !destructor
     procedure,pass :: put          => put_operators_list      !put sparse operator
     procedure,pass :: update       => update_operators_list   !update matrix, preserve metadata
     procedure,pass :: append       => append_operators_list   !put sparse operator
     procedure,pass :: show         => show_operators_list     !show operators_list to screen
     procedure,pass :: load         => load_operators_list     !load dense matrix operator
     procedure,pass :: dump         => dump_op_operators_list  !dump dense matrix operator
     procedure,pass :: get          => get_all_operators_list  !get {key,operator,type}
     procedure,pass :: op           => get_op_operators_list   !return operator given:key,indx
     procedure,pass :: key          => get_key_operators_list  !return key given: indx
     procedure,pass :: type         => get_type_operators_list  !return type given: key,indx
     procedure,pass :: dq           => get_dq_operators_list
     procedure,pass :: keys         => keys_operators_list     !return all the keys
     procedure,pass :: types        => types_operators_list     !return all the types
     procedure,pass :: has_key      => has_key_operators_list  !True if key exists
     procedure,pass :: is_valid     => is_valid_operators_list !True if operators_list is valid
     procedure,pass :: has_valid_dq => has_valid_dq_operators_list
     procedure,pass :: shape        => shape_operators_list
     procedure,pass :: write        => write_operators_list  !write operators list
     procedure,pass :: read         => read_operators_list  !read operators list
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
  function construct_from_single_operator(key,op,type,dq) result(self)
    type(operators_list)                     :: self
    character(len=*),intent(in)              :: key
    type(sparse_matrix),intent(in)           :: op
    character(len=*),intent(in),optional     :: type
    real(8),dimension(:),intent(in)          :: dq
    character(len=16)                        :: type_
    type_=''       ;if(present(type))type_=str(type)
    if(size(dq)==0)stop "construct_from_single_operator: size(dq) must be > 0"
    call self%free()
    allocate(self%root)
    call self%put(key,op,type_,dq)
  end function construct_from_single_operator


  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor: given a list of keys+operators
  !+------------------------------------------------------------------+
  function construct_from_array_operator(keys,ops,types,dqs) result(self)
    type(operators_list)                                       :: self
    character(len=*),intent(in),dimension(:)                   :: keys
    type(sparse_matrix),intent(in),dimension(size(keys))       :: ops
    character(len=*),intent(in),dimension(size(keys)),optional :: types
    real(8),dimension(:,:),intent(in)                          :: dqs !shape=(Nqn,size(keys))
    character(len=16),dimension(size(keys))                    :: types_
    integer                                                    :: i
    !
    types_=''
    if(present(types))then
       do i=1,size(keys)
          types_(i)=str(types(i))
       enddo
    endif
    !
    if(size(dqs,1)==0)stop "construct_from_array_operator: size(dqs,1) must be > 0"
    if(size(dqs,2)/=size(keys))stop "construct_from_array_operator: size(dqs,2) must equal size(keys)"
    !
    call self%free()
    allocate(self%root)
    do i=1,size(keys)
       call self%put(keys(i),ops(i),types_(i),dq=dqs(:,i))
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
       if(allocated(c%dq))deallocate(c%dq)
       deallocate(c)
    enddo
    self%size=0
    deallocate(self%root)
    nullify(self%root,p,c)
  end subroutine free_operators_list





  !##################################################################
  !##################################################################
  !       PUT/LOAD/APPEND  - GET/DUMP OPERATORS IN/FROM A LIST
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Put a sparse matrix as operator in the operators_list
  !+------------------------------------------------------------------+
  subroutine put_operators_list(self,key,op,type,dq)
    class(operators_list),intent(inout)      :: self
    character(len=*),intent(in)              :: key
    type(sparse_matrix),intent(in)           :: op
    character(len=*),intent(in),optional     :: type
    real(8),dimension(:),intent(in)          :: dq
    type(optype),pointer                     :: p,c
    logical                                  :: iadd
    character(len=16)                        :: type_
    !
    type_=''       ;if(present(type))type_=str(type)
    if(size(dq)==0)stop "put_operators_list: dq must have size > 0"
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
       c%ctype        = str(type_)
       c%ope          = op
       c%dq           = dq
    else                        !QN does not exist: create a new element
       allocate(p%next)
       p%next%ckey    = str(key)
       p%next%ctype   = str(type_)
       p%next%index   = p%index+1
       p%next%ope     = op
       p%next%dq      = dq
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
  !PURPOSE: Update the matrix of an existing operator while preserving
  !         its key, type and dq metadata.
  !+------------------------------------------------------------------+
  subroutine update_operators_list(self,key,op)
    class(operators_list),intent(inout) :: self
    character(len=*),intent(in)         :: key
    type(sparse_matrix),intent(in)      :: op
    type(optype),pointer                :: c
    !
    if(.not.associated(self%root))&
         stop "update_operators_list: empty operators list"
    c=>self%root%next
    do
       if(.not.associated(c))&
            stop "update_operators_list: key not found: "//str(key)
       if(str(c%ckey)==str(key))then
          c%ope=op
          exit
       endif
       c=>c%next
    enddo
    nullify(c)
  end subroutine update_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE:  Append == Put a sparse matrix as operator in the operators_list
  !+------------------------------------------------------------------+
  subroutine append_operators_list(self,key,op,type,dq)
    class(operators_list),intent(inout)      :: self
    character(len=*),intent(in)              :: key
    type(sparse_matrix),intent(in)           :: op
    character(len=*),intent(in),optional     :: type
    real(8),intent(in),dimension(:)          :: dq
    character(len=16)                        :: type_
    type_='';if(present(type))type_=str(type)
    if(size(dq)==0)stop "append_operators_list: dq must have size > 0"
    call self%put(str(key),op,str(type_),dq)
  end subroutine append_operators_list




  !+------------------------------------------------------------------+
  !PURPOSE:  Load a dense matrix as operator in the operators_list
  !+------------------------------------------------------------------+
  subroutine load_operators_list(self,key,op,type,dq)
    class(operators_list),intent(inout)      :: self
    character(len=*),intent(in)              :: key
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in)     :: op
#else
    real(8),dimension(:,:),intent(in)        :: op
#endif
    character(len=*),intent(in),optional     :: type
    real(8),intent(in),dimension(:)          :: dq
    character(len=16)                        :: type_
    type_='';if(present(type))type_=str(type)    
    if(.not.associated(self%root))allocate(self%root)
    if(size(dq)==0)stop "load_operators_list: dq must have size > 0"
    call self%put(key,as_sparse(op),type_,dq)
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
  !              RETRIEVE CONTENT: OP, KEY, DQ
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Get {key,operator,type} of the list given an index
  !+------------------------------------------------------------------+
  subroutine get_all_operators_list(self,index,key,op,type,dq)
    class(operators_list),intent(inout)                   :: self
    integer,intent(in)                                    :: index
    character(len=*),intent(out)                          :: key
    type(sparse_matrix),intent(out)                       :: op
    character(len=*),intent(out),optional                 :: type
    real(8),intent(out),dimension(:),allocatable,optional :: dq
    integer                                               :: index_
    type(optype),pointer                                  :: c
    logical                                               :: ifound
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
    if(present(dq))dq=c%dq
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
    if(.not.present(index).AND..not.present(key))&
    stop "get_op_operators_list: no input given: use index=i OR key=str"
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
    !
    if(.not.ifound)then
       if(present(key))then
          stop "get_op_operators_list error: key not found: "//str(key)
       else
          stop "get_op_operators_list error: index not found"
       endif
    endif
    !
    op = c%ope    
    !
    c=>null()
  end function get_op_operators_list



  !+------------------------------------------------------------------+
  !PURPOSE: Return dq of the operator in the list as sparse matrix given:
  ! + key: the operator corresponding to the key value
  ! + indx: the operator corresponding to  the indx value
  !+------------------------------------------------------------------+
  function get_dq_operators_list(self,key,index) result(dq)
    class(operators_list)                :: self
    character(len=*),intent(in),optional :: key
    integer,intent(in),optional          :: index
    real(8),dimension(:),allocatable     :: dq
    integer                              :: index_
    type(optype),pointer                 :: c
    logical                              :: ifound
    !
    index_=self%size;if(present(index))index_=index
    if(index_>self%size.OR.index_<=0)stop "get_dq_operators_list: index !in [1,self.size]"
    if(.not.present(index).AND..not.present(key))&
    stop "get_dq_operators_list: no input given: use index=i OR key=str"
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
    !
    if(.not.ifound)then
       if(present(key))then
          stop "get_dq_operators_list error: key not found: "//str(key)
       else
          stop "get_dq_operators_list error: index not found"
       endif
    else
      if(.not.allocated(c%dq))stop "get_dq_operators_list error: dq not allocated for key: "//str(key)
    endif
    !
    allocate(dq, source=c%dq)
    !
    c=>null()
  end function get_dq_operators_list





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
    class(operators_list)                :: self
    integer,intent(in),optional          :: index
    character(len=*),intent(in),optional :: key
    character(len=:),allocatable :: type
    integer                      :: index_
    type(optype),pointer         :: c
    logical                      :: ifound
    !
    index_=self%size;if(present(index))index_=index
    if(index_>self%size.OR.index_<=0)stop "get_op_operators_list: index !in [1,self.size]"
    if(.not.present(index).AND..not.present(key))&
         stop "get_type_operators_list: no input given: use index=i OR key=str"
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
    if(.not.ifound)then
       if(present(key))then
          stop "get_type_operators_list error: key not found: "//str(key)
       else
          stop "get_type_operators_list error: index not found"
       endif
    endif
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
  !PURPOSE: Return all the types in the operators_list
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
  function is_valid_operators_list(self,dim,qdim) result(bool)
    class(operators_list),intent(in) :: self
    integer,intent(in),optional      :: dim
    integer,intent(in),optional      :: qdim
    integer                          :: dim_
    logical                          :: bool
    type(optype),pointer             :: c
    bool = .true.
    if(.not.associated(self%root))return
    dim_ = 0;if(present(dim))dim_=dim
    if(present(qdim))then
       if(qdim<=0)then
          bool=.false.
          return
       endif
    endif
    c => self%root%next
    do 
       if(.not.associated(c))exit
       if(dim_==0)then
          dim_ = c%ope%Nrow
       endif
       bool = bool.AND.(all([c%ope%Nrow,c%ope%Ncol] == [dim_,dim_]))
       if(.not.allocated(c%dq))then
          bool=.false.
       else
          bool=bool.AND.size(c%dq)>0
          if(present(qdim))bool=bool.AND.size(c%dq)==qdim
       endif
       c => c%next
    enddo
    c=>null()
  end function is_valid_operators_list


  !+------------------------------------------------------------------+
  !PURPOSE: Return True if every operator has a well-formed dq. If qdim
  ! is present, every dq must also have that size.
  !+------------------------------------------------------------------+
  function has_valid_dq_operators_list(self,qdim) result(bool)
    class(operators_list),intent(in) :: self
    integer,intent(in),optional      :: qdim
    logical                          :: bool
    integer                          :: qdim_
    type(optype),pointer             :: c
    !
    bool=.false.
    if(.not.associated(self%root))return
    if(self%size==0)return
    !
    qdim_=0
    c=>self%root%next
    do
       if(.not.associated(c))exit
       if(.not.allocated(c%dq))return
       if(size(c%dq)==0)return
       if(present(qdim))then
          if(size(c%dq)/=qdim)return
       else
          if(qdim_==0)qdim_=size(c%dq)
          if(size(c%dq)/=qdim_)return
       endif
       c=>c%next
    enddo
    bool=.true.
    c=>null()
  end function has_valid_dq_operators_list








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
    type(optype),pointer                :: c
    character(len=*),optional           :: file
    integer                             :: unit_
    unit_=6
    if(present(unit))unit_=unit
    if(present(file))open(free_unit(unit_),file=str(file))
    !
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    !
    write(unit_,"(A7,I12)")"Size :",self%size
    write(unit_,"(A18)")"------------------"
    c => self%root%next
    do
       if(.not.associated(c))exit
       write(unit_,"(A7,I12)")"Index: ",c%index
       write(unit_,"(A7,A)")"Key  : ",str(c%ckey)
       write(unit_,"(A7,A)")"Type : ",str(c%ctype)
       write(unit_,*)"dq   : ",c%dq
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
    type(optype),pointer                :: c
    integer                             :: unit_
    !
    unit_=-1
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    if(unit_==-1)stop "write_operators_list error: no input +file or +unit given"
    !
    if(.not.self%is_valid())stop "write_operators_list error: invalid operators list"
    !
    ! Format 3 requires dq metadata for every operator.
    write(unit_,"(A,1X,I0)")operators_list_header,operators_list_format
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
       write(unit_,*)size(c%dq)
       write(unit_,*)c%dq
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
    integer                             :: i,ListSize,format_version,ndq,ios
    integer                             :: unit_
    character(len=512)                  :: header_line,header
    character(len=512)                  :: key
    character(len=512)                  :: type
    type(sparse_matrix)                 :: ope
    real(8),dimension(:),allocatable    :: dq
    !
    unit_=-1
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    if(unit_==-1)stop "read_operators_list error: no input +file or +unit given"
    !
    read(unit_,"(A)",iostat=ios)header_line
    if(ios/=0)stop "read_operators_list error: unable to read format header"
    read(header_line,*,iostat=ios)header,format_version
    if(ios/=0)then
       write(*,"(A)")"read_operators_list error: legacy operators-list format detected."
       write(*,"(A)")"This file is readable only with nssDMRG version <= 5.1.2."
       stop "read_operators_list error: incompatible file format"
    endif
    if(str(header)/=operators_list_header)then
       write(*,"(A)")"read_operators_list error: legacy operators-list format detected."
       write(*,"(A)")"This file is readable only with nssDMRG version <= 5.1.2."
       stop "read_operators_list error: incompatible file format"
    endif
    if(format_version/=operators_list_format)then
       write(*,"(A,I0,A,I0)")"read_operators_list error: file format ",format_version,&
            " is incompatible with supported format ",operators_list_format
       stop "read_operators_list error: unsupported operators-list format"
    endif
    !
    call self%free()
    read(unit_,*,iostat=ios)ListSize
    if(ios/=0)stop "read_operators_list error: unable to read list size"
    if(ListSize<0)stop "read_operators_list error: invalid list size"
    do i=1,ListSize
       read(unit_,*,iostat=ios)key
       if(ios/=0)stop "read_operators_list error: unable to read operator key"
       read(unit_,*,iostat=ios)type
       if(ios/=0)stop "read_operators_list error: unable to read operator type"
       if(str(type)=="none")type=""
       if(allocated(dq))deallocate(dq)
       read(unit_,*,iostat=ios)ndq
       if(ios/=0)stop "read_operators_list error: unable to read dq size"
       if(ndq<=0)stop "read_operators_list error: invalid dq size"
       allocate(dq(ndq))
       read(unit_,*,iostat=ios)dq
       if(ios/=0)stop "read_operators_list error: unable to read dq"
       call ope%read(unit=unit_)
       call self%append(key=str(key),op=ope,type=str(type),dq=dq)
    end do
    if(allocated(dq))deallocate(dq)
    call ope%free()
    if(present(file))close(unit_)
  end subroutine read_operators_list





  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Equality between two operators_lists (full copy)
  !+------------------------------------------------------------------+
  subroutine equality_operators_list(A,B)
    type(operators_list),intent(inout) :: A
    type(operators_list),intent(in)    :: B
    integer                            :: i
    call A%free()
    do i=1,size(B)
       call A%put(B%key(index=i),B%op(index=i),B%type(index=i),dq=B%dq(index=i))
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
  real(8),dimension(2,5)                 :: dqs
  real(8),dimension(:),allocatable       :: dq_get
  integer,parameter                     :: sec=1


  Gamma13=kron(Sx,Sz)
  Gamma03=kron(S0,Sz)


  print*,"TEST DQ CONSTRUCTOR, ACCESSORS, GET AND OVERWRITE"
  dqs(:,1)=[ 0d0, 0d0]
  dqs(:,2)=[ 0d0, 0d0]
  dqs(:,3)=[ 1d0,-1d0]
  dqs(:,4)=[ 0d0, 0d0]
  dqs(:,5)=[-1d0, 0d0]
  my_list = operators_list(&
       ['H0','Sz','Sp','P ','C '],&
       [sparse(Hzero),sparse(Sz),sparse(Splus),sparse(Sz),sparse(Sz)],&
       ['bose ','bose ','bose ','sign ','fermi'],&
       dqs=dqs)
  call assert_true(my_list%is_valid(dim=2),"dq constructor: algebraic validity")
  call assert_true(my_list%has_valid_dq(qdim=2),"dq constructor: complete metadata")
  call assert_true(my_list%is_valid(dim=2,qdim=2),"dq constructor: complete validity")
  call assert_dq(my_list%dq(key="Sp"),[1d0,-1d0],"dq accessor by key")
  call assert_dq(my_list%dq(index=5),[-1d0,0d0],"dq accessor by index")
  call my_list%get(index=3,key=key,op=a,type=type,dq=dq_get)
  call assert_true(allocated(dq_get),"get: dq is allocated")
  call assert_dq(dq_get,[1d0,-1d0],"get: dq value")
  call my_list%get(index=4,key=key,op=a,type=type,dq=dq_get)
  call assert_dq(dq_get,[0d0,0d0],"get: P dq")
  call my_list%put("Sp",sparse(Splus),"bose",dq=[1d0,-1d0])
  call assert_dq(my_list%dq(key="Sp"),[1d0,-1d0],"put overwrite sets dq value")
  call my_list%update("Sp",sparse(Splus))
  call assert_dq(my_list%dq(key="Sp"),[1d0,-1d0],"update preserves dq")
  call assert_true(str(my_list%type(key="Sp"))=="bose","update preserves type")
  copy_list=my_list
  call assert_dq(copy_list%dq(key="Sp"),[1d0,-1d0],"deep copy: dq value")
  call assert_dq(copy_list%dq(key="P"),[0d0,0d0],"deep copy: P dq")
  call my_list%put("Sp",sparse(Splus),"bose",dq=[1d0,-1d0])
  call assert_dq(copy_list%dq(key="Sp"),[1d0,-1d0],"deep copy: independent metadata")
  a_list=operators_list("Sp",sparse(Splus),"bose",dq=[1d0,-1d0])
  call assert_true(a_list%has_valid_dq(qdim=2),"single constructor: valid dq")
  call assert_dq(a_list%dq(key="Sp"),[1d0,-1d0],"single constructor: dq value")
  call my_list%free()
  call copy_list%free()
  call a_list%free()
  print*,"DQ CORE TESTS: PASS"
  print*,""


  print*,"TEST CONSTRUCTOR, PUT, SHOW, FREE"
  my_list = operators_list(&
       ['H0','Sz','Sp','P ','C '],&
       [sparse(Hzero),sparse(Sz),sparse(Splus),sparse(Sz),sparse(Sz)],&
       ['bose ','bose ','bose ','sign ','fermi'],dqs=dqs)
  call my_list%show()
  call my_list%free()
  call wait(sec)


  print*,"TEST LOAD matrices"
  call my_list%load("H0",Hzero,'b',dq=[0d0,0d0])
  call my_list%load("Sz",Sz,'s',dq=[0d0,0d0])
  call my_list%load("Sp",Splus,'fermi',dq=[1d0,-1d0])
  print*,"TEST SHOW"
  call my_list%show()
  call my_list%free()
  call wait(sec)



  print*,"TEST (CONSTRUCT + )APPEND matrices"
  call my_list%append("H0",as_sparse(Hzero),'b',dq=[0d0,0d0])
  call my_list%append("Sz",as_sparse(Sz),'b',dq=[0d0,0d0])
  call my_list%append("Sp",as_sparse(Splus),'Bosonic',dq=[1d0,-1d0])
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
  call my_list%append("W",as_sparse(Gamma03),'b',dq=[0d0,0d0])
  print*,my_list%is_valid()
  call my_list%free
  print*,""
  call wait(sec)



  call my_list%append("H0",as_sparse(Hzero),'b',dq=[0d0,0d0])
  call my_list%append("Sz",as_sparse(Sz),'b',dq=[0d0,0d0])
  call my_list%load("Sp",Splus,'b',dq=[1d0,-1d0])



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
  call a_list%append("gamma13",as_sparse(Gamma13),'b',dq=[1d0,-1d0])
  call a_list%append("gamma03",as_sparse(Gamma03),'b',dq=[0d0,0d0])
  call a_list%append("Gamma33",as_sparse(kron(Sz,Sz)),'b',dq=[0d0,0d0])

  clist(1) = my_list
  clist(2) = a_list

  call clist(1)%show()
  call clist(2)%show()
  print*,""
  call wait(sec)

  print*,"TEST WRITE/READ"
  call a_list%free()
  call a_list%append("gamma13",as_sparse(Gamma13),'b',dq=[1d0,-1d0])
  call a_list%append("gamma03",as_sparse(Gamma03),'b',dq=[0d0,0d0])
  call a_list%append("Gamma33",as_sparse(kron(Sz,Sz)),'b',dq=[0d0,0d0])

  print*,"write:"
  call a_list%write(file="a_list_write.dat")

  print*,"show:"
  call a_list%show()

  print*,"free:"
  call a_list%free()

  print*,"read:"
  call a_list%read(file="a_list_write.dat")

  call assert_true(size(a_list)==3,"write/read: list size")
  call assert_dq(a_list%dq(key="gamma13"),[1d0,-1d0],"write/read: gamma13 dq")
  call assert_dq(a_list%dq(key="gamma03"),[0d0,0d0],"write/read: gamma03 dq")
  call assert_dq(a_list%dq(key="Gamma33"),[0d0,0d0],"write/read: Gamma33 dq")
  call assert_true(a_list%has_valid_dq(qdim=2),"write/read: complete dq metadata")

  print*,"show again:"
  call a_list%show()

  print*,"DQ WRITE/READ TESTS: PASS"

contains

  subroutine assert_true(condition,message)
    logical,intent(in)          :: condition
    character(len=*),intent(in) :: message
    if(.not.condition)then
       write(*,"(A)")"FAILED: "//trim(message)
       error stop 1
    endif
  end subroutine assert_true


  subroutine assert_dq(actual,expected,message)
    real(8),dimension(:),intent(in) :: actual,expected
    character(len=*),intent(in)     :: message
    call assert_true(size(actual)==size(expected),trim(message)//": size")
    call assert_true(all(abs(actual-expected)<1d-12),trim(message)//": value")
  end subroutine assert_dq

end program testOPERATORS_TUPLE
#endif
