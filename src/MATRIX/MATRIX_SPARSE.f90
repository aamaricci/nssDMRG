MODULE MATRIX_SPARSE  
  USE SCIFOR, only: str,free_unit,assert_shape,zeye,eye
  USE AUX_FUNCS, only: show_fmt,append,fopen
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private


#ifdef _CMPLX
  complex(8)         :: zero=cmplx(0d0,0d0)
#else
  real(8)            :: zero=0d0
#endif
#ifdef _MPI
#ifdef _CMPLX  
  integer, parameter :: MPI_VAL_TYPE = MPI_DOUBLE_COMPLEX
#else
  integer, parameter :: MPI_VAL_TYPE = MPI_DOUBLE_PRECISION
#endif
#endif
  integer            :: i,j



  !SPARSE ROW OF THE SPARSE MATRIX: note this is dynamic array
  type sparse_row
     sequence
     integer                             :: size
     integer,dimension(:),allocatable    :: cols
#ifdef _CMPLX
     complex(8),dimension(:),allocatable :: vals
#else
     real(8),dimension(:),allocatable    :: vals
#endif
  end type sparse_row


  !SPARSE MATRIX STRUCTURE
  type sparse_matrix
     type(sparse_row),dimension(:),allocatable :: row
     integer                                   :: Nrow=0
     integer                                   :: Ncol=0
     logical                                   :: status=.false.
   contains
     procedure,pass :: init       => sp_init_matrix
     procedure,pass :: free       => sp_free_matrix
     procedure,pass :: load       => sp_load_matrix
     procedure,pass :: copy       => sp_copy_matrix
     procedure,pass :: dump       => sp_dump_matrix
     procedure,pass :: fast_insert=> sp_fast_insert_element
     procedure,pass :: insert     => sp_insert_element
     procedure,pass :: get        => sp_get_element
     procedure,pass :: show       => sp_show_matrix
     procedure,pass :: display    => sp_display_matrix
     procedure,pass :: spy        => sp_spy_matrix
     procedure,pass :: as_matrix  => sp_as_matrix
     procedure,pass :: dgr        => sp_dgr_matrix
     procedure,pass :: t          => sp_transpose_matrix
     procedure,pass :: nnz        => sp_nnz_matrix
     procedure,pass :: dot        => sp_matmul_vector
     procedure,pass :: write      => sp_write_matrix
     procedure,pass :: read       => sp_read_matrix
#ifdef _MPI
     procedure,pass :: bcast      => sp_bcast_matrix
     procedure,pass :: pdot       => sp_p_matmul_matrix
     procedure,pass :: bytes      => sp_get_bytes_matrix
#endif
  end type sparse_matrix


  interface sparse_matrix
     module procedure :: sp_construct_matrix
  end interface sparse_matrix

  interface as_sparse
     module procedure :: sp_construct_matrix
  end interface as_sparse

  interface as_matrix
     module procedure :: sp_as_matrix
  end interface as_matrix

  interface sparse
     module procedure :: sp_construct_matrix
  end interface sparse

  !EQUALITY with scalar and function (A=B, A=cmplx)
  interface assignment(=)
     module procedure :: sp_matrix_equal_scalar
     module procedure :: sp_matrix_equal_matrix
  end interface assignment(=)

  !ADDITION
  interface operator (+)
     module procedure :: sp_plus_matrix
  end interface operator (+)

  !SUBTRACTION
  interface operator (-)
     module procedure :: sp_minus_matrix
  end interface operator (-)

  !SCALAR PRODUCT
  interface operator(*)
     module procedure :: sp_left_product_matrix_i
     module procedure :: sp_left_product_matrix_d
#ifdef _CMPLX
     module procedure :: sp_left_product_matrix_c
#endif
     !
     module procedure :: sp_right_product_matrix_i
     module procedure :: sp_right_product_matrix_d
#ifdef _CMPLX
     module procedure :: sp_right_product_matrix_c
#endif
  end interface operator(*)

  !SCALAR DIVISION
  interface operator(/)
     module procedure :: sp_right_division_matrix_i
     module procedure :: sp_right_division_matrix_d
#ifdef _CMPLX
     module procedure :: sp_right_division_matrix_c
#endif
  end interface operator(/)


  !Matrix-Matrix PRODUCT
  interface operator(.m.)
     module procedure :: sp_matmul_matrix
  end interface operator(.m.)

  !Parallel Matrix-Matric product
#ifdef _MPI
  interface operator(.pm.)
     module procedure :: sp_p_matmul_matrix
  end interface operator(.pm.)
#endif

  !KRONECKER PRODUCT
  interface operator(.x.)
     module procedure :: sp_kron_matrix
  end interface operator(.x.)

  interface sp_kron
     module procedure :: sp_restricted_kron_matrix
  end interface sp_kron


  !RETURN SHAPE OF THE SPARSE MATRIX
  intrinsic :: shape
  interface shape
     module procedure :: sp_shape_matrix
  end interface shape

  !EXTEND TRANSPOSE TO SPARSE MATRIX
  intrinsic :: transpose
  interface transpose
     module procedure :: sp_transpose_matrix
  end interface transpose

  !EXTEND TRANSPOSE TO SPARSE MATRIX
  interface hconjg
     module procedure :: sp_hconjg_matrix
  end interface hconjg

  intrinsic :: matmul
  interface matmul
     module procedure :: sp_matmul_matrix
     module procedure :: sp_matmul_vector
  end interface matmul

  interface sp_filter
     module procedure :: sp_filter_matrix_1
     module procedure :: sp_filter_matrix_2
  end interface sp_filter

  interface sp_add3
     module procedure :: sp_plus3_matrix
  end interface sp_add3

#ifdef _MPI
  interface  AllGather_MPI
     module procedure ::  sp_array_all_gather_1
     module procedure ::  sp_array_all_gather_2
  end interface AllGather_MPI
#endif


  public :: sparse_matrix
  public :: as_sparse
  public :: as_matrix
  public :: sparse
  public :: assignment(=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.m.)
  public :: operator(.x.)
  public :: sp_kron
  public :: shape
  public :: transpose
  public :: hconjg
  public :: matmul
  public :: sp_eye
  public :: sp_filter
  public :: sp_add3
#ifdef _MPI
  public :: operator(.pm.)
  public :: AllGather_MPI
#endif




contains       







  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  elemental subroutine sp_init_matrix(sparse,Nrow,Ncol)
    class(sparse_matrix),intent(inout) :: sparse
    integer,intent(in)                 :: Nrow,Ncol
    integer                            :: i
    !
    if(sparse%status)call sparse%free()
    ! 
    sparse%Nrow=Nrow
    sparse%Ncol=Ncol
    !
    if(allocated(sparse%row))deallocate(sparse%row)
    allocate(sparse%row(Nrow))
    do i=1,Nrow
       sparse%row(i)%size=0
       allocate(sparse%row(i)%vals(0)) !empty array
       allocate(sparse%row(i)%cols(0)) !empty array
    end do
    sparse%status=.true.
  end subroutine sp_init_matrix


  function sp_construct_matrix(matrix) result(self)
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in) :: matrix
#else
    real(8),dimension(:,:),intent(in)    :: matrix
#endif
    type(sparse_matrix)                  :: self
    call self%load(matrix)
  end function sp_construct_matrix




  subroutine sp_copy_matrix(sparse,b)
    class(sparse_matrix),intent(inout) :: sparse
    type(sparse_matrix),intent(in)     :: b
    integer                            :: col
#ifdef _CMPLX
    complex(8)                         :: val
#else
    real(8)                            :: val
#endif
    if(.not.b%status)stop "sp_matrix_equal_matrix: B.status=F"
    call sparse%free()
    call sparse%init(b%Nrow,b%Ncol)
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val = b%row(i)%vals(j)
          call sparse%insert(val,i,col)          
       enddo
    enddo
  end subroutine sp_copy_matrix



  !+------------------------------------------------------------------+
  !PURPOSE: free an entire sparse matrix
  !+------------------------------------------------------------------+
  elemental subroutine sp_free_matrix(sparse)    
    class(sparse_matrix),intent(inout) :: sparse
    integer                            :: i
    !
    do i=1,sparse%Nrow
       sparse%row(i)%Size  = 0
       if(allocated(sparse%row(i)%vals))deallocate(sparse%row(i)%vals)
       if(allocated(sparse%row(i)%cols))deallocate(sparse%row(i)%cols)
    enddo
    if(allocated(sparse%row))deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
  end subroutine sp_free_matrix




  !+------------------------------------------------------------------+
  !PURPOSE: load a dense matrix into a sparse one
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix(sparse,matrix)
    class(sparse_matrix),intent(inout)   :: sparse
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in) :: matrix
#else
    real(8),dimension(:,:),intent(in)    :: matrix
#endif
    integer                              :: Ndim1,Ndim2
    !
    call sparse%free()
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    call sparse%init(Ndim1,Ndim2) !set status=T
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=zero)call sp_insert_element(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix(sparse,matrix)
    class(sparse_matrix),intent(in)         :: sparse
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(inout) :: matrix
#else
    real(8),dimension(:,:),intent(inout)    :: matrix
#endif
    integer                                 :: Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    call assert_shape(matrix,[sparse%Nrow,sparse%Ncol],"sp_dump_matrix","Matrix")
    Matrix = zero
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%vals(j)
       enddo
    enddo
  end subroutine sp_dump_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  function sp_as_matrix(sparse) result(matrix)
    class(sparse_matrix),intent(in)               :: sparse
#ifdef _CMPLX
    complex(8),dimension(sparse%Nrow,sparse%Ncol) :: matrix
#else
    real(8),dimension(sparse%Nrow,sparse%Ncol)    :: matrix
#endif
    matrix = zero
    if(.not.sparse%status)return
    do i=1,sparse%Nrow
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = sparse%row(i)%vals(j)
       enddo
    enddo
  end function sp_as_matrix


  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element(sparse,value,i,j)
    class(sparse_matrix),intent(inout) :: sparse
#ifdef _CMPLX
    complex(8),intent(in)              :: value
#else
    real(8),intent(in)                 :: value
#endif
    integer,intent(in)                 :: i,j
    integer                            :: column,pos
    logical                            :: iadd
    !
    if(.not.sparse%status)stop "sp_insert_element: sparse.status=F"
    column = j
    !
    iadd = .false.                          !check if column already exist
    if(any(sparse%row(i)%cols == column))then         !
       pos = binary_search(sparse%row(i)%cols,column) !find the position  column in %cols
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       sparse%row(i)%vals(pos)=sparse%row(i)%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it 
       call append(sparse%row(i)%vals,value)
       call append(sparse%row(i)%cols,column)
       sparse%row(i)%Size = sparse%row(i)%Size + 1
    endif
    !
    if(sparse%row(i)%Size > sparse%Ncol)stop "sp_insert_element ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_insert_element

  !no addition, no check... wild cow-boy 
  subroutine sp_fast_insert_element(sparse,value,i,j)
    class(sparse_matrix),intent(inout) :: sparse
#ifdef _CMPLX
    complex(8),intent(in)              :: value
#else
    real(8),intent(in)                 :: value
#endif
    integer,intent(in)                 :: i,j
    integer                            :: column,pos
    !
    column = j
    !
    call append(sparse%row(i)%vals,value)
    call append(sparse%row(i)%cols,column)
    sparse%row(i)%Size = sparse%row(i)%Size + 1
    !
    if(sparse%row(i)%Size > sparse%Ncol)stop "sp_insert_element ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_fast_insert_element


  !+------------------------------------------------------------------+
  !PURPOSE: get the element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  function sp_get_element(sparse,i,j) result(val)
    class(sparse_matrix),intent(in) :: sparse    
    integer,intent(in)              :: i,j
#ifdef _CMPLX
    complex(8)                      :: val
#else
    real(8)                         :: val
#endif
    integer                         :: pos
    if(.not.sparse%status)stop "sp_get_element: sparse.status=F"
    val=zero
    if(.not.any(sparse%row(i)%cols==j))return
    pos=binary_search(sparse%row(i)%cols,j)
    val=sparse%row(i)%vals(pos)
  end function sp_get_element


  !+------------------------------------------------------------------+
  !PURPOSE: return the number of non-zero elements
  !+------------------------------------------------------------------+
  function sp_nnz_matrix(sparse) result(nnz)
    class(sparse_matrix),intent(in) :: sparse    
    integer                         :: nnz
    integer                         :: i
    nnz=0
    if(.not.sparse%status)return
    do i=1,sparse%Nrow
       nnz=nnz + sparse%row(i)%size
    enddo
  end function sp_nnz_matrix




  !+------------------------------------------------------------------+
  !PURPOSE: Filter the sparse matrix along given Row-states or Row- and
  ! Col-states
  !+------------------------------------------------------------------+


  function sp_filter_matrix_1(A,states) result(Ak)
    class(sparse_matrix), intent(in)    :: A
    integer,dimension(:),intent(in)     :: states
    type(sparse_matrix)                 :: Ak
    integer                             :: Nstates, Nrow, M
    integer, dimension(:), allocatable  :: local_map
    integer, dimension(:), allocatable  :: nnz_count
    integer                             :: i,j,jj,istate,jstate
#ifdef _CMPLX
    complex(8)                          :: val
#else
    real(8)                             :: val
#endif
    !
    !
    if(.not.A%status)stop "sp_filter_matrix_1: A.status=F"
    !
    Nstates = size(states)
    Nrow = A%Nrow
    !
    call Ak%free()
    call Ak%init(Nstates,Nstates)
    !
    !
    ! A Map to invert the states arrray.
    ! This is used later to lookup with O(1) access non-zero elements on a
    ! given row
    allocate(local_map(Nrow))
    local_map = 0
    do istate = 1, Nstates
       local_map(states(istate)) = istate
    enddo
    !
    ! Count how many nnz elements of A will be present in A(k) (a block of A)
    allocate(nnz_count(Nstates))
    nnz_count = 0
    do istate = 1, Nstates
       i = states(istate)
       do jj = 1, A%row(i)%Size
          j  = A%row(i)%cols(jj)
          ! Controlla se anche la colonna appartiene agli stati filtrati
          if (local_map(j) > 0) nnz_count(istate) = nnz_count(istate) + 1
       enddo
    enddo
    !
    !Re-allocate vals and cols of Ak to nnz value so we avoid call to append:
    do istate = 1, Nstates
       M = nnz_count(istate)
       if (M == 0) cycle
       Ak%row(istate)%Size = M
       if(allocated(Ak%row(istate)%vals))deallocate(Ak%row(istate)%vals)
       if(allocated(Ak%row(istate)%cols))deallocate(Ak%row(istate)%cols)
       allocate(Ak%row(istate)%vals(M))
       allocate(Ak%row(istate)%cols(M))
    enddo
    !
    ! Build  Ak 
    nnz_count = 0 ! Re-use nnz_count
    do istate = 1, Nstates
       i = states(istate)
       do jj = 1, A%row(i)%Size
          j  = A%row(i)%cols(jj)
          jstate = local_map(j)
          if(jstate==0)cycle
          val = A%row(i)%vals(jj)
          nnz_count(istate) = nnz_count(istate) + 1
          Ak%row(istate)%vals(nnz_count(istate)) = val
          Ak%row(istate)%cols(nnz_count(istate)) = jstate
       enddo
    enddo
  end function sp_filter_matrix_1


  function sp_filter_matrix_2(A,Istates,Jstates) result(Ak)
    class(sparse_matrix), intent(in)    :: A
    integer,dimension(:),intent(in)     :: Istates,Jstates
    type(sparse_matrix)                 :: Ak
    integer                             :: NIstates, NJstates, Nrow, Ncol, M
    integer, dimension(:), allocatable  :: local_map !this is now for Jstates
    integer, dimension(:), allocatable  :: nnz_count
    integer                             :: i,j,ii,jj,istate,jstate
#ifdef _CMPLX
    complex(8)                          :: val
#else
    real(8)                             :: val
#endif
    !
    ! if(.not.A%status)stop "sp_filter_matrix_2: A.status=F"
    !
    NIstates = size(Istates)
    NJstates = size(Jstates)
    Ncol     = A%Ncol
    !
    call Ak%free()
    call Ak%init(NIstates,NJstates)
    !
    ! A Map to invert the states arrray.
    ! This is used later to lookup with O(1) access non-zero elements on a
    ! given row
    allocate(local_map(Ncol))
    local_map = 0
    do jstate = 1, NJstates
       local_map(Jstates(jstate)) = jstate
    enddo
    !
    ! Count how many nnz elements of A will be present in A(k) (a block of A)
    allocate(nnz_count(NIstates))
    nnz_count = 0
    do istate = 1, NIstates
       i = Istates(istate)
       do jj = 1, A%row(i)%Size
          j  = A%row(i)%cols(jj)
          ! Controlla se l'indice di colonna è valido e se appartiene a Jstates
          if (j > Ncol) cycle
          if (local_map(j) > 0)nnz_count(istate) = nnz_count(istate) + 1
       enddo
    enddo
    !
    !Re-allocate vals and cols of Ak to nnz value so we avoid call to append:
    do istate = 1, NIstates
       M = nnz_count(istate)
       if (M == 0) cycle
       Ak%row(istate)%Size = M
       if(allocated(Ak%row(istate)%vals))deallocate(Ak%row(istate)%vals)
       if(allocated(Ak%row(istate)%cols))deallocate(Ak%row(istate)%cols)
       allocate(Ak%row(istate)%vals(M))
       allocate(Ak%row(istate)%cols(M))
    enddo
    !
    ! Build  Ak 
    nnz_count = 0 ! Re-use nnz_count
    do istate = 1, NIstates
       i = Istates(istate)
       do jj = 1, A%row(i)%Size
          j     = A%row(i)%cols(jj)
          if(j > Ncol) cycle        
          jstate = local_map(j)
          if(jstate==0)cycle
          val = A%row(i)%vals(jj)
          nnz_count(istate) = nnz_count(istate) + 1
          Ak%row(istate)%vals(nnz_count(istate)) = val
          Ak%row(istate)%cols(nnz_count(istate)) = jstate
       enddo
    enddo
  end function sp_filter_matrix_2






  !##################################################################
  !##################################################################
  !              SHOW  / SPY
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: show matrix
  !+------------------------------------------------------------------+
  subroutine sp_show_matrix(sparse,fmt,file)
    class(sparse_matrix)      :: sparse
    character(len=*),optional :: fmt,file
    character(len=12)         :: fmt_
    integer                   :: unit_
    integer                   :: Ns,NN
    character(len=64)         :: format
#ifdef _CMPLX
    complex(8)                :: val
#else
    real(8)                   :: val
#endif
    unit_=6
    if(present(file))open(free_unit(unit_),file=str(file))
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
#ifdef _CMPLX
    format='(A1,'//str(fmt_)//',A1,'//str(fmt_)//',A1,1x)'
#else
    format='(A1,'//str(fmt_)//',1x)'
#endif
    if(.not.sparse%status)then
       write(*,*)"sparse.status=F: nothing to show"
       return
    endif
    Ns=sparse%Nrow
    do i=1,sparse%Nrow
       do j=1,sparse%Ncol
          val = sp_get_element(sparse,i,j)
#ifdef _CMPLX
          write(unit_,"("//str(sparse%Ncol)//str(format)//")",advance='no')"(",dreal(val),",",dimag(val),")"
#else
          write(unit_,"("//str(sparse%Ncol)//"F5.1)",advance='no')val
#endif
       enddo
       write(unit_,*)
    enddo
    if(present(file))close(unit_)
  end subroutine sp_show_matrix

  subroutine sp_display_matrix(sparse)
    class(sparse_matrix) :: sparse
    integer              :: unit_
    integer              :: Ns
    character(len=20)    :: fmtR,fmtI
#ifdef _CMPLX
    complex(8)           :: val
#else
    real(8)              :: val
#endif
    unit_=6
    fmtI='(I10)'
    fmtR='(ES10.3)'
    if(.not.sparse%status)return
    do i=1,sparse%Nrow
       Ns = size(sparse%row(i)%cols)
       if(Ns==0)cycle
       do j=1,Ns
          write(unit_,"(A1,2I5,A1,2F8.3)",advance='no')"[",i,sparse%row(i)%cols(j),"]",sparse%row(i)%vals(j)
       enddo
       write(unit_,"(A1)",advance='yes')""
       write(unit_,*)
    enddo
  end subroutine sp_display_matrix



  subroutine sp_write_matrix(sparse,file,unit)
    class(sparse_matrix)      :: sparse
    character(len=*),optional :: file
    integer,optional          :: unit
    integer                   :: unit_
    integer                   :: Ns
    !
    unit_=-1
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    if(unit_==-1)stop "sp_write_matrix error: no input +file or +unit given"
    !
    if(.not.sparse%status)return
    write(unit_,*)sparse%Nrow,sparse%Ncol
    do i=1,sparse%Nrow
       write(unit_,*)sparse%row(i)%size
       do j=1,sparse%row(i)%size
          write(unit_,*)sparse%row(i)%cols(j),sparse%row(i)%vals(j)
       enddo
    enddo
    if(present(file))close(unit_)
  end subroutine sp_write_matrix




  subroutine sp_read_matrix(sparse,file,unit)
    class(sparse_matrix),intent(inout) :: sparse
    character(len=*),optional          :: file
    integer,optional                   :: unit
    integer                            :: unit_
    integer                            :: Ns,Nrow,Ncol
    integer                            :: irow,i,nsize
    integer                            :: row,col
    integer,allocatable                :: cols(:)
#ifdef _CMPLX
    complex(8),allocatable             :: vals(:)
#else
    real(8),allocatable                :: vals(:)
#endif
    !
    unit_=-1
    if(present(file))open(free_unit(unit_),file=str(file))
    if(present(unit))unit_=unit
    if(unit_==-1)stop "sp_read_matrix error: no input +file or +unit given"
    !
    if(sparse%status)call sparse%free()
    read(unit_,*)Nrow,Ncol
    call sparse%init(Nrow,Ncol)
    allocate(vals(Ncol));vals=0
    allocate(cols(Ncol));cols=0
    do irow=1,Nrow
       read(unit_,*)nsize
       if(nsize==0)cycle
       do i=1,nsize
          read(unit_,*)cols(i),vals(i)
       enddo
       sparse%row(irow)%size=nsize
       call append(sparse%row(irow)%cols,cols(1:nsize))
       call append(sparse%row(irow)%vals,vals(1:nsize))
    enddo
    if(present(file))close(unit_)
    deallocate(vals,cols)
  end subroutine sp_read_matrix



  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+  
  subroutine sp_spy_matrix(sparse,header)
    class(sparse_matrix)    :: sparse
    character ( len = * )   :: header
    integer                 :: N1,N2
    character ( len = 255 ) :: command_filename
    integer                 :: command_unit
    character ( len = 255 ) :: data_filename
    integer                 :: data_unit
    integer                 :: i, j
    integer                 :: nz_num
    !
    !  Create data file.
    !
    if(.not.sparse%status)return
    !
    N1 = sparse%Nrow
    N2 = sparse%Ncol
    data_filename = trim ( header ) // '_data.dat'
    open (unit=free_unit(data_unit), file = data_filename, status = 'replace' )
    nz_num = 0
    do i=1,N1
       do j=1,sparse%row(i)%size
          write(data_unit,'(2x,i6,2x,i6)') sparse%row(i)%cols(j),i
          nz_num = nz_num + 1
       enddo
    enddo
    close(data_unit)
    !
    !  Create command file.
    !
    command_filename = "plot_"//str(header)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)')'set title "',nz_num,' nonzeros for '//str(header)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine sp_spy_matrix





  !##################################################################
  !##################################################################
  !               SPARSE MATRIX KRON PRODUCT 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Perform simple Kroenecker product of two sparse matrices
  !AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  !+------------------------------------------------------------------+
  function sp_kron_matrix(A,B) result(AxB)
    type(sparse_matrix), intent(in)     :: A,B
    type(sparse_matrix)                 :: AxB
    integer                             :: i,icol,j,k,kcol,l,ii
    integer                             :: indx_row,indx_col
    integer,dimension(:),allocatable    :: cols
#ifdef _CMPLX
    complex(8),dimension(:),allocatable :: vals
    complex(8)                          :: value
#else
    real(8),dimension(:),allocatable    :: vals
    real(8)                             :: value
#endif
    if(.not.A%status)stop "sp_kron_matrix: A.status=F"
    if(.not.B%status)stop "sp_kron_matrix: B.status=F"
    call AxB%free()
    call AxB%init(a%Nrow*b%Nrow,a%Ncol*b%Ncol)
    allocate(cols(a%Ncol*b%Ncol))
    allocate(vals(a%Ncol*b%Ncol))
    do indx_row = 1,A%Nrow*B%Nrow
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       i = (indx_row-1)/B%Nrow+1
       ii=0
       vals=0
       cols=0
       do icol=1,A%row(i)%size
          j = A%row(i)%cols(icol)
          do kcol=1,B%row(k)%size
             l = B%row(k)%cols(kcol)
             indx_col = l + (j-1)*B%Ncol
             ii = ii+1            
             ! value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
             ! !
             ! call append(AxB%row(indx_row)%vals,value)
             ! call append(AxB%row(indx_row)%cols,indx_col)
             vals(ii) = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
             cols(ii) = indx_col
             AxB%row(indx_row)%Size = AxB%row(indx_row)%Size + 1
          enddo
       enddo
       call append(AxB%row(indx_row)%vals,vals(1:ii))
       call append(AxB%row(indx_row)%cols,cols(1:ii))
    enddo
    deallocate(vals,cols)
  end function sp_kron_matrix

  function sp_restricted_kron_matrix(A,B,states) result(AxB)
    type(sparse_matrix), intent(in)     :: A,B
    integer,dimension(:),intent(in)     :: states
    type(sparse_matrix)                 :: AxB,Ap,Bp
    integer                             :: i,icol,j,k,kcol,l,istate,jstate,ii
    integer                             :: indx_row,indx_col
    integer,dimension(:),allocatable    :: cols
#ifdef _CMPLX
    complex(8),dimension(:),allocatable :: vals
    complex(8)                          :: val,Aval,Bval
#else
    real(8),dimension(:),allocatable    :: vals
    real(8)                             :: val,Aval,Bval
#endif
    !
    call AxB%free()
    call AxB%init(size(states),size(states))
    !
    allocate(cols(size(states)))
    allocate(vals(size(states)))
    do istate = 1,size(states)
       indx_row=states(istate)
       i = (indx_row-1)/B%Nrow+1
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       ii=0
       vals=0
       cols=0
       do jstate=1,size(states)
          indx_col=states(jstate)
          j = (indx_col-1)/B%Ncol+1
          l = mod(indx_col,B%Ncol);if(l==0)l=B%Ncol
          if(&
               (.not.any(A%row(i)%cols==j))&
               .OR. &
               (.not.any(B%row(k)%cols==l)) )cycle
          Aval = A%get(i,j)
          Bval = B%get(k,l)
          ii   = ii+1
          vals(ii) = Aval*Bval
          cols(ii) = jstate
          ! call append(AxB%row(istate)%vals,Aval*Bval)
          ! call append(AxB%row(istate)%cols,jstate)
          AxB%row(istate)%Size = AxB%row(istate)%Size + 1
       enddo
       call append(AxB%row(istate)%vals,vals(1:ii))
       call append(AxB%row(istate)%cols,cols(1:ii))
    enddo
    deallocate(vals,cols)    
  end function sp_restricted_kron_matrix




  function sp_matmul_matrix(A,B) result(AxB)
    type(sparse_matrix), intent(in)     :: A,B
    type(sparse_matrix)                 :: Bt,AxB
    integer                             :: i,icol,j,jcol,k,ii
    integer,dimension(:),allocatable    :: cols
#ifdef _CMPLX
    complex(8),dimension(:),allocatable :: vals
    complex(8)                          :: value
#else
    real(8),dimension(:),allocatable    :: vals
    real(8)                             :: value
#endif
    integer, dimension(:), allocatable  :: indx_A, indx_B
    integer                             :: count
    !
    if(.not.A%status)stop "sp_matmul_matrix: A.status=F"
    if(.not.B%status)stop "sp_matmul_matrix: B.status=F"
    !
    !Assume A & B are known to all NODES
    call AxB%free
    call AxB%init(a%Nrow,b%Ncol)
    Bt = B%t()
    !
    allocate(cols(b%Ncol))
    allocate(vals(b%Ncol))
    !
    do i=1,AxB%Nrow    !==A.Nrow
       ii  = 0
       cols= 0
       vals= 0
       do j=1,AxB%Ncol !==Bt.Nrow=B.Ncol
          ! if(.NOT.check_intersection(A%row(i)%cols, Bt%row(j)%cols))cycle
          call get_intersection(A%row(i)%cols, Bt%row(j)%cols) !count = 0 cycle 
          value = zero
          do k=1,count
             icol = indx_A(k)
             jcol = indx_B(k)
             value    = value + A%row(i)%vals(icol)*Bt%row(j)%vals(jcol)
          enddo
          if(value==zero)cycle
          ii = ii + 1
          vals(ii) = value
          cols(ii) = j
       enddo
       AxB%row(i)%Size = ii
       call append(AxB%row(i)%vals,vals(1:ii))
       call append(AxB%row(i)%cols,cols(1:ii))
    enddo
    call Bt%free
    deallocate(cols,vals)
  contains
    !
    logical function check_intersection(A, B)
      integer, dimension(:), intent(in) :: A, B
      integer                           :: i !local copy, no interference
      do i=1,size(A)
         check_intersection=any(B==A(i))
         if(check_intersection)exit
      enddo
    end function check_intersection
    !
    !Order O(size(A)+size(B)) << O(size(A)*size(B))
    subroutine get_intersection(A, B)
      integer, dimension(:), intent(in) :: A, B
      integer                           :: max_intersections
      integer                           :: i, j !local copies no interference
      !
      max_intersections = min(size(A), size(B))
      if (allocated(indx_A)) deallocate(indx_A)
      if (allocated(indx_B)) deallocate(indx_B)
      allocate(indx_A(max_intersections))
      allocate(indx_B(max_intersections))
      !
      i = 1
      j = 1
      count = 0
      do while (i <= size(A) .and. j <= size(B))
         if (A(i) < B(j)) then
            i = i + 1
         else if (A(i) > B(j)) then
            j = j + 1
         else ! A(i) == B(j)
            count = count + 1
            indx_A(count) = i
            indx_B(count) = j
            i = i + 1
            j = j + 1
         end if
      end do
    end subroutine get_intersection
    !
  end function sp_matmul_matrix



  function sp_matmul_vector(H,v) result(Hv)
    class(sparse_matrix), intent(in)    :: H
    integer                             :: Nloc
#ifdef _CMPLX
    complex(8),dimension(:)             :: v
    complex(8),dimension(:),allocatable :: Hv
    complex(8)                          :: val
#else
    real(8),dimension(:)                :: v
    real(8),dimension(:),allocatable    :: Hv
    real(8)                             :: val
#endif
    integer                             :: jcol
    !
    if(.not.H%status)stop "sp_matmul_vector: H.status=F"
    !
    if(allocated(Hv))deallocate(Hv)
    allocate(Hv(size(v)))
    Hv=zero
    Nloc=size(v)
    do i=1,Nloc
       matmul: do jcol=1, H%row(i)%Size
          val = H%row(i)%vals(jcol)
          j   = H%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
  end function sp_matmul_vector




  !##################################################################
  !##################################################################
  !               SPARSE MATRIX BASIC ALGEBRA 
  !##################################################################
  !##################################################################
  function sp_dgr_matrix(a) result(c)
    class(sparse_matrix), intent(in) :: a
    type(sparse_matrix)              :: c
    integer                          :: i, j, k, col, M, pos
    integer, allocatable             :: nnz_per_row(:), idx(:)
#ifdef _CMPLX
    complex(8)                       :: val
#else
    real(8)                          :: val
#endif
    if(.not.a%status)stop "sp_transpose_matrix: A.status=F"
    call c%init(a%Ncol,a%Nrow)       !tranpose
    !
    !Count how many elements are per row with a given col-index
    allocate(nnz_per_row(c%Nrow)) ! c%Nrow == a%Ncol
    nnz_per_row = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          nnz_per_row(col) = nnz_per_row(col) + 1
       enddo
    enddo
    !Allocate C manually to correct size
    do k = 1, c%Nrow
       M = nnz_per_row(k)
       c%row(k)%size = M
       if(allocated(c%row(k)%vals)) deallocate(c%row(k)%vals)
       if(allocated(c%row(k)%cols)) deallocate(c%row(k)%cols)
       allocate(c%row(k)%vals(M))
       allocate(c%row(k)%cols(M))
    enddo
    !Build C using a counter to get the incremental col position:
    allocate(idx(c%Nrow))
    idx = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          k      = a%row(i)%cols(j)  ! A-col => new C-row index
          val    = a%row(i)%vals(j)
          idx(k) = idx(k) + 1      !advance counter 
          pos    = idx(k)
          c%row(k)%cols(pos) = i
#ifdef _CMPLX
          c%row(k)%vals(pos) = conjg(val)
#else
          c%row(k)%vals(pos) = val
#endif
       enddo
    enddo
    deallocate(nnz_per_row, idx)
  end function sp_dgr_matrix


  function sp_transpose_matrix(a) result(c)
    class(sparse_matrix), intent(in) :: a
    type(sparse_matrix)              :: c
    integer                          :: i, j, k, col, M, pos
    integer, allocatable             :: nnz_per_row(:), idx(:)
#ifdef _CMPLX
    complex(8)                       :: val
#else
    real(8)                          :: val
#endif
    if(.not.a%status)stop "sp_transpose_matrix: A.status=F"
    call c%init(a%Ncol,a%Nrow)       !tranpose
    !
    !Count how many elements are per row with a given col-index
    allocate(nnz_per_row(c%Nrow)) ! c%Nrow == a%Ncol
    nnz_per_row = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          nnz_per_row(col) = nnz_per_row(col) + 1
       enddo
    enddo
    !Allocate C manually to correct size
    do k = 1, c%Nrow
       M = nnz_per_row(k)
       c%row(k)%size = M
       if(allocated(c%row(k)%vals)) deallocate(c%row(k)%vals)
       if(allocated(c%row(k)%cols)) deallocate(c%row(k)%cols)
       allocate(c%row(k)%vals(M))
       allocate(c%row(k)%cols(M))
    enddo
    !Build C using a counter to get the incremental col position:
    allocate(idx(c%Nrow))
    idx = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          k      = a%row(i)%cols(j)  ! A-col => new C-row index
          val    = a%row(i)%vals(j)
          idx(k) = idx(k) + 1      !advance counter 
          pos    = idx(k)
          c%row(k)%cols(pos) = i
          c%row(k)%vals(pos) = val
       enddo
    enddo
    deallocate(nnz_per_row, idx)
  end function sp_transpose_matrix


  function sp_hconjg_matrix(a) result(c)
    class(sparse_matrix), intent(in) :: a
    type(sparse_matrix)              :: c
    integer                          :: i, j, k, col, M, pos
    integer, allocatable             :: nnz_per_row(:), idx(:)
#ifdef _CMPLX
    complex(8)                       :: val
#else
    real(8)                          :: val
#endif
    if(.not.a%status)stop "sp_transpose_matrix: A.status=F"
    call c%init(a%Ncol,a%Nrow)       !tranpose
    !
    !Count how many elements are per row with a given col-index
    allocate(nnz_per_row(c%Nrow)) ! c%Nrow == a%Ncol
    nnz_per_row = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          nnz_per_row(col) = nnz_per_row(col) + 1
       enddo
    enddo
    !Allocate C manually to correct size
    do k = 1, c%Nrow
       M = nnz_per_row(k)
       c%row(k)%size = M
       if(allocated(c%row(k)%vals)) deallocate(c%row(k)%vals)
       if(allocated(c%row(k)%cols)) deallocate(c%row(k)%cols)
       allocate(c%row(k)%vals(M))
       allocate(c%row(k)%cols(M))
    enddo
    !Build C using a counter to get the incremental col position:
    allocate(idx(c%Nrow))
    idx = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          k      = a%row(i)%cols(j)  ! A-col => new C-row index
          val    = a%row(i)%vals(j)
          idx(k) = idx(k) + 1      !advance counter 
          pos    = idx(k)
          c%row(k)%cols(pos) = i
#ifdef _CMPLX
          c%row(k)%vals(pos) = conjg(val)
#else
          c%row(k)%vals(pos) = val
#endif
       enddo
    enddo
    deallocate(nnz_per_row, idx)
  end function sp_hconjg_matrix


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix equality spA = spB. Deep copy
  !+------------------------------------------------------------------+
  subroutine sp_matrix_equal_matrix(a,b)
    type(sparse_matrix),intent(inout) :: a
    type(sparse_matrix),intent(in)    :: b
    if(.not.b%status)stop "sp_matrix_equal_matrix: B.status=F"
    call a%free()
    call a%init(b%Nrow,b%Ncol)
    do i=1,b%Nrow
       a%row(i)%size = b%row(i)%size
       call append(a%row(i)%cols,b%row(i)%cols)
       call append(a%row(i)%vals,b%row(i)%vals)
    enddo
  end subroutine sp_matrix_equal_matrix



  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix scalar equality spA = const. 
  !+------------------------------------------------------------------+
  subroutine sp_matrix_equal_scalar(a,c)
    type(sparse_matrix),intent(inout) :: a
#ifdef _CMPLX
    complex(8),intent(in)             :: c
#else
    real(8),intent(in)                :: c
#endif
    if(.not.a%status)stop "sp_matrix_equal_matrix: A.status=F"
    do i=1,a%Nrow
       a%row(i)%vals = c
    enddo
  end subroutine sp_matrix_equal_scalar



  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix addition spA + spB = spC
  !+------------------------------------------------------------------+
  !   function sp_plus_matrix(a,b) result(c)
  !     type(sparse_matrix), intent(in) :: a,b
  !     type(sparse_matrix)             :: c
  !     integer                         :: col
  ! #ifdef _CMPLX
  !     complex(8)                      :: val
  ! #else
  !     real(8)                         :: val
  ! #endif
  !     if(.not.a%status)stop "sp_plus_matrix error: a.status=F"
  !     if(.not.b%status)stop "sp_plus_matrix error: b.status=F"
  !     if(a%Nrow/=b%Nrow)stop "sp_plus_matrix error: a.Nrow != b.Nrow"
  !     if(a%Ncol/=b%Ncol)stop "sp_plus_matrix error: a.Ncol != b.Ncol"
  !     c=a                         !copy a into c
  !     do i=1,b%Nrow
  !        do j=1,b%row(i)%size
  !           col = b%row(i)%cols(j)
  !           val = b%row(i)%vals(j)
  !           call c%insert(val,i,col)
  !        enddo
  !     enddo
  !   end function sp_plus_matrix

  function sp_plus_matrix(a,b) result(c)
    type(sparse_matrix), intent(in) :: a,b
    type(sparse_matrix)             :: c
    integer                         :: i, na, nb, p, q, nn, col_a, col_b
    integer,allocatable             :: cols(:)
#ifdef _CMPLX
    complex(8), allocatable         :: vals(:)
    complex(8)                      :: va, vb
#else
    real(8), allocatable            :: vals(:)
    real(8)                         :: va, vb
#endif
    if(.not.a%status) stop "sp_plus_matrix_merge error: a.status=F"
    if(.not.b%status) stop "sp_plus_matrix_merge error: b.status=F"
    if(a%Nrow /= b%Nrow) stop "sp_plus_matrix_merge error: Nrow mismatch"
    if(a%Ncol /= b%Ncol) stop "sp_plus_matrix_merge error: Ncol mismatch"
    !
    call c%free()
    call c%init(a%Nrow,a%Ncol)
    !
    do i=1,a%Nrow
       na = a%row(i)%size
       nb = b%row(i)%size
       if(na+nb == 0) cycle
       ! allocate worst-case storage
       allocate(cols(na+nb))
       allocate(vals(na+nb))
       p = 1
       q = 1
       nn = 0
       do while(p <= na .or. q <= nb)
          if(p > na) then
             ! copy remaining from b
             nn = nn + 1
             cols(nn) = b%row(i)%cols(q)
             vals(nn) = b%row(i)%vals(q)
             q = q + 1
          elseif(q > nb) then
             nn = nn + 1
             cols(nn) = a%row(i)%cols(p)
             vals(nn) = a%row(i)%vals(p)
             p = p + 1
          else
             col_a = a%row(i)%cols(p)
             col_b = b%row(i)%cols(q)
             if(col_a == col_b) then
                va = a%row(i)%vals(p)
                vb = b%row(i)%vals(q)
                nn = nn + 1
                cols(nn) = col_a
                vals(nn) = va + vb
                p = p + 1; q = q + 1
             elseif(col_a < col_b) then
                nn = nn + 1
                cols(nn) = col_a
                vals(nn) = a%row(i)%vals(p)
                p = p + 1
             else
                nn = nn + 1
                cols(nn) = col_b
                vals(nn) = b%row(i)%vals(q)
                q = q + 1
             endif
          endif
       end do
       !
       c%row(i)%size = nn
       call append(c%row(i)%cols,cols(1:nn))
       call append(c%row(i)%vals,vals(1:nn))
       deallocate(cols,vals)
    end do
  end function sp_plus_matrix



  ! Merge three sparse matrices a, b, c into d in one rowwise pass.
  function sp_plus3_matrix(a,b,c) result(d)
    type(sparse_matrix), intent(in) :: a,b,c
    type(sparse_matrix)             :: d
    integer                         :: i, na, nb, nc, p, q, r, nn
    integer                         :: ca, cb, cc

    integer, allocatable            :: cols(:)
#ifdef _CMPLX
    complex(8), allocatable         :: vals(:)
    complex(8)                      :: va, vb, vc
    complex(8)                      :: sumv
#else
    real(8), allocatable            :: vals(:)
    real(8)                         :: va, vb, vc
    real(8)                         :: sumv
#endif
    if(.not.a%status .or. .not.b%status .or. .not.c%status) stop "sp_plus3_matrix: status error"
    if(a%Nrow /= b%Nrow .or. a%Nrow /= c%Nrow) stop "sp_plus3_matrix: Nrow mismatch"
    if(a%Ncol /= b%Ncol .or. a%Ncol /= c%Ncol) stop "sp_plus3_matrix: Ncol mismatch"
    !
    call d%free()
    call d%init(a%Nrow, a%Ncol)
    !
    do i = 1, a%Nrow
       na = a%row(i)%size; nb = b%row(i)%size; nc = c%row(i)%size
       if(na+nb+nc == 0) cycle
       ! allocate worst-case once for this row
       allocate(cols(na+nb+nc))
       allocate(vals(na+nb+nc))
       !
       p = 1
       q = 1
       r = 1
       nn = 0
       !
       do while(p<=na .or. q<=nb .or. r<=nc)
          !
          ca = huge(1)
          cb = huge(1)
          cc = huge(1)
          if(p<=na) ca = a%row(i)%cols(p)
          if(q<=nb) cb = b%row(i)%cols(q)
          if(r<=nc) cc = c%row(i)%cols(r)
          !
          if(ca <= cb .AND. ca <= cc) then
             va   = a%row(i)%vals(p)
             sumv = va
             if(cb == ca) then
                sumv = sumv + b%row(i)%vals(q)
                q = q+1
             endif
             if(cc == ca) then
                sumv = sumv + c%row(i)%vals(r)
                r = r+1
             endif
             p = p + 1
             if(abs(sumv)>1d-16)then
                nn = nn + 1
                cols(nn) = ca
                vals(nn) = sumv
             endif
          elseif(cb <= ca .AND. cb <= cc) then
             vb   = b%row(i)%vals(q)
             sumv = vb
             if(cc == cb) then
                sumv = sumv + c%row(i)%vals(r)
                r = r+1
             end if
             q = q + 1
             if(abs(sumv)>1d-16) then
                nn = nn + 1
                cols(nn) = cb
                vals(nn) = sumv
             endif
          else
             vc   = c%row(i)%vals(r)
             sumv = vc
             r = r + 1
             if(abs(sumv)>1d-16) then
                nn = nn + 1
                cols(nn) = cc
                vals(nn) = sumv
             end if
          end if
       end do
       !
       d%row(i)%size = nn
       if(nn > 0) then
          call append(d%row(i)%cols, cols(1:nn))
          call append(d%row(i)%vals, vals(1:nn))
       end if
       deallocate(cols, vals)
       !
    end do
  end function sp_plus3_matrix



  
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix difference spA - spB = spC
  !+------------------------------------------------------------------+
  function sp_minus_matrix(a,b) result(c)
    type(sparse_matrix), intent(in) :: a,b
    type(sparse_matrix)             :: c
    integer                         :: col
#ifdef _CMPLX
    complex(8)                      :: val
#else
    real(8)                         :: val
#endif
    if(.not.a%status)stop "sp_minus_matrix error: a.status=F"
    if(.not.b%status)stop "sp_minus_matrix error: b.status=F"
    if(a%Nrow/=b%Nrow)stop "sp_minus_matrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "sp_minus_matrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val =-b%row(i)%vals(j)
          call c%insert(val,i,col)
       enddo
    enddo
  end function sp_minus_matrix




  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix left scalar product const*spA = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_left_product_matrix_i(C,A) result(B)
    integer,intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
#ifdef _CMPLX
    complex(8)                     :: val
#else
    real(8)                        :: val
#endif
    if(.not.A%status)stop "sp_left_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo
  end function sp_left_product_matrix_i

  function sp_left_product_matrix_d(C,A) result(B)
    real(8),intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
#ifdef _CMPLX
    complex(8)                     :: val
#else
    real(8)                        :: val
#endif
    if(.not.A%status)stop "sp_left_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo
  end function sp_left_product_matrix_d

#ifdef _CMPLX
  function sp_left_product_matrix_c(C,A) result(B)
    complex(8),intent(in)          :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    complex(8)                     :: val
    if(.not.A%status)stop "sp_left_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo
  end function sp_left_product_matrix_c
#endif



  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar product spA*const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_right_product_matrix_i(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
#ifdef _CMPLX
    complex(8)                     :: val
#else
    real(8)                        :: val
#endif
    if(.not.A%status)stop "sp_right_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo
  end function sp_right_product_matrix_i

  function sp_right_product_matrix_d(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
#ifdef _CMPLX
    complex(8)                     :: val
#else
    real(8)                        :: val
#endif
    if(.not.A%status)stop "sp_right_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo
  end function sp_right_product_matrix_d

#ifdef _CMPLX
  function sp_right_product_matrix_c(A,C) result(B)
    complex(8),intent(in)          :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    complex(8)                     :: val
    if(.not.A%status)stop "sp_right_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo
  end function sp_right_product_matrix_c
#endif



  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar division spA/const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_right_division_matrix_i(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
#ifdef _CMPLX
    complex(8)                     :: val
#else
    real(8)                        :: val
#endif
    if(.not.A%status)stop "sp_right_division_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals/C)
    enddo
  end function sp_right_division_matrix_i

  function sp_right_division_matrix_d(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
#ifdef _CMPLX
    complex(8)                     :: val
#else
    real(8)                        :: val
#endif
    if(.not.A%status)stop "sp_right_division_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals/C)
    enddo
  end function sp_right_division_matrix_d

#ifdef _CMPLX
  function sp_right_division_matrix_c(A,C) result(B)
    complex(8),intent(in)          :: C
    type(sparse_matrix),intent(in) :: A
    type(sparse_matrix)            :: B
    integer                        :: col
    complex(8)                     :: val
    if(.not.A%status)stop "sp_right_division_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals/C)
    enddo
  end function sp_right_division_matrix_c
#endif


  !+------------------------------------------------------------------+
  !PURPOSE:  Return the shape of a sparse matrix
  !+------------------------------------------------------------------+
  function sp_shape_matrix(self) result(shape)
    class(sparse_matrix),intent(in) :: self
    integer,dimension(2)            :: shape
    shape = [self%Nrow,self%Ncol]
  end function sp_shape_matrix



  !+------------------------------------------------------------------+
  !PURPOSE:  Return the identiy sparse matrix of given dimension
  !+------------------------------------------------------------------+
  function sp_eye(ndim) result(self)
    type(sparse_matrix) :: self
    integer             :: ndim
#ifdef _CMPLX
    call self%load(zeye(ndim))
#else
    call self%load(eye(ndim))
#endif
  end function sp_eye







#ifdef _MPI
  subroutine sp_bcast_matrix(self,comm,root)
    class(sparse_matrix), intent(inout) :: self
    integer,intent(in),optional         :: comm
    integer,intent(in),optional         :: root
    integer                             :: comm_
    integer                             :: root_
    integer                             :: rank, ierr, i
    integer                             :: Nrow, Ncol, Nsize
    logical                             :: master
    !
    if(.not.check_MPI())stop "sp_bcast error: check_MPI=F"
    comm_  = MPI_COMM_WORLD;if(present(comm))comm_ = comm
    root_  = 0 ;if(present(root))root_=root
    rank   = get_Rank_MPI(comm_)
    ! master = get_Master_MPI(comm_)
    master = (rank == root_)
    !
    if(master)Nrow = self%Nrow ; call Bcast_MPI(comm_,Nrow, root=root_)
    if(master)Ncol = self%Ncol ; call Bcast_MPI(comm_,Ncol, root=root_)
    !
    if(.not.master)call self%init(Nrow, Ncol)
    !
    do i=1,Nrow
       if(master)Nsize = self%row(i)%size
       call Bcast_MPI(comm_,Nsize, root=root_)
       !
       self%row(i)%size = Nsize !tautology for master
       if(Nsize==0)cycle
       if (.NOT.master) then
          if(allocated(self%row(i)%cols)) deallocate(self%row(i)%cols)
          if(allocated(self%row(i)%vals)) deallocate(self%row(i)%vals)
          allocate(self%row(i)%cols(Nsize))
          allocate(self%row(i)%vals(Nsize))
       endif
       !
       call Bcast_MPI(comm_,self%row(i)%cols, root=root_)       
       call Bcast_MPI(comm_,self%row(i)%vals, root=root_)
    end do
    call Barrier_MPI
  end subroutine sp_bcast_matrix



  !+------------------------------------------------------------------+
  ! PURPOSE:  Performs an "All-Gather" operation on a distributed
  !           array of sparse_matrix objects. After the call, every
  !           process will have the complete, assembled array.
  ! USAGE:    Assumes the array is distributed, where for each index `i`,
  !           the matrix A(i) is initialized (i.e., A(i)%status == .true.)
  !           on exactly ONE process.
  !+------------------------------------------------------------------+
  subroutine sp_array_all_gather_1(comm,A)
    class(sparse_matrix),dimension(:),intent(inout) :: A
    integer,intent(in)                              :: comm
    integer,allocatable,dimension(:)                :: local_ownership, global_ownership
    integer                                         :: i, N, rank, owner_rank, ierr
    logical                                         :: master
    !
    if (.not. check_MPI()) stop "sp_array_all_gather: check_MPI=F"
    !
    N = size(A)
    rank   = get_Rank_MPI(comm)
    master = get_Master_MPI(comm)
    !
    ! Step 1: Determine the owner of each array element.
    allocate(local_ownership(N))
    allocate(global_ownership(N))
    !
    ! Each process creates a local map. If it owns A(i), it marks it with its rank.
    ! Otherwise, it marks it with -1 (or any value lower than a valid rank).
    do i=1,N
       if (A(i)%status) then
          local_ownership(i) = rank
       else
          local_ownership(i) = -1
       endif
    enddo
    !
    ! Use MPI_Allreduce with MPI_MAX. For each index `i`, the process with the
    ! valid rank will "win", and every process will get the correct owner rank.
    call MPI_Allreduce(local_ownership, global_ownership, N, MPI_INTEGER, MPI_MAX, comm, ierr)
    !
    ! Step 2: Iterate through the array and broadcast each element from its owner.
    do i = 1, N
       owner_rank = global_ownership(i)
       !
       ! If owner_rank is -1, it means NO process initialized this element.
       ! Skip the broadcast, leave A(i) uninitialized state on all processes.
       if (owner_rank == -1) then
          cycle
       else
          !Owner exist, bcast to all nodes
          call A(i)%bcast(comm,root=owner_rank)
       endif
       !
    enddo
    !
    deallocate(global_ownership)
    deallocate(local_ownership)
    call Barrier_MPI(comm)
  end subroutine sp_array_all_gather_1


  subroutine sp_array_all_gather_2(comm,A)
    class(sparse_matrix),dimension(:,:),intent(inout) :: A
    integer,intent(in)                                :: comm
    integer,allocatable,dimension(:,:)                :: local_ownership,global_ownership
    integer                                           :: i,j,N,M,rank,owner_rank,ierr
    logical                                           :: master
    !
    if (.not. check_MPI()) stop "sp_array_all_gather: check_MPI=F"
    !
    N = size(A,1)
    M = size(A,2)
    !
    rank   = get_Rank_MPI(comm)
    master = get_Master_MPI(comm)
    !
    ! Step 1: Determine the owner of each array element.
    allocate(local_ownership(N,M))
    allocate(global_ownership(N,M))
    global_ownership=0
    ! Each process creates a local map. If it owns A(i), it marks it with its rank.
    ! Otherwise, it marks it with -1 (or any value lower than a valid rank).
    do i=1,N
       do j=1,M
          if (A(i,j)%status) then
             local_ownership(i,j) = rank
          else
             local_ownership(i,j) = -1
          endif
       enddo
    enddo
    !
    ! Use MPI_Allreduce with MPI_MAX. For each index `i`, the process with the
    ! valid rank will "win", and every process will get the correct owner rank.
    call MPI_Allreduce(local_ownership, global_ownership, size(local_ownership), MPI_INTEGER, MPI_MAX, comm, ierr)
    !
    ! Step 2: Iterate through the array and broadcast each element from its owner.
    do i=1,N
       do j=1,M
          owner_rank = global_ownership(i,j)
          ! If owner_rank is -1, it means NO process initialized this element.
          ! Skip the broadcast, leave A(i) uninitialized state on all processes.
          if (owner_rank == -1) then
             cycle
          else
             !Owner exist, bcast to all nodes
             call A(i,j)%bcast(comm,root=owner_rank)
          endif
       enddo
    enddo
    !
    deallocate(global_ownership)
    deallocate(local_ownership)
    call Barrier_MPI(comm)
  end subroutine sp_array_all_gather_2






  function sp_get_bytes_matrix(self) result(kbytes)
    class(sparse_matrix), intent(inout) :: self
    real(8)                             :: kbytes
    integer                             :: i,ierr
    real(8)                             :: I_kb_size  
    real(8)                             :: DATA_kb_size
    integer                             :: typesize
    logical                             :: master
    if(.not.check_MPI())stop "sp_get_kbytes error: check_MPI=F"
    call MPI_Type_size(MPI_INTEGER, typesize, ierr)
#ifdef _CMPLX
    call MPI_Type_size(MPI_DOUBLE_COMPLEX, typesize, ierr)
#else
    call MPI_Type_size(MPI_DOUBLE_PRECISION, typesize, ierr)
#endif
    I_kb_size    = dble(typesize)/1000d0
    DATA_kb_size = dble(typesize)/1000d0
    master = get_Master_MPI(MPI_COMM_WORLD)
    kbytes=0d0
    if(master)then
       kbytes = 2*I_kb_size                            !Nrow + Ncol
       do i=1,self%Nrow
          kbytes=kbytes+I_kb_size                      !row(i).Nsize
          kbytes=kbytes+2*(self%row(i)%size)*I_kb_size !row(i).cols + row(i).vals
       enddo
    endif
  end function sp_get_bytes_matrix

#endif


  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  recursive function binary_search(Ain,value) result(bsresult)
    integer,intent(in)           :: Ain(:), value
    integer                      :: bsresult, mid
    integer,dimension(size(Ain)) :: A,Order
    !
    a = ain
    call sort_array(a,Order)
    !
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
    !
    bsresult = Order(bsresult)
    !
  end function binary_search



  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  subroutine sort_array(array,order)
    implicit none
    integer,dimension(:)                    :: array
    integer,dimension(size(array))          :: order
    integer,dimension(size(array))          :: backup
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort(array, order,1, size(array))
    do i=1,size(array)
       backup(i)=array(order(i))
    enddo
    array=backup
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:) :: array
      integer, dimension(:) :: order
      integer               :: left
      integer               :: right
      integer               :: i
      integer               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:) :: order
      integer               :: first, second
      integer               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    integer function qsort_rand( lower, upper )
      integer               :: lower, upper
      real(8)               :: r
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    !---------------------------------------------!
    function compare(f,g)
      implicit none
      integer               :: f,g
      integer               :: compare
      if(f<g) then
         compare=-1
      else
         compare=1
      endif
    end function compare
  end subroutine sort_array












#ifdef _MPI
  function sp_p_matmul_matrix(A,B) result(AxB)
    class(sparse_matrix), intent(in)    :: A,B
    type(sparse_matrix)                 :: Bt,AxB,AxB_local
    integer                             :: i,icol,j,jcol,k, ierr
    integer                             :: comm_
    integer                             :: rank, ncpu
    integer                             :: i_start, i_end, Nrows_per_cpu, remainder
    integer                             :: p_i_start, p_i_end
    integer                             :: p_rank, Nsize,ii
    integer,dimension(:),allocatable    :: cols
#ifdef _CMPLX
    complex(8),dimension(:),allocatable :: vals
    complex(8)                          :: value
#else
    real(8),dimension(:),allocatable    :: vals
    real(8)                             :: value
#endif
    integer, dimension(:), allocatable  :: indx_A, indx_B
    integer                             :: count
    !
    if(.not.check_MPI())stop "sp_p_matmul_matrix error: check_MPI=F"
    comm_ = MPI_COMM_WORLD!;if(present(comm))comm_=comm
    rank = get_Rank_MPI(comm_)
    ncpu = get_Size_MPI(comm_)
    !
    if(.not.A%status)stop "sp_p_matmul_matrix: A.status=F"
    if(.not.B%status)stop "sp_p_matmul_matrix: B.status=F"
    !
    Bt = B%t()
    !
    ! Local Computation
    ! Each process computes its assigned rows and stores them in a local matrix.
    call AxB_local%init(A%Nrow, B%Ncol)
    allocate(cols(B%Ncol))
    allocate(vals(B%Ncol))
    do i=1+rank, A%Nrow, ncpu
       !
       ii  = 0
       vals= 0
       cols= 0
       do j=1,Bt%Nrow
          call get_intersection(A%row(i)%cols, Bt%row(j)%cols)
          value = zero
          do k=1,count
             icol = indx_A(k)
             jcol = indx_B(k)
             value = value + A%row(i)%vals(icol)*Bt%row(j)%vals(jcol)
          enddo
          if(value==zero)cycle
          ii = ii+1
          vals(ii) = value
          cols(ii) = j
          AxB_local%row(i)%Size = AxB_local%row(i)%Size + 1
       enddo
       !
       call append(AxB_local%row(i)%vals,vals(1:ii))!value)
       call append(AxB_local%row(i)%cols,cols(1:ii))!j)
    enddo
    call Bt%free
    deallocate(vals,cols)
    !
    ! All-Gather Results
    ! The AxB_local are now gathered onto all processes.
    ! We do this by iterating through each process and having it broadcast its
    ! computed rows to everyone else.
    call AxB%init(A%Nrow, B%Ncol)
    !
    do p_rank = 0, ncpu - 1
       ! Broadcast each computed row from process 'p_rank' to all others.
       do i = 1+p_rank,A%Nrow,ncpu
          Nsize=0
          if (rank == p_rank) Nsize = AxB_local%row(i)%size
          call MPI_Bcast(Nsize, 1, MPI_INTEGER, p_rank, comm_, ierr)
          !
          ! All processes now know the size of row 'i'.
          AxB%row(i)%size = Nsize
          if(allocated(AxB%row(i)%cols)) deallocate(AxB%row(i)%cols)
          if(allocated(AxB%row(i)%vals)) deallocate(AxB%row(i)%vals)
          if (Nsize > 0) then
             if (rank /= p_rank) then
                allocate(AxB%row(i)%cols(Nsize))
                allocate(AxB%row(i)%vals(Nsize))
             else
                allocate(AxB%row(i)%cols(Nsize))
                allocate(AxB%row(i)%vals(Nsize))
                AxB%row(i)%cols = AxB_local%row(i)%cols
                AxB%row(i)%vals = AxB_local%row(i)%vals
             endif
             call MPI_Bcast(AxB%row(i)%cols, Nsize, MPI_INTEGER, p_rank, comm_, ierr)
             call MPI_Bcast(AxB%row(i)%vals, Nsize, MPI_VAL_TYPE, p_rank, comm_, ierr)
          else
             allocate(AxB%row(i)%cols(0))
             allocate(AxB%row(i)%vals(0))
          endif
       enddo
    enddo
    !
    call AxB_local%free()
    !
  contains
    !
    subroutine get_intersection(A, B)
      integer, dimension(:), intent(in) :: A, B
      integer                           :: max_intersections
      integer                           :: i, j
      max_intersections = min(size(A), size(B))
      if (allocated(indx_A)) deallocate(indx_A)
      if (allocated(indx_B)) deallocate(indx_B)
      allocate(indx_A(max_intersections))
      allocate(indx_B(max_intersections))
      i = 1
      j = 1
      count = 0
      do while (i <= size(A) .and. j <= size(B))
         if (A(i) < B(j)) then
            i = i + 1
         else if (A(i) > B(j)) then
            j = j + 1
         else ! A(i) == B(j)
            count = count + 1
            indx_A(count) = i
            indx_B(count) = j
            i = i + 1
            j = j + 1
         end if
      end do
    end subroutine get_intersection
    !
  end function sp_p_matmul_matrix
#endif






end module MATRIX_SPARSE








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
program testSPARSE_MATRICES
  USE MATRIX_SPARSE
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none


  integer                                          :: i,j
  type(sparse_matrix)                              :: spH,spK,a,b,c,avec(2)


  type(sparse_matrix), dimension(:), allocatable :: vecA
  type(sparse_matrix), dimension(:,:), allocatable :: arrA

#ifdef _CMPLX
  complex(8),dimension(4,4)                        :: GammaX
  complex(8),dimension(:,:),allocatable            :: Amat,Bmat,Cmat
  complex(8),dimension(2,2),parameter              :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  complex(8),dimension(2,2),parameter              :: S0=pauli_0
  complex(8),dimension(2,2),parameter              :: Sz=pauli_z
  complex(8),dimension(2,2),parameter              :: Sx=pauli_x
  complex(8),dimension(2,2),parameter              :: Splus=reshape([zero,zero,one,zero],[2,2])
  complex(8),dimension(4,4)                        :: Gamma13,Gamma03
  complex(8)                                       :: myone=dcmplx(1d0,0d0),myzero=dcmplx(0d0,0d0)
#else
  real(8),dimension(4,4)                           :: GammaX
  real(8),dimension(:,:),allocatable               :: Amat,Bmat,Cmat
  real(8),dimension(2,2),parameter                 :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  real(8),dimension(2,2),parameter                 :: S0=pauli_0
  real(8),dimension(2,2),parameter                 :: Sz=pauli_z
  real(8),dimension(2,2),parameter                 :: Sx=pauli_x
  real(8),dimension(2,2),parameter                 :: Splus=reshape([zero,zero,one,zero],[2,2])
  real(8),dimension(4,4)                           :: Gamma13,Gamma03
  real(8)                                          :: myone=1d0,myzero=0d0
#endif
  type(sparse_matrix),dimension(:),allocatable     :: Olist
  integer                                          :: irank,comm,rank,ierr,ncpu
  logical                                          :: master=.true.


  integer                                          :: Nvals,Ncols
  !Variables used in MPI_derived_type:
  integer,parameter                                :: mpiBlockNum=3         !# of derived type components
  integer,dimension(mpiBlockNum)                   :: mpiBlockLen !Size of the components  for a single derived_type_array
  integer,dimension(mpiBlockNum)                   :: mpiBlockType     !Type of the components (e,M) of a single
  !element of the derived_type_array
  integer(MPI_ADDRESS_KIND),dimension(mpiBlockNum) :: mpiBlockDisp
  integer(MPI_ADDRESS_KIND)                        :: base
  !from the MAN page for MPI_TYPE_CREATE_STRUCT
  integer                                          :: mpiSparse_Row



#ifdef _MPI
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  ncpu = get_Size_MPI(comm)
  master = get_Master_MPI(comm)
#endif

  Gamma13=kron(Sx,Sz)
  Gamma03=kron(S0,Sz)


  if(master)print*,"test INIT"
  call spH%init(2,2)
  if(master)print*,shape(spH)
  if(master)print*,all(shape(spH) == [2,2])
  if(master)print*,all(shape(spH) == [3,3])
  if(master)print*,""


  if(master)print*,"test CONSTRUCTOR 1: sparse_matrix(matrix)"
  a = sparse_matrix(Sz)
  if(master)call a%show()
  if(master)print*,shape(a)
  if(master)print*,all(shape(a) == [2,2])
  if(master)print*,all(shape(a) == [3,3])
  if(master)print*,"a.NNZ=",a%nnz()
  call a%free()
  if(master)print*,""


  if(master)print*,"test CONSTRUCTOR 2: as_sparse(pauli_x*pauli_z)"
  a = as_sparse(Gamma13)
  if(master)call a%show()
  if(master)print*,shape(a)
  if(master)print*,all(shape(a) == [2,2])
  if(master)print*,all(shape(a) == [4,4])
  if(master)print*,"a.NNZ=",a%nnz()
  if(master)print*,""
  call a%free()



  if(master)print*,"test CONSTRUCTOR 3: avec(1:2)= [as_sparse(pauli_x),as_sparse(pauli_z)]"  
  avec = [sparse(Sx),as_sparse(Sz)]
  if(master)call avec(1)%show()
  if(master)call avec(2)%show()
  if(master)print*,shape(avec(1))
  call avec%free()
  if(master)print*,""


  if(master)print*,"test FREE"
  call spH%free()
  if(master)print*,""



  if(master)print*,"test LOAD and PRINT"
  call spH%load(kron(S0,Sz))
  if(master)call spH%show()
  if(master)print*,"spH.NNZ=",spH%nnz()
  if(master)print*,""

  if(master)print*,"test GET ELEMENT"
  write(*,*)"spH(2,2)=",spH%get(2,2)
  write(*,*)"spH(3,4)=",spH%get(3,4)
  if(master)print*,""


  if(master)print*,"test INSERT ELEMENT"
  call spH%insert(myone,1,4)
  call spH%insert(-myone,4,4)
  if(master)call spH%show()  
  if(master)print*,""


  if(master)print*,"test SPY"
  if(master)call spH%spy("spH")
  if(master)print*,""


  if(master)print*,"test DUMP"
  gammaX=myzero
  do i=1,4
     if(master)write(*,*)(gammaX(i,j),j=1,4)
  enddo
  if(master)print*,""
  if(master)call spH%dump(gammaX)
  do i=1,4
     if(master)write(*,*)(gammaX(i,j),j=1,4)
  enddo

  gammaX=myzero
  do i=1,4
     if(master)write(*,*)(gammaX(i,j),j=1,4)
  enddo
  if(master)print*,""
  if(master)gammaX = spH%as_matrix()
  if(master)call spH%show()
  do i=1,4
     if(master)write(*,*)(gammaX(i,j),j=1,4)
  enddo
  if(master)print*,""


  if(master)print*,"test spK=spH"  
  spK=spH
  if(master)call spK%show()


  if(master)print*,"test spH=zero"  
  spH=myzero
  if(master)call spH%show()
  spH=spK



  if(master)print*,"test ADDITION a+b=c"
  if(master)print*,"a=sigma_0"
  call a%init(2,2)
  call a%load(S0)
  if(master)call a%show()

  if(master)print*,"b=sigma_X"  
  call b%init(2,2)
  call b%load(Sx)
  if(master)call b%show()

  if(master)print*,"c=sigma_0 + sigma_X"
  c = a+b
  if(master)call c%show()

  call a%free()
  call b%free()
  call c%free()



  if(master)print*,"test SUBTRACTION a-b=c"
  if(master)print*,"a=sigma_0"
  call a%init(2,2)
  call a%load(S0)
  if(master)call a%show()

  if(master)print*,"b=sigma_Z"  
  call b%init(2,2)
  call b%load(Sz)
  if(master)call b%show()


  if(master)print*,"c=sigma_0 - sigma_Z"  
  c = a-b
  if(master)call c%show()


  call a%free()
  call b%free()
  call c%free()



  if(master)print*,"test LEFT SCALAR PRODUCT b=a*const"
  if(master)print*,"a=sigma_0"
  call a%init(2,2)
  call a%load(S0)
  if(master)call a%show()

  if(master)print*,"b=2*a"  
  b = 2*a
  if(master)call b%show()

  if(master)print*,"b=2d0*a"  
  b = 2d0*a
  if(master)call b%show()

  if(master)print*,"test RIGHT SCALAR PRODUCT b=const*a"
  if(master)print*,"b=a*2"  
  b = a*2
  if(master)call b%show()

  if(master)print*,"b=a*2d0"  
  b = a*2d0
  if(master)call b%show()


  if(master)print*,"test RIGHT SCALAR DIVISDION b=a/const"
  if(master)print*,"b=a/2"  
  b = a/2
  if(master)call b%show()

  if(master)print*,"b=a/2d0"  
  b = a/2d0
  if(master)call b%show()




  if(master)print*,"test KRON PRODUCT 1"
  call a%free()
  call b%free()
  call c%free()
  call spH%free()
  !
  call spH%load(gamma13)
  call a%load(Sz)
  call b%load(Sx)

  if(master)print*,"a=sigma_0"  
  if(master)call a%show()  
  if(master)print*,"b=sigma_Z"  
  if(master)call b%show()  
  if(master)print*,"c=a.x.b"  
  c = a.x.b
  if(master)call c%show()
  if(master)print*,"spH=sigma_0xsigma_Z"  
  if(master)call spH%show()
  if(master)print*,""


  if(master)print*,""
  if(master)print*,"test KRON PRODUCT 2"
  allocate(Amat(2,2),Bmat(2,2))
  allocate(Cmat(4,4))
  Amat = dble(transpose(reshape([1,2,3,4],[2,2])))
  Bmat = dble(transpose(reshape([0,5,6,7],[2,2])))
  Cmat = dble(transpose(reshape([0,5,0,10,5,7,12,14,0,15,0,20,18,21,24,28],[4,4])))
  call a%load(Amat)
  call b%load(Bmat)
  if(master)print*,"A = 1 2    B = 0 5"
  if(master)print*,"    3 4        6 7"
  if(master)call a%show()  
  if(master)call b%show()  

  if(master)print*,"c=a.x.b"  
  c = a.x.b
  if(master)call c%show()
  if(master)print*,"C = 0  5  0  10"
  if(master)print*,"    6  7  12 14"
  if(master)print*,"    0  15 0  20"
  if(master)print*,"    18 21 24 28"

  call a%free()
  call b%free()
  call c%free()
  call spH%free()
  deallocate(Amat,Bmat,Cmat)
  if(master)print*,""




  if(master)print*,""
  if(master)print*,"test KRON PRODUCT 3"
  allocate(Amat(3,2),Bmat(2,3))
  Amat = dble(transpose(reshape([1,2,3,4,1,0],[2,3])))
  Bmat = dble(transpose(reshape([0,5,2,6,7,3],[3,2])))
  call a%load(Amat)
  call b%load(Bmat)
  !
  if(master)print*," A = 1 2    B = 0 5 2"
  if(master)print*,"     3 4        6 7 3"
  if(master)print*,"     1 0             "
  if(master)call a%show()  
  if(master)call b%show()
  !

  allocate(Cmat(6,6))
  Cmat = dble(transpose(reshape([&
       0,5,2,0,10,4, &
       6,7,3,12,14,6,&
       0,15,6,0,20,8,&
       18,21,9,24,28,12,&
       0,5,2,0,0,0,&
       6,7,3,0,0,0],&
       [4,4])))

  if(master)print*,"c=a.x.b"  
  c = a.x.b
  if(master)call c%show()
  if(master)print*,"C = 0      5    2    0     10    4"
  if(master)print*,"    6      7    3   12     14    6"
  if(master)print*,"    0     15    6    0     20    8"    
  if(master)print*,"    18     21    9   24     28   12"    
  if(master)print*,"    0      5    2    0      0    0"    
  if(master)print*,"    6      7    3    0      0    0"




  call a%free()
  call b%free()
  call c%free()
  call spH%free()

  if(master)print*, "TEST TRANSPOSE CONJUGATE"
  call a%load(Sx+Sz)
  if(master)call a%show()
  b = hconjg(a)
  if(master)call b%show()

  if(any( a%as_matrix()-b%as_matrix() /= zero) )then
     if(master)write(*,*)"Wrong TRANSPOSE"
  else
     if(master)write(*,*)"Good TRANSPOSE"
  endif





  if(master)print*, "TEST TRANSPOSE CONJUGATE"
  call a%load(Splus)
  if(master)call a%show()
  b = a%t()
  if(master)call b%show()

  if(any( a%as_matrix()-b%as_matrix() /= zero) )then
     if(master)write(*,*)"Wrong TRANSPOSE"
  else
     if(master)write(*,*)"Good TRANSPOSE"
  endif


  deallocate(Amat,Bmat,Cmat)

  if(master)print*,""
  if(master)print*,"test KRON PRODUCT 3"

  allocate(Amat(5,5));Amat=zero
  Amat(1,2) = 1d0
  do i=2,5-1
     Amat(i,i-1) = 1d0
     Amat(i,i+1) = 1d0    
  enddo
  Amat(5,5-1) = 1d0

  allocate(Bmat(5,5))
  Bmat = dble((reshape([1,0,1,0,1,1,0,1,0,1,1,0,1,0,1,1,0,1,0,1,1,0,1,0,1],[5,5])))

  call a%load(Amat)
  call b%load(Bmat)
  !
  if(master)print*,"A"
  if(master)call a%show()
  if(master)print*,"B"
  if(master)call b%show()
  if(master)print*,""





  allocate(Cmat(5,5))
  Cmat = matmul(Amat,Bmat)
  do i=1,5
     if(master)write(*,"(5F9.3,1x)")(Cmat(i,j),j=1,5)
  enddo

  if(master)print*,""
  if(master)print*,"c=a.m.b"
  c = a.m.b
  if(master)call c%show()
  if(master)print*,c%nnz()

  if(master)print*,""
  if(master)print*,"c=a.pm.b"
  c = a.pm.b
  if(master)call c%show()
  if(master)print*,c%nnz()
  if(master)print*,""
  if(master)print*,""


  if(master)print*,"TEST APPEND TO SPARSE VECTOR"
  a = sparse(Gamma03)

  do i=1,12
     if(mod(i,2)==0)then        
        call append_sparse(Olist,a)
     else
        call append_matrix(Olist,Gamma13)
     endif
  enddo

  do i=1,size(Olist)
     if(master)print*,i,mod(i,2)==0
     if(master)call Olist(i)%show()
  enddo

  call Olist%free()









  if(master)print*,"CHECK COPY B <= A"  
  call a%free()
  call b%free()

  if(master)print*,"A:"
  a = as_sparse(kron(Gamma13,Sz))
  if(master)call a%show()

  if(master)print*,"B.copy(A)"
  call b%copy(a)
  if(master)call b%show()
  call b%free()

  if(master)print*,"B=A"
  b = a
  if(master)call b%show()





  call a%free()
  call b%free()


  if(master)print*,"CHECK WRITE/READ A"
  a = as_sparse(kron(Splus,Gamma13))
  b = as_sparse(0d0*Sz)
  if(master)print*,"A:"
  if(master)call a%show()
  if(master)print*,"B:"
  if(master)call b%show()

  if(master)print*,"write A&B:"
  if(master)call a%write(file="sp_a.dat")
  if(master)call b%write(file="sp_b.dat")

  if(master)print*,"free:"
  call a%free()
  call b%free()


  if(master)print*,"read A&B:"
  if(master)call a%read(file="sp_a.dat")
  if(master)call b%read(file="sp_b.dat")

  if(master)print*,"show A:"
  if(master)call a%show()
  if(master)print*,"show B:"
  if(master)call b%show()












  if(master)print*,"CHECK A.Bcast"
  call a%free()
  call Barrier_MPI(comm)
  if(master)then
     a = as_sparse(kron(Sz,Gamma13))
     if(master)print*,"Master a=Sz.x.G_13"
     if(master)call a%show()
     if(master)print*,""
  endif
  call Barrier_MPI(comm)

  if(rank==1)then
     print*,"Node a"
     call a%show()
  endif
  call Barrier_MPI(comm)
  if(master)print*,""
  if(master)print*,"Bcast a: 0 --> node"
  call a%bcast()
  call Barrier_MPI(comm)

  if(rank==1)then
     print*,"Node a"
     call a%show()
  endif
  call Barrier_MPI(comm)
  if(master)print*,""






  if(master)print*,"CHECK AllGather_MPI(comm,A) case 1 (Ncpu</>N[3])"
  allocate(vecA(3))
  !
  do i=1+rank,size(vecA),ncpu
     if(mod(i,2)==0)then
        print*,rank,i,"EVEN"
        vecA(i) = as_sparse(Sz*i)
     else
        print*,rank,i,"ODD"
        vecA(i) = as_sparse(S0*i)
     endif
  enddo
  call Barrier_MPI(comm)
  !
  do irank=0,ncpu-1
     if(irank==rank)then
        print *, "--- BEFORE All-Gather on rank:",rank,", ---"
        do i = 1, size(vecA)
           print *, "A(", i, ") status: ", vecA(i)%status
           if (vecA(i)%status) call vecA(i)%show()
        enddo
     endif
     call Barrier_MPI(comm)
  enddo
  !
  call AllGather_MPI(comm,vecA)
  !
  do irank=0,ncpu-1
     if(irank==rank)then
        print *, "--- AFTER All-Gather on rank:",rank,", ---"
        do i = 1, size(vecA)
           print *, "A(", i, ") status: ", vecA(i)%status
           if (vecA(i)%status) call vecA(i)%show()
        enddo
     endif
     call Barrier_MPI(comm)
  enddo
  !
  do i=1,size(vecA)
     call vecA(i)%free()
  enddo
  deallocate(vecA)
  print*,""





  if(master)print*,"CHECK AllGather_MPI(comm,A) Case 2: if()cycle"
  allocate(vecA(5))
  !
  do i=1+rank,size(vecA),ncpu
     if(i==3)then
        if(master)print*, "INFO: Index 3 is being skipped by all processes."
        cycle
     endif
     vecA(i) = as_sparse(S0*i)
  enddo
  !
  call Barrier_MPI(comm)
  !
  if(master)then
     print *, "--- BEFORE All-Gather on rank:",rank,", ---"
     do i = 1, size(vecA)
        print *, "A(", i, ") status: ", vecA(i)%status
        if (vecA(i)%status) call vecA(i)%show()
     enddo
  endif
  call Barrier_MPI(comm)
  !
  call AllGather_MPI(comm,vecA)
  !
  call Barrier_MPI(comm)
  do irank=0,ncpu-1
     if(irank==rank)then
        print *, "--- AFTER All-Gather on rank:",rank,", ---"
        do i = 1, size(vecA)
           print *, "A(", i, ") status: ", vecA(i)%status
           if (vecA(i)%status) then
              print*, "  -> Matrix A(",i,") received and initialized:"
              call vecA(i)%show()
           else
              print*, "  -> Matrix A(",i,") is correctly uninitialized."
           endif
        enddo
     endif
     call Barrier_MPI(comm)
  enddo
  !
  do i=1,size(vecA)
     call vecA(i)%free()
  enddo
  deallocate(vecA)



  call finalize_MPI()








contains


#ifdef _MPI
  subroutine sp_Bcast(comm,matrix)
    integer, intent(in)                 :: comm
    class(sparse_matrix), intent(inout) :: matrix
    integer                             :: rank, ierr, i
    integer                             :: Nrow, Ncol, Nsize
    logical                             :: master
    !
    if(.not.check_MPI())stop "sp_bcast error: check_MPI=F"
    rank   = get_Rank_MPI(comm)
    master = get_Master_MPI(comm)
    !
    if(master)Nrow = matrix%Nrow ; call Bcast_MPI(comm,Nrow)
    if(master)Ncol = matrix%Ncol ; call Bcast_MPI(comm,Ncol)
    !
    if(.not.master)call matrix%init(Nrow, Ncol)
    !
    do i=1,Nrow
       if(master)Nsize = matrix%row(i)%size
       call Bcast_MPI(comm,Nsize)
       !
       matrix%row(i)%size = Nsize !tautology for master
       if(Nsize==0)cycle
       if (.not.master) then
          if(allocated(matrix%row(i)%cols)) deallocate(matrix%row(i)%cols)
          if(allocated(matrix%row(i)%vals)) deallocate(matrix%row(i)%vals)
          allocate(matrix%row(i)%cols(Nsize))
          allocate(matrix%row(i)%vals(Nsize))
       endif
       !
       call Bcast_MPI(comm,matrix%row(i)%cols)
       call Bcast_MPI(comm,matrix%row(i)%vals)
    end do
  end subroutine sp_Bcast
#endif


  subroutine append_sparse(self,sparse)
    type(sparse_matrix),dimension(:),allocatable,intent(inout) :: self
    type(sparse_matrix)                                         :: sparse
    type(sparse_matrix),dimension(:),allocatable                :: tmp
    integer                                                     :: N
    if(allocated(self))then
       N = size(self)
       allocate(tmp(N+1))
       tmp(1:N) = self
       !This is quite unbelieavable: provided = operation among object is defined
       call move_alloc(tmp,self)
    else
       N = 1
       allocate(self(N))
    endif
    !
    self(N) = sparse
    !
  end subroutine append_sparse

  subroutine append_matrix(self,matrix)
    type(sparse_matrix),dimension(:),allocatable,intent(inout) :: self
#ifdef _CMPLX
    complex(8),dimension(:,:)                                  :: matrix
#else
    real(8),dimension(:,:)                                  :: matrix
#endif
    type(sparse_matrix),dimension(:),allocatable               :: tmp
    integer                                                    :: N
    if(allocated(self))then
       N = size(self)
       allocate(tmp(N+1))
       tmp(1:N) = self
       call move_alloc(tmp,self)
    else
       N = 1
       allocate(self(N))
    endif
    !
    self(N) = as_sparse(matrix)
    !
  end subroutine append_matrix
end program testSPARSE_MATRICES
#endif






!   function sp_filter_matrix_1(A,states) result(Ak)
!     class(sparse_matrix), intent(in)    :: A
!     integer,dimension(:),intent(in)     :: states
!     type(sparse_matrix)                 :: Ak
!     integer                             :: ii,istate,jstate
!     integer,dimension(:),allocatable    :: cols
! #ifdef _CMPLX
!     complex(8),dimension(:),allocatable :: vals
!     complex(8)                          :: val
! #else
!     real(8),dimension(:),allocatable    :: vals
!     real(8)                             :: val
! #endif
!     !
!     if(.not.A%status)stop "sp_filter_matrix_1: A.status=F"
!     !
!     call Ak%free()
!     call Ak%init(size(states),size(states))
!     !
!     allocate(vals(size(states)))
!     allocate(cols(size(states)))
!     !
!     do istate = 1,size(states)
!        i   = states(istate)
!        ii  = 0
!        vals= 0
!        cols= 0
!        do jstate=1,size(states)
!           j  = states(jstate)
!           val= A%get(i,j)
!           if(val==zero)cycle
!           ii = ii+1
!           vals(ii)=val
!           cols(ii)=jstate
!           Ak%row(istate)%Size = Ak%row(istate)%Size + 1
!        enddo
!        call append(Ak%row(istate)%vals,vals(1:ii))
!        call append(Ak%row(istate)%cols,cols(1:ii))
!     enddo
!     deallocate(vals,cols)
!   end function sp_filter_matrix_1

!   function sp_filter_matrix_2(A,Istates,Jstates) result(Ak)
!     class(sparse_matrix), intent(in) :: A
!     integer,dimension(:),intent(in)  :: Istates,Jstates
!     type(sparse_matrix)              :: Ak
!     integer                             :: ii,istate,jstate
!     integer,dimension(:),allocatable    :: cols
! #ifdef _CMPLX
!     complex(8),dimension(:),allocatable :: vals
!     complex(8)                          :: val
! #else
!     real(8),dimension(:),allocatable    :: vals
!     real(8)                             :: val
! #endif
!     !
!     call Ak%free()
!     call Ak%init(size(Istates),size(Jstates))
!     !
!     allocate(vals(size(Jstates)))
!     allocate(cols(size(Jstates)))
!     !
!     do istate = 1,size(Istates)
!        i   = Istates(istate)
!        ii  = 0
!        vals= 0
!        cols= 0
!        do jstate=1,size(Jstates)
!           j=Jstates(jstate)
!           val = A%get(i,j)
!           if(val==zero)cycle
!           ii = ii+1
!           vals(ii)=val
!           cols(ii)=jstate
!           Ak%row(istate)%Size = Ak%row(istate)%Size + 1
!           !
!        enddo
!        call append(Ak%row(istate)%vals,vals(1:ii))
!        call append(Ak%row(istate)%cols,cols(1:ii))
!     enddo
!     deallocate(vals,cols)
!   end function sp_filter_matrix_2










! !+------------------------------------------------------------------+
! !PURPOSE:  Perform simple Kroenecker product of two sparse matrices
! !AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
! !+------------------------------------------------------------------+
! function sp_kron_matrix(A,B) result(AxB)
!   type(sparse_matrix), intent(in) :: A,B
!   type(sparse_matrix)             :: AxB
!   integer                         :: i,icol,j,k,kcol,l
!   integer                         :: indx_row,indx_col
!   real(8)                      :: value
!   call AxB%free()
!   call AxB%init(a%Nrow*b%Nrow,a%Ncol*b%Ncol)
!   do i=1,A%Nrow
!      do k=1,B%Nrow

!         do icol=1,A%row(i)%size
!            j = A%row(i)%cols(icol)
!            do kcol=1,B%row(k)%size
!               l = B%row(k)%cols(kcol)
!               !
!               indx_row = k + (i-1)*B%Nrow
!               indx_col = l + (j-1)*B%Ncol
!               value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
!               call AxB%insert(value,indx_row,indx_col)
!               !
!            enddo
!         enddo
!      enddo
!   enddo
! end function sp_kron_matrix











! ! Sul master
! position = 0
! call MPI_Pack(matrix%Nrow, 1, MPI_INTEGER, buffer, buf_size, position, comm, ierr)
! call MPI_Pack(matrix%Ncol, 1, MPI_INTEGER, buffer, buf_size, position, comm, ierr)
! do i = 1, matrix%Nrow
!     call MPI_Pack(matrix%row(i)%size, 1, MPI_INTEGER, buffer, buf_size, position, comm, ierr)
!     if (matrix%row(i)%size > 0) then
!         call MPI_Pack(matrix%row(i)%cols, matrix%row(i)%size, MPI_INTEGER, ...)
!         call MPI_Pack(matrix%row(i)%vals, matrix%row(i)%size, MPI_DOUBLE_COMPLEX, ...)
!     endif
! end do
! call MPI_Bcast(buffer, position, MPI_PACKED, ...)

! ! Sui worker
! call MPI_Bcast(buffer, ..., comm, ...)
! position = 0
! call MPI_Unpack(buffer, buf_size, position, Nrow_b, 1, MPI_INTEGER, comm, ierr)
! call MPI_Unpack(buffer, buf_size, position, Ncol_b, 1, MPI_INTEGER, comm, ierr)
! call matrix%init(Nrow_b, Ncol_b)
! do i = 1, Nrow_b
!     call MPI_Unpack(buffer, buf_size, position, row_size, 1, MPI_INTEGER, ...)
!     matrix%row(i)%size = row_size
!     if (row_size > 0) then
!         allocate(...)
!         call MPI_Unpack(buffer, buf_size, position, matrix%row(i)%cols, row_size, ...)
!         call MPI_Unpack(buffer, buf_size, position, matrix%row(i)%vals, row_size, ...)
!     endif
! end do








!   call Barrier_MPI(Comm)








!   call a%free()
!   call a%init(4,4)

!   if(allocated(a%row(1)%vals))deallocate(a%row(1)%vals)
!   if(allocated(a%row(1)%cols))deallocate(a%row(1)%cols)
!   allocate(a%row(1)%vals(1))
!   allocate(a%row(1)%cols(1))

!   if(master)then
!      a = as_sparse(Gamma13)
!      call a%show()
!      print*,""
!   endif

!   call MPI_GET_ADDRESS(a%row(1)%size, MpiBlockDisp(1), ierr)
!   call MPI_GET_ADDRESS(a%row(1)%cols, MpiBlockDisp(2), ierr)
!   call MPI_GET_ADDRESS(a%row(1)%vals, MpiBlockDisp(3), ierr)


!   base=MpiBlockDisp(1)
!   MpiBlockDisp=MpiBlockDisp-base


!   mpiBlockLen(1)=1
!   mpiBlockLen(2)=1
!   mpiBlockLen(3)=1

!   mpiBlockType(1)=MPI_INTEGER
!   mpiBlockType(2)=MPI_INTEGER
! #ifdef _CMPLX
!   mpiBlockType(3)=MPI_DOUBLE_COMPLEX
! #else
!   mpiBlockType(3)=MPI_DOUBLE_PRECISION
! #endif

!   call MPI_TYPE_CREATE_STRUCT(mpiBlockNum,mpiBlockLen,MpiBlockDisp,MpiBlockType,MpiSparse_Row,ierr)
!   call MPI_TYPE_COMMIT(mpiSparse_Row,ierr)


!   call MPI_BCAST(a%row(1),1,mpiSparse_Row,0,MPI_COMM_WORLD,ierr)

!   ! a%row(1)%size=1
!   ! a%row(1)%vals(1)=1d0
!   ! a%row(1)%cols(1)=3
!   print*,""
!   call a%show()
!   call a%free()
!   print*,""
!   print*,""
!   call Barrier_MPI(comm)

!   call a%free()
!   call Barrier_MPI(comm)
!   if(master)then
!      a = as_sparse(kron(Gamma13,Sz))
!      print*,"Master a=\G_13.x.S_z:"
!      call a%show()
!      print*,""
!   endif
!   call Barrier_MPI(comm)

!   if(rank==1)then
!      print*,"Node a"
!      call a%show()
!   endif
!   call Barrier_MPI(comm)
!   if(master)print*,""
!   if(master)print*,"Bcast a: 0 --> node"
!   call sp_Bcast(comm,a)
!   call Barrier_MPI(comm)

!   if(rank==1)then
!      print*,"Node a"
!      call a%show()
!   endif
!   call Barrier_MPI(comm)
!   if(master)print*,""
