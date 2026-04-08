MODULE SITES
  USE SCIFOR, only: str,diag,zeros,kron,pauli_x,pauli_y,pauli_z,xi,operator(.kx.),zero,one,print_matrix,diag,to_lower
  USE INPUT_VARS
  USE AUX_FUNCS
  USE MATRIX_SPARSE
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  USE HLOCAL
#ifdef _MPI
  USE SCIFOR,only: check_MPI,get_master_MPI   
#endif  
  implicit none
  private

  type site
     integer                                     :: Dim=1
     type(sectors_list),dimension(:),allocatable :: sectors
     type(operators_list)                        :: operators
     character(len=:),allocatable                :: OpName
     character(len=:),allocatable                :: SiteType
   contains
     procedure,pass     :: free        => free_site
     procedure,pass     :: put         => put_op_site
     procedure,pass     :: load        => load_op_site
     procedure,pass     :: get_basis   => get_basis_site
     procedure,pass     :: set_basis   => set_basis_site
     procedure,pass     :: show        => show_site
     procedure,pass     :: is_valid    => is_valid_site
     procedure,pass     :: okey        => okey_site
     procedure,pass     :: name        => OpName_site
     procedure,pass     :: type        => SiteType_site
  end type site


  !GENERIC CONSTRUCTOR
  interface site
     module procedure :: constructor_site
  end interface site

  !EQUALITY 
  interface assignment(=)
     module procedure :: equality_site
  end interface assignment(=)


  public :: site
  public :: assignment(=)
  !TEMPLATES:
  public :: spin_site
  public :: electron_site


contains




  !##################################################################
  !##################################################################
  !       LIST CONSTRUCTOR/DESTRUCTOR
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Destructor
  !+------------------------------------------------------------------+
  subroutine free_site(self)
    class(site) :: self
    self%dim   = 1
    call self%operators%free()
    if(allocated(self%sectors))then
       call self%sectors%free()
       deallocate(self%sectors)
    endif
    if(allocated(self%OpName))deallocate(self%OpName)
    if(allocated(self%SiteType))deallocate(self%SiteType)
  end subroutine free_site



  !+------------------------------------------------------------------+
  !PURPOSE:  Intrinsic constructor
  !+------------------------------------------------------------------+
  function constructor_site(Dim,sectors,operators,opname,sitetype) result(self)
    integer,intent(in)              :: Dim
    type(sectors_list),intent(in)   :: sectors(:)
    type(operators_list),intent(in) :: operators
    character(len=:),allocatable    :: opname,sitetype
    type(site)                      :: self
    integer                         :: i
    self%Dim       = Dim
    self%operators = operators
    allocate(self%sectors(size(sectors)))
    do i=1,size(self%sectors)
       self%sectors(i)   = sectors(i)
    enddo
    allocate(self%OpName, source=opname)
    allocate(self%SiteType, source=sitetype)
  end function constructor_site








  !##################################################################
  !##################################################################
  !       PUT/LOAD
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Put a sparse operator in the site dictionary
  !+------------------------------------------------------------------+
  subroutine put_op_site(self,key,op,type)
    class(site)                    :: self
    character(len=*),intent(in)    :: key
    type(sparse_matrix),intent(in) :: op
    character(len=*),intent(in)    :: type
    call self%operators%put(str(key),op,type)
  end subroutine put_op_site



  !+------------------------------------------------------------------+
  !PURPOSE:  Load a dense operator in the site dictionary
  !+------------------------------------------------------------------+
  subroutine load_op_site(self,key,op,type)
    class(site)                          :: self
    character(len=*),intent(in)          :: key
    character(len=*),intent(in)          :: type
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in) :: op
#else
    real(8),dimension(:,:),intent(in)    :: op
#endif
    call self%operators%load(str(key),op,type)
  end subroutine load_op_site


  !+------------------------------------------------------------------+
  !PURPOSE:  Get Basis of the sector
  !+------------------------------------------------------------------+
  subroutine get_basis_site(self,basis,indx)
    class(site)      :: self
    type(tbasis)     :: basis
    integer,optional :: indx
    integer          :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))&
         stop "SET_SECTORS_BLOCK ERROR: indx out of range"
    call basis%free()
    basis  = self%sectors(indx_)%basis()
  end subroutine get_basis_site



  !+------------------------------------------------------------------+
  !PURPOSE:  Put a QN array in the site
  !+------------------------------------------------------------------+
  subroutine set_basis_site(self,basis,indx)
    class(site)      :: self
    type(tbasis)     :: basis
    integer,optional :: indx
    integer          :: indx_
    indx_=1;if(present(indx))indx_=indx
    if(indx_<1.OR.indx_>size(self%sectors))&
         stop "SET_SECTORS_SITE ERROR: indx out of range"
    self%sectors(indx) = sectors_list( basis )
  end subroutine set_basis_site


  function okey_site(self,iorb,ispin,isite,ilink) result(string)
    class(site)                  :: self
    integer,optional             :: iorb,isite,ispin
    character(len=1),optional    :: ilink
    integer                      :: iorb_,isite_,ispin_
    character(len=1)             :: ilink_
    character(len=:),allocatable :: string
    iorb_ =0;if(present(iorb))iorb_=iorb
    ispin_=0;if(present(ispin))ispin_=ispin
    isite_=0;if(present(isite))isite_=isite
    ilink_="n";if(present(ilink))ilink_=ilink
    !
    !if(iorb_==0.AND.ispin_==0)stop "Okey_Site ERROR: iorb == ispin == 0"
    string = okey(iorb_,ispin_,isite_,ilink_)
    !
  end function okey_site


  function OpName_site(self) result(string)
    class(site)                  :: self
    character(len=:),allocatable :: string
    allocate(string, source=self%OpName)
  end function OpName_site


  function SiteType_site(self) result(string)
    class(site)                  :: self
    character(len=:),allocatable :: string
    ! allocate(string, source=self%SiteType)
    string = to_lower(str(self%SiteType))
  end function SiteType_site



  !##################################################################
  !##################################################################
  !              OPERATIONS / ASSIGNEMENTS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Equality between two sites
  !+------------------------------------------------------------------+
  subroutine equality_site(A,B)
    type(site),intent(inout) :: A
    type(site),intent(in)    :: B
    integer                  :: i
    call A%free
    A%Dim       = B%Dim
    A%operators = B%operators
    allocate(A%sectors(size(B%sectors)))
    do i=1,size(A%sectors)
       A%sectors(i)   = B%sectors(i)
    enddo
    allocate(A%OpName, source=B%OpName)
    allocate(A%SiteType, source=B%SiteType)
  end subroutine equality_site


  function is_valid_site(self) result(bool)
    class(site)                           :: self
    logical                               :: bool
    integer,dimension(size(self%sectors)) :: Dims
    integer                               :: i
    bool = self%operators%is_valid(self%Dim)
    do i=1,size(self%sectors)
       Dims = dim(self%sectors(i))
    enddo
    bool=bool.AND.(self%dim==product(Dims))
  end function is_valid_site




  !##################################################################
  !##################################################################
  !              SHOW 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Pretty print the site structure
  !+------------------------------------------------------------------+
  subroutine show_site(self,fmt)
    class(site)               :: self
    character(len=*),optional :: fmt
    character(len=32)         :: fmt_
    integer                   :: i
    fmt_=str(show_fmt);if(present(fmt))fmt_=str(fmt)
    write(*,"(A14,I6)")"Site Dim     =",self%Dim
    write(*,"(A15,A)")"Site Type    = ",self%SiteType
    do i=1,size(self%sectors)
       write(*,"(A14,I6)")"Site Sectors =",i
       call self%sectors(i)%show()
    enddo
    write(*,"(A15,A)")"Op Name    = ",self%OpName
    write(*,"(A14)")"Site Operators:"
    call self%operators%show(fmt=fmt_)
  end subroutine show_site





  !##################################################################
  !##################################################################
  !           SPIN PRESET
  !##################################################################
  !##################################################################
  function spin_site(sun,hvec) result(self)
    integer                               :: sun
    real(8),dimension(2),optional         :: hvec
    type(site)                            :: self
    character(len=:),allocatable          :: key
    integer                               :: ispin
#ifdef _CMPLX
    complex(8),dimension(2)               :: h_
    complex(8),dimension(:,:),allocatable :: H,Sz,Sp,Sx
#else
    real(8),dimension(2)                  :: h_
    real(8),dimension(:,:),allocatable    :: H,Sz,Sp,Sx
#endif
logical :: master=.true.
    !
#ifdef _MPI
    if(check_MPI())master  = get_master_MPI()
#endif
    !
#ifdef _DEBUG
    if(master)write(LOGfile,*)"DEBUG: Spin SITE with SU"//str(sun)
#endif
    h_ = zero; if(present(hvec))h_=hvec
    !
    call self%free()
    self%SiteType="SPIN"
    self%OpName="S"
    !
    select case(sun)
    case default;stop "spin_site ERROR: SU(N) value not supported"
    case (2)
       self%Dim = 2
       !
       allocate(Sz(2,2));Sz=reshape([one,zero,zero,-one]/2,[2,2])
       allocate(Sp(2,2));Sp=reshape([zero,zero,one,zero],[2,2])
       allocate(H(2,2));H=h_(1)*pauli_x+h_(2)*pauli_z
       !
       !> Build Sectors:
       allocate(self%sectors(1))
       self%sectors(1) = sectors_list( tbasis([0.5d0,-0.5d0],qdim=1) )
    case(3)
       self%Dim = 3
       !
       allocate(Sz(3,3));Sz=diag([one,zero,-one])
       allocate(Sp(3,3));Sp=zero;Sp(1,2)=one*sqrt(2d0);Sp(2,3)=one*sqrt(2d0)
       allocate(Sx(3,3));Sx=reshape([zero,one,zero,one,zero,one,zero,one,zero],[3,3])/sqrt(2d0)
       allocate(H(3,3));H=h_(1)*Sx+h_(2)*Sz
       !
       !> Build Sectors:       
       allocate(self%sectors(1))
       self%sectors(1) = sectors_list( tbasis([1d0,0d0,-1d0],qdim=1) )
    end select
    !
    !> Build the H operator (it contains local magnetic field)
    call self%put("H",sparse(H),'bosonic')
    !> Build all the S operators (Sz=S(spin=1), S+=S(spin=2), S-=H.c. S+ )
    call self%put(self%OpName//self%okey(0,1,ilink="n"),sparse(Sz),"bosonic") !S_z_next
    call self%put(self%OpName//self%okey(0,2,ilink="n"),sparse(Sp),"bosonic") !S_+_next
    if(PBCdmrg)then
       call self%put(self%OpName//self%okey(0,1,ilink="p"),sparse(Sz),"bosonic") !S_z_prev
       call self%put(self%OpName//self%okey(0,2,ilink="p"),sparse(Sp),"bosonic") !S_+_prev          
    endif
    !
  end function spin_site







  !##################################################################
  !##################################################################
  !           ELECTRON = FERMION SU(2) PRESET
  !##################################################################
  !##################################################################
  !Fock = [|0,0>,|1,0>,|0,1>,|1,1>] <- |up,dw> <- cycle UP first 1_dw x 1_up
  function electron_site(hloc) result(self)
#ifdef _CMPLX
    complex(8),dimension(Nspin*Norb,Nspin*Norb),optional :: hloc  
    complex(8),dimension(Nspin*Norb,Nspin*Norb)          :: hloc_
    complex(8),dimension(:,:),allocatable                :: H,P
    complex(8),dimension(:,:),allocatable                :: Op
#else
    real(8),dimension(Nspin*Norb,Nspin*Norb),optional    :: hloc  
    real(8),dimension(Nspin*Norb,Nspin*Norb)             :: hloc_
    real(8),dimension(:,:),allocatable                   :: H,P
    real(8),dimension(:,:),allocatable                   :: Op
#endif
    type(site)                                           :: self
    integer,dimension(:),allocatable                     :: Basis
    integer                                              :: iorb,ispin
    character(len=:),allocatable                         :: key
    logical :: master=.true.
    !
#ifdef _MPI
    if(check_MPI())master  = get_master_MPI()
#endif
  
#ifdef _DEBUG
    if(master)write(LOGfile,*)"DEBUG: ELECTRON SITE"
#endif
    call self%free()
    self%SiteType="FERMION"
    self%OpName="C"
    !
    !
    hloc_ = zeros(Nspin*Norb,Nspin*Norb);if(present(hloc))hloc_=hloc
    call Init_LocalFock_Space()
    !
    self%Dim = 4**Norb
    !
    !Set local H operator
    if(master)write(LOGfile,"(A,I3,A1,I3,A)")&
         "Using Hloc, shape=[",Nspin*Norb,",",Nspin*Norb,"]"
    if(master)call print_matrix(Hloc_)
    if(master)write(LOGfile,"(A)")""
    !
    !> Build local Hamiltonian:
    H = build_Hlocal_operator(hloc_)
    call self%put("H",sparse(H),"bosonic")
    !> Build all the C operators and sign P operators: OBC/PBC_next
    do ispin=1,2
       do iorb=1,Norb
          Op = build_C_operator(iorb,ispin)
          Key= self%OpName//self%okey(iorb,ispin,ilink="n")
          call self%put(Key,sparse(Op),"fermionic")
       enddo
    enddo
    P = Build_FermionicSign()
    Key= "P"//self%okey(0,0,ilink="n")
    call self%put(key,sparse(P),"psign")
    !> Build all the C operators and sign P operators: PBC_prev
    if(PBCdmrg)then
       do ispin=1,2
          do iorb=1,Norb
             Op = build_C_operator(iorb,ispin)
             Key= self%OpName//self%okey(iorb,ispin,ilink="p")
             call self%put(Key,sparse(Op),"fermionic")
          enddo
       enddo
       P = Build_FermionicSign()
       Key= "P"//self%okey(0,0,ilink="p")
       call self%put(key,sparse(P),"psign")
    endif
    !
    !> Build QN for the local basis:
    allocate(self%sectors(1))
    Basis = Build_BasisStates()
    self%sectors(1) = sectors_list( tbasis(Basis,Qdim=Get_QNdimension()) )
    !
    call Delete_LocalFock_Space()
    !
  end function electron_site




END MODULE SITES








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
program testSITES
  USE SCIFOR
  USE INPUT_VARS
  USE MATRIX_SPARSE
  USE LIST_OPERATORS
  USE TUPLE_BASIS
  USE LIST_SECTORS
  USE SITES
  implicit none
  type(site)                          :: my_site,a,b
  type(tbasis)                        :: sz_basis
#ifdef _CMPLX
  complex(8),dimension(2,2),parameter :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  complex(8),dimension(2,2),parameter :: S0=pauli_0
  complex(8),dimension(2,2),parameter :: Sz=pauli_z
  complex(8),dimension(2,2),parameter :: Sx=pauli_x
  complex(8),dimension(2,2),parameter :: Splus=reshape([zero,zero,one,zero],[2,2])
  complex(8),dimension(4,4)           :: Gamma13,Gamma03
#else
  real(8),dimension(2,2),parameter    :: Hzero=reshape([zero,zero,zero,zero],[2,2])
  real(8),dimension(2,2),parameter    :: S0=pauli_0
  real(8),dimension(2,2),parameter    :: Sz=pauli_z
  real(8),dimension(2,2),parameter    :: Sx=pauli_x
  real(8),dimension(2,2),parameter    :: Splus=reshape([zero,zero,one,zero],[2,2])
  real(8),dimension(4,4)              :: Gamma13,Gamma03
#endif
  Gamma13=kron(Sx,Sz)
  Gamma03=kron(S0,Sz)

  sz_basis = tbasis([0.5d0,-0.5d0],Qdim=1)

  my_site = site(&
       dim      = 2, &
       sectors  = [sectors_list(sz_basis)],&
       operators= operators_list(['H0','Sz','Sp'],&
       [sparse(Hzero),sparse(Sz),sparse(Splus)],['b   ','sign','bose']),&
       opname='S',&
       sitetype='spin')
  print*,"Is site valid:",my_site%is_valid()
  call my_site%show()



  a = my_site

  print*,"Test =: a=my_site"
  call a%show()
  print*,a%is_valid()
  print*,"modify a, check a and my_site"
  call a%put("G5",sparse(Gamma03),'b')
  print*,"Is A site valid:",a%is_valid()
  print*,"Is My_site site valid:",my_site%is_valid()


  print*,"Test PAULI SITE"
  b = spin_site(2)
  call b%show()
  call b%free()


  call read_input("DMRG.conf")
  b = electron_site()
  call b%show()
  print*,"Is site valid:",b%is_valid()
  call b%free()
end program testSITES
#endif
