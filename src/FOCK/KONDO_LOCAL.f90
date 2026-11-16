MODULE KONDO_LOCAL
  USE SCIFOR, only: zero,one,eye,kron,diag,to_lower
  implicit none
  private

  integer,parameter :: KondoFixedDim = 8
  integer,parameter :: KondoMobileDim=12

  ! Build the local Hamiltonian of one conduction orbital coupled to a
  ! permanently occupied spin-1/2 impurity.  The basis is
  ! |S_imp> tensor |n_up,n_down>, with the electronic states ordered as
  ! |0>, |up>, |down>, |up down>.  Jz multiplies s_z.S_z, Jxy multiplies
  ! (s_+.S_- + s_-.S_+)/2, and optional Vk multiplies the conduction
  ! density.  The returned matrix has shape [8,8].
  public :: Build_KondoFixed_Hamiltonian

  ! Return the conduction-electron annihilation operator for spin=1 (up)
  ! or spin=2 (down), including the identity on the impurity spin.  The
  ! returned matrix has shape [8,8]; invalid spin values are fatal.
  public :: Build_KondoFixed_C

  ! Return the local fermion-parity operator (-1)^(n_up+n_down).  The
  ! localized spin does not contribute to fermionic parity.
  public :: Build_KondoFixed_Parity

  ! Return the conduction-electron or impurity-spin operator selected by
  ! component='z', '+' or '-'.  The returned matrix has shape [8,8]; an
  ! unsupported component is fatal.
  public :: Build_KondoFixed_ElectronSpin
  public :: Build_KondoFixed_ImpuritySpin

  ! Return the impurity occupation operator.  For the fixed-spin basis it
  ! is the [8,8] identity because the impurity is present in every state.
  public :: Build_KondoFixed_Number

  ! Return the flattened QN table consumed by tbasis(...,Qdim=2).  The two
  ! entries per state are total (N_up,N_down), including the impurity spin.
  ! The result has size 16 and follows the local basis ordering above.
  public :: Build_KondoFixed_BasisQN

  ! Build the local Hamiltonian of one conduction orbital coupled to an
  ! empty or singly occupied mobile impurity.  The basis is
  ! |0,up,down>_imp tensor |0,up,down,up down>_c.  Jz, Jxy and optional Vk
  ! use the fixed-site conventions; optional Eimp multiplies Nimp.  The
  ! empty-impurity block is decoupled and the occupied block reproduces the
  ! fixed eight-state Hamiltonian up to Eimp.  The result has shape [12,12].
  public :: Build_KondoMobile_Hamiltonian

  ! Return the conduction-electron or projected-impurity annihilation
  ! operator for spin=1 (up) or spin=2 (down).  The impurity operator
  ! includes the electronic parity required to anticommute with C.  Both
  ! results have shape [12,12]; invalid spin values are fatal.
  public :: Build_KondoMobile_C
  public :: Build_KondoMobile_D

  ! Return the total local fermion parity (-1)^(Nc+Nimp).  This parity is
  ! shared by conduction and impurity hopping and has shape [12,12].
  public :: Build_KondoMobile_Parity

  ! Return the conduction-electron or projected-impurity spin component
  ! selected by component='z', '+' or '-'.  The impurity spin vanishes in
  ! the empty sector.  Invalid components terminate with error stop.
  public :: Build_KondoMobile_ElectronSpin
  public :: Build_KondoMobile_ImpuritySpin

  ! Return the impurity occupation projector with eigenvalues zero on the
  ! four empty states and one on the eight occupied states.
  public :: Build_KondoMobile_Number

  ! Return the flattened QN table consumed by tbasis(...,Qdim=3).  The
  ! charges are (Nup_total,Ndown_total,Nimp), with total spin populations
  ! including both conduction electrons and the projected impurity.
  public :: Build_KondoMobile_BasisQN

contains


  function Build_KondoFixed_Hamiltonian(Jz,Jxy,Vk) result(H)
    real(8),intent(in)                    :: Jz,Jxy
    real(8),intent(in),optional           :: Vk
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: H
#else
    real(8),dimension(:,:),allocatable    :: H
#endif
    real(8)                               :: Vk_
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Sez,Sep,Sem,Siz,Sip,Sim,Nc
#else
    real(8),dimension(:,:),allocatable    :: Sez,Sep,Sem,Siz,Sip,Sim,Nc
#endif
    !
    Vk_=0d0
    if(present(Vk))Vk_=Vk
    !
    Sez=Build_KondoFixed_ElectronSpin('z')
    Sep=Build_KondoFixed_ElectronSpin('+')
    Sem=Build_KondoFixed_ElectronSpin('-')
    Siz=Build_KondoFixed_ImpuritySpin('z')
    Sip=Build_KondoFixed_ImpuritySpin('+')
    Sim=Build_KondoFixed_ImpuritySpin('-')
    Nc =Build_KondoFixed_Density()
    !
    H=Jz*matmul(Sez,Siz)+0.5d0*Jxy*&
         (matmul(Sep,Sim)+matmul(Sem,Sip))+Vk_*Nc
    !
  end function Build_KondoFixed_Hamiltonian





  function Build_KondoFixed_C(spin) result(C)
    integer,intent(in)                    :: spin
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: C
    complex(8),dimension(4,4)             :: Ce
#else
    real(8),dimension(:,:),allocatable    :: C
    real(8),dimension(4,4)                :: Ce
#endif
    !
    if(spin<1.OR.spin>2)error stop 'Build_KondoFixed_C: spin must be 1 (up) or 2 (down)'
    !
    Ce=zero
    select case(spin)
    case(1)
       Ce(1,2)=one
       Ce(3,4)=one
    case(2)
       Ce(1,3)=one
       Ce(2,4)=-one
    end select
    C=kron(eye(2),Ce)
    !
  end function Build_KondoFixed_C







  function Build_KondoFixed_Parity() result(P)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: P
    complex(8),dimension(4,4)             :: Pe
#else
    real(8),dimension(:,:),allocatable    :: P
    real(8),dimension(4,4)                :: Pe
#endif
    !
    Pe=diag([one,-one,-one,one])
    P =kron(eye(2),Pe)
    !
  end function Build_KondoFixed_Parity







  function Build_KondoFixed_ElectronSpin(component) result(S)
    character(len=*),intent(in)           :: component
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: S
    complex(8),dimension(:,:),allocatable :: Cup,Cdw,Se
#else
    real(8),dimension(:,:),allocatable    :: S
    real(8),dimension(:,:),allocatable    :: Cup,Cdw,Se
#endif
    !
    !
    Cup=Build_Electron_C(1)
    Cdw=Build_Electron_C(2)
    select case(to_lower(component(1:1)))
    case('z');Se=0.5d0*(matmul(Dagger(Cup),Cup)-matmul(Dagger(Cdw),Cdw))
    case('+');Se=matmul(Dagger(Cup),Cdw)
    case('-');Se=matmul(Dagger(Cdw),Cup)
    case default
       error stop "Build_KondoFixed_ElectronSpin: component must be 'z', '+' or '-'"
    end select
    !
    S=kron(eye(2),Se)
    !
  end function Build_KondoFixed_ElectronSpin




  function Build_KondoFixed_ImpuritySpin(component) result(S)
    character(len=*),intent(in)           :: component
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: S
    complex(8),dimension(2,2)             :: Si
#else
    real(8),dimension(:,:),allocatable    :: S
    real(8),dimension(2,2)                :: Si
#endif
    !
    Si=zero
    select case(to_lower(component(1:1)))
    case('z');Si=diag([0.5d0*one,-0.5d0*one])
    case('+');Si(1,2)=one
    case('-');Si(2,1)=one
    case default
       error stop "Build_KondoFixed_ImpuritySpin: component must be 'z', '+' or '-'"
    end select
    !
    S=kron(Si,eye(4))
    !
  end function Build_KondoFixed_ImpuritySpin






  function Build_KondoFixed_Number() result(Nimp)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Nimp
#else
    real(8),dimension(:,:),allocatable    :: Nimp
#endif
    !
    Nimp=eye(KondoFixedDim)
    !
  end function Build_KondoFixed_Number






  function Build_KondoFixed_BasisQN() result(QN)
    integer,dimension(:),allocatable :: QN
    integer                          :: istate,impurity,electron
    integer                          :: Nup,Ndown
    !
    allocate(QN(2*KondoFixedDim))
    do istate=1,KondoFixedDim
       impurity=(istate-1)/4+1
       electron=mod(istate-1,4)
       Nup=merge(1,0,btest(electron,0))+merge(1,0,impurity==1)
       Ndown=merge(1,0,btest(electron,1))+merge(1,0,impurity==2)
       QN(2*istate-1)=Nup
       QN(2*istate)=Ndown
    enddo
    !
  end function Build_KondoFixed_BasisQN






  function Build_KondoFixed_Density() result(Nc)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Nc
    complex(8),dimension(:,:),allocatable :: Cup,Cdw
#else
    real(8),dimension(:,:),allocatable    :: Nc
    real(8),dimension(:,:),allocatable    :: Cup,Cdw
#endif
    !
    Cup=Build_KondoFixed_C(1)
    Cdw=Build_KondoFixed_C(2)
    Nc=matmul(Dagger(Cup),Cup)+matmul(Dagger(Cdw),Cdw)
    !
  end function Build_KondoFixed_Density










  function Build_KondoMobile_Hamiltonian(Jz,Jxy,Vk,Eimp) result(H)
    real(8),intent(in)                    :: Jz,Jxy
    real(8),intent(in),optional           :: Vk,Eimp
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: H
    complex(8),dimension(:,:),allocatable :: Sez,Sep,Sem
    complex(8),dimension(:,:),allocatable :: Siz,Sip,Sim,Nc,Nimp
#else
    real(8),dimension(:,:),allocatable    :: H
    real(8),dimension(:,:),allocatable    :: Sez,Sep,Sem
    real(8),dimension(:,:),allocatable    :: Siz,Sip,Sim,Nc,Nimp
#endif
    real(8)                               :: Vk_,Eimp_
    !
    ! Select optional local energies.
    Vk_  =0d0
    Eimp_=0d0
    if(present(Vk))  Vk_  =Vk
    if(present(Eimp))Eimp_=Eimp
    !
    ! Construct the Kondo exchange and local potential terms.
    Sez =Build_KondoMobile_ElectronSpin('z')
    Sep =Build_KondoMobile_ElectronSpin('+')
    Sem =Build_KondoMobile_ElectronSpin('-')
    Siz =Build_KondoMobile_ImpuritySpin('z')
    Sip =Build_KondoMobile_ImpuritySpin('+')
    Sim =Build_KondoMobile_ImpuritySpin('-')
    Nc  =Build_KondoMobile_Density()
    Nimp=Build_KondoMobile_Number()
    H=Jz*matmul(Sez,Siz)+0.5d0*Jxy*&
         (matmul(Sep,Sim)+matmul(Sem,Sip))+Vk_*Nc+Eimp_*Nimp
    !
  end function Build_KondoMobile_Hamiltonian





  function Build_KondoMobile_C(spin) result(C)
    integer,intent(in)                    :: spin
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: C
    complex(8),dimension(4,4)             :: Ce
#else
    real(8),dimension(:,:),allocatable    :: C
    real(8),dimension(4,4)                :: Ce
#endif
    !
    if(spin<1.or.spin>2)&
         error stop 'Build_KondoMobile_C: spin must be 1 (up) or 2 (down)'
    !
    Ce=Build_Electron_C(spin)
    C =kron(eye(3),Ce)
    !
  end function Build_KondoMobile_C





  function Build_KondoMobile_D(spin) result(D)
    integer,intent(in)                    :: spin
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: D
    complex(8),dimension(3,3)             :: Di
    complex(8),dimension(4,4)             :: Pe
#else
    real(8),dimension(:,:),allocatable    :: D
    real(8),dimension(3,3)                :: Di
    real(8),dimension(4,4)                :: Pe
#endif
    !
    if(spin<1.or.spin>2)&
         error stop 'Build_KondoMobile_D: spin must be 1 (up) or 2 (down)'
    !
    ! The electronic modes precede the projected impurity modes in the
    ! local Fock ordering.  Pe supplies their Jordan-Wigner sign.
    Di=zero
    Di(1,spin+1)=one
    Pe=diag([one,-one,-one,one])
    D =kron(Di,Pe)
    !
  end function Build_KondoMobile_D





  function Build_KondoMobile_Parity() result(P)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: P
    complex(8),dimension(3,3)             :: Pi
    complex(8),dimension(4,4)             :: Pe
#else
    real(8),dimension(:,:),allocatable    :: P
    real(8),dimension(3,3)                :: Pi
    real(8),dimension(4,4)                :: Pe
#endif
    !
    Pi=diag([one,-one,-one])
    Pe=diag([one,-one,-one,one])
    P =kron(Pi,Pe)
    !
  end function Build_KondoMobile_Parity





  function Build_KondoMobile_ElectronSpin(component) result(S)
    character(len=*),intent(in)           :: component
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: S
    complex(8),dimension(:,:),allocatable :: Cup,Cdw,Se
#else
    real(8),dimension(:,:),allocatable    :: S
    real(8),dimension(:,:),allocatable    :: Cup,Cdw,Se
#endif
    !
    Cup=Build_Electron_C(1)
    Cdw=Build_Electron_C(2)
    select case(to_lower(component(1:1)))
    case('z');Se=0.5d0*(matmul(Dagger(Cup),Cup)-&
         matmul(Dagger(Cdw),Cdw))
    case('+');Se=matmul(Dagger(Cup),Cdw)
    case('-');Se=matmul(Dagger(Cdw),Cup)
    case default
       error stop &
            "Build_KondoMobile_ElectronSpin: component must be 'z', '+' or '-'"
    end select
    S=kron(eye(3),Se)
    !
  end function Build_KondoMobile_ElectronSpin





  function Build_KondoMobile_ImpuritySpin(component) result(S)
    character(len=*),intent(in)           :: component
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: S
    complex(8),dimension(3,3)             :: Si
#else
    real(8),dimension(:,:),allocatable    :: S
    real(8),dimension(3,3)                :: Si
#endif
    !
    Si=zero
    select case(to_lower(component(1:1)))
    case('z');Si=diag([zero,0.5d0*one,-0.5d0*one])
    case('+');Si(2,3)=one
    case('-');Si(3,2)=one
    case default
       error stop &
            "Build_KondoMobile_ImpuritySpin: component must be 'z', '+' or '-'"
    end select
    S=kron(Si,eye(4))
    !
  end function Build_KondoMobile_ImpuritySpin





  function Build_KondoMobile_Number() result(Nimp)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Nimp
    complex(8),dimension(3,3)             :: Ni
#else
    real(8),dimension(:,:),allocatable    :: Nimp
    real(8),dimension(3,3)                :: Ni
#endif
    !
    Ni  =diag([zero,one,one])
    Nimp=kron(Ni,eye(4))
    !
  end function Build_KondoMobile_Number





  function Build_KondoMobile_BasisQN() result(QN)
    integer,dimension(:),allocatable :: QN
    integer                          :: istate,impurity,electron
    integer                          :: Nup,Ndown,Nimp
    !
    allocate(QN(3*KondoMobileDim))
    do istate=1,KondoMobileDim
       impurity=(istate-1)/4
       electron=mod(istate-1,4)
       Nup=merge(1,0,btest(electron,0))+merge(1,0,impurity==1)
       Ndown=merge(1,0,btest(electron,1))+merge(1,0,impurity==2)
       Nimp=merge(0,1,impurity==0)
       QN(3*istate-2)=Nup
       QN(3*istate-1)=Ndown
       QN(3*istate)=Nimp
    enddo
    !
  end function Build_KondoMobile_BasisQN





  function Build_KondoMobile_Density() result(Nc)
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Nc
    complex(8),dimension(:,:),allocatable :: Cup,Cdw
#else
    real(8),dimension(:,:),allocatable    :: Nc
    real(8),dimension(:,:),allocatable    :: Cup,Cdw
#endif
    !
    Cup=Build_KondoMobile_C(1)
    Cdw=Build_KondoMobile_C(2)
    Nc=matmul(Dagger(Cup),Cup)+matmul(Dagger(Cdw),Cdw)
    !
  end function Build_KondoMobile_Density






  function Build_Electron_C(spin) result(C)
    integer,intent(in)        :: spin
#ifdef _CMPLX
    complex(8),dimension(4,4) :: C
#else
    real(8),dimension(4,4)    :: C
#endif
    !
    if(spin<1.OR.spin>2)&
         error stop 'Build_Electron_C: spin must be 1 (up) or 2 (down)'
    !
    C=zero
    select case(spin)
    case(1)
       C(1,2)=one
       C(3,4)=one
    case(2)
       C(1,3)=one
       C(2,4)=-one
    end select
    !
  end function Build_Electron_C




  
  function Dagger(A) result(Adgr)
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in) :: A
    complex(8),dimension(size(A,2),size(A,1)) :: Adgr
#else
    real(8),dimension(:,:),intent(in)    :: A
    real(8),dimension(size(A,2),size(A,1))    :: Adgr
#endif
    !
#ifdef _CMPLX
    Adgr=conjg(transpose(A))
#else
    Adgr=transpose(A)
#endif
    !
  end function Dagger


END MODULE KONDO_LOCAL




!##################################################################
!##################################################################
!                          TEST
!##################################################################
!##################################################################
#ifdef _TEST
program testKONDO_LOCAL
  USE SCIFOR, only: zero,one,eye
  USE KONDO_LOCAL
  implicit none

  real(8),parameter                 :: tolerance=1d-12
  integer,dimension(:),allocatable :: basis
  integer,dimension(2,8)           :: qn
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable :: H,Cup,Cdown,P
  complex(8),dimension(:,:),allocatable :: Sez,Sep,Sem,Siz,Sip,Sim
  complex(8),dimension(:,:),allocatable :: expected_H,identity8,zero8
  complex(8),dimension(8)               :: singlet,triplet
#else
  real(8),dimension(:,:),allocatable    :: H,Cup,Cdown,P
  real(8),dimension(:,:),allocatable    :: Sez,Sep,Sem,Siz,Sip,Sim
  real(8),dimension(:,:),allocatable    :: expected_H,identity8,zero8
  real(8),dimension(8)                  :: singlet,triplet
#endif
  !
  ! Build the fixed-impurity local problem.
  H     =Build_KondoFixed_Hamiltonian(1d0,1d0)
  Cup   =Build_KondoFixed_C(1)
  Cdown =Build_KondoFixed_C(2)
  P     =Build_KondoFixed_Parity()
  Sez   =Build_KondoFixed_ElectronSpin('z')
  Sep   =Build_KondoFixed_ElectronSpin('+')
  Sem   =Build_KondoFixed_ElectronSpin('-')
  Siz   =Build_KondoFixed_ImpuritySpin('z')
  Sip   =Build_KondoFixed_ImpuritySpin('+')
  Sim   =Build_KondoFixed_ImpuritySpin('-')
  basis =Build_KondoFixed_BasisQN()
  qn    =reshape(basis,[2,8])
  identity8=one*eye(8)
  allocate(expected_H(8,8),zero8(8,8))
  expected_H=zero
  zero8=zero
  expected_H(2,2)=0.25d0
  expected_H(3,3)=-0.25d0
  expected_H(3,6)=0.5d0
  expected_H(6,3)=0.5d0
  expected_H(6,6)=-0.25d0
  expected_H(7,7)=0.25d0
  !
  ! Check dimensions, Hermiticity and the singlet/triplet spectrum.
  call assert_true(all(shape(H)==[8,8]),'Hamiltonian dimension')
  call assert_matrix_close(H,Dagger(H),'Hamiltonian Hermiticity')
  call assert_matrix_close(H,expected_H,'Kondo Hamiltonian matrix')
  singlet=zero
  singlet(3)=one/sqrt(2d0)
  singlet(6)=-one/sqrt(2d0)
  triplet=zero
  triplet(3)=one/sqrt(2d0)
  triplet(6)=one/sqrt(2d0)
  call assert_state_close(matmul(H,singlet),-0.75d0*singlet,&
       'singlet energy')
  call assert_state_close(matmul(H,triplet),0.25d0*triplet,&
       'triplet energy')
  !
  ! Check the fermion and parity algebras.
  call assert_matrix_close(matmul(Cup,Dagger(Cup))+&
       matmul(Dagger(Cup),Cup),identity8,'Cup anticommutator')
  call assert_matrix_close(matmul(Cdown,Dagger(Cdown))+&
       matmul(Dagger(Cdown),Cdown),identity8,'Cdown anticommutator')
  call assert_matrix_close(matmul(Cup,Cdown)+matmul(Cdown,Cup),zero8,&
       'cross-spin anticommutator')
  call assert_matrix_close(matmul(P,P),identity8,'parity square')
  call assert_matrix_close(matmul(P,Cup)+matmul(Cup,P),zero8,'parity Cup')
  call assert_matrix_close(matmul(P,Cdown)+matmul(Cdown,P),zero8,&
       'parity Cdown')
  !
  ! Check the two commuting SU(2) algebras.
  call assert_matrix_close(matmul(Sep,Sem)-matmul(Sem,Sep),2d0*Sez,&
       'electron spin algebra')
  call assert_matrix_close(matmul(Sip,Sim)-matmul(Sim,Sip),2d0*Siz,&
       'impurity spin algebra')
  call assert_matrix_close(matmul(Sep,Sip)-matmul(Sip,Sep),zero8,&
       'electron/impurity commutation')
  !
  ! Check QN shifts from every nonzero matrix element.
  call assert_operator_dq(Cup,qn,[-1,0],'Cup dq')
  call assert_operator_dq(Cdown,qn,[0,-1],'Cdown dq')
  call assert_operator_dq(Sip,qn,[1,-1],'Simp+ dq')
  call assert_operator_dq(Sep,qn,[1,-1],'Sele+ dq')
  call assert_operator_dq(H,qn,[0,0],'Hamiltonian dq')
  !
  ! Check the mobile-impurity local representation independently.  The
  ! explicit comparison with the fixed occupied block is a separate test.
  call test_mobile_algebra()
  call test_mobile_fixed_embedding()
  !
  print *,'KONDO_LOCAL TEST: PASS'

contains


  subroutine test_mobile_algebra()
    integer,dimension(:),allocatable      :: mobile_basis
    integer,dimension(3,12)               :: mobile_qn
    integer                               :: istate
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: H0,Himp
    complex(8),dimension(:,:),allocatable :: Cupm,Cdownm,Dup,Ddown
    complex(8),dimension(:,:),allocatable :: Pm,Nimp
    complex(8),dimension(:,:),allocatable :: Sezm,Sepm,Semm
    complex(8),dimension(:,:),allocatable :: Sizm,Sipm,Simm
    complex(8),dimension(:,:),allocatable :: identity12,zero12
    complex(8),dimension(:,:),allocatable :: projected_up,projected_down
#else
    real(8),dimension(:,:),allocatable    :: H0,Himp
    real(8),dimension(:,:),allocatable    :: Cupm,Cdownm,Dup,Ddown
    real(8),dimension(:,:),allocatable    :: Pm,Nimp
    real(8),dimension(:,:),allocatable    :: Sezm,Sepm,Semm
    real(8),dimension(:,:),allocatable    :: Sizm,Sipm,Simm
    real(8),dimension(:,:),allocatable    :: identity12,zero12
    real(8),dimension(:,:),allocatable    :: projected_up,projected_down
#endif
    !
    ! Construct every operator entering the mobile local contract.
    H0      =Build_KondoMobile_Hamiltonian(1d0,1d0)
    Himp    =Build_KondoMobile_Hamiltonian(1d0,1d0,Eimp=0.37d0)
    Cupm    =Build_KondoMobile_C(1)
    Cdownm  =Build_KondoMobile_C(2)
    Dup     =Build_KondoMobile_D(1)
    Ddown   =Build_KondoMobile_D(2)
    Pm      =Build_KondoMobile_Parity()
    Nimp    =Build_KondoMobile_Number()
    Sezm    =Build_KondoMobile_ElectronSpin('z')
    Sepm    =Build_KondoMobile_ElectronSpin('+')
    Semm    =Build_KondoMobile_ElectronSpin('-')
    Sizm    =Build_KondoMobile_ImpuritySpin('z')
    Sipm    =Build_KondoMobile_ImpuritySpin('+')
    Simm    =Build_KondoMobile_ImpuritySpin('-')
    mobile_basis=Build_KondoMobile_BasisQN()
    mobile_qn=reshape(mobile_basis,[3,12])
    identity12=one*eye(12)
    allocate(zero12(12,12),projected_up(12,12),projected_down(12,12))
    zero12=zero
    projected_up=zero
    projected_down=zero
    do istate=1,8
       projected_up(istate,istate)=one
    enddo
    do istate=1,4
       projected_down(istate,istate)=one
    enddo
    do istate=9,12
       projected_down(istate,istate)=one
    enddo
    !
    ! Check dimensions, Hermiticity and the optional impurity energy.
    call assert_true(all(shape(H0)==[12,12]),&
         'mobile Hamiltonian dimension')
    call assert_matrix_close(H0,Dagger(H0),&
         'mobile Hamiltonian Hermiticity')
    call assert_matrix_close(Himp-H0,0.37d0*Nimp,&
         'mobile impurity energy')
    call assert_matrix_close(H0(1:4,:),zero12(1:4,:),&
         'mobile empty rows decoupled')
    call assert_matrix_close(H0(:,1:4),zero12(:,1:4),&
         'mobile empty columns decoupled')
    !
    ! Check the canonical conduction-electron algebra.
    call assert_matrix_close(matmul(Cupm,Dagger(Cupm))+&
         matmul(Dagger(Cupm),Cupm),identity12,&
         'mobile Cup anticommutator')
    call assert_matrix_close(matmul(Cdownm,Dagger(Cdownm))+&
         matmul(Dagger(Cdownm),Cdownm),identity12,&
         'mobile Cdown anticommutator')
    call assert_matrix_close(matmul(Cupm,Cdownm)+&
         matmul(Cdownm,Cupm),zero12,&
         'mobile cross-spin C anticommutator')
    !
    ! Check the projected impurity algebra.  The missing state in each
    ! anticommutator is the forbidden opposite-spin double occupation.
    call assert_matrix_close(matmul(Dup,Dagger(Dup))+&
         matmul(Dagger(Dup),Dup),projected_up,&
         'mobile Dup projected anticommutator')
    call assert_matrix_close(matmul(Ddown,Dagger(Ddown))+&
         matmul(Dagger(Ddown),Ddown),projected_down,&
         'mobile Ddown projected anticommutator')
    call assert_matrix_close(matmul(Dup,Ddown)+&
         matmul(Ddown,Dup),zero12,&
         'mobile cross-spin D anticommutator')
    !
    ! Check that conduction and impurity fermions anticommute locally.
    call assert_matrix_close(matmul(Cupm,Dup)+matmul(Dup,Cupm),&
         zero12,'mobile Cup Dup anticommutator')
    call assert_matrix_close(matmul(Cupm,Ddown)+matmul(Ddown,Cupm),&
         zero12,'mobile Cup Ddown anticommutator')
    call assert_matrix_close(matmul(Cdownm,Dup)+matmul(Dup,Cdownm),&
         zero12,'mobile Cdown Dup anticommutator')
    call assert_matrix_close(matmul(Cdownm,Ddown)+&
         matmul(Ddown,Cdownm),zero12,&
         'mobile Cdown Ddown anticommutator')
    !
    ! Check total parity for both fermionic species.
    call assert_matrix_close(matmul(Pm,Pm),identity12,&
         'mobile parity square')
    call assert_matrix_close(matmul(Pm,Cupm)+matmul(Cupm,Pm),&
         zero12,'mobile parity Cup')
    call assert_matrix_close(matmul(Pm,Cdownm)+matmul(Cdownm,Pm),&
         zero12,'mobile parity Cdown')
    call assert_matrix_close(matmul(Pm,Dup)+matmul(Dup,Pm),&
         zero12,'mobile parity Dup')
    call assert_matrix_close(matmul(Pm,Ddown)+matmul(Ddown,Pm),&
         zero12,'mobile parity Ddown')
    !
    ! Check the spin algebras and their mutual commutation.
    call assert_matrix_close(matmul(Sepm,Semm)-matmul(Semm,Sepm),&
         2d0*Sezm,'mobile electron spin algebra')
    call assert_matrix_close(matmul(Sipm,Simm)-matmul(Simm,Sipm),&
         2d0*Sizm,'mobile impurity spin algebra')
    call assert_matrix_close(matmul(Sepm,Sipm)-matmul(Sipm,Sepm),&
         zero12,'mobile electron/impurity spin commutation')
    !
    ! Check the impurity-number spectrum in the selected basis ordering.
    do istate=1,4
       call assert_true(abs(Nimp(istate,istate))<tolerance,&
            'mobile empty Nimp eigenvalue')
    enddo
    do istate=5,12
       call assert_true(abs(Nimp(istate,istate)-one)<tolerance,&
            'mobile occupied Nimp eigenvalue')
    enddo
    call assert_matrix_close(Nimp(1:4,5:12),zero12(1:4,5:12),&
         'mobile Nimp off-diagonal block')
    !
    ! Check QN shifts from every nonzero matrix element.
    call assert_operator_dq(Cupm,mobile_qn,[-1,0,0],&
         'mobile Cup dq')
    call assert_operator_dq(Cdownm,mobile_qn,[0,-1,0],&
         'mobile Cdown dq')
    call assert_operator_dq(Dup,mobile_qn,[-1,0,-1],&
         'mobile Dup dq')
    call assert_operator_dq(Ddown,mobile_qn,[0,-1,-1],&
         'mobile Ddown dq')
    call assert_operator_dq(Sipm,mobile_qn,[1,-1,0],&
         'mobile Simp+ dq')
    call assert_operator_dq(Sepm,mobile_qn,[1,-1,0],&
         'mobile Sele+ dq')
    call assert_operator_dq(Nimp,mobile_qn,[0,0,0],&
         'mobile Nimp dq')
    call assert_operator_dq(H0,mobile_qn,[0,0,0],&
         'mobile Hamiltonian dq')
    !
  end subroutine test_mobile_algebra



  subroutine test_mobile_fixed_embedding()
    integer,dimension(:),allocatable      :: fixed_basis,mobile_basis
    integer,dimension(2,8)                :: fixed_qn
    integer,dimension(3,12)               :: mobile_qn
    integer                               :: istate
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hfixed,Hmobile
    complex(8),dimension(:,:),allocatable :: Cfixed,Cmobile
    complex(8),dimension(:,:),allocatable :: Sfixed,Smobile
    complex(8),dimension(:,:),allocatable :: Nfixed,Nmobile
#else
    real(8),dimension(:,:),allocatable    :: Hfixed,Hmobile
    real(8),dimension(:,:),allocatable    :: Cfixed,Cmobile
    real(8),dimension(:,:),allocatable    :: Sfixed,Smobile
    real(8),dimension(:,:),allocatable    :: Nfixed,Nmobile
#endif
    !
    ! The occupied mobile sector must reproduce the fixed-spin model.
    Hfixed =Build_KondoFixed_Hamiltonian(0.73d0,1.17d0,0.23d0)
    Hmobile=Build_KondoMobile_Hamiltonian(0.73d0,1.17d0,0.23d0)
    call assert_matrix_close(Hmobile(5:12,5:12),Hfixed,&
         'mobile/fixed Hamiltonian embedding')
    !
    Cfixed =Build_KondoFixed_C(1)
    Cmobile=Build_KondoMobile_C(1)
    call assert_matrix_close(Cmobile(5:12,5:12),Cfixed,&
         'mobile/fixed Cup embedding')
    Cfixed =Build_KondoFixed_C(2)
    Cmobile=Build_KondoMobile_C(2)
    call assert_matrix_close(Cmobile(5:12,5:12),Cfixed,&
         'mobile/fixed Cdown embedding')
    !
    Sfixed =Build_KondoFixed_ElectronSpin('z')
    Smobile=Build_KondoMobile_ElectronSpin('z')
    call assert_matrix_close(Smobile(5:12,5:12),Sfixed,&
         'mobile/fixed electron Sz embedding')
    Sfixed =Build_KondoFixed_ImpuritySpin('+')
    Smobile=Build_KondoMobile_ImpuritySpin('+')
    call assert_matrix_close(Smobile(5:12,5:12),Sfixed,&
         'mobile/fixed impurity S+ embedding')
    !
    Nfixed =Build_KondoFixed_Number()
    Nmobile=Build_KondoMobile_Number()
    call assert_matrix_close(Nmobile(5:12,5:12),Nfixed,&
         'mobile/fixed Nimp embedding')
    !
    ! The first two mobile charges reduce to the fixed charges, while the
    ! third charge identifies the occupied impurity sector.
    fixed_basis =Build_KondoFixed_BasisQN()
    mobile_basis=Build_KondoMobile_BasisQN()
    fixed_qn =reshape(fixed_basis,[2,8])
    mobile_qn=reshape(mobile_basis,[3,12])
    do istate=1,8
       call assert_true(all(mobile_qn(1:2,istate+4)==&
            fixed_qn(:,istate)),&
            'mobile/fixed QN embedding')
       call assert_true(mobile_qn(3,istate+4)==1,&
            'mobile occupied Nimp charge')
    enddo
    !
  end subroutine test_mobile_fixed_embedding


  subroutine assert_true(condition,message)
    logical,intent(in)          :: condition
    character(len=*),intent(in) :: message
    !
    if(.not.condition)error stop 'testKONDO_LOCAL: '//trim(message)
    !
  end subroutine assert_true


  subroutine assert_state_close(actual,expected,message)
#ifdef _CMPLX
    complex(8),dimension(:),intent(in) :: actual,expected
#else
    real(8),dimension(:),intent(in)    :: actual,expected
#endif
    character(len=*),intent(in)        :: message
    !
    call assert_true(size(actual)==size(expected),trim(message)//' size')
    call assert_true(maxval(abs(actual-expected))<tolerance,&
         trim(message)//' value')
    !
  end subroutine assert_state_close


  subroutine assert_matrix_close(actual,expected,message)
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in) :: actual,expected
#else
    real(8),dimension(:,:),intent(in)    :: actual,expected
#endif
    character(len=*),intent(in)          :: message
    !
    call assert_true(all(shape(actual)==shape(expected)),&
         trim(message)//' shape')
    call assert_true(maxval(abs(actual-expected))<tolerance,&
         trim(message)//' value')
    !
  end subroutine assert_matrix_close


  subroutine assert_operator_dq(operator,charges,expected,message)
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in) :: operator
#else
    real(8),dimension(:,:),intent(in)    :: operator
#endif
    integer,dimension(:,:),intent(in)    :: charges
    integer,dimension(:),intent(in)      :: expected
    character(len=*),intent(in)          :: message
    integer                              :: irow,icol
    !
    do icol=1,size(operator,2)
       do irow=1,size(operator,1)
          if(abs(operator(irow,icol))<tolerance)cycle
          call assert_true(all(charges(:,irow)-charges(:,icol)==expected),&
               message)
       enddo
    enddo
    !
  end subroutine assert_operator_dq


  function Dagger(operator) result(adjoint)
#ifdef _CMPLX
    complex(8),dimension(:,:),intent(in)             :: operator
    complex(8),dimension(size(operator,2),size(operator,1)) :: adjoint
#else
    real(8),dimension(:,:),intent(in)                 :: operator
    real(8),dimension(size(operator,2),size(operator,1))    :: adjoint
#endif
    !
#ifdef _CMPLX
    adjoint=conjg(transpose(operator))
#else
    adjoint=transpose(operator)
#endif
    !
  end function Dagger


end program testKONDO_LOCAL
#endif
