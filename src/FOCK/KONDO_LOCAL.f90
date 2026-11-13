MODULE KONDO_LOCAL
  USE SCIFOR, only: zero,one,eye,kron,diag,to_lower
  implicit none
  private

  integer,parameter :: KondoFixedDim=8

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
  print *,'KONDO_LOCAL TEST: PASS'

contains


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
