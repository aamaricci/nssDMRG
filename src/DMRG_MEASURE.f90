module DMRG_MEASURE
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK
  implicit none
  private


  !Measuring:
  public :: Init_Measure_DMRG
  public :: End_measure_DMRG
  public :: Measure_DMRG
  public :: Measure_Op_DMRG
  public :: Build_Op_DMRG
  public :: Advance_Op_DMRG
  public :: Advance_Corr_DMRG
  public :: Average_Op_DMRG
  public :: Measure_Corr_DMRG
  public :: Measure_Product_DMRG
  !Predefined procedures:
  public :: Measure_SpinSpin_DMRG          !get the spin-spin correlation ("s")
  public :: Measure_DensityDensity_DMRG    !get the density-density correlation ("f")
  public :: Measure_FermionBond_DMRG       !get the fermion bond energy E_ij=\sum_ab <t_ijab c_ai.c_bj>+h.c.
  !public :: Measure_KineticEnergy_DMRG     !get the kinetic energy sum_ij E_ij
  public :: Measure_SpinBond_DMRG          !get the spin bond energy E_ij= H.<S_i.S_j>
  public :: Measure_SpinExchangeEnergy_DMRG!get the spin-exchange energy sum_ij E_ij
  public :: Measure_LocalEnergy_DMRG       !get the local energy <H_i> (contains interaction and local terms: crystal field, external fields, etc.)
  public :: Measure_Energy_DMRG            !a convenience wrapper returning Etotal and partial Ebond+Eloc
  public :: Measure_Polarization_DMRG
  public :: Measure_Structure_DMRG
  public :: Measure_Structure_Aq_DMRG
  public :: Measure_String_DMRG
  public :: Measure_Parity_DMRG
  !Helper procedures:    
  public :: Structure_Factor_DMRG
  public :: Correlation_Distance_DMRG
  ! 
  public :: Write_DMRG

  interface Measure_Corr_DMRG
     module procedure :: Measure_Corr_keys_DMRG
     module procedure :: Measure_Corr_ops_DMRG
  end interface Measure_Corr_DMRG

  interface Measure_Product_DMRG
     module procedure :: Measure_Product_keys_DMRG
     module procedure :: Measure_Product_ops_DMRG
  end interface Measure_Product_DMRG

  interface Measure_Polarization_DMRG
     module procedure :: Measure_Polarization_key_DMRG
     module procedure :: Measure_Polarization_op_DMRG
  end interface

  interface Measure_Structure_DMRG
     module procedure :: Measure_Structure_keys_DMRG
     module procedure :: Measure_Structure_ops_DMRG
  end interface

  interface Measure_Structure_Aq_DMRG
     module procedure :: Measure_Structure_Aq_key_DMRG
     module procedure :: Measure_Structure_Aq_op_DMRG
  end interface

  interface Correlation_Distance_DMRG
     module procedure :: Correlation_Distance_DMRG_d
     module procedure :: Correlation_Distance_DMRG_c
  end interface

  interface Structure_Factor_DMRG
     module procedure :: Structure_Factor_DMRG_d 
     module procedure :: Structure_Factor_DMRG_c
  end interface

  interface Write_DMRG
     module procedure :: write_user_scalar
     module procedure :: write_user_array
     module procedure :: write_user_matrix
  end interface Write_DMRG


  interface  Measure_DMRG
     module procedure :: Measure_DMRG_scalar
     module procedure :: Measure_DMRG_vector
  end interface Measure_DMRG
!
  integer                                      :: Nsb,isb,vecDim
  real(8),dimension(:),allocatable             :: qn,qm
  real(8),dimension(:),allocatable             :: dq
  type(sparse_matrix),allocatable,dimension(:) :: Olist
  type(tstates),dimension(:),allocatable       :: Li,Ri
  type(tstates),dimension(:),allocatable       :: Lmap,Rmap
  logical                                      :: measure_status=.false.
  character(:),allocatable                     :: string
!
!
!Measure State: a meta-object containing all the SB info needed for a measurement
! * sb_states
! * gs_vector
! * sb_sector
! * quantum numbers
!
contains



  !Measure "polarization operator" Z. 
  !see: R. Resta and S. Sorella, `Electron Localization in the Insulating State`, Phys. Rev. Lett. 82, 370 (1999).    
  ! and also in https://arxiv.org/abs/2412.05975:
  !
  ! Z=< product_j exp[i 2*pi*x_j*charge_j/length] >.
  !
  ! Charge must be real diagonal in the local basis and parity even.
  ! Computes the polarization marker using an operator saved under key on each site.
  ! Uses all sites unless positions is supplied. By default, the coordinates are
  ! the physical site numbers and length is the full chain length.
  ! Checks that the operator conserves the quantum numbers and is parity even.
  ! Reuses an open measurement session, or opens and closes its own session.
  function Measure_Polarization_key_DMRG(key,positions,x,length) result(z)
    character(len=*),intent(in)     :: key
    integer,optional,intent(in)     :: positions(:)
    real(8),optional,intent(in)     :: x(:),length
    complex(8)                      :: z
    type(sparse_matrix),allocatable :: Ops(:)
    integer,allocatable             :: p(:)
    integer                         :: i
    logical                         :: owned
    !
    owned=.not.measure_status
    !
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)error stop 'DMRG measurement state unavailable'
    !
    call measurement_positions(p,positions)
    !
    allocate(Ops(size(p)))
    do i=1,size(p)
       if(.NOT.dot(p(i))%operators%has_key(key))error stop 'Measure_Polarization_DMRG: missing charge key'
       if(any(abs(dot(p(i))%operators%dq(key=key))>100d0*epsilon(1d0)))&
            error stop 'Measure_Polarization_DMRG: charge must conserve QNs'
       if(is_odd_fermion_type(dot(p(i))%operators%type(key=key)))&
            error stop 'Measure_Polarization_DMRG: charge must be parity even'
       Ops(i)=dot(p(i))%operators%op(key=key)
    enddo
    !
    z=polarization_product(Ops,p,x,length)
    !
    do i=1,size(p)
       call Ops(i)%free()
    enddo
    if(owned)call End_Measure_DMRG()
  end function Measure_Polarization_key_DMRG

  ! Computes the same polarization marker, but takes the local operator directly.
  ! Applies this operator on each selected site. It must be real and diagonal
  ! in the physical basis, conserve the quantum numbers, and be parity even.
  ! The result keeps both the real and imaginary parts, even for a real state.
  function Measure_Polarization_op_DMRG(Op,positions,x,length) result(z)
    type(sparse_matrix),intent(in)  :: Op
    integer,optional,intent(in)     :: positions(:)
    real(8),optional,intent(in)     :: x(:),length
    complex(8)                      :: z
    type(sparse_matrix),allocatable :: Ops(:)
    integer,allocatable             :: p(:)
    integer                         :: i
    logical                         :: owned
    !
    owned=.not.measure_status
    !
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)error stop 'DMRG measurement state unavailable'
    !
    call measurement_positions(p,positions)
    !
    allocate(Ops(size(p)))
    do i=1,size(p)
       Ops(i)=Op
    enddo
    !
    z=polarization_product(Ops,p,x,length)
    !
    do i=1,size(p)
       call Ops(i)%free()
    enddo
    if(owned)call End_Measure_DMRG()
  end function

  ! String-order parameter: 
  ! see: M. den Nijs and K. Rommelse, `Preroughening transitions
  !      in crystal surfaces and valence-bond phases in quantum spin chains`, Phys. Rev. B 40, 4709 (1989).
  !      and https://arxiv.org/abs/2412.05975. 
  ! Endpoints included as Q, exponential on the strict interior.
  ! Q must be a real diagonal, parity-even, QN-conserving local operator.
  ! Measures a string between sites i and j: Op at the two ends, with
  ! exp(i*theta*Op) on every site strictly between them. theta defaults to pi.
  ! With total on-site Sz as Op, this probes hidden spin correlations.
  ! For adjacent sites there is no interior string, only the two end operators.
  ! Returns the correlation for this interval, without averaging over its origin.
  function Measure_String_DMRG(Op,i,j,theta) result(value)
    type(sparse_matrix),intent(in)  :: Op
    integer,intent(in)              :: i,j
    real(8),optional,intent(in)     :: theta
    complex(8)                      :: value
    type(sparse_matrix),allocatable :: Ops(:)
    integer,allocatable             :: p(:)
    real(8),allocatable             :: x(:)
    real(8)                         :: angle
    integer                         :: k
    logical                         :: owned
    !
    angle=pi;if(present(theta))angle=theta
    owned=.not.measure_status
    !
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)error stop 'DMRG measurement state unavailable'
    !
    if(i<1.or.j>left%length+right%length.or.j<=i)&
         error stop 'Measure_String_DMRG: require 1 <= i < j <= L'
    p=[(k,k=i,j)]
    allocate(Ops(size(p)),x(size(p)))
    !
    x         = angle/pi2
    x(1)      =0d0
    x(size(p))=0d0
    do k=1,size(p)
       Ops(k)=Op
    enddo
    !
    value=polarization_product(Ops,p,x,1d0,endpoints=.true.)
    !
    !
    do k=1,size(p)
       call Ops(k)%free()
    enddo
    if(owned)call End_Measure_DMRG()
  end function Measure_String_DMRG


  ! Parity String-order parameter: 
  ! see: M. den Nijs and K. Rommelse, `Preroughening transitions
  !      in crystal surfaces and valence-bond phases in quantum spin chains`, Phys. Rev. B 40, 4709 (1989).
  !      and https://arxiv.org/abs/2412.05975. 
  ! Interval (i,j], default reference density 2 and angle pi.
  ! Measures exp(i*theta*sum(Op-reference)) from site i+1 through site j.
  ! The interval excludes i and includes j; theta defaults to pi and reference to 2.
  ! For integer density and theta=pi, this measures whether the excess charge
  ! in the interval is even or odd. It does not add operators at the two ends.
  function Measure_Parity_DMRG(Op,i,j,reference,theta) result(value)
    type(sparse_matrix),intent(in)  :: Op
    integer,intent(in)              :: i,j
    real(8),optional,intent(in)     :: reference,theta
    complex(8)                      :: value
    type(sparse_matrix),allocatable :: Ops(:)
    integer,allocatable             :: p(:)
    real(8),allocatable             :: x(:)
    real(8)                         :: angle,nref
    integer                         :: k
    logical                         :: owned
    !
    angle=pi;if(present(theta))angle=theta
    nref=2d0;if(present(reference))nref=reference
    owned=.not.measure_status
    !
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)error stop 'DMRG measurement state unavailable'
    !
    if(i<1.or.j>left%length+right%length.or.j<=i)&
         error stop 'Measure_Parity_DMRG: require 1 <= i < j <= L'
    !
    p=[(k,k=i+1,j)]
    allocate(Ops(size(p)),x(size(p)))
    x=angle/(pi2)
    do k=1,size(p)
       Ops(k)=Op
    enddo
    !
    value=polarization_product(Ops,p,x,1d0)*exp(cmplx(0d0,-angle*nref*(j-i),8))
    !
    do k=1,size(p)
       call Ops(k)%free()
    enddo
    if(owned)call End_Measure_DMRG()
  end function

  ! Builds the local exponential factors and measures their product in the state.
  ! Forms each exponential in the physical basis before any DMRG truncation.
  ! If endpoints is true, uses the bare operators at the first and last selected
  ! sites instead of their exponentials; this also serves the spin string.
  ! In a real build, carries the product as complex matrices to retain its phase.
  function polarization_product(Ops,positions,x,length,endpoints) result(z)
    type(sparse_matrix),intent(in)  :: Ops(:)
    integer,intent(in)              :: positions(:)
    real(8),optional,intent(in)     :: x(:),length
    logical,optional,intent(in)     :: endpoints
    logical                         :: with_endpoints
    complex(8)                      :: z,value
    real(8),allocatable             :: xp(:)
    real(8)                         :: period,theta
    real(8),parameter               :: tol=100d0*epsilon(1d0)
    integer                         :: i,j,k,N,M,d,p
#ifdef _CMPLX
    type(sparse_matrix),allocatable :: factors(:)
    real(8),allocatable             :: dqs(:,:)
    character(len=4),allocatable    :: types(:)
#else
    complex(8),allocatable          :: KL(:,:),KR(:,:),psi(:,:),transformed(:,:)
    real(8),allocatable             :: full(:,:)
    integer                         :: il,ir,mpi_error
#endif
    with_endpoints=.false.;if(present(endpoints))with_endpoints=endpoints
    !
    N =left%length+right%length
    M =size(positions)
    xp=dble(positions)
    if(present(x))then
       if(size(x)/=M)error stop 'Measure_Polarization_DMRG: incompatible coordinates'
       xp=x
    endif
    period=dble(N);if(present(length))period=length
    if(.not.(period>0d0))error stop 'Measure_Polarization_DMRG: length must be positive'
#ifdef _CMPLX
    allocate(factors(M),dqs(size(current_target_qn),M),types(M))
    dqs=0d0
    types='none'
#endif
    do k=1,M
       p=positions(k)
       d=dot(p)%Dim
       if(Ops(k)%Nrow/=d.or.Ops(k)%Ncol/=d)error stop 'Measure_Polarization_DMRG: charge dimension'
#ifdef _CMPLX
       call factors(k)%init(d,d)
#endif
       theta=pi2*xp(k)/period
       do i=1,d
          ! Exponentiate BEFORE truncation, in the physical local basis.
          do j=1,Ops(k)%row(i)%size
             value=Ops(k)%row(i)%vals(j)
             if(Ops(k)%row(i)%cols(j)/=i.and.abs(value)>tol)&
                  error stop 'Measure_Polarization_DMRG: charge must be diagonal'
             if(abs(aimag(value))>tol)error stop 'Measure_Polarization_DMRG: charge must be real'
          enddo
#ifdef _CMPLX
          value=exp(cmplx(0d0,theta*real(Ops(k)%get(i,i),8),8))
          if(with_endpoints.and.(k==1.or.k==M))value=Ops(k)%get(i,i)
          call factors(k)%insert(value,i,i)
#endif
       enddo
    enddo
#ifdef _CMPLX
    z=Measure_Product_DMRG(factors,dqs,types,positions)
    do k=1,M
       call factors(k)%free()
    enddo
#else
    ! A real ground state can have complex z. Keep the propagated product
    ! complex rather than approximating it by cos of the mean position.
    KL=block_product('l')
    KR=block_product('r')
    allocate(full(size(sb_states),1));full=0d0
#ifdef _MPI
    if(MpiStatus)then
       call gather_vector_MPI(MpiComm,gs_vector(:,1:1),full)
       call MPI_Bcast(full,size(full),MPI_DOUBLE_PRECISION,0,MpiComm,mpi_error)
    else
       full=gs_vector(:,1:1)
    endif
#else
    full=gs_vector(:,1:1)
#endif
    allocate(psi(right%Dim,left%Dim));psi=0d0
    do i=1,size(sb_states)
       ir=mod(sb_states(i)-1,right%Dim)+1
       il=(sb_states(i)-1)/right%Dim+1
       psi(ir,il)=full(i,1)
    enddo
    transformed=matmul(matmul(KR,psi),transpose(KL))
    z=sum(conjg(psi)*transformed)
    !
  contains
    !
    ! Builds the whole product on the left or right block, one site at a time.
    ! Uses the saved rotations to move it into the bases kept by DMRG, preserving
    ! the site ordering used when the block was grown, including for PBC.
    function block_product(label) result(Kprod)
      character(len=1),intent(in) :: label
      complex(8),allocatable      :: Kprod(:,:),Xsite(:,:),Udense(:,:)
      type(sparse_matrix)         :: U
      integer                     :: iter,nstep,physical
      nstep=left%length;if(label=='r')nstep=right%length
      physical=physical_position(1,label)
      Kprod=site_phase(physical)
      do iter=1,nstep-1
         if(MpiMaster)then
            if(label=='l')U=left%omatrices%op(key=str(iter))
            if(label=='r')U=right%omatrices%op(key=str(iter))
         endif
#ifdef _MPI
         if(MpiStatus)call U%bcast(comm=MpiComm)
#endif
         Udense  =cmplx(U%as_matrix(),0d0,8)
         Kprod   =matmul(conjg(transpose(Udense)),matmul(Kprod,Udense))
         physical=physical_position(iter+1,label)
         Xsite   =site_phase(physical)
         if(label=='l')then
            if(PBCdmrg.and.mod(iter,2)==0)then
               Kprod=kron(Xsite,Kprod)
            else
               Kprod=kron(Kprod,Xsite)
            endif
         else
            if(PBCdmrg.and.mod(iter,2)==0)then
               Kprod=kron(Kprod,Xsite)
            else
               Kprod=kron(Xsite,Kprod)
            endif
         endif
         call U%free()
      enddo
    end function
    !
    !
    ! Builds the diagonal factor for one physical site. An unselected site gets
    ! the identity; a selected site gets its exponential, or its bare operator
    ! when it is an endpoint of a spin string.
    function site_phase(pos) result(Xsite)
      integer,intent(in)     :: pos
      complex(8),allocatable :: Xsite(:,:)
      integer                :: a,b,dim
      real(8)                :: angle
      dim=dot(pos)%Dim
      allocate(Xsite(dim,dim));Xsite=0d0
      do a=1,dim
         Xsite(a,a)=1d0
      enddo
      do b=1,M
         if(positions(b)/=pos)cycle
         angle=2d0*acos(-1d0)*xp(b)/period
         do a=1,dim
            Xsite(a,a)=exp(cmplx(0d0,angle*Ops(b)%get(a,a),8))
            if(with_endpoints.and.(b==1.or.b==M))Xsite(a,a)=Ops(b)%get(a,a)
         enddo
      enddo
    end function
    !
    !
    ! Finds the physical site number for a given step in the growth of this block.
    ! This lets the product follow the block ordering while using the correct
    ! physical coordinates for the phases.
    function physical_position(index,label) result(pos)
      integer,intent(in)          :: index
      character(len=1),intent(in) :: label
      integer                     :: pos,j
      pos=0
      if(label=='l')then
         do j=1,left%length
            if(b2gMap(j)==index)pos=j
         enddo
      else
         do j=left%length+1,N
            if(b2gMap(N+1-j)==index)pos=j
         enddo
      endif
      if(pos==0)error stop 'Measure_Polarization_DMRG: growth index not mapped'
    end function
#endif
  end function


  ! Measures the structure factor of two named operators. By default, exact
  ! mode measures every ordered pair. The three symmetry flags together allow
  ! it to measure just one triangle. Bulk mode samples central pairs at each
  ! distance instead; in OBC this assumes the selected bulk is homogeneous.
  ! See measure_structure_pairs below for normalization and distance options.
  subroutine Measure_Structure_keys_DMRG(keyA,keyB,q,Fq,positions,connected,Cij,mode,same_operator,hermitian,commuting,rmax,norigins,window_weights,Cr)
    character(len=*),intent(in) :: keyA,keyB
    real(8),intent(in) :: q(:)
    complex(8),allocatable,intent(out) :: Fq(:)
    integer,optional,intent(in) :: positions(:)
    logical,optional,intent(in) :: connected
    complex(8),allocatable,optional,intent(out) :: Cij(:,:)
    character(len=*),optional,intent(in) :: mode
    logical,optional,intent(in) :: same_operator,hermitian,commuting,window_weights
    integer,optional,intent(in) :: rmax,norigins
    complex(8),allocatable,optional,intent(out) :: Cr(:)
    integer,allocatable :: p(:)
    logical :: owned
    owned=.not.measure_status
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)error stop 'DMRG measurement state unavailable'
    call measurement_positions(p,positions)
    if(present(same_operator))then
       if(same_operator.and.keyA/=keyB)error stop 'Measure_Structure_DMRG: different operator keys'
    endif
    call measure_structure_pairs(pair,p,q,Fq,Cij,mode,same_operator,hermitian,commuting,rmax,norigins,window_weights,Cr)
    if(owned)call End_Measure_DMRG()
  contains
    function pair(i,j) result(c)
      integer,intent(in) :: i,j
      complex(8) :: c
      c=Measure_Corr_DMRG(keyA,keyB,i,j,connected)
    end function
  end subroutine Measure_Structure_keys_DMRG

  ! Same measurement with local matrices supplied directly. dqA/dqB describe
  ! their quantum-number changes; typeA/typeB describe their operator types.
  ! The symmetry flags are promises by the caller, not tests of the state.
  subroutine Measure_Structure_ops_DMRG(OpA,dqA,OpB,dqB,q,Fq,positions,typeA,typeB,connected,Cij,mode,same_operator,hermitian,commuting,rmax,norigins,window_weights,Cr)
    type(sparse_matrix),intent(in) :: OpA,OpB
    real(8),intent(in) :: dqA(:),dqB(:),q(:)
    complex(8),allocatable,intent(out) :: Fq(:)
    integer,optional,intent(in) :: positions(:)
    character(len=*),optional,intent(in) :: typeA,typeB
    logical,optional,intent(in) :: connected
    complex(8),allocatable,optional,intent(out) :: Cij(:,:)
    character(len=*),optional,intent(in) :: mode
    logical,optional,intent(in) :: same_operator,hermitian,commuting,window_weights
    integer,optional,intent(in) :: rmax,norigins
    complex(8),allocatable,optional,intent(out) :: Cr(:)
    integer,allocatable :: p(:)
    logical :: owned
    owned=.not.measure_status
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)error stop 'DMRG measurement state unavailable'
    call measurement_positions(p,positions)
    if(present(same_operator))then
       if(same_operator)then
          if(OpA%Nrow/=OpB%Nrow.or.OpA%Ncol/=OpB%Ncol)&
               error stop 'Measure_Structure_DMRG: different operator dimensions'
          if(any(abs(OpA%as_matrix()-OpB%as_matrix())>100d0*epsilon(1d0)))&
               error stop 'Measure_Structure_DMRG: different operators'
          if(size(dqA)/=size(dqB))error stop 'Measure_Structure_DMRG: different QN dimensions'
          if(any(abs(dqA-dqB)>100d0*epsilon(1d0)))&
               error stop 'Measure_Structure_DMRG: different QN shifts'
       endif
    endif
    call measure_structure_pairs(pair,p,q,Fq,Cij,mode,same_operator,hermitian,commuting,rmax,norigins,window_weights,Cr)
    if(owned)call End_Measure_DMRG()
  contains
    function pair(i,j) result(c)
      integer,intent(in) :: i,j
      complex(8) :: c
      c=Measure_Corr_DMRG(OpA,dqA,OpB,dqB,i,j,typeA,typeB,connected)
    end function
  end subroutine Measure_Structure_ops_DMRG

  ! Shared measurement loop. Exact mode sums all pairs and divides by the
  ! number of selected sites; it stores Cij only when requested. If all three
  ! symmetry flags are true, Cij=Cji is real and only i<=j is measured.
  ! Bulk mode measures up to norigins central pairs at each distance, including
  ! r=0. It requires consecutive ascending positions and does not return Cij.
  ! Cr(-rmax:rmax) holds the sampled distance averages. By default Fq is their
  ! unweighted Fourier sum. window_weights=true adds the finite-window factor
  ! 1-|r|/M. Truncating at rmax omits all longer-distance contributions.
  subroutine measure_structure_pairs(pair,p,q,Fq,Cij,mode,same_operator,hermitian,commuting,rmax,norigins,window_weights,Cr)
    interface
       function pair(i,j) result(c)
         integer,intent(in) :: i,j
         complex(8) :: c
       end function
    end interface
    integer,intent(in) :: p(:)
    real(8),intent(in) :: q(:)
    complex(8),allocatable,intent(out) :: Fq(:)
    complex(8),allocatable,optional,intent(out) :: Cij(:,:),Cr(:)
    character(len=*),optional,intent(in) :: mode
    logical,optional,intent(in) :: same_operator,hermitian,commuting,window_weights
    integer,optional,intent(in) :: rmax,norigins
    character(len=:),allocatable :: choice
    complex(8),allocatable :: distance(:)
    complex(8) :: c
    logical :: same,herm,comm,weights,symmetric
    integer :: M,i,j,first,r,Rcut,origins,nsample,available
    real(8) :: weight
    choice='exact';if(present(mode))choice=trim(to_lower(mode))
    same=.false.;if(present(same_operator))same=same_operator
    herm=.false.;if(present(hermitian))herm=hermitian
    comm=.false.;if(present(commuting))comm=commuting
    weights=.false.;if(present(window_weights))weights=window_weights
    symmetric=same.and.herm.and.comm
    M=size(p);allocate(Fq(size(q)));Fq=0d0
    select case(choice)
    case('exact')
       if(present(rmax).or.present(norigins).or.present(window_weights).or.present(Cr))&
            error stop 'Measure_Structure_DMRG: distance options require bulk mode'
       if(present(Cij))allocate(Cij(M,M))
       do i=1,M
          first=1;if(symmetric)first=i
          do j=first,M
             c=pair(p(i),p(j))
             if(symmetric)c=cmplx(real(c,8),0d0,8)
             if(present(Cij))Cij(i,j)=c
             if(symmetric.and.i/=j)then
                Fq=Fq+2d0*cos(q*dble(p(j)-p(i)))*c
                if(present(Cij))Cij(j,i)=c
             else
                Fq=Fq+exp(cmplx(0d0,1d0,8)*q*dble(p(j)-p(i)))*c
             endif
          enddo
       enddo
       Fq=Fq/dble(M)
    case('bulk')
       if(present(Cij))error stop 'Measure_Structure_DMRG: bulk mode returns Cr, not Cij'
       if(any(p(2:)-p(:M-1)/=1))&
            error stop 'Measure_Structure_DMRG: bulk positions must be consecutive and ascending'
       Rcut=min(M-1,M/2);if(present(rmax))Rcut=rmax
       if(Rcut<0.or.Rcut>=M)error stop 'Measure_Structure_DMRG: require 0 <= rmax < window size'
       origins=1;if(present(norigins))origins=norigins
       if(origins<1)error stop 'Measure_Structure_DMRG: norigins must be positive'
       allocate(distance(-Rcut:Rcut));distance=0d0
       do r=0,Rcut
          available=M-r;nsample=min(origins,available)
          first=(available-nsample)/2+1
          do i=first,first+nsample-1
             j=i+r
             distance(r)=distance(r)+pair(p(i),p(j))
             if(r>0.and..not.symmetric)distance(-r)=distance(-r)+pair(p(j),p(i))
          enddo
          distance(r)=distance(r)/dble(nsample)
          if(symmetric)distance(r)=cmplx(real(distance(r),8),0d0,8)
          if(r>0)then
             if(symmetric)then
                distance(-r)=distance(r)
             else
                distance(-r)=distance(-r)/dble(nsample)
             endif
          endif
       enddo
       do r=-Rcut,Rcut
          weight=1d0;if(weights)weight=1d0-dble(abs(r))/dble(M)
          Fq=Fq+weight*exp(cmplx(0d0,1d0,8)*q*dble(r))*distance(r)
       enddo
       if(present(Cr))call move_alloc(distance,Cr)
    case default
       error stop 'Measure_Structure_DMRG: mode must be exact or bulk'
    end select
  end subroutine measure_structure_pairs

  ! Measures <Aq^dagger Aq>/M directly, with Aq=sum_j exp(i*q*j) A_j.
  ! This route needs no translation symmetry and does not measure site pairs.
  ! A_j must be Hermitian, parity even, and conserve all quantum numbers.
  ! The key version can use different local matrices on different sites.
  subroutine Measure_Structure_Aq_key_DMRG(key,q,Fq,positions,connected)
    character(len=*),intent(in) :: key
    real(8),intent(in) :: q(:)
    complex(8),allocatable,intent(out) :: Fq(:)
    integer,optional,intent(in) :: positions(:)
    logical,optional,intent(in) :: connected
    type(sparse_matrix),allocatable :: Ops(:)
    integer,allocatable :: p(:)
    integer :: i
    logical :: owned
    owned=.not.measure_status
    if(owned)call Init_Measure_DMRG()
    if(.not.measure_status)error stop 'Measure_Structure_Aq_DMRG: measurement state unavailable'
    call measurement_positions(p,positions)
    allocate(Ops(size(p)))
    do i=1,size(p)
       if(.not.dot(p(i))%operators%has_key(key))error stop 'Measure_Structure_Aq_DMRG: missing key'
       if(any(abs(dot(p(i))%operators%dq(key=key))>100d0*epsilon(1d0)))&
            error stop 'Measure_Structure_Aq_DMRG: operator must conserve QNs'
       if(is_odd_fermion_type(dot(p(i))%operators%type(key=key)))&
            error stop 'Measure_Structure_Aq_DMRG: operator must be parity even'
       Ops(i)=dot(p(i))%operators%op(key)
    enddo
    call measure_structure_aq_product(Ops,p,q,Fq,connected)
    do i=1,size(p)
       call Ops(i)%free()
    enddo
    if(owned)call End_Measure_DMRG()
  end subroutine

  ! Same direct Fourier measurement with a uniform local matrix supplied by
  ! the caller. dq must be zero; type, if supplied, must be parity even.
  ! positions may be any distinct physical sites, not just a consecutive window.
  subroutine Measure_Structure_Aq_op_DMRG(Op,dq,q,Fq,positions,connected,type)
    type(sparse_matrix),intent(in) :: Op
    real(8),intent(in) :: dq(:),q(:)
    complex(8),allocatable,intent(out) :: Fq(:)
    integer,optional,intent(in) :: positions(:)
    logical,optional,intent(in) :: connected
    character(len=*),optional,intent(in) :: type
    type(sparse_matrix),allocatable :: Ops(:)
    integer,allocatable :: p(:)
    integer :: i
    logical :: owned
    owned=.not.measure_status
    if(owned)call Init_Measure_DMRG()
    if(.not.measure_status)error stop 'Measure_Structure_Aq_DMRG: measurement state unavailable'
    if(size(dq)/=size(current_target_qn))error stop 'Measure_Structure_Aq_DMRG: QN dimension'
    if(any(abs(dq)>100d0*epsilon(1d0)))error stop 'Measure_Structure_Aq_DMRG: operator must conserve QNs'
    if(is_odd_fermion_type(type))error stop 'Measure_Structure_Aq_DMRG: operator must be parity even'
    call measurement_positions(p,positions)
    allocate(Ops(size(p)))
    do i=1,size(p)
       Ops(i)=Op
    enddo
    call measure_structure_aq_product(Ops,p,q,Fq,connected)
    do i=1,size(p)
       call Ops(i)%free()
    enddo
    if(owned)call End_Measure_DMRG()
  end subroutine

  ! Split Aq into cosine and sine sums C+i*S. Since the physical Hermitian
  ! operators commute on distinct sites, Aq^dagger Aq=C*C+S*S.
  ! Build that second moment independently: multiplying the final truncated
  ! C and S matrices would insert unwanted projections and lose contributions.
  ! Connected correlations subtract |<Aq>|^2, even in an inhomogeneous state.
  subroutine measure_structure_aq_product(Ops,p,q,Fq,connected)
    type(sparse_matrix),intent(in) :: Ops(:)
    integer,intent(in) :: p(:)
    real(8),intent(in) :: q(:)
    complex(8),allocatable,intent(out) :: Fq(:)
    logical,optional,intent(in) :: connected
    type(sparse_matrix) :: CL,SL,TL,CR,SR,TR
    real(8),allocatable :: dq0(:)
    real(8) :: value,mean_c,mean_s
    integer :: i,k,N,M
    logical :: subtract_mean
    N=left%length+right%length;M=size(p)
    do i=1,M
       if(Ops(i)%Nrow/=dot(p(i))%Dim.or.Ops(i)%Ncol/=dot(p(i))%Dim)&
            error stop 'Measure_Structure_Aq_DMRG: local dimension mismatch'
       if(any(abs(Ops(i)%as_matrix()-transpose(conjg(cmplx(Ops(i)%as_matrix(),kind=8))))&
            >100d0*epsilon(1d0)))error stop 'Measure_Structure_Aq_DMRG: operator must be Hermitian'
    enddo
    subtract_mean=.false.;if(present(connected))subtract_mean=connected
    allocate(dq0(size(current_target_qn)),Fq(size(q)));dq0=0d0
    do k=1,size(q)
       call build_fourier_moments_block(Ops,p,q(k),'l',CL,SL,TL)
       call build_fourier_moments_block(Ops,p,q(k),'r',CR,SR,TR)
       value=Average_Op_DMRG(TL,1)+Average_Op_DMRG(TR,N)&
            +2d0*real(Average_Corr_LR_DMRG(CL,dq0,CR,dq0),8)&
            +2d0*real(Average_Corr_LR_DMRG(SL,dq0,SR,dq0),8)
       if(subtract_mean)then
          mean_c=Average_Op_DMRG(CL,1)+Average_Op_DMRG(CR,N)
          mean_s=Average_Op_DMRG(SL,1)+Average_Op_DMRG(SR,N)
          value=value-mean_c**2-mean_s**2
       endif
       Fq(k)=cmplx(value/dble(M),0d0,8)
       call CL%free();call SL%free();call TL%free()
       call CR%free();call SR%free();call TR%free()
    enddo
  end subroutine

  ! Reconstruct a block's cosine sum C, sine sum S, and second moment T.
  ! Rotate all three with each saved U, then add the new physical site:
  ! T_new=T_old + 2*cos(q*j)*C_old*A_j + 2*sin(q*j)*S_old*A_j + A_j^2.
  ! Products here act on old-block and new-site spaces before the next rotation.
  ! This preserves the result of the real-space pair construction despite
  ! truncation. A site outside positions contributes zero, not an identity.
  subroutine build_fourier_moments_block(Ops,p,q,side,C,S,T)
    type(sparse_matrix),intent(in) :: Ops(:)
    integer,intent(in) :: p(:)
    real(8),intent(in) :: q
    character(len=1),intent(in) :: side
    type(sparse_matrix),intent(out) :: C,S,T
    type(sparse_matrix) :: U,A,A2,Ib,Is,Cnew,Snew,Tnew
    integer,allocatable :: physical(:),selected(:)
    integer :: nstep,N,j,i,it,pos,d
    real(8) :: wc,ws
    logical :: site_first
    N=left%length+right%length
    nstep=left%length;if(side=='r')nstep=right%length
    allocate(physical(nstep),selected(N));physical=0;selected=0
    do i=1,size(p)
       selected(p(i))=i
    enddo
    if(side=='l')then
       do j=1,left%length
          physical(b2gMap(j))=j
       enddo
    else
       do j=left%length+1,N
          physical(b2gMap(N+1-j))=j
       enddo
    endif
    if(any(physical==0))error stop 'Measure_Structure_Aq_DMRG: incomplete growth map'
    pos=physical(1)
    call local_factor(pos)
    C=wc*A;S=ws*A;T=A2
    do it=1,nstep-1
       if(MpiMaster)then
          if(side=='l')U=left%omatrices%op(key=str(it))
          if(side=='r')U=right%omatrices%op(key=str(it))
       endif
#ifdef _MPI
       if(MpiStatus)call U%bcast(comm=MpiComm)
#endif
       call rotate(C);call rotate(S);call rotate(T)
       pos=physical(it+1)
       call local_factor(pos)
       Ib=Id(C%Nrow);Is=Id(A%Nrow)
       site_first=side=='r'
       if(PBCdmrg.and.mod(it,2)==0)site_first=.not.site_first
       Cnew=tensor(C,Is)+wc*tensor(Ib,A)
       Snew=tensor(S,Is)+ws*tensor(Ib,A)
       Tnew=tensor(T,Is)+2d0*wc*tensor(C,A)+2d0*ws*tensor(S,A)+tensor(Ib,A2)
       C=Cnew;S=Snew;T=Tnew
       call U%free();call Ib%free();call Is%free()
       call Cnew%free();call Snew%free();call Tnew%free()
    enddo
    call A%free();call A2%free()
  contains
    ! Use the physical coordinate in the phase, including in PBC growth order.
    subroutine local_factor(pos)
      integer,intent(in) :: pos
      call A%free();call A2%free()
      d=dot(pos)%Dim;i=selected(pos)
      if(i>0)then
         A=Ops(i)
      else
         call A%init(d,d)
      endif
      A2=matmul(A,A)
      wc=cos(q*dble(pos));ws=sin(q*dble(pos))
    end subroutine
    ! Keep the sum and its second moment as separate projected operators.
    subroutine rotate(X)
      type(sparse_matrix),intent(inout) :: X
#ifdef _MPI
      if(MpiStatus)then
         X=(U%dgr().pm.X).pm.U
      else
         X=matmul(matmul(U%dgr(),X),U)
      endif
#else
      X=matmul(matmul(U%dgr(),X),U)
#endif
    end subroutine
    ! Match the tensor ordering of the saved left/right block growth.
    function tensor(block_op,site_op) result(X)
      type(sparse_matrix),intent(in) :: block_op,site_op
      type(sparse_matrix) :: X
      if(site_first)then
         X=site_op.x.block_op
      else
         X=block_op.x.site_op
      endif
    end function
  end subroutine build_fourier_moments_block

  ! Prepares the list of sites to measure, using the full chain by default.
  ! Rejects an empty list, site numbers outside the chain, and repeated sites.
  ! Keeps the supplied order; it does not sort the sites or make them consecutive.
  subroutine measurement_positions(p,positions)
    integer,allocatable,intent(out) :: p(:)
    integer,optional,intent(in) :: positions(:)
    integer :: i,N
    N=left%length+right%length
    p=[(i,i=1,N)]
    if(present(positions))p=positions
    if(size(p)==0)error stop 'Measure_Structure_DMRG: empty window'
    if(any(p<1).or.any(p>N))error stop 'Measure_Structure_DMRG: invalid position'
    do i=1,size(p)
       if(count(p==p(i))/=1)error stop 'Measure_Structure_DMRG: repeated position'
    enddo
  end subroutine






  ! F_AB(q) = sum_ij exp[i q (x_j-x_i)] C_AB(i,j) / M.
  ! q is in radians (lattice spacing one); arbitrary grids are accepted.
  ! Fourier transforms an existing complex correlation matrix without measuring
  ! any new correlations. Sums every pair with its phase and divides by M,
  ! the number of sites represented by the matrix.
  ! Uses coordinates x if given, otherwise 1,...,M. q is given in radians
  ! per coordinate unit and need not lie on a particular momentum grid.
  function Structure_Factor_DMRG_c(C,q,x) result(Fq)
    complex(8),intent(in)       :: C(:,:)
    real(8),intent(in)          :: q(:)
    real(8),optional,intent(in) :: x(:)
    complex(8),allocatable      :: Fq(:),phase(:)
    real(8),allocatable         :: xp(:)
    integer                     :: M,i,k
    M=size(C,1)
    if(M==0.or.size(C,2)/=M)error stop 'Structure_Factor_DMRG: nonempty square C required'
    xp=[(dble(i),i=1,M)]
    if(present(x))then
       if(size(x)/=M)error stop 'Structure_Factor_DMRG: incompatible coordinates'
       xp=x
    endif
    allocate(Fq(size(q)))
    do k=1,size(q)
       phase=exp(cmplx(0d0,q(k),8)*xp)
       Fq(k)=dot_product(phase,matmul(C,phase))/dble(M)
    enddo
  end function

  ! Accepts a real correlation matrix and passes it to the complex transform.
  ! The result is still complex, since a real matrix need not be symmetric.
  function Structure_Factor_DMRG_d(C,q,x) result(Fq)
    real(8),intent(in)          :: C(:,:),q(:)
    real(8),optional,intent(in) :: x(:)
    complex(8),allocatable      :: Fq(:)
    Fq=Structure_Factor_DMRG_c(cmplx(C,0d0,8),q,x)
  end function






  
  ! Signed separation bins, including r=0. Counts retain the OBC weights.
  ! Groups the correlation matrix entries by the signed distance x_j-x_i.
  ! Returns each distance r, its mean correlation Cr, and its pair count.
  ! Includes zero and negative distances; unused distances have zero count.
  ! The counts are needed when rebuilding the Fourier sum, since different
  ! distances generally have different numbers of pairs in an open chain.
  subroutine Correlation_Distance_DMRG_c(C,positions,r,Cr,counts)
    complex(8),intent(in)              :: C(:,:)
    integer,intent(in)                 :: positions(:)
    integer,allocatable,intent(out)    :: r(:),counts(:)
    complex(8),allocatable,intent(out) :: Cr(:)
    integer                            :: d,i,j,k,M,span
    M=size(positions)
    if(M==0.or.any(shape(C)/=[M,M]))error stop 'Correlation_Distance_DMRG: incompatible input'
    span=maxval(positions)-minval(positions)
    r=[(d,d=-span,span)]
    allocate(Cr(size(r)),counts(size(r)));Cr=0d0;counts=0
    do j=1,M
       do i=1,M
          k=positions(j)-positions(i)+span+1
          Cr(k)=Cr(k)+C(i,j)
          counts(k)=counts(k)+1
       enddo
    enddo
    do k=1,size(r)
       if(counts(k)>0)Cr(k)=Cr(k)/dble(counts(k))
    enddo
  end subroutine

  ! Accepts a real correlation matrix and passes it to the complex distance
  ! averaging routine, returning the same distances, averages, and pair counts.
  subroutine Correlation_Distance_DMRG_d(C,positions,r,Cr,counts)
    real(8),intent(in)                 :: C(:,:)
    integer,intent(in)                 :: positions(:)
    integer,allocatable,intent(out)    :: r(:),counts(:)
    complex(8),allocatable,intent(out) :: Cr(:)
    call Correlation_Distance_DMRG_c(cmplx(C,0d0,8),positions,r,Cr,counts)
  end subroutine







  
  !##################################################################
  !          INIT / END MEASUREMENT: allocate/deallocate
  !##################################################################
  subroutine Init_Measure_dmrg(msg)
    real(8),dimension(:),allocatable :: qn
    character(len=*),optional        :: msg
    integer,dimension(2)             :: omat_dims
    integer                          :: ilat,i,f,m,istate,mpi_error,map_length
    logical                          :: found_measure_state
    logical                          :: need_measure_state
    character(len=:),allocatable     :: default_suffix,current_suffix
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: init measure"
#endif
    !
    !Check if BLOCK/MEASURE states is available.
    need_measure_state = (.not.allocated(sb_states)).or.(.not.allocated(gs_vector)).OR.(size(sb_sector)==0)
    !
    !Load Measure Blocks (if required)
    if(need_measure_state)call sb_load_measure_blocks()
    !
    !Load Umatrices
    if(MpiMaster)then 
      if(size(left%omatrices)<=1)then
         default_suffix = suffix_dmrg('left',type='i')//".restart"
         current_suffix = suffix_dmrg('left')//".restart"
         call left%load_umat(str(default_suffix),left%length)
         if(size(left%omatrices)<=1.and.str(current_suffix)/=str(default_suffix))&
              call left%load_umat(str(current_suffix),left%length)
      endif
      !
      if(size(right%omatrices)<=1)then
         default_suffix = suffix_dmrg('right',type='i')//".restart"
         current_suffix = suffix_dmrg('right')//".restart"
         call right%load_umat(str(default_suffix),right%length)
         if(size(right%omatrices)<=1.and.str(current_suffix)/=str(default_suffix))&
              call right%load_umat(str(current_suffix),right%length)
      endif
    endif
    !
    !A missing state must not allow measurements on incomplete data.
    if(need_measure_state)then
       call sb_load_measure_state(found_measure_state)
       if(.not.found_measure_state)then
          if(MpiMaster)write(LOGfile,*)"Init_Measure_DMRG: no saved SuperBlock measurement state found."
          measure_status=.false.
          error stop 1
       endif
    endif
    !
    !
    if(MpiMaster)omat_dims = [size(left%omatrices),size(right%omatrices)]
#ifdef _MPI
    ! Use the native inout collective before testing the received dimensions.
    if(MpiStatus)call MPI_Bcast(omat_dims,2,MPI_INTEGER,0,MpiComm,mpi_error)
#endif
    if(any(omat_dims==1))then
       measure_status=.false.
       return
    endif
    !
    if(measure_status)call End_Measure_DMRG()
    !
    string="";if(present(msg))string=msg
    !
    if(MpiMaster)call start_timer("Start measuring..."//str(string))
    !
    call sb_build_dims(quiet=.true.)
    !
    Nsb  = size(sb_sector)
    !
    allocate(LI(Nsb))
    allocate(RI(Nsb))
    allocate(Lmap(Nsb))
    allocate(Rmap(Nsb))
    do isb=1,Nsb
       qn  = sb_sector%qn(index=isb)
       LI(isb)%states = sb2block_states(qn,'left')
       RI(isb)%states = sb2block_states(qn,'right')
       !Inverse maps are required to filter rectangular dq-changing
       !operator blocks without repeatedly searching the sector states.
       allocate(Lmap(isb)%states(left%Dim));Lmap(isb)%states=0
       allocate(Rmap(isb)%states(right%Dim));Rmap(isb)%states=0
       do istate=1,size(LI(isb)%states)
          Lmap(isb)%states(LI(isb)%states(istate))=istate
       enddo
       do istate=1,size(RI(isb)%states)
          Rmap(isb)%states(RI(isb)%states(istate))=istate
       enddo
    enddo
    suffix=label_DMRG('u')
    !
    measure_status=.true.
    !
    !    
    !add setup the map from local to global index here:
    if(allocated(b2gMap))deallocate(b2gMap)
    ! Final enlarged blocks can exceed the input Ldmrg after finite sweeps.
    ! Map their actual lengths, including when loading a measurement checkpoint.
    map_length=max(left%length,right%length)
    allocate(b2gMap(map_length))
    if(PBCdmrg)then
       !The symmetric PBC blocks have odd length, map_length = 2*m+1
       f = (map_length+1)/2          !mid-point
       m = (map_length-1)/2
       !
       b2gMap(f) = 1
       do i=1,m
          b2gMap(f+i)  = 2*i
          b2gMap(f-i)  = 2*i+1
       enddo
    else
       b2gMap = (/(i,i=1,map_length)/)
    endif
  end subroutine Init_Measure_dmrg





  subroutine End_measure_DMRG()
    type(sparse_matrix) :: Ileft,Iright
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: end measure"
#endif
    !Measure_DMRG finalizes its own measurement session.  Make an
    !additional explicit call harmless instead of accessing state that
    !has already been released (in particular the allocatable string).
    if(.not.measure_status)return
    if(MpiMaster)call stop_timer("Done "//str(string))
    if(allocated(string))deallocate(string)
    if(allocated(Olist))deallocate(Olist)
    if(allocated(Li))deallocate(Li)
    if(allocated(Ri))deallocate(Ri)
    if(allocated(Lmap))deallocate(Lmap)
    if(allocated(Rmap))deallocate(Rmap)
    call sb_delete_dims()
    if(.not.block_umat_cache)then   
      if(MpiMaster)then
        !U(1) is the identity in the initial one-site block basis.  Do
        !not use left/right%Dim here: at the end of a run those are the
        !dimensions of the final enlarged blocks.  A subsequent call to
        !Init_Measure_DMRG reloads U(2),...,U(L-1) and relies on U(1) to
        !reconstruct operators on the first two growth sites.
        Ileft  = id(init_left%Dim)
        Iright = id(init_right%Dim)
        call left%omatrices%free()
        call right%omatrices%free()
        call left%put_omat("1",Ileft)
        call right%put_omat("1",Iright)
        call Ileft%free()
        call Iright%free()
      endif
    endif
    measure_status=.false.
  end subroutine End_measure_DMRG


  



  subroutine Error_measure_DMRG
    if(MpiMaster)then
       write(LOGfile,*)"Init_Measure_DMRG: either Block=L,R block have size(Block.omatrices)==1."
       write(LOGfile,*)"Init_Measure_DMRG: No rotation matrices are present.                    "
       write(LOGfile,*)"Init_Measure_DMRG: No measurements are possible.                        "
    endif
  end subroutine Error_measure_DMRG








  !##################################################################
  !          Measure local operators on a given set of
  !     positions, write results on file and return average value
  !##################################################################
  subroutine Measure_DMRG_scalar(Op,pos,file,avOp)
    type(sparse_matrix),intent(in)            :: Op
    integer,dimension(:),optional             :: pos
    character(len=*),optional                 :: file
    real(8),dimension(:),allocatable,optional :: avOp
    real(8),dimension(:),allocatable          :: vals
    integer,dimension(:),allocatable          :: pos_
    character(:),allocatable                  :: msg
    character(len=1)                          :: label
    character(len=128)                        :: file_
    integer                                   :: i,ipos,L,R,Np
    type(sparse_matrix)                       :: Oi
    integer                                   :: it,j,dims(2)
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: measure scalar"
#endif
    !
    msg="";if(present(file))msg=file
    !
    call Init_measure_dmrg(msg)
    if(.not.measure_status)then
       call Error_measure_DMRG
       return
    endif
    !
    L = left%length ; R = right%length
    Np = L+R;if(present(pos))Np= size(pos)
    !
    allocate(pos_(Np))
    pos_=arange(1,Np);if(present(pos))pos_=pos
    allocate(vals(Np))
    !
    do i=1,Np
       ipos    = pos_(i)
       vals(i) = Measure_Op_DMRG(Op,ipos)
#ifdef _DEBUG
       if(MpiMaster)write(LOGfile,*)ipos,vals(i)
#endif
       if(MpiMaster)call eta(i,Np)
    enddo
    if(MpiMaster.AND.present(file))call Write_DMRG(trim(file),vals,pos_)
    call End_measure_dmrg()
    !
    if(present(avOp))then
       if(allocated(avOp))deallocate(avOp)
       allocate(avOp, source=vals)
    endif
  end subroutine Measure_dmrg_scalar



  subroutine Measure_DMRG_vector(Op,pos,file,avOp)
    type(sparse_matrix),dimension(:),intent(in) :: Op
    integer,dimension(:),optional               :: pos
    character(len=*),optional                   :: file
    real(8),dimension(:,:),allocatable,optional :: avOp
    real(8),dimension(:,:),allocatable          :: vals
    integer,dimension(:),allocatable            :: pos_
    character(:),allocatable                    :: msg
    character(len=1)                            :: label
    character(len=128)                          :: file_
    integer                                     :: i,ipos,L,R,Np,M
    type(sparse_matrix)                         :: Oi
    integer                                     :: it,j,dims(2)
    !
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: measure vector"
#endif
    !
    msg="";if(present(file))msg=file
    !
    call Init_measure_dmrg(msg)
    if(.not.measure_status)then
       call Error_measure_DMRG
       return
    endif
    !
    M  = size(Op)
    L  = left%length
    R  = right%length
    Np = L+R;if(present(pos))Np= size(pos)
    !
    allocate(pos_(Np))
    pos_=arange(1,Np);if(present(pos))pos_=pos
    !
    allocate(vals(M,Np))
    !
    do i=1,Np
       ipos = pos_(i)
       do j=1,M
          vals(j,i) = Measure_Op_DMRG(Op(j),ipos)
       enddo
#ifdef _DEBUG
       if(MpiMaster)write(LOGfile,*)ipos,(vals(j,i),j=1,M)
#endif
       if(MpiMaster)call eta(i,Np)
    enddo
    if(MpiMaster.AND.present(file))call Write_DMRG(trim(file),vals,pos_)
    call End_measure_dmrg()
    !
    if(present(avOp))then
       if(allocated(avOp))deallocate(avOp)
       allocate(avOp, source=vals)
    endif
  end subroutine Measure_DMRG_vector






  
  !##################################################################
  !              Measure local Operator Op
  !Purpose: return the average value <gs|Op|gs> for a given Op
  !##################################################################
  function Measure_Op_DMRG(Op,pos) result(avOp)
    type(sparse_matrix),intent(in) :: Op
    integer                        :: pos
    type(sparse_matrix)            :: Oi
    real(8)                        :: avOp
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: measure Op",pos
#endif
    !
    avOp = zero
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    Oi   = Build_Op_dmrg(Op,pos)
    Oi   = Advance_Op_dmrg(Oi,pos)
    avOp = Average_Op_dmrg(Oi,pos)
    call Oi%free()
  end function Measure_Op_DMRG





  !##################################################################
  !              MEASURE A GENERIC STATIC CORRELATION
  !##################################################################
  ! Measure a static two-point function from two operator keys:
  ! $ C_{AB}(i,j)=\langle\Psi|O_A(i)O_B(j)|\Psi\rangle $.
  ! If connected=.true. return instead
  ! $ C^c_{AB}(i,j)=C_{AB}(i,j)
  !                    -\langle O_A(i)\rangle\langle O_B(j)\rangle$.
  ! Matrix, quantum-number shift and grading are read together from
  ! LIST_OPERATORS; this is therefore the preferred public interface.
  !
  ! Example:
  ! corr = Measure_Corr_DMRG(keyA,keyB,i,j,connected=.true.)
  !
  function Measure_Corr_keys_DMRG(keyA,keyB,posA,posB,connected) result(corr)
    character(len=*),intent(in)       :: keyA,keyB
    integer,intent(in)                :: posA,posB
    logical,optional,intent(in)        :: connected
#ifdef _CMPLX
    complex(8)                        :: corr
#else
    real(8)                           :: corr
#endif
    type(sparse_matrix)               :: OpA,OpB
    real(8),dimension(:),allocatable  :: dqA,dqB
    character(len=:),allocatable      :: typeA,typeB
    integer                           :: N
    !
    corr=zero
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    !
    N=left%length+right%length
    if(posA<1.OR.posA>N)stop "Measure_Corr_DMRG ERROR: posA not in [1,Ldmrg]"
    if(posB<1.OR.posB>N)stop "Measure_Corr_DMRG ERROR: posB not in [1,Ldmrg]"
    if(.not.dot(posA)%operators%has_key(keyA))&
         stop "Measure_Corr_DMRG ERROR: keyA missing in site operator list"
    if(.not.dot(posB)%operators%has_key(keyB))&
         stop "Measure_Corr_DMRG ERROR: keyB missing in site operator list"
    !
    OpA   = dot(posA)%operators%op(key=keyA)
    OpB   = dot(posB)%operators%op(key=keyB)
    dqA   = dot(posA)%operators%dq(key=keyA)
    dqB   = dot(posB)%operators%dq(key=keyB)
    typeA = dot(posA)%operators%type(key=keyA)
    typeB = dot(posB)%operators%type(key=keyB)
    corr  = Measure_Corr_ops_DMRG(OpA,dqA,OpB,dqB,posA,posB,typeA,typeB,connected)
    !
    call OpA%free()
    call OpB%free()
  end function Measure_Corr_keys_DMRG




  ! Low-level interface for composite or conjugated operators.
  ! The caller supplies both shifts, defined by
  ! $ O|q\rangle\longmapsto|q+dq(O)\rangle$, 
  ! and the optional types used to determine the fermionic grading.
  !
  ! Example:
  !
  ! corr = Measure_Corr_DMRG(OpA,dqA,OpB,dqB,i,j,typeA,typeB,connected=.true.)
  !
  function Measure_Corr_ops_DMRG(OpA,dqA,OpB,dqB,posA,posB,typeA,typeB,connected) result(corr)
    type(sparse_matrix),intent(in) :: OpA,OpB
    real(8),dimension(:),intent(in):: dqA,dqB
    integer,intent(in)             :: posA,posB
    character(len=*),optional      :: typeA,typeB
    logical,optional,intent(in)    :: connected
#ifdef _CMPLX
    complex(8)                     :: corr
#else
    real(8)                        :: corr
#endif
    type(sparse_matrix)            :: Oi,Oj,Oij
    character(len=:),allocatable   :: typeA_,typeB_
    integer                        :: L,N
    logical                        :: oddA,oddB,connected_
    real(8)                        :: avA,avB
    real(8),parameter              :: dq_tol=100d0*epsilon(1d0)
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: measure Corr",posA,posB
#endif
    !
    typeA_="";if(present(typeA))typeA_=to_lower(str(typeA))
    typeB_="";if(present(typeB))typeB_=to_lower(str(typeB))
    !
    corr=zero
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    !
    L=left%length
    N=L+right%length
    if(posA<1.OR.posA>N)stop "Measure_Corr_DMRG ERROR: posA not in [1,Ldmrg]"
    if(posB<1.OR.posB>N)stop "Measure_Corr_DMRG ERROR: posB not in [1,Ldmrg]"
    if(size(dqA)/=size(current_target_qn))&
         stop "Measure_Corr_DMRG ERROR: size(dqA) != QN dimension"
    if(size(dqB)/=size(current_target_qn))&
         stop "Measure_Corr_DMRG ERROR: size(dqB) != QN dimension"
    !
    !Check for fermion parity of the operator. 
    !If odd then it is a fermion-like operator then 
    !Wigner-Jordan strings shoule be included in the measurement.
    oddA=is_odd_fermion_type(typeA_)
    oddB=is_odd_fermion_type(typeB_)
    connected_=.false.;if(present(connected))connected_=connected
    !
    !Check for forbidden expectation values
    !A forbidden expectation value is zero, not an invalid input:
    !  dq(O_A O_B)=dq_A+dq_B=0,
    !  p(O_A O_B)=p_A+p_B=0 mod 2.
    corr=zero
    if(any(abs(dqA+dqB)>dq_tol).OR.(oddA.neqv.oddB))return
    !
    if(posA==posB)then
       !At equal positions preserve the requested product order O_A.O_B.
       Oij  = matmul(OpA,OpB)
       Oi   = Build_Op_DMRG(Oij,posA)
       Oj   = Advance_Op_DMRG(Oi,posA)
       corr = Average_Op_DMRG(Oj,posA)
    elseif((posA<=L.AND.posB<=L).OR.(posA>L.AND.posB>L))then
       if(oddA)then
          Oij=Build_Fermion_Corr_Block_DMRG(OpA,OpB,posA,posB)
       else
          Oij=Build_Corr_Block_DMRG(OpA,OpB,posA,posB)
       endif
       corr = Average_Op_DMRG(Oij,posA)
    else
       !Put the operator belonging to the left block first. For two odd
       !operators this canonical reordering contributes one minus sign if
       !the requested product was originally right-operator times left.
       if(posA<=L)then
          if(oddA)then
             Oi=Build_Fermion_LR_End_DMRG(OpA,posA,'l')
             Oj=Build_Fermion_LR_End_DMRG(OpB,posB,'r')
          else
             Oi=Build_Op_DMRG(OpA,posA)
             Oi=Advance_Op_DMRG(Oi,posA)
             Oj=Build_Op_DMRG(OpB,posB)
             Oj=Advance_Op_DMRG(Oj,posB)
          endif
          corr=Average_Corr_LR_DMRG(Oi,dqA,Oj,dqB)
       else
          if(oddA)then
             Oi=Build_Fermion_LR_End_DMRG(OpB,posB,'l')
             Oj=Build_Fermion_LR_End_DMRG(OpA,posA,'r')
          else
             Oi=Build_Op_DMRG(OpB,posB)
             Oi=Advance_Op_DMRG(Oi,posB)
             Oj=Build_Op_DMRG(OpA,posA)
             Oj=Advance_Op_DMRG(Oj,posA)
          endif
          corr=Average_Corr_LR_DMRG(Oi,dqB,Oj,dqA)
          if(oddA)corr=-corr
       endif
    endif
    !The disconnected term can be non-zero only for parity-even,
    !dq=0 operators. Avoid measurements known to vanish by symmetry.
    if(connected_)then
       avA=zero;avB=zero
       if(all(abs(dqA)<=dq_tol).AND..not.oddA)avA=Measure_Op_DMRG(OpA,posA)
       if(all(abs(dqB)<=dq_tol).AND..not.oddB)avB=Measure_Op_DMRG(OpB,posB)
       corr=corr-avA*avB
    endif
    !
    call Oi%free()
    call Oj%free()
    call Oij%free()
  end function Measure_Corr_ops_DMRG




  !##################################################################
  !         MEASURE AN ORDERED PRODUCT OF LOCAL OPERATORS
  !##################################################################
  !> Read every factor from its site's LIST_OPERATORS entry.  The order
  !> of keys and positions is the order of the operator product; the
  !> same position may appear more than once.
  function Measure_Product_keys_DMRG(keys,positions) result(value)
    character(len=*),intent(in)          :: keys(:)
    integer,intent(in)                   :: positions(:)
#ifdef _CMPLX
    complex(8)                            :: value
#else
    real(8)                               :: value
#endif
    type(sparse_matrix),allocatable       :: Ops(:)
    real(8),allocatable                   :: dqs(:,:)
    character(len=:),allocatable          :: types(:)
    integer                               :: a,M,N,Qdim
    !
    value=zero
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    M=size(keys)
    N=left%length+right%length
    if(M==0.OR.size(positions)/=M)&
         stop "Measure_Product_DMRG ERROR: incompatible factors and positions"
    if(any(positions<1).OR.any(positions>N))&
         stop "Measure_Product_DMRG ERROR: position not in [1,Nsites]"
    !
    Qdim=size(current_target_qn)
    allocate(Ops(M),dqs(Qdim,M))
    allocate(character(len=64)::types(M))
    do a=1,M
       if(.not.dot(positions(a))%operators%has_key(trim(keys(a))))&
            stop "Measure_Product_DMRG ERROR: missing site operator key"
       Ops(a)=dot(positions(a))%operators%op(trim(keys(a)))
       dqs(:,a)=dot(positions(a))%operators%dq(trim(keys(a)))
       types(a)=dot(positions(a))%operators%type(key=trim(keys(a)))
    enddo
    value=Measure_Product_ops_DMRG(Ops,dqs,types,positions)
    do a=1,M
       call Ops(a)%free()
    enddo
  end function Measure_Product_keys_DMRG




  !> Evaluate <O_1(p_1)...O_M(p_M)> in exactly the supplied order.
  !> dqs(:,a) and types(a) describe Ops(a); type="fermionic" marks an
  !> odd factor.  This interface also accepts conjugated/composite local
  !> operators that have no key in LIST_OPERATORS.
  function Measure_Product_ops_DMRG(Ops,dqs,types,positions) result(value)
    type(sparse_matrix),intent(in)        :: Ops(:)
    real(8),intent(in)                   :: dqs(:,:)
    character(len=*),intent(in)          :: types(:)
    integer,intent(in)                   :: positions(:)
#ifdef _CMPLX
    complex(8)                            :: value
#else
    real(8)                               :: value
#endif
    type(sparse_matrix),allocatable       :: SiteOps(:),BlockOps(:)
    type(sparse_matrix)                   :: Oleft,Oright,Psite,Tmp
    integer,allocatable                   :: BlockPos(:)
    logical,allocatable                   :: active(:)
    real(8),allocatable                   :: dqLeft(:),dqRight(:)
    integer                               :: a,k,M,N,L,Qdim,Np,ib,odd_count
    real(8),parameter                     :: dq_tol=100d0*epsilon(1d0)
    !
    value=zero
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    M=size(Ops)
    L=left%length
    N=L+right%length
    Qdim=size(current_target_qn)
    if(M==0.OR.size(types)/=M.OR.size(positions)/=M)&
         stop "Measure_Product_DMRG ERROR: incompatible factors and positions"
    if(size(dqs,1)/=Qdim.OR.size(dqs,2)/=M)&
         stop "Measure_Product_DMRG ERROR: shape(dqs) != [Qdim,M]"
    if(any(positions<1).OR.any(positions>N))&
         stop "Measure_Product_DMRG ERROR: position not in [1,Nsites]"
    !
    !A fixed-sector expectation value vanishes unless the total shift
    !and fermionic grading are both zero.  Keep the two block shifts for
    !the rectangular L/R sector contraction below.
    allocate(dqLeft(Qdim),dqRight(Qdim))
    dqLeft=0d0;dqRight=0d0;odd_count=0
    do a=1,M
       if(positions(a)<=L)then
          dqLeft=dqLeft+dqs(:,a)
       else
          dqRight=dqRight+dqs(:,a)
       endif
       if(is_odd_fermion_type(types(a)))odd_count=odd_count+1
    enddo
    if(any(abs(dqLeft+dqRight)>dq_tol).OR.mod(odd_count,2)/=0)return
    !
    !Build the ordinary tensor-product representation site by site.
    !For each odd factor at p, its Jordan--Wigner string contributes
    !P_k on every k<p; then append the factor itself at p.  Traversing
    !the factors in input order preserves signs without sorting them.
    !In particular, repeated positions are multiplied locally before
    !any later DMRG truncation.
    allocate(SiteOps(N),active(N));active=.false.
    do a=1,M
       if(is_odd_fermion_type(types(a)))then
          do k=1,positions(a)-1
             Psite=local_parity_operator(k)
             call append_at_site(k,Psite)
             call Psite%free()
          enddo
       endif
       call append_at_site(positions(a),Ops(a))
    enddo
    !
    !Each side is built in its final renormalized basis.  An empty side
    !contributes the identity.  The final contraction handles operators
    !with nonzero but compensating shifts on opposite sides of the cut.
    do ib=1,2
       if(ib==1)then
          Np=count(active(1:L))
       else
          Np=count(active(L+1:N))
       endif
       if(Np==0)then
          if(ib==1)Oleft=id(left%Dim)
          if(ib==2)Oright=id(right%Dim)
          cycle
       endif
       allocate(BlockOps(Np),BlockPos(Np))
       Np=0
       do k=1,N
          if(ib==1.AND.k>L)cycle
          if(ib==2.AND.k<=L)cycle
          if(.not.active(k))cycle
          Np=Np+1
          BlockOps(Np)=SiteOps(k)
          BlockPos(Np)=k
       enddo
       if(ib==1)Oleft=Build_Product_Block_DMRG(BlockOps,BlockPos)
       if(ib==2)Oright=Build_Product_Block_DMRG(BlockOps,BlockPos)
       do k=1,Np
          call BlockOps(k)%free()
       enddo
       deallocate(BlockOps,BlockPos)
    enddo
    value=Average_Corr_LR_DMRG(Oleft,dqLeft,Oright,dqRight)
    !
    do k=1,N
       if(active(k))call SiteOps(k)%free()
    enddo
    call Oleft%free()
    call Oright%free()
  contains
    !Multiply a new local factor on the right of the factors already
    !assigned to this site, preserving the original product order.
    subroutine append_at_site(pos,Op)
      integer,intent(in)             :: pos
      type(sparse_matrix),intent(in) :: Op
      if(active(pos))then
         Tmp=matmul(SiteOps(pos),Op)
         SiteOps(pos)=Tmp
         call Tmp%free()
      else
         SiteOps(pos)=Op
         active(pos)=.true.
      endif
    end subroutine append_at_site
  end function Measure_Product_ops_DMRG




  !##################################################################
  !                 PHYSICAL CORRELATION WRAPPERS
  !##################################################################
  !Return <S_i.S_j>=<Sz_i Sz_j>+1/2 <S+_i S-_j>+1/2 <S-_i S+_j>.
  function Measure_SpinSpin_DMRG(posA,posB) result(corr)
    integer,intent(in)       :: posA,posB
#ifdef _CMPLX
    complex(8)               :: corr
#else
    real(8)                  :: corr
#endif
    type(sparse_matrix)      :: SzA,SzB,SpA,SpB,SmA,SmB
    real(8),allocatable      :: dqzA(:),dqzB(:),dqpA(:),dqpB(:)
    character(:),allocatable :: key
    !
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)then
       corr=zero
       return
    endif
    key="S"//dot(posA)%okey(0,1,ilink="n")
    SzA=dot(posA)%operators%op(key);dqzA=dot(posA)%operators%dq(key)
    key="S"//dot(posB)%okey(0,1,ilink="n")
    SzB=dot(posB)%operators%op(key);dqzB=dot(posB)%operators%dq(key)
    key="S"//dot(posA)%okey(0,2,ilink="n")
    SpA=dot(posA)%operators%op(key);dqpA=dot(posA)%operators%dq(key)
    key="S"//dot(posB)%okey(0,2,ilink="n")
    SpB=dot(posB)%operators%op(key);dqpB=dot(posB)%operators%dq(key)
    SmA=hconjg(SpA)
    SmB=hconjg(SpB)
    !
    corr = Measure_Corr_ops_DMRG(SzA,dqzA,SzB,dqzB,posA,posB,"bosonic","bosonic")
    corr = corr + 0.5d0*Measure_Corr_ops_DMRG(SpA,dqpA,SmB,-dqpB,posA,posB,"bosonic","bosonic")
    corr = corr + 0.5d0*Measure_Corr_ops_DMRG(SmA,-dqpA,SpB,dqpB,posA,posB,"bosonic","bosonic")
    !
    call SzA%free();call SzB%free()
    call SpA%free();call SpB%free()
    call SmA%free();call SmB%free()
  end function Measure_SpinSpin_DMRG


  !> Return all spin-orbital resolved density correlations
  !> \f[ C_{ab}(i,j)=\langle n_{ia}n_{jb}\rangle,
  !>     \qquad n_a=c_a^\dagger c_a, \quad a,b=1,\ldots,N_{so}. \f]
  !> The compound index follows the site convention
  !> \f$a=i_{orb}+(i_{spin}-1)N_{orb}\f$.
  function Measure_DensityDensity_DMRG(posA,posB,connected) result(corr)
    integer,intent(in)       :: posA,posB
    logical,optional,intent(in):: connected
#ifdef _CMPLX
    complex(8)               :: corr(Nspin*Norb,Nspin*Norb)
#else
    real(8)                  :: corr(Nspin*Norb,Nspin*Norb)
#endif
    type(sparse_matrix)      :: C,NopA(Nspin*Norb),NopB(Nspin*Norb)
    real(8),allocatable      :: dq0(:)
    character(:),allocatable :: key
    integer                  :: io,jo,iorb,ispin,Nso
    !
    corr=zero
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    !
    Nso=Nspin*Norb
    allocate(dq0(size(current_target_qn)));dq0=0d0
    do io=1,Nso
       iorb=mod(io-1,Norb)+1
       ispin=(io-1)/Norb+1
       !
       key="C"//dot(posA)%okey(iorb,ispin,ilink="n")
       C=dot(posA)%operators%op(key)
       NopA(io)=matmul(hconjg(C),C)
       call C%free()
       !
       key="C"//dot(posB)%okey(iorb,ispin,ilink="n")
       C=dot(posB)%operators%op(key)
       NopB(io)=matmul(hconjg(C),C)
       call C%free()
    enddo
    !
    do io=1,Nso
       do jo=1,Nso
          corr(io,jo)=Measure_Corr_ops_DMRG(NopA(io),dq0,NopB(jo),dq0,&
               posA,posB,"bosonic","bosonic",connected)
       enddo
    enddo
    !Every NopB(jo) is used by every row io.  Freeing NopB(io) inside
    !the preceding loop would invalidate it before the next row.
    do io=1,Nso
       call NopA(io)%free()
       call NopB(io)%free()
    enddo
  end function Measure_DensityDensity_DMRG




  !##################################################################
  !                 MEASURE ENERGY COMPONENTS
  !##################################################################  
  !> Return the energy of the effective DMRG Hamiltonian and split it as
  !>
  !>   Etotal = Ebond + Eloc.
  !>
  !> Etotal is evaluated from the Hamiltonians already stored in the two
  !> final blocks plus the bond(s) connecting them.  This is both faster
  !> and, after truncation, more accurate than reconstructing every bond
  !> from physical two-point functions.  In fact, inside a renormalized
  !> block the projected product (U^+ A U)(U^+ B U) is not generally
  !> equal to U^+ A B U, whereas the stored block Hamiltonian is exactly
  !> the operator used in the final superblock diagonalization.
  !>
  !> Eloc and its optional scalar components are measured independently.
  !> Site-resolved local energies are deliberately left to
  !> Measure_LocalEnergy_DMRG.  Ebond is then obtained as Etotal-Eloc and
  !> includes all non-local terms represented by Hij.
  !>
  !> The optional scalar outputs precede the optional operators so that
  !> an obsolete positional call passing Kij or Hi fails at compile time.
  subroutine Measure_Energy_DMRG(Hij,Ebond,Eloc,Etotal,E0loc,Eint,Eshift,H0loc,Hint,Hshift)
#ifdef _CMPLX
    complex(8),intent(in)                    :: Hij(:,:)
#else
    real(8),intent(in)                       :: Hij(:,:)
#endif
    real(8),intent(out)                      :: Ebond,Eloc,Etotal
    real(8),optional,intent(out)             :: E0loc,Eint,Eshift
    type(sparse_matrix),optional,intent(in)  :: H0loc,Hint,Hshift
    type(sparse_matrix)                      :: Hleft,Hright
    character(len=1)                         :: site_type
    real(8)                                  :: Econnect
    integer                                  :: L,N
    !
    if(.not.allocated(dot))stop "Measure_Energy_DMRG ERROR: DMRG sites are not initialized"
    if(.not.left%operators%has_key("H"))&
         stop "Measure_Energy_DMRG ERROR: missing left-block Hamiltonian"
    if(.not.right%operators%has_key("H"))&
         stop "Measure_Energy_DMRG ERROR: missing right-block Hamiltonian"
    !Initialize the measurement and obtain the local contribution first.
    !The same state/maps are then reused for all remaining contractions.
    Eloc=Measure_LocalEnergy_DMRG(H0loc=H0loc,Hint=Hint,Hshift=Hshift,&
         E0loc=E0loc,Eint=Eint,Eshift=Eshift)
    if(.not.measure_status)return
    L=left%length
    N=L+right%length
    site_type=to_lower(dot(1)%SiteType(1:1))
    !Only the bond(s) joining the final blocks have to be contracted
    !explicitly; every internal bond is already contained in Hleft/Hright.
    Econnect=0d0
    select case(site_type)
    case("f","e")
       Econnect=Measure_FermionBond_DMRG(L,L+1,Hij)
       if(PBCdmrg)Econnect=Econnect+Measure_FermionBond_DMRG(1,N,Hij)
    case("s")
       Econnect=Measure_SpinBond_DMRG(L,L+1,Hij)
       if(PBCdmrg)Econnect=Econnect+Measure_SpinBond_DMRG(1,N,Hij)
    case default
       stop "Measure_Energy_DMRG ERROR: unsupported site type"
    end select
    !In an MPI run the full block operators are owned by the master.
    !Broadcast them before the distributed expectation-value kernels.
#ifdef _MPI
    if(MpiStatus)then
       if(MpiMaster)Hleft =left%operators%op("H")
       if(MpiMaster)Hright=right%operators%op("H")
       call Hleft%bcast()
       call Hright%bcast()
    else
       Hleft =left%operators%op("H")
       Hright=right%operators%op("H")
    endif
#else
    Hleft =left%operators%op("H")
    Hright=right%operators%op("H")
#endif
    Etotal=Average_Op_DMRG(Hleft,1)+Average_Op_DMRG(Hright,N)+Econnect
    Ebond =Etotal-Eloc
    call Hleft%free()
    call Hright%free()
  end subroutine Measure_Energy_DMRG






  !> Reconstruct the total kinetic energy as a sum of physical nearest-
  !> neighbour correlators for the uniform hopping matrix Hij.  This is
  !> useful as an observable, but it need not coincide exactly with the
  !> kinetic part of the truncated effective DMRG Hamiltonian.  Use
  !> Measure_Energy_DMRG when an energy consistent with the final DMRG
  !> diagonalization is required.
  !> USE WITH CAUTION: This function is not returning the correct estimate of the
  !> kinetic energy...
  function Measure_KineticEnergy_DMRG(Hij) result(Ekin)
#ifdef _CMPLX
    complex(8),intent(in)                    :: Hij(:,:)
#else
    real(8),intent(in)                       :: Hij(:,:)
#endif
    real(8)                                  :: Ekin,Eij
    integer                                  :: i,N
    !
    Ekin=0d0
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    N=left%length+right%length
    !
    !Open-chain bonds: (1,2),...,(N-1,N).
    if(MpiMaster)call start_timer("get Ekin: open-chain")
    do i=1,N-1
       Eij=Measure_FermionBond_DMRG(i,i+1,Hij)
       Ekin=Ekin+Eij
       if(MpiMaster)call eta(i,N-1)
    enddo
    if(MpiMaster)call stop_timer()
    !The boundary bond is stored as (1,N), preserving the same
    !upper-triangular convention used for all open-chain bonds.
    if(PBCdmrg)then
      if(MpiMaster)call start_timer("get Ekin: PBC terms")
       Eij=Measure_FermionBond_DMRG(1,N,Hij)
       Ekin=Ekin+Eij
       if(MpiMaster)call stop_timer()
    endif
  end function Measure_KineticEnergy_DMRG




  




  !> Measure the complete local Hamiltonian stored with key "H" on
  !> every site.  If local Fock-space operators H0loc, Hint or Hshift
  !> are supplied, their extensive expectation values are evaluated
  !> independently and returned through the matching optional scalar.
  !>
  !> "H" is the operator constructed by spin_site/electron_site and may
  !> contain fields, crystal-field terms and local interactions.  The
  !> optional component operators never enter LIST_OPERATORS or BLOCK:
  !> they are propagated only here, along the stored DMRG growth path.
  !>
  !> H0loc and Hshift are one-body operators and could alternatively be
  !> obtained from a local 1RDM.  Propagating the already contracted local
  !> operators is cheaper when that complete 1RDM is not otherwise needed.
  function Measure_LocalEnergy_DMRG(Hi,H0loc,Hint,Hshift,E0loc,Eint,Eshift) result(Eloc)
    type(sparse_matrix),optional,intent(out) :: Hi
    type(sparse_matrix),optional,intent(in)  :: H0loc,Hint,Hshift
    real(8),optional,intent(out)             :: E0loc,Eint,Eshift
    real(8)                                  :: Eloc,Ei
    type(sparse_matrix)                      :: Hsite
    integer                                  :: i,N
    !
    Eloc=0d0
    if(present(E0loc))E0loc=0d0
    if(present(Eint))Eint=0d0
    if(present(Eshift))Eshift=0d0
    if(present(E0loc).neqv.present(H0loc))&
         stop "Measure_LocalEnergy_DMRG ERROR: H0loc and E0loc must be supplied together"
    if(present(Eint).neqv.present(Hint))&
         stop "Measure_LocalEnergy_DMRG ERROR: Hint and Eint must be supplied together"
    if(present(Eshift).neqv.present(Hshift))&
         stop "Measure_LocalEnergy_DMRG ERROR: Hshift and Eshift must be supplied together"
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    N=left%length+right%length
    if(present(Hi))call Hi%init(N,N)
    !
    !When the complete decomposition is supplied and no site-resolved
    !output is requested, its three O(N) block sums also give Eloc.
    !Avoid the more expensive propagation of H_i from every position.
    if(present(H0loc)) E0loc =Measure_LocalTerm_DMRG(H0loc)
    if(present(Hint))  Eint  =Measure_LocalTerm_DMRG(Hint)
    if(present(Hshift))Eshift=Measure_LocalTerm_DMRG(Hshift)
    if(.not.present(Hi).AND.present(H0loc).AND.present(Hint).AND.present(Hshift))then
       Eloc=E0loc+Eint+Eshift
       return
    endif
    !
    !
    !Measure every local Hamiltonian in its physical position.  Hi is
    !diagonal because it is a site-resolved container, not an operator
    !acting in the many-body Hilbert space.
    if(MpiMaster)call start_timer("get Eloc")
    do i=1,N
       if(.not.dot(i)%operators%has_key("H"))&
            stop "Measure_LocalEnergy_DMRG ERROR: missing local H operator"
       Hsite=dot(i)%operators%op("H")
       Ei=Measure_Op_DMRG(Hsite,i)
       Eloc=Eloc+Ei
       if(present(Hi))then
          if(Ei/=0d0)then
#ifdef _CMPLX
             call Hi%insert(cmplx(Ei,0d0,8),i,i)
#else
             call Hi%insert(Ei,i,i)
#endif
          endif
       endif
       call Hsite%free()
       if(MpiMaster)call eta(i,N)
    enddo
    if(MpiMaster)call stop_timer()
    !Each scalar component is accumulated as sum_i O_i with one rotation
    !per DMRG growth step.  This is O(N), unlike measuring every site
    !operator independently, which would repeatedly traverse the blocks.
  end function Measure_LocalEnergy_DMRG




  !Measure the extensive sum of one uniform local operator.  The input
  !matrix acts in the original site Fock basis and is never stored in a
  !site or block operator list.
  function Measure_LocalTerm_DMRG(Olocal) result(Eterm)
    type(sparse_matrix),intent(in) :: Olocal
    type(sparse_matrix)            :: Oleft,Oright
    real(8)                        :: Eterm
    integer                        :: i,N
    N=left%length+right%length
    do i=1,N
       if(Olocal%Nrow/=dot(i)%Dim.OR.Olocal%Ncol/=dot(i)%Dim)&
            stop "Measure_LocalTerm_DMRG ERROR: incompatible local-operator dimension"
    enddo
    Oleft =Build_LocalSum_Block_DMRG(Olocal,'l')
    Oright=Build_LocalSum_Block_DMRG(Olocal,'r')
    Eterm=Average_Op_DMRG(Oleft,1)+Average_Op_DMRG(Oright,N)
    call Oleft%free()
    call Oright%free()
  end function Measure_LocalTerm_DMRG




  !Build sum_i O_i directly in the final basis of one DMRG block.  At
  !each growth step the accumulated operator is first rotated with the
  !stored U matrix and then enlarged by adding O on the new local site.
  function Build_LocalSum_Block_DMRG(Olocal,side) result(Osum)
    type(sparse_matrix),intent(in) :: Olocal
    character(len=1),intent(in)    :: side
    type(sparse_matrix)            :: Osum,Osite,U,Oold,Onew
    integer                        :: L,R,nstep,it,D,Dsite
    L=left%length
    R=right%length
    select case(side)
    case('l');nstep=L
    case('r');nstep=R
    case default;stop "Build_LocalSum_Block_DMRG ERROR: side not in [l,r]"
    end select
    Osum=Olocal
    do it=1,nstep-1
       select case(side)
       case('l');if(MpiMaster)U=left%omatrices%op(key=str(it))
       case('r');if(MpiMaster)U=right%omatrices%op(key=str(it))
       end select
#ifdef _MPI
       if(MpiStatus)then
          call U%bcast()
          Osum=(U%dgr().pm.Osum).pm.U
       else
          Osum=matmul(matmul(U%dgr(),Osum),U)
       endif
#else
       Osum=matmul(matmul(U%dgr(),Osum),U)
#endif
       D=Osum%Nrow
       Dsite=Olocal%Nrow
       Osite=Olocal
       select case(side)
       case('l')
          if(PBCdmrg.AND.mod(it,2)==0)then
             Oold=Id(Dsite).x.Osum
             Onew=Osite.x.Id(D)
          else
             Oold=Osum.x.Id(Dsite)
             Onew=Id(D).x.Osite
          endif
       case('r')
          if(PBCdmrg.AND.mod(it,2)==0)then
             Oold=Osum.x.Id(Dsite)
             Onew=Id(D).x.Osite
          else
             Oold=Id(Dsite).x.Osum
             Onew=Osite.x.Id(D)
          endif
       end select
       Osum=Oold+Onew
       call Osite%free()
       call Oold%free()
       call Onew%free()
    enddo
    call U%free()
  end function Build_LocalSum_Block_DMRG





  !> Return the expectation value of one fermionic hopping bond,
  !> using exactly the convention employed by connect_fermion_blocks:
  !> \f[ K_{ij}=\sum_{ab}\left[
  !> H_{ab}\langle c^\dagger_{ia}c_{jb}\rangle+
  !> H_{ab}^*\langle c^\dagger_{jb}c_{ia}\rangle\right]. \f]
  !> Measure_Corr_ops_DMRG supplies all Jordan--Wigner strings, also
  !> when the two endpoints belong to different DMRG blocks.
  !>
  !> Hij acts in the compound spin-orbital space
  !> \f$a=i_{orb}+(i_{spin}-1)N_{orb}\f$.  Only the directed correlator
  !> \f$G_{ab}(i,j)=\langle c^\dagger_{ia}c_{jb}\rangle\f$ is evaluated;
  !> its Hermitian conjugate is added analytically at the end.
  function Measure_FermionBond_DMRG(posA,posB,Hij) result(Eij)
    integer,intent(in)                    :: posA,posB
#ifdef _CMPLX
    complex(8),intent(in)                 :: Hij(:,:)
    complex(8)                            :: corr,Ebond
#else
    real(8),intent(in)                    :: Hij(:,:)
    real(8)                               :: corr,Ebond
#endif
    real(8)                               :: Eij
    type(sparse_matrix)                   :: Ca,Cb,Cdag
    real(8),allocatable                   :: dqA(:),dqB(:)
    character(len=:),allocatable          :: key
    integer                               :: io,jo,iorb,jorb,ispin,jspin,N,Nso
    real(8),parameter                     :: imag_tol=1d-10
    !
    Eij=0d0
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    !
    N=left%length+right%length
    Nso=Nspin*Norb
    if(posA<1.OR.posA>N)stop "Measure_FermionBond_DMRG ERROR: posA not in [1,Nsites]"
    if(posB<1.OR.posB>N)stop "Measure_FermionBond_DMRG ERROR: posB not in [1,Nsites]"
    if(posA==posB)stop "Measure_FermionBond_DMRG ERROR: equal positions"
    if(size(Hij,1)/=Nso.OR.size(Hij,2)/=Nso)&
         stop "Measure_FermionBond_DMRG ERROR: shape(Hij) != [Nso,Nso]"
    if(dot(posA)%SiteType(1:1)/="F".AND.dot(posA)%SiteType(1:1)/="f")&
         stop "Measure_FermionBond_DMRG ERROR: posA is not a fermion site"
    if(dot(posB)%SiteType(1:1)/="F".AND.dot(posB)%SiteType(1:1)/="f")&
         stop "Measure_FermionBond_DMRG ERROR: posB is not a fermion site"
    !
    !Accumulate one oriented half of the bond Hamiltonian:
    !  Ebond = sum_ab Hij(a,b) <C^+_a(posA) C_b(posB)>.
    !The conjugate half is not measured separately.
    Ebond=zero
    do io=1,Nso
       !Convert the flattened index io to the site key convention.
       iorb=mod(io-1,Norb)+1
       ispin=(io-1)/Norb+1
       key="C"//dot(posA)%okey(iorb,ispin,ilink="n")
       if(.not.dot(posA)%operators%has_key(key))&
            stop "Measure_FermionBond_DMRG ERROR: missing C operator at posA"
       Ca=dot(posA)%operators%op(key)
       dqA=dot(posA)%operators%dq(key)
       !If dq(C) is the annihilation shift, dq(C^+)=-dq(C).
       Cdag=hconjg(Ca)
       !
       do jo=1,Nso
          !A zero hopping does not contribute and requires no operator
          !construction or many-body contraction.
          if(Hij(io,jo)==zero)cycle
          jorb=mod(jo-1,Norb)+1
          jspin=(jo-1)/Norb+1
          key="C"//dot(posB)%okey(jorb,jspin,ilink="n")
          if(.not.dot(posB)%operators%has_key(key))&
               stop "Measure_FermionBond_DMRG ERROR: missing C operator at posB"
          Cb=dot(posB)%operators%op(key)
          dqB=dot(posB)%operators%dq(key)
          !
          !This call handles same-block/LR cases and inserts the full
          !Jordan--Wigner string required by the two odd operators.
          corr=Measure_Corr_ops_DMRG(Cdag,-dqA,Cb,dqB,posA,posB,&
               "fermionic","fermionic")
          Ebond=Ebond+Hij(io,jo)*corr
          !
          call Cb%free()
       enddo
       call Ca%free()
       call Cdag%free()
    enddo
    !The second hopping direction is the Hermitian conjugate of the
    !first one.  Forming 2 Re[...] avoids a redundant DMRG contraction.
#ifdef _CMPLX
    if(abs(aimag(Ebond))>imag_tol*max(1d0,abs(real(Ebond,8))))then
       if(MpiMaster)write(LOGfile,*)"Measure_FermionBond_DMRG WARNING: finite imaginary directed bond",aimag(Ebond)
    endif
    Eij=2d0*real(Ebond,8)
#else
    Eij=2d0*Ebond
#endif
  end function Measure_FermionBond_DMRG





  !> Return the spin-exchange energy of one bond, with the same
  !> convention used by connect_spin_blocks:
  !> \f[ J_{ij}=H_{11}\langle S_i^zS_j^z\rangle+
  !> H_{22}\langle S_i^+S_j^-\rangle+
  !> H_{22}^*\langle S_i^-S_j^+\rangle. \f]
  !>
  !> In the site convention, component 1 is Sz and component 2 is S+.
  !> The current spin Hamiltonian accepts only diagonal Hij: Hij(1,1)
  !> is the longitudinal coupling and Hij(2,2) the transverse one.
  function Measure_SpinBond_DMRG(posA,posB,Hij) result(Eij)
    integer,intent(in)                    :: posA,posB
#ifdef _CMPLX
    complex(8),intent(in)                 :: Hij(:,:)
    complex(8)                            :: corrzz,corrpm,Ediag,Etrans
#else
    real(8),intent(in)                    :: Hij(:,:)
    real(8)                               :: corrzz,corrpm,Ediag,Etrans
#endif
    real(8)                               :: Eij
    type(sparse_matrix)                   :: SzA,SzB,SpA,SpB,SmB
    real(8),allocatable                   :: dqzA(:),dqzB(:),dqpA(:),dqpB(:)
    character(len=:),allocatable          :: key
    integer                               :: N
    real(8),parameter                     :: imag_tol=1d-10
    !
    Eij=0d0
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    !
    N=left%length+right%length
    if(posA<1.OR.posA>N)stop "Measure_SpinBond_DMRG ERROR: posA not in [1,Nsites]"
    if(posB<1.OR.posB>N)stop "Measure_SpinBond_DMRG ERROR: posB not in [1,Nsites]"
    if(posA==posB)stop "Measure_SpinBond_DMRG ERROR: equal positions"
    if(size(Hij,1)/=2.OR.size(Hij,2)/=2)&
         stop "Measure_SpinBond_DMRG ERROR: shape(Hij) != [2,2]"
    if(dot(posA)%SiteType(1:1)/="S".AND.dot(posA)%SiteType(1:1)/="s")&
         stop "Measure_SpinBond_DMRG ERROR: posA is not a spin site"
    if(dot(posB)%SiteType(1:1)/="S".AND.dot(posB)%SiteType(1:1)/="s")&
         stop "Measure_SpinBond_DMRG ERROR: posB is not a spin site"
    if(Hij(1,2)/=zero.OR.Hij(2,1)/=zero)&
         stop "Measure_SpinBond_DMRG ERROR: off-diagonal spin Hij is not implemented"
    !
    !Load Sz and S+ together with their physical quantum-number shifts.
    !S- is obtained by Hermitian conjugation and has shift -dq(S+).
    key="S"//dot(posA)%okey(0,1,ilink="n")
    SzA=dot(posA)%operators%op(key);dqzA=dot(posA)%operators%dq(key)
    key="S"//dot(posB)%okey(0,1,ilink="n")
    SzB=dot(posB)%operators%op(key);dqzB=dot(posB)%operators%dq(key)
    key="S"//dot(posA)%okey(0,2,ilink="n")
    SpA=dot(posA)%operators%op(key);dqpA=dot(posA)%operators%dq(key)
    key="S"//dot(posB)%okey(0,2,ilink="n")
    SpB=dot(posB)%operators%op(key);dqpB=dot(posB)%operators%dq(key)
    SmB=hconjg(SpB)
    !
    !Measure one longitudinal and one directed transverse correlator.
    !The S-_i S+_j contribution is the Hermitian conjugate of the
    !transverse term and is added below as twice its real part.
    corrzz=Measure_Corr_ops_DMRG(SzA,dqzA,SzB,dqzB,posA,posB,&
         "bosonic","bosonic")
    corrpm=Measure_Corr_ops_DMRG(SpA,dqpA,SmB,-dqpB,posA,posB,&
         "bosonic","bosonic")
    Ediag=Hij(1,1)*corrzz
    Etrans=Hij(2,2)*corrpm
#ifdef _CMPLX
    if(abs(aimag(Ediag))>imag_tol*max(1d0,abs(real(Ediag,8))))then
       if(MpiMaster)write(LOGfile,*)"Measure_SpinBond_DMRG WARNING: finite imaginary Sz.Sz bond",aimag(Ediag)
    endif
    Eij=real(Ediag,8)+2d0*real(Etrans,8)
#else
    Eij=Ediag+2d0*Etrans
#endif
    !
    call SzA%free();call SzB%free()
    call SpA%free();call SpB%free();call SmB%free()
  end function Measure_SpinBond_DMRG




  !> Sum the uniform spin-exchange Hamiltonian over all physical bonds.
  !> If requested, Jij stores each bond energy once in the upper triangle.
  !> Thus Espin is the sum of the stored Jij values, without a factor 1/2.
  function Measure_SpinExchangeEnergy_DMRG(Hij,Jij) result(Espin)
#ifdef _CMPLX
    complex(8),intent(in)                    :: Hij(:,:)
#else
    real(8),intent(in)                       :: Hij(:,:)
#endif
    type(sparse_matrix),optional,intent(out) :: Jij
    real(8)                                  :: Espin,Eij
    integer                                  :: i,N
    !
    Espin=0d0
    if(.not.measure_status)call Init_Measure_DMRG()
    if(.not.measure_status)return
    N=left%length+right%length
    if(present(Jij))call Jij%init(N,N)
    !
    !Use the same bond enumeration and sparse storage convention as
    !Measure_KineticEnergy_DMRG.
    do i=1,N-1
       Eij=Measure_SpinBond_DMRG(i,i+1,Hij)
       Espin=Espin+Eij
       if(present(Jij).AND.Eij/=0d0)then
#ifdef _CMPLX
          call Jij%insert(cmplx(Eij,0d0,8),i,i+1)
#else
          call Jij%insert(Eij,i,i+1)
#endif
       endif
    enddo
    if(PBCdmrg)then
       Eij=Measure_SpinBond_DMRG(1,N,Hij)
       Espin=Espin+Eij
       if(present(Jij).AND.Eij/=0d0)then
#ifdef _CMPLX
          call Jij%insert(cmplx(Eij,0d0,8),1,N)
#else
          call Jij%insert(Eij,1,N)
#endif
       endif
    endif
  end function Measure_SpinExchangeEnergy_DMRG
  
  




  !> Return the fermionic grading stored in LIST_OPERATORS:
  !> \f[ p(O)=N_f(O)\pmod 2. \f]
  !> type="fermionic" denotes p=1; P=(-1)^N has type="psign" and p=0.
  function is_odd_fermion_type(otype) result(odd)
    character(len=*),optional :: otype
    character(:),allocatable  :: otype_
    logical                   :: odd
    otype_="";if(present(otype))otype_=to_lower(str(otype))
    odd=.false.
    if(len(otype_)>0)odd=otype_(1:1)=="f"
  end function is_odd_fermion_type




  !##################################################################
  !       BUILD A PRODUCT OF LOCAL FACTORS INSIDE ONE DMRG BLOCK
  !##################################################################
  !> Build \f$\prod_a O_a(i_a)\f$ in one final block basis. Physical
  !> positions are converted to growth indices through b2gMap. Starting
  !> from the first non-trivial factor, each subsequent growth step is
  !> \f[ O\longmapsto U_n^\dagger O U_n
  !>                 \longmapsto (U_n^\dagger O U_n)\otimes X_{n+1}, \f]
  !> where \f$X_{n+1}\f$ is either a requested local factor or identity.
  function Build_Product_Block_DMRG(Ops,positions) result(Oprod)
    type(sparse_matrix),intent(in) :: Ops(:)
    integer,intent(in)             :: positions(:)
    type(sparse_matrix)            :: Oprod
    type(sparse_matrix)            :: U,X
    character(len=1)               :: label
    integer,allocatable            :: growth(:)
    integer                        :: L,R,N,Np,ibeg,iend,it,ipos,j,ifactor
    !
    L=left%length
    R=right%length
    N=L+R
    Np=size(positions)
    if(Np==0.OR.size(Ops)/=Np)&
         stop "Build_Product_Block_DMRG ERROR: incompatible factors and positions"
    if(any(positions<1).OR.any(positions>N))&
         stop "Build_Product_Block_DMRG ERROR: position not in [1,Ldmrg]"
    if(any(positions<=L).AND.any(positions>L))&
         stop "Build_Product_Block_DMRG ERROR: factors belong to different blocks"
    do j=1,Np
       if(count(positions==positions(j))/=1)&
            stop "Build_Product_Block_DMRG ERROR: repeated position"
    enddo
    !
    label='l';if(all(positions>L))label='r'
    allocate(growth(Np))
    do j=1,Np
       growth(j)=growth_index(positions(j))
    enddo
    ifactor=minloc(growth,dim=1)
    !
    ibeg=growth(ifactor)
    iend=L;if(label=='r')iend=R
    !
    Oprod=Build_Op_DMRG(Ops(ifactor),positions(ifactor),set_basis=.false.)
    !
    !At growth step n, rotate the old block and append the local factor
    !at the site entering at n+1 (or the identity when no factor acts):
    !  O_n -> U_n^\dagger O_n U_n -> O_{n+1}=O_n\otimes X_{n+1}.
    do it=ibeg,iend-1
       call get_U_and_rotate(it)
       ipos=physical_position(it+1)
       ifactor=0
       do j=1,Np
          if(positions(j)==ipos)ifactor=j
       enddo
       if(ifactor>0)then
          X=Ops(ifactor)
       else
          X=Id(dot(ipos)%Dim)
       endif
       call enlarge_operator(X,it)
       call X%free()
    enddo
    !
    call U%free()
    deallocate(growth)
    !
  contains
    !Map a physical position to the corresponding block-growth index.
    function growth_index(pos) result(index)
      integer,intent(in) :: pos
      integer            :: index
      if(label=='l')then
         index=b2gMap(pos)
      else
         index=b2gMap(N+1-pos)
      endif
    end function growth_index
    !
    !Inverse map: physical site introduced at a given growth step.
    function physical_position(index) result(pos)
      integer,intent(in) :: index
      integer            :: pos,p
      pos=0
      select case(label)
      case('l')
         do p=1,L
            if(b2gMap(p)==index)pos=p
         enddo
      case('r')
         do p=L+1,N
            if(b2gMap(N+1-p)==index)pos=p
         enddo
      end select
      if(pos==0)stop "Build_Product_Block_DMRG ERROR: growth index not mapped"
    end function physical_position

    subroutine get_U_and_rotate(iter)
      integer,intent(in) :: iter
      select case(label)
      case('l');if(MpiMaster)U=left%omatrices%op(key=str(iter))
      case('r');if(MpiMaster)U=right%omatrices%op(key=str(iter))
      end select
#ifdef _MPI
      if(MpiStatus)then
         call U%bcast()
         Oprod=(U%dgr().pm.Oprod).pm.U
      else
         Oprod=matmul(matmul(U%dgr(),Oprod),U)
      endif
#else
      Oprod=matmul(matmul(U%dgr(),Oprod),U)
#endif
    end subroutine get_U_and_rotate

    subroutine enlarge_operator(Op,iter)
      type(sparse_matrix),intent(in) :: Op
      integer,intent(in)             :: iter
      select case(label)
      case('l')
         if(PBCdmrg)then
            if(mod(iter,2)==0)then
               Oprod=Op.x.Oprod
            else
               Oprod=Oprod.x.Op
            endif
         else
            Oprod=Oprod.x.Op
         endif
      case('r')
         if(PBCdmrg)then
            if(mod(iter,2)==0)then
               Oprod=Oprod.x.Op
            else
               Oprod=Op.x.Oprod
            endif
         else
            Oprod=Op.x.Oprod
         endif
      end select
    end subroutine enlarge_operator
  end function Build_Product_Block_DMRG




  !> Build two odd endpoints in the same block. For i<j:
  !> \f[ O_iO_j=(O_iP_i)P_{i+1}\cdots P_{j-1}O_j,
  !>     \qquad P_k=(-1)^{N_k}. \f]
  !> If the requested order is j,i, one fermionic exchange adds -1.
  function Build_Fermion_Corr_Block_DMRG(OpA,OpB,posA,posB) result(Oij)
    type(sparse_matrix),intent(in) :: OpA,OpB
    integer,intent(in)             :: posA,posB
    type(sparse_matrix)            :: Oij
    type(sparse_matrix),allocatable:: Ops(:)
    type(sparse_matrix)            :: Psite
    integer,allocatable            :: positions(:)
    integer                        :: isite,pmin,pmax,j,Np
    !
    if(posA==posB)stop "Build_Fermion_Corr_Block_DMRG ERROR: equal positions"
    pmin=min(posA,posB)
    pmax=max(posA,posB)
    Np=pmax-pmin+1
    allocate(Ops(Np),positions(Np))
    do j=1,Np
       isite=pmin+j-1
       positions(j)=isite
       if(isite==pmin)then
          if(posA==pmin)then
             Ops(j)=OpA
          else
             Ops(j)=OpB
          endif
          Psite=local_parity_operator(isite)
          Ops(j)=matmul(Ops(j),Psite)
          call Psite%free()
       elseif(isite==pmax)then
          if(posA==pmax)then
             Ops(j)=OpA
          else
             Ops(j)=OpB
          endif
       else
          Ops(j)=local_parity_operator(isite)
       endif
    enddo
    Oij=Build_Product_Block_DMRG(Ops,positions)
    !Changing O_A(posA) O_B(posB) into physical left-to-right order is
    !one exchange of odd operators when A is the right endpoint.
    if(posA>posB)then
       do j=1,Oij%Nrow
          Oij%row(j)%vals=-one*Oij%row(j)%vals
       enddo
    endif
    do j=1,Np
       call Ops(j)%free()
    enddo
    deallocate(Ops,positions)
  end function Build_Fermion_Corr_Block_DMRG




  !> Build one half of a Jordan--Wigner string crossing the L/R cut:
  !> \f[ O_iP_i\cdots P_L \f] for side='l', and
  !> \f[ P_{L+1}\cdots P_{j-1}O_j \f] for side='r'.
  function Build_Fermion_LR_End_DMRG(Op,pos,side) result(Ostring)
    type(sparse_matrix),intent(in) :: Op
    integer,intent(in)             :: pos
    character(len=1),intent(in)    :: side
    type(sparse_matrix)            :: Ostring
    type(sparse_matrix),allocatable:: Ops(:)
    type(sparse_matrix)            :: Psite
    integer,allocatable            :: positions(:)
    integer                        :: L,R,N,isite,pmin,pmax,j,Np
    !
    L=left%length;R=right%length;N=L+R
    select case(side)
    case('l')
       if(pos<1.OR.pos>L)stop "Build_Fermion_LR_End_DMRG ERROR: left endpoint not in L"
       pmin=pos;pmax=L
    case('r')
       if(pos<=L.OR.pos>N)stop "Build_Fermion_LR_End_DMRG ERROR: right endpoint not in R"
       pmin=L+1;pmax=pos
    case default
       stop "Build_Fermion_LR_End_DMRG ERROR: side not in [l,r]"
    end select
    Np=pmax-pmin+1
    allocate(Ops(Np),positions(Np))
    do j=1,Np
       isite=pmin+j-1
       positions(j)=isite
       select case(side)
       case('l')
          if(isite==pos)then
             Psite=local_parity_operator(isite)
             Ops(j)=matmul(Op,Psite)
             call Psite%free()
          else
             Ops(j)=local_parity_operator(isite)
          endif
       case('r')
          if(isite==pos)then
             Ops(j)=Op
          else
             Ops(j)=local_parity_operator(isite)
          endif
       end select
    enddo
    Ostring=Build_Product_Block_DMRG(Ops,positions)
    do j=1,Np
       call Ops(j)%free()
    enddo
    deallocate(Ops,positions)
  end function Build_Fermion_LR_End_DMRG




  !> Return the local parity \f$P_i=(-1)^{N_i}\f$ using the same key
  !> convention employed by enlarge_block and connect_fermion_blocks.
  function local_parity_operator(pos) result(P)
    integer,intent(in)       :: pos
    type(sparse_matrix)      :: P
    character(:),allocatable :: key
    key="P"//dot(pos)%okey(0,0,ilink="n")
    if(.not.dot(pos)%operators%has_key(key))&
         stop "Measure_Corr_DMRG ERROR: missing local fermionic sign operator"
    P=dot(pos)%operators%op(key)
  end function local_parity_operator




  !##################################################################
  !       BUILD A TWO-SITE OPERATOR INSIDE THE SAME DMRG BLOCK
  !##################################################################
  !> Build two parity-even operators before every later truncation:
  !> \f[ \widetilde O_{AB}=U^\dagger(O_AO_B)U. \f]
  !> Projecting them separately would instead produce
  !> \f$(U^\dagger O_AU)(U^\dagger O_BU)\f$, inserting \f$UU^\dagger\f$.
  function Build_Corr_Block_DMRG(OpA,OpB,posA,posB) result(Oij)
    type(sparse_matrix),intent(in) :: OpA,OpB
    integer,intent(in)             :: posA,posB
    type(sparse_matrix)            :: Oij
    type(sparse_matrix)            :: Ops(2)
    integer                        :: positions(2),L
    !
    L=left%length
    if(posA==posB)stop "Build_Corr_Block_DMRG ERROR: equal positions"
    if((posA<=L).neqv.(posB<=L))&
         stop "Build_Corr_Block_DMRG ERROR: positions belong to different blocks"
    !
    Ops=[OpA,OpB]
    positions=[posA,posB]
    Oij=Build_Product_Block_DMRG(Ops,positions)
    call Ops(1)%free();call Ops(2)%free()
  end function Build_Corr_Block_DMRG




  !##################################################################
  !          AVERAGE OF O_LEFT x O_RIGHT ON THE SUPERBLOCK
  !##################################################################
  !> Contract operators on opposite sides of the superblock cut:
  !> \f[ C=\langle\Psi|O_L\otimes O_R|\Psi\rangle. \f]
  !> For an output sector \f$q\f$, the input sector is
  !> \f[ q'=q-dq_L, \qquad dq_L+dq_R=0. \f]
  !> sp_filter therefore constructs the rectangular maps
  !> \f$O_{L,R}:\mathcal H(q')\rightarrow\mathcal H(q)\f$ before the
  !> serial or distributed tensor-product contraction.
  function Average_Corr_LR_DMRG(Oleft,dqLeft,Oright,dqRight) result(corr)
    type(sparse_matrix),intent(in) :: Oleft,Oright
    real(8),intent(in)             :: dqLeft(:),dqRight(:)
#ifdef _CMPLX
    complex(8)                     :: corr,Otmp
    complex(8),allocatable         :: Ov(:)
#else
    real(8)                        :: corr,Otmp
    real(8),allocatable            :: Ov(:)
#endif
    type(sparse_matrix)            :: Al,Br
    real(8),allocatable            :: qrow(:),qcol(:)
    integer                        :: irow,icol
    real(8),parameter              :: dq_tol=100d0*epsilon(1d0)
    !
    if(any(abs(dqLeft+dqRight)>dq_tol))then
       corr=zero
       return
    endif
    allocate(Ov(size(gs_vector,1)));Ov=zero
    !
    do irow=1,Nsb
       qrow=sb_sector%qn(index=irow)
       qcol=qrow-dqLeft
       if(.not.sb_sector%has_qn(qcol))cycle
       icol=sb_sector%index(qn=qcol)
       !Rows belong to the output sector irow and columns to the input
       !sector icol. The inverse maps make both blocks rectangular.
       Al=sp_filter(Oleft,LI(irow)%states,Lmap(icol)%states,&
            size(LI(icol)%states))
       Br=sp_filter(Oright,RI(irow)%states,Rmap(icol)%states,&
            size(RI(icol)%states))
#ifdef _MPI
       if(MpiStatus)then
          call Apply_AxB_Measure_MPI(Al,Br,irow,icol,gs_vector(:,1),Ov)
       else
          call Apply_AxB_Measure(Al,Br,Offset(irow),Offset(icol),gs_vector(:,1),Ov)
       endif
#else
       call Apply_AxB_Measure(Al,Br,Offset(irow),Offset(icol),gs_vector(:,1),Ov)
#endif
       call Al%free()
       call Br%free()
    enddo
    !
#ifdef _MPI
    if(MpiStatus)then
       Otmp=dot_product(gs_vector(:,1),Ov)
       corr=zero
       call AllReduce_MPI(MpiComm,Otmp,corr)
    else
       corr=dot_product(gs_vector(:,1),Ov)
    endif
#else
    corr=dot_product(gs_vector(:,1),Ov)
#endif
    deallocate(Ov)
  end function Average_Corr_LR_DMRG




  !Apply a rectangular tensor product A x B from col_offset to
  !row_offset. Superblock vectors use the right index as the fast index.
  subroutine Apply_AxB_Measure(Aop,Bop,row_offset,col_offset,v,Ov)
    type(sparse_matrix),intent(in) :: Aop,Bop
    integer,intent(in)             :: row_offset,col_offset
#ifdef _CMPLX
    complex(8),intent(in)          :: v(:)
    complex(8),intent(inout)       :: Ov(:)
    complex(8),allocatable         :: C(:,:)
    complex(8)                     :: val
#else
    real(8),intent(in)             :: v(:)
    real(8),intent(inout)          :: Ov(:)
    real(8),allocatable            :: C(:,:)
    real(8)                        :: val
#endif
    integer                        :: ai,aj,bi,bj,ja,jb,jc,i,j
    !
    if(.not.Aop%status.OR..not.Bop%status)return
    allocate(C(Bop%Nrow,Aop%Ncol));C=zero
    do aj=1,Aop%Ncol
       do bi=1,Bop%Nrow
          do jb=1,Bop%row(bi)%Size
             bj=Bop%row(bi)%cols(jb)
             val=Bop%row(bi)%vals(jb)
             jc=bj+(aj-1)*Bop%Ncol
             j=jc+col_offset
             C(bi,aj)=C(bi,aj)+val*v(j)
          enddo
       enddo
    enddo
    do bi=1,Bop%Nrow
       do ai=1,Aop%Nrow
          i=bi+(ai-1)*Bop%Nrow+row_offset
          do ja=1,Aop%row(ai)%Size
             aj=Aop%row(ai)%cols(ja)
             val=Aop%row(ai)%vals(ja)
             Ov(i)=Ov(i)+val*C(bi,aj)
          enddo
       enddo
    enddo
    deallocate(C)
  end subroutine Apply_AxB_Measure



#ifdef _MPI
  !Distributed version of Apply_AxB_Measure. The intermediate matrix is
  !transposed between the left- and right-distributed layouts exactly as
  !in the direct superblock H*v implementation.
  subroutine Apply_AxB_Measure_MPI(Aop,Bop,k,q,v,Ov)
    type(sparse_matrix),intent(in) :: Aop,Bop
    integer,intent(in)             :: k,q
#ifdef _CMPLX
    complex(8),intent(in)          :: v(:)
    complex(8),intent(inout)       :: Ov(:)
    complex(8),allocatable         :: C(:,:),Ct(:,:),vt(:),Ovt(:)
    complex(8)                     :: val
#else
    real(8),intent(in)             :: v(:)
    real(8),intent(inout)          :: Ov(:)
    real(8),allocatable            :: C(:,:),Ct(:,:),vt(:),Ovt(:)
    real(8)                        :: val
#endif
    integer                        :: ai,aj,bi,bj,ja,jb,jc,i,j
    integer                        :: mpiArow,mpiAcol,mpiBrow
    integer                        :: abcomm,i_start,i_end
    !
    if(.not.Aop%status.OR..not.Bop%status)return
    mpiAcol=mpiDls(q)
    mpiArow=mpiDls(k)
    mpiBrow=mpiDrs(k)
    allocate(C(Bop%Nrow,mpiAcol));C=zero
    do aj=1,mpiAcol
       do bi=1,Bop%Nrow
          do jb=1,Bop%row(bi)%Size
             bj=Bop%row(bi)%cols(jb)
             val=Bop%row(bi)%vals(jb)
             jc=bj+(aj-1)*Bop%Ncol
             j=jc+mpiOffset(q)
             C(bi,aj)=C(bi,aj)+val*v(j)
          enddo
       enddo
    enddo
    !Use the larger active communicator while changing distribution.
    if(mpiNactive(q)>mpiNactive(k))then
       abcomm=mpiSBCOMM(q)
    else
       abcomm=mpiSBCOMM(k)
    endif
    allocate(Ct(Aop%Ncol,mpiBrow));Ct=zero
    call vector_transpose_MPI(Bop%Nrow,mpiAcol,C,Aop%Ncol,mpiBrow,Ct,abcomm)
    allocate(Ovt(Aop%Nrow*mpiBrow));Ovt=zero
    do bi=1,mpiBrow
       do ai=1,Aop%Nrow
          i=ai+(bi-1)*Aop%Nrow
          do ja=1,Aop%row(ai)%Size
             aj=Aop%row(ai)%cols(ja)
             val=Aop%row(ai)%vals(ja)
             Ovt(i)=Ovt(i)+val*Ct(aj,bi)
          enddo
       enddo
    enddo
    allocate(vt(Bop%Nrow*mpiArow));vt=zero
    call vector_transpose_MPI(Aop%Nrow,mpiBrow,Ovt,Bop%Nrow,mpiArow,vt,mpiSBCOMM(k))
    i_start=1+mpiOffset(k)
    i_end=Bop%Nrow*mpiArow+mpiOffset(k)
    Ov(i_start:i_end)=Ov(i_start:i_end)+vt
    deallocate(C,Ct,vt,Ovt)
  end subroutine Apply_AxB_Measure_MPI
#endif


  


  



 
  !##################################################################
  !              BUILD LOCAL OPERATOR 
  !Purpose: return the O(i) at a site I of the chain given an 
  !         operator O in the local dot basis:
  !##################################################################
  function Build_Op_dmrg(Op,pos,set_basis) result(Oi)
    type(sparse_matrix),intent(in) :: Op
    integer                        :: pos
    type(sparse_matrix)            :: Oi
    logical,optional               :: set_basis   
    !
    character(len=1)               :: label
    type(sparse_matrix)            :: U
    integer                        :: L,R,N
    integer                        :: i,dB(2),d,dOp(2)
    logical                        :: set_basis_
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Build Op"
#endif
    !
    set_basis_ = .false. ;if(present(set_basis))set_basis_=set_basis
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length
    R = right%length
    N = L+R
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Build_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    !Get index in the block from the position pos in the chain:
    !i = M(pos)
    !recall that M: OBC: 1+2+...Ldmrg-2+Ldmrg-1+Ldmrg
    !               PBC: Ldmrg+Ldmrg-2..+1+..+Ldmrg-1
    if(pos<=L)then
       i=b2gMap(pos)
    else
       i=b2gMap(N+1-pos)
    endif
    !
    !Build Operator on the chain at position pos:   
    if(i==1)then
       Oi = Op
    else
       dOp = shape(Op)
       select case(label)
       case('l')
          if(MpiMaster)then
             if(left%omatrices%has_key(str(i-1)))then
                dB = shape(left%omatrices%op(key=str(i-1)));D=dB(2)
             else
                dB = shape(left%omatrices%op(key=str(i)));D=dB(1)/dOp(1)
             endif
          endif
#ifdef _MPI
          if(MpiStatus)call Bcast_MPI(MpiComm,D)
#endif
          if(PBCdmrg)then
             if(mod(i,2)==0)then
                Oi = Id(D).x.Op !o--x
             else
                Oi = Op.x.Id(D) !x--o
             endif
          else
             Oi = Id(D).x.Op    !o--x
          endif
       case('r')
          if(MpiMaster)then
             if(right%omatrices%has_key(str(i-1)))then
                dB = shape(right%omatrices%op(key=str(i-1)));D=dB(2)
             else
                dB = shape(right%omatrices%op(key=str(i)));D=dB(1)/dOp(1)
             endif
          endif
#ifdef _MPI
          if(MpiStatus)call Bcast_MPI(MpiComm,D)
#endif
          if(PBCdmrg)then
             if(mod(i,2)==0)then
                Oi = Op.x.Id(D) !x--o
             else
                Oi = Id(D).x.Op !o--x
             endif
          else
             Oi = Op.x.Id(D)    !x--o
          endif
       end select
       !
       !Set local Basis if required
       if(set_basis_)call Urotate()
       !
    endif
    !
    call U%free()
    !
  contains
    !
    subroutine Urotate()
      select case(label)
      case("l")
         if(MpiMaster)U = left%omatrices%op(key=str(i))
      case("r")
         if(MpiMaster)U = right%omatrices%op(key=str(i))
      end select
#ifdef _MPI
      if(MpiStatus)then
         call U%bcast(MpiComm)
         Oi = U%dgr().pm.(Oi.pm.U)
      else
         Oi = matmul(U%dgr(),matmul(Oi,U))
      endif
#else
      Oi = matmul(U%dgr(),matmul(Oi,U))
#endif
    end subroutine Urotate

  end function Build_Op_dmrg







  !##################################################################
  !                   ADVANCE OPERATOR 
  !Purpose: advance the operator O(i) Nstep from site I 
  !##################################################################
  function Advance_Op_dmrg(Op,pos,nstep) result(Oi)
    type(sparse_matrix),intent(in)   :: Op
    integer                          :: pos
    integer,optional                 :: nstep
    type(sparse_matrix)              :: Oi,U
    character(len=1)                 :: label
    integer                          :: L,R,N
    integer                          :: i,istart,iend,it
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Advance Op"
#endif
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length
    R = right%length
    N = L+R
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Advance_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    !Get index in the block from the position pos in the chain:
    if(pos<=L)then
       i=b2gMap(pos)
    else
       i=b2gMap(N+1-pos)
    endif
    !
    istart  = i
    select case(label)
    case ("l")
       istart = i ; iend   = L-1 ; if(present(nstep))iend=istart+nstep
       if(iend>L-1)stop "Advance_Op_DMRG ERROR: iend > L-1"
    case ("r") 
       istart = i ; iend   = R-1 ; if(present(nstep))iend=istart+nstep
       if(iend>R-1)stop "Advance_Op_DMRG ERROR: iend > R-1"
    end select
    !
    !Evolve to SB basis
    Oi = Op
    select case(label)
    case ("l") 
       do it=istart,iend
          if(MpiMaster) U  = left%omatrices%op(key=str(it))
          call Urotate()
          if(PBCdmrg)then
             if(mod(it,2)==0)then
                Oi = Id(dot(it)%dim).x.Oi
             else
                Oi = Oi.x.Id(dot(it)%dim)
             endif
          else
             Oi = Oi.x.Id(dot(it)%dim)
          endif
       enddo
    case ("r") 
       do it=istart,iend
          if(MpiMaster) U  = right%omatrices%op(key=str(it))
          call Urotate()
          if(PBCdmrg)then
             if(mod(it,2)==0)then
                Oi = Oi.x.Id(dot(it)%dim)
             else
                Oi = Id(dot(it)%dim).x.Oi
             endif
          else
             Oi = Id(dot(it)%dim).x.Oi
          endif
       enddo
    end select
    !
    call U%free()
    !
    !
  contains
    !
    subroutine Urotate()
#ifdef _MPI
      if(MpiStatus)then
         call U%bcast()
         Oi = (U%dgr().pm.Oi).pm.U
      else
         Oi = matmul(matmul(U%dgr(),Oi),U)
      endif
#else
      Oi = matmul(matmul(U%dgr(),Oi),U)
#endif
    end subroutine Urotate
    !
  end function Advance_Op_dmrg







  !##################################################################
  !                   ADVANCE CORRELATION FUNCTION 
  !Purpose: advance the correlation O(i) Nstep from site I 
  !##################################################################
  function Advance_Corr_dmrg(Op,pos,nstep) result(Oi)
    type(sparse_matrix),intent(in)   :: Op
    integer                          :: pos
    integer,optional                 :: nstep
    type(sparse_matrix)              :: Oi,U
    character(len=1)                 :: label
    integer                          :: L,R,N
    integer                          :: i,istart,iend,it
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Advance Correlator"
#endif
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length             !
    R = right%length            !
    N = L+R                       !== Ldmrg
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Advance_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    !Get index in the block from the position pos in the chain:
    if(pos<=L)then
       i=b2gMap(pos)
    else
       i=b2gMap(N+1-pos)
    endif
    !
    istart  = i
    select case(label)
    case ("l")
       istart = i ; iend   = L-1 ; if(present(nstep))iend=istart+nstep
       if(iend>L-1)stop "Advance_Op_DMRG ERROR: iend > L-1"
    case ("r") 
       istart = i ; iend   = R-1 ; if(present(nstep))iend=istart+nstep
       if(iend>R-1)stop "Advance_Op_DMRG ERROR: iend > R-1"
    end select
    !
    !
    Oi = Op
    select case(label)
    case ("l")
       do it=istart+1,iend
          if(MpiMaster)U  = left%omatrices%op(key=str(it))
          call Urotate
          if(PBCdmrg)then
             if(mod(it,2)==0)then
                Oi = Id(dot(it)%dim).x.Oi
             else
                Oi = Oi.x.Id(dot(it)%dim)
             endif
          else
             Oi = Oi.x.Id(dot(it)%dim)
          endif
       enddo
    case ("r")
       do it=istart+1,iend
          if(MpiMaster)U  = right%omatrices%op(key=str(it))
          call Urotate()
          if(PBCdmrg)then
             if(mod(it,2)==0)then
                Oi = Oi.x.Id(dot(it)%dim)
             else
                Oi = Id(dot(it)%dim).x.Oi
             endif
          else
             Oi = Id(dot(it)%dim).x.Oi
          endif
       enddo
    end select
    !
    call U%free()
    !
  contains
    !
    subroutine Urotate()
#ifdef _MPI
      if(MpiStatus)then
         call U%bcast()
         Oi = (U%dgr().pm.Oi).pm.U
      else
         Oi = matmul(matmul(U%dgr(),Oi),U)
      endif
#else
      Oi = matmul(matmul(U%dgr(),Oi),U)
#endif
    end subroutine Urotate

  end function Advance_Corr_dmrg





  !##################################################################
  !                   AVERAGE OPERATOR 
  !Purpose: take the average of an operator O on the last step basis 
  !##################################################################
  function Average_Op_dmrg(Oi,pos) result(Oval)
    type(sparse_matrix),intent(in)   :: Oi
    integer                          :: pos
    character(len=1)                 :: label
    real(8)                          :: Oval,Otmp
    type(sparse_matrix)              :: Psi
    integer                          :: L,R,N
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Average Op"
#endif
    !
    !The lenght of the last block contributing to the SB construction-> \psi
    L = left%length
    R = right%length
    N = L+R
    !
    !Check:
    if(pos<1.OR.pos>N)stop "Average_op_dmrg error: Pos not in [1,Ldmrg]"
    !
    !Get label of the block holding the site at position pos:
    label='l'; if(pos>L)label='r'
    !
    allocate(Olist(Nsb))
    !
    !Measure using PSI matrix:
    do isb=1,Nsb
       select case(label)
       case default;stop "Average_op_dmrg error: label not [l,r]"
       case ("l");Olist(isb) = sp_filter(Oi,LI(isb)%states)
       case ("r");Olist(isb) = sp_filter(Oi,RI(isb)%states)
       end select
    enddo
    !
#ifdef _MPI
    if(MpiStatus)then
       Otmp = dot_product(gs_vector(:,1), OdotV_MPI_direct(Olist,gs_vector(:,1),label))
       Oval = zero
       call AllReduce_MPI(MpiComm,Otmp,Oval)
    else
       Oval = dot_product(gs_vector(:,1), OdotV_direct(Olist,gs_vector(:,1),label))
    endif
#else
    Oval = dot_product(gs_vector(:,1), OdotV_direct(Olist,gs_vector(:,1),label))
#endif
    !
    do isb=1,Nsb
       call Olist(isb)%free()
    enddo
    deallocate(Olist)
    !
  end function Average_Op_dmrg


  !#################################
  !#################################


  function OdotV_direct(Op,v,direction) result(Ov)
    integer                          :: Nsb,Nloc
    type(sparse_matrix),dimension(:) :: Op
    character(len=*)                 :: direction
#ifdef _CMPLX
    complex(8),dimension(:)             :: v
    complex(8),dimension(size(v))       :: Ov
    complex(8)                          :: val
#else
    real(8),dimension(:)             :: v
    real(8),dimension(size(v))       :: Ov
    real(8)                          :: val
#endif
    integer                          :: i,j,k,n
    integer                          :: ir,il,jr,jl,it
    integer                          :: ia,ib,ic,ja,jb,jc,jcol
    !
    Ov=zero
    !> loop over all the SB sectors:
    sector: do  k=1,size(sb_sector)
       select case(to_lower(direction))
       case("l","left","sys","s")
          !> apply the H^L x 1^r: need to T v and Ov
          do ir=1,Drs(k)
             do il=1,Dls(k)
                i = ir + (il-1)*Drs(k) + offset(k)
                do jcol=1,Op(k)%row(il)%Size
                   val = Op(k)%row(il)%vals(jcol)
                   jl  = Op(k)%row(il)%cols(jcol)
                   j   = ir + (jl-1)*Drs(k) + offset(k)
                   Ov(i) = Ov(i) + val*v(j)
                end do
             enddo
          enddo
          !
       case("r","right","env","e")
          !> apply the 1^L x H^r
          do il=1,Dls(k)
             do ir=1,Drs(k)
                i = ir + (il-1)*Drs(k) + offset(k)           
                do jcol=1,Op(k)%row(ir)%Size
                   val = Op(k)%row(ir)%vals(jcol)
                   jr  = Op(k)%row(ir)%cols(jcol)
                   j   = jr + (il-1)*Drs(k) + offset(k)
                   Ov(i) = Ov(i) + val*v(j)
                end do
             enddo
          enddo
          !
       end select
       !
    enddo sector
  end function OdotV_direct


#ifdef _MPI
  function OdotV_MPI_direct(Op,v,direction) result(Ov)
    integer                             :: Nsb,Nloc
    type(sparse_matrix),dimension(:)    :: Op
    character(len=*)                    :: direction
#ifdef _CMPLX
    complex(8),dimension(:)             :: v
    complex(8),dimension(size(v))       :: Ov
    complex(8)                          :: val
    complex(8),dimension(:),allocatable :: vt,Hvt
#else
    real(8),dimension(:)                :: v
    real(8),dimension(size(v))          :: Ov
    real(8)                             :: val
    real(8),dimension(:),allocatable    :: vt,Hvt
#endif
    integer                             :: i,j,k,n
    integer                             :: ir,il,jr,jl,it
    integer                             :: ia,ib,ic,ja,jb,jc,jcol
    integer                             :: i_start,i_end
    !
    Ov=zero
    !> loop over all the SB sectors:
    sector: do  k=1,size(sb_sector)
       select case(to_lower(direction))
       case("l","left","sys","s")
          !> apply the H^L x 1^r: need to T v and Ov
          allocate(vt(mpiDrs(k)*Dls(k))) ;vt=zero
          allocate(Hvt(mpiDrs(k)*Dls(k)));Hvt=zero
          i_start = 1 + mpiOffset(k)
          i_end   = mpiDl(k)+mpiOffSet(k)
          call vector_transpose_MPI(Drs(k),mpiDls(k),v(i_start:i_end),Dls(k),mpiDrs(k),vt, mpiSBCOMM(k))
          do il=1,mpiDrs(k)
             do ir=1,Dls(k)
                i = ir + (il-1)*Dls(k)
                do jcol=1,Op(k)%row(ir)%Size
                   val = Op(k)%row(ir)%vals(jcol)
                   jr  = Op(k)%row(ir)%cols(jcol)
                   j   = jr + (il-1)*Dls(k)
                   Hvt(i) = Hvt(i) + val*vt(j)
                end do
             enddo
          enddo
          deallocate(vt) ; allocate(vt(Drs(k)*mpiDls(k))) ; vt=zero
          call vector_transpose_MPI(Dls(k),mpiDrs(k),Hvt,Drs(k),mpiDls(k),vt, mpiSBCOMM(k))
          Ov(i_start:i_end) = Ov(i_start:i_end) + Vt
          deallocate(vt,Hvt)
          !
       case("r","right","env","e")
          !> apply the 1^L x H^r
          do il=1,mpiDls(k)
             do ir=1,Drs(k)
                i = ir + (il-1)*Drs(k) + mpiOffset(k)
                do jcol=1,Op(k)%row(ir)%Size
                   val = Op(k)%row(ir)%vals(jcol)
                   jr  = Op(k)%row(ir)%cols(jcol)
                   j   = jr + (il-1)*Drs(k) + mpiOffset(k)
                   Ov(i) = Ov(i) + val*v(j)
                end do
             enddo
          enddo
          !
       end select
       !
    enddo sector
  end function OdotV_MPI_direct
#endif







!##################################################################
!                        WRITE OUTPUT
!##################################################################
  subroutine write_user_scalar(file,val,x)
    character(len=*) :: file
    real(8)          :: val
    integer,optional :: x
    integer          :: x_
    integer          :: i,Eunit
    x_ = left%length;if(present(x))x_=x
    Eunit     = fopen(str(file)//str(suffix),append=.true.)
    write(Eunit,*)x_,val
    close(Eunit)
  end subroutine write_user_scalar

  subroutine write_user_array(file,vals,x)
    character(len=*) :: file
    real(8)          :: vals(:)
    integer,optional :: x(size(vals))
    integer          :: x_(size(vals))
    integer          :: i,Eunit
    x_=arange(1,size(vals));if(present(x))x_=x
    Eunit     = fopen(str(file)//str(suffix),append=.true.)
    do i=1,size(vals)
       write(Eunit,*)x_(i),vals(i)
    enddo
    close(Eunit)
  end subroutine write_user_array

  subroutine write_user_matrix(file,vals,x)
    character(len=*) :: file
    real(8)          :: vals(:,:) !M,N
    integer,optional :: x(size(vals,2))
    integer          :: x_(size(vals,2))
    integer          :: i,j,Eunit
    x_=arange(1,size(vals,2));if(present(x))x_=x
    Eunit     = fopen(str(file)//str(suffix),append=.true.)
    do i=1,size(vals,2)
       write(Eunit,*)x_(i),(vals(j,i),j=1,size(vals,1))
    enddo
    close(Eunit)
  end subroutine write_user_matrix


















END MODULE DMRG_MEASURE
