module DMRG_MEASURE
  USE SCIFOR, only: to_lower
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
  public :: Measure_SpinSpin_DMRG
  public :: Measure_DensityDensity_DMRG
  public :: Write_DMRG

  interface Measure_Corr_DMRG
     module procedure :: Measure_Corr_keys_DMRG
     module procedure :: Measure_Corr_ops_DMRG
  end interface Measure_Corr_DMRG

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


  !##################################################################
  !          INIT / END MEASUREMENT: allocate/deallocate
  !##################################################################
  subroutine Init_Measure_dmrg(msg)
    real(8),dimension(:),allocatable :: qn
    character(len=*),optional        :: msg
    integer,dimension(2)             :: omat_dims
    integer                          :: ilat,i,f,m,istate
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
    !Load Measure State (if required) (soft check)
    if(need_measure_state)then
       call sb_load_measure_state(found_measure_state)
       if(.not.found_measure_state)then
          if(MpiMaster)write(LOGfile,*)"Init_Measure_DMRG: no saved SuperBlock measurement state found."
       endif
    endif
    !
    !
    if(MpiMaster)omat_dims = [size(left%omatrices),size(right%omatrices)]
#ifdef _MPI
    call Bcast_MPI(MpiComm,omat_dims)
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
    allocate(b2gMap(Ldmrg))
    if(PBCdmrg)then
       !Ldmrg = 2*m+1
       f = (Ldmrg+1)/2          !mid-point
       m = (Ldmrg-1)/2
       !
       b2gMap(f) = 1
       do i=1,m
          b2gMap(f+i)  = 2*i
          b2gMap(f-i)  = 2*i+1
       enddo
    else
       b2gMap = (/(i,i=1,Ldmrg)/)
    endif
  end subroutine Init_Measure_dmrg





  subroutine End_measure_DMRG()
    type(sparse_matrix) :: Ileft,Iright
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: end measure"
#endif
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











!##################################################################
!                   AVERAGE OPERATOR 
!Purpose: take the average of an operator O on the last step basis 
!##################################################################
! function Average_Op_dmrg(Oi,pos) result(Oval)
!   type(sparse_matrix),intent(in)   :: Oi
!   integer                          :: pos
!   character(len=1)                 :: label
!   real(8)                          :: Oval
!   type(sparse_matrix)              :: Psi
!   integer                          :: L,R,N
!   !
!   !The lenght of the last block contributing to the SB construction-> \psi
!   L = left%length-1
!   R = right%length-1
!   N = L+R
!   !
!   !Check:
!   if(pos<1.OR.pos>N)stop "Average_op_dmrg error: Pos not in [1,Ldmrg]"
!   !
!   !Get label of the block holding the site at position pos:
!   label='l'; if(pos>L)label='r'
!   !
!   !Measure using PSI matrix:
!   select case(label)
!   case ("l")
!      if(any(shape(psi_left)/=shape(Oi)))stop "average_op_dmrg ERROR: shape(psi_left) != shape(Oi)"
!      Psi  = as_sparse(psi_left)
!      Oval = trace(as_matrix(matmul(matmul(Psi%dgr(),Oi),Psi)))
!   case ("r")
!      if(any(shape(psi_right)/=shape(Oi)))stop "average_op_dmrg ERROR: shape(psi_right) != shape(Oi)"
!      Psi  = as_sparse(psi_right)
!      Oval = trace(as_matrix(matmul(matmul(Psi%dgr(),Oi),Psi)))
!   end select
!   call Psi%free()
! end function Average_Op_dmrg



! subroutine Measure_Op_dmrg(Op,file,ref,avOp)
!   type(sparse_matrix),intent(in)            :: Op
!   character(len=*)                          :: file
!   character(len=1)                          :: label
!   real(8) :: ref
!   real(8),dimension(:),allocatable,optional :: avOp
!   real(8)                                   :: val
!   integer                                   :: it,i,L,R,N,j,pos,dims(2)
!   type(sparse_matrix)                       :: Oi,U,Psi,Ok,I_R,I_L
!   !
!   suffix=label_DMRG('u')
!   !
!   L = left%length-1           !the len of the last block used to create the SB->\psi
!   R = right%length-1
!   N = L+R
!   !
!   if(present(avOp))then
!      if(allocated(avOp))deallocate(avOp)
!      allocate(avOp(N))
!   endif
!   !
!   call start_timer()
!   do pos=1,N
!      Oi = Build_Op_dmrg(Op,pos)
!      Oi = Advance_Op_dmrg(Oi,pos)
!      val= Average_Op_dmrg(Oi,pos)
!      if(present(avOp))avOp(pos)=val
!      call write_user_scalar(trim(file),val,x=pos)
!      call Oi%free()
!      call progress(pos,N)
!   enddo
!   call stop_timer("Done "//str(file))




!   U =  right%omatrices%op(index=R)
!   dims=shape(U)
!   ! I_R = Id(dot%dim).x.Id(dims(2))!(matmul(U%t(),U))
!   I_R = Id(dot%dim*dims(2))
!   U =  left%omatrices%op(index=R)
!   dims=shape(U)
!   ! I_L = Id(dims(2)).x.Id(dot%dim)!(matmul(U%t(),U)).x.Id(dot%dim)
!   I_L = Id(dims(2)*dot%dim)!(matmul(U%t(),U)).x.Id(dot%dim)


!   print*,"size(GSpsi)",size(gs_vector(:,1))



!   print*,""
!   print*,"- - - - - - - - - - - - - - - - -"
!   print*," METHOD 1: O.x.I - I.x.O full"
!   print*,"- - - - - - - - - - - - - - - - -"
!   print*,""


!   print*,""
!   print*,"pos=1"
!   print*,""
!   print*,"Method <psi|O.x.I_R|psi>"
!   print*,shape(I_R)
!   pos=1
!   Oi = Build_Op_Dmrg(Op,pos)
!   do it=1,L
!      U  = left%omatrices%op(index=it)
!      Oi = (matmul(matmul(U%t(),Oi),U)).x.Id(dot%dim)
!   enddo
!   Oi = sp_kron(Oi,I_R,sb_states)
!   print*,shape(Oi)
!   val = dot_product(gs_vector(:,1),Oi%dot(gs_vector(:,1)))
!   print*,pos,val    


!   print*,""
!   print*,"pos=N"
!   print*,""
!   print*,"Method <psi|I_L.x.O|psi>"
!   pos=N
!   print*,shape(I_L)
!   Oi = Build_Op_Dmrg(Op,pos)
!   do it=1,R
!      U  = right%omatrices%op(index=it)
!      Oi = Id(dot%dim).x.(matmul(matmul(U%t(),Oi),U))
!   enddo
!   Oi = sp_kron(I_L,Oi,sb_states)
!   print*,shape(Oi)
!   val = dot_product(gs_vector(:,1),OI%dot(gs_vector(:,1)))
!   print*,pos,val    



!   print*,""
!   print*,"- - - - - - - - - - - - - - - - -"
!   print*," METHOD 3: direct O_L, O_R"
!   print*,"- - - - - - - - - - - - - - - - -"
!   print*,""


!   Nsb  = size(sb_sector)
!   allocate(Dls(Nsb),Drs(Nsb),Offset(Nsb),Oleft(Nsb),Oright(Nsb))
!   allocate(LI(Nsb),RI(Nsb))
!   Offset=0
!   do isb=1,Nsb
!      qn   = sb_sector%qn(index=isb)
!      Dls(isb)= sector_qn_dim(left%sectors(1),qn)
!      Drs(isb)= sector_qn_dim(right%sectors(1),current_target_qn - qn)
!      if(isb>1)Offset(isb)=Offset(isb-1)+Dls(isb-1)*Drs(isb-1)
!      LI(isb)%states = sb2block_states(qn,'left')
!      RI(isb)%states = sb2block_states(qn,'right')
!   enddo

!   pos=1
!   Oi = Build_Op_Dmrg(Op,pos)
!   do it=1,L
!      U  = left%omatrices%op(index=it)
!      Oi = (matmul(matmul(U%t(),Oi),U)).x.Id(dot%dim)
!   enddo
!   do isb=1,Nsb
!      !> get: Oi*^L 
!      Oleft(isb) = sp_filter(Oi,LI(isb)%states)
!   enddo
!   val = dot_product(gs_vector(:,1), OdotV_direct(Oleft,gs_vector(:,1),'left'))
!   print*,pos,val


!   pos=N
!   Oi = Build_Op_Dmrg(Op,pos)
!   do it=1,R
!      U  = right%omatrices%op(index=it)
!      Oi = Id(dot%dim).x.(matmul(matmul(U%t(),Oi),U))
!   enddo
!   do isb=1,Nsb
!      !> get: Oi*^R
!      Oright(isb) = sp_filter(Oi,RI(isb)%states)
!   enddo
!   val = dot_product(gs_vector(:,1), OdotV_direct(Oright,gs_vector(:,1),'right'))
!   print*,pos,val


!   if(allocated(Dls))deallocate(Dls)
!   if(allocated(Drs))deallocate(Drs)
!   if(allocated(Offset))deallocate(Offset)
!   if(allocated(Oleft))deallocate(Oleft)
!   if(allocated(Oright))deallocate(Oright)
!   if(allocated(Li))deallocate(Li)
!   if(allocated(Ri))deallocate(Ri)


! contains



!   function OdotV_direct(Op,v,direction) result(Ov)
!     integer                          :: Nsb,Nloc
!     type(sparse_matrix),dimension(:) :: Op
!     real(8),dimension(:)             :: v
!     character(len=*)                 :: direction
!     real(8),dimension(size(v))       :: Ov
!     real(8)                          :: val
!     integer                          :: i,j,k,n
!     integer                          :: ir,il,jr,jl,it
!     integer                          :: ia,ib,ic,ja,jb,jc,jcol
!     real(8)                          :: aval,bval
!     !
!     Ov=zero
!     !> loop over all the SB sectors:
!     select case(to_lower(direction))
!     case("l","left","sys","s")
!        do k=1,size(sb_sector)
!           !> apply the H^L x 1^r: need to T v and Ov
!           do ir=1,Drs(k)
!              do il=1,Dls(k)
!                 i = ir + (il-1)*Drs(k) + offset(k)
!                 do jcol=1,Op(k)%row(il)%Size
!                    val = Op(k)%row(il)%vals(jcol)
!                    jl  = Op(k)%row(il)%cols(jcol)
!                    j   = ir + (jl-1)*Drs(k) + offset(k)
!                    Ov(i) = Ov(i) + val*v(j)
!                 end do
!              enddo
!           enddo
!        enddo
!     case("r","right","env","e")
!        do k=1,size(sb_sector)
!           !> apply the 1^L x H^r
!           do il=1,Drs(k)
!              do ir=1,Dls(k)
!                 i = il + (ir-1)*Drs(k) + offset(k)           
!                 do jcol=1,Op(k)%row(il)%Size
!                    val = Op(k)%row(il)%vals(jcol)
!                    jl  = Op(k)%row(il)%cols(jcol)
!                    j   = jl + (ir-1)*Drs(k) + offset(k)
!                    Ov(i) = Ov(i) + val*v(j)
!                 end do
!              enddo
!           enddo
!        enddo
!     end select
!     !
!   end function OdotV_direct

! end subroutine Measure_Op_dmrg
