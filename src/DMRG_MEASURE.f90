module DMRG_MEASURE
  USE SCIFOR, only: to_lower
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
  USE DMRG_SUPERBLOCK
  USE DMRG_SUPERBLOCK_SETUP
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
  public :: Write_DMRG

  interface Write_DMRG
     module procedure :: write_user_scalar
     module procedure :: write_user_array
     module procedure :: write_user_matrix
  end interface Write_DMRG


  interface  Measure_DMRG
     module procedure :: Measure_DMRG_scalar
     module procedure :: Measure_DMRG_vector
  end interface Measure_DMRG


  integer                                      :: Nsb,isb,vecDim
  real(8),dimension(:),allocatable             :: qn,qm
  real(8),dimension(:),allocatable             :: dq
  type(sparse_matrix),allocatable,dimension(:) :: Olist
  type(tstates),dimension(:),allocatable       :: Li,Ri
  logical                                      :: measure_status=.false.

  character(:),allocatable                     :: string

contains


  !##################################################################
  !          INIT / END MEASUREMENT: allocate/deallocate
  !##################################################################
  subroutine Init_Measure_dmrg(msg)
    real(8),dimension(:),allocatable :: qn
    character(len=*),optional        :: msg
    integer,dimension(2)             :: omat_dims
    integer                          :: ilat,i,f,m
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: init measure"
#ifdef _MPI
    if(MpiStatus.AND.MpiMaster)write(LOGfile,*)"DEBUG: using MPI"
#endif
#endif
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
    call sb_build_dims()
    !
    Nsb  = size(sb_sector)
    !
    allocate(LI(Nsb))
    allocate(RI(Nsb))
    do isb=1,Nsb
       qn  = sb_sector%qn(index=isb)
       LI(isb)%states = sb2block_states(qn,'left')
       RI(isb)%states = sb2block_states(qn,'right')
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
       ! g2bMap(1) = f
       do i=1,m
          b2gMap(f+i)  = 2*i
          b2gMap(f-i)  = 2*i+1
          ! g2bMap(2*i)  = f+i
          ! g2bMap(2*i+1)= f-i
       enddo
    else
       b2gMap = (/(i,i=1,Ldmrg)/)
       ! g2bMap = b2gMap
    endif
  end subroutine Init_Measure_dmrg



  subroutine End_measure_DMRG()
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: end measure"
#endif
    if(MpiMaster)call stop_timer("Done "//str(string))
    if(allocated(string))deallocate(string)
    if(allocated(Olist))deallocate(Olist)
    if(allocated(Li))deallocate(Li)
    if(allocated(Ri))deallocate(Ri)
    call sb_delete_dims()
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
    if(.not.measure_status)call Init_Measure_DMRG()
    Oi   = Build_Op_dmrg(Op,pos)
    Oi   = Advance_Op_dmrg(Oi,pos)
    avOp = Average_Op_dmrg(Oi,pos)
    call Oi%free()
  end function Measure_Op_DMRG






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
    L = left%length ; R = right%length
    Np = L+R;if(present(pos))Np= size(pos)
    !
    allocate(vals(Np))
    allocate(pos_(Np))
    pos_=arange(1,Np);if(present(pos))pos_=pos
    !
    call Init_measure_dmrg(msg)
    if(.not.measure_status)then
       call Error_measure_DMRG
       return
    endif
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
    call Init_measure_dmrg(msg)
    if(.not.measure_status)then
       call Error_measure_DMRG
       return
    endif
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
    integer                        :: i,dB(2),d
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
    i=b2gMap(pos)    ; if(pos>L)i=b2gMap(N+1-pos)
    !
    !Build Operator on the chain at position pos:   
    if(i==1)then
       Oi = Op
    else
       select case(label)
       case('l')
          if(MpiMaster)dB = shape(left%omatrices%op(index=i-1));D=dB(2)
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
          if(MpiMaster)dB = shape(right%omatrices%op(index=i-1));D=dB(2)
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
         if(MpiMaster)U = left%omatrices%op(index=i)
      case("r")
         if(MpiMaster)U = right%omatrices%op(index=i)
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
    i=b2gMap(pos)    ; if(pos>L)i=b2gMap(N+1-pos)
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
          if(MpiMaster) U  = left%omatrices%op(index=it)
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
          if(MpiMaster) U  = right%omatrices%op(index=it)
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
    i=b2gMap(pos)    ; if(pos>L)i=b2gMap(N+1-pos)
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
          if(MpiMaster)U  = left%omatrices%op(index=it)
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
          if(MpiMaster)U  = right%omatrices%op(index=it)          
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






  !#################################
  !#################################



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
