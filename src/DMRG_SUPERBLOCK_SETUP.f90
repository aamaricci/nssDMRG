MODULE DMRG_SUPERBLOCK_SETUP
  USE DMRG_GLOBAL
  USE DMRG_CONNECT
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private

  !Local variables
  integer                                        :: tNso
  integer                                        :: isb,jsb
  integer                                        :: i,j
  integer                                        :: ispin
  integer                                        :: iorb,jorb
  integer                                        :: io,jo


  public :: Setup_SuperBlock_Sparse
  public :: Setup_SuperBlock_Direct
  public :: spMatVec_sparse_main
  public :: spMatVec_direct_main
  public :: spMatVec_MPI_direct_main
  public :: spMatVec_direct_lazy_main
#ifdef _MPI
  public :: spMatVec_MPI_direct_lazy_main
#endif
  !
  !-> used in DMRG_MEASURE to perform H|gs>
  public :: sb2block_states
  !
contains

	


  
  !##################################################################
  !              SETUP THE SUPERBLOCK HAMILTONIAN
  !                      SPARSE MODE
  !    H^SB = H^L x 1^R  + 1^L x H^R + H^LR
  !    H^LR = sum_p A_p x B_p
  ! 
  ! * sparse: get the sparse global SB Hamiltonian spHsb
  !##################################################################
  !POSSIBLY INCLUDE MPI HERE... this is probably a less efficient version
  !note that only the SB Hamiltonian needs to be constructed in parallel form
  !the other operators (small) are stored by each cpu.
  !In principle one could store any sparse matrix in parallel and build H^SB as MPI too.
  subroutine Setup_SuperBlock_Sparse()
    integer                      :: m_left,m_right
    character(len=:),allocatable :: type
    type(sparse_matrix)          :: H2
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Setup SB Sparse"
#endif
    !
    if(MpiMaster)call start_timer("Setup SB Sparse")
    t0=t_start()
    !
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Sparse ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Sparse ERROR: Missing right.H operator in the list"
    !
    type=str(left%type())
    if(type/=str(right%type()))&
         stop "Setup_SuperBlock_Sparse ERROR: left.Type != right.Type"
    !
    m_left = left%dim
    m_right= right%dim
    !
    select case(to_lower(type(1:1)))
    case default;stop "Setup_SuperBlock_Sparse ERROR: wrong left/right.Type"
    case ("s")
       H2 = connect_spin_blocks(left,right,sb_states,link="n")
       if(PBCdmrg)H2 = H2 + connect_spin_blocks(left,right,sb_states,link="p")
    case ("f","e")
       H2 = connect_fermion_blocks(left,right,sb_states,link="n")
       if(PBCdmrg)H2 = H2 + connect_fermion_blocks(left,right,sb_states,link="p")
    end select
    !
    spHsb= H2 & 
         + sp_kron(left%operators%op("H"),id(m_right),sb_states) &
         + sp_kron(id(m_left),right%operators%op("H"),sb_states) 
    !
    if(MpiMaster)call stop_timer("Setup SB Sparse")
    t_setup_sb_sparse=t_stop()
    t_connect_blocks=0d0
    !
  end subroutine Setup_SuperBlock_Sparse










  !##################################################################
  !              SETUP THE SUPERBLOCK HAMILTONIAN
  !                       DIRECT MODE
  !    H^SB = H^L x 1^R  + 1^L x H^R + H^LR
  !    H^LR = sum_k A_k x B_k; k={q,p}, Q.N. + elements in H^LR
  ! 
  ! * direct: get vec{A},vec{B},vec{Hleft},vec{Hright}
  !   the sparse matrices which reconstruct H^SB above in terms of
  !   conserved QN blocks:
  !   + H^a=L,R = [H^a[q1], H^a[q2], ... , [H^a[qN]]
  !     with each H^a[qi] a given block of conserved QN qi
  !   + A or B  = [A[q1,q1'], ... , A[qN,qN']]
  !     with each A[q,q'] connectin the blocks of differend QN q and q'
  !##################################################################
  subroutine Setup_SuperBlock_Direct()
    character(len=:),allocatable                 :: type
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Setup SB Direct"
#endif
    !
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
    !
    type=str(left%type())
    if(type/=str(right%type()))&
         stop "Setup_SuperBlock_Direct ERROR: left.Type != right.Type"
    !
    t0=t_start()
    select case(to_lower(type(1:1)))
    case default;stop "Setup_SuperBlock_Direct ERROR: wrong left/right.Type"
    case ("s")    ;call Setup_SuperBlock_Spin_Direct()
    case ("f","e");call Setup_SuperBlock_Fermion_Direct()
    end select
    t_setup_sb_direct=t_stop()
  end subroutine Setup_SuperBlock_Direct












  !##################################################################
  !                          SPIN CASE
  !##################################################################
  subroutine Setup_SuperBlock_Spin_Direct()
    integer                                      :: Nso,Nsb
    integer                                      :: it,isb,jsb,ierr,sizeA,sizeB
    real(8),dimension(:),allocatable             :: qn,qm
    type(tstates),dimension(:),allocatable       :: Aj,Bj
    real(8),dimension(:),allocatable             :: dq
    integer,dimension(:,:,:),allocatable         :: tMap
    type(sparse_matrix),allocatable,dimension(:) :: Sl_n,Sl_p !Left
    type(sparse_matrix),allocatable,dimension(:) :: Sr_n,Sr_p !Right
    type(sparse_matrix)                          :: Hl,Hr
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable        :: Hij
#else
    real(8),dimension(:,:),allocatable           :: Hij
#endif
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Setup SB Direct - spin"
#endif
    !
    if(MpiMaster)call start_timer("Setup SB Direct, Nsb: "//str(size(sb_sector)))
    !
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
    !
    !
    !> GET THE USER DEFINED MODEL HAMILTONIAN PARAMETERS:
    ! Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)
    !
    !
    Nso  = Nspin
    tNso = 3                    !Sz.Sz + S+.S- + S-.S+ ([...]<->[...]) 
    if(PBCdmrg)tNso=6           !Sz.Sz + S+.S- + S-.S+ (->[...][...]<-)
    nsb  = size(sb_sector)
    !
    !Massive allocation
    if(allocated(tMap))deallocate(tMap)
    allocate(tMap(tNso,1,1))
    !Creating the sequence of operators A*_q, B*_q
    ! which decompose the term H^LR of the
    ! super-block Hamiltonian.
    it = 0
    do i=1,tNso
       it = it+1
       tMap(i,1,1)=it
    enddo
    !
    !
    allocate(RowOffset(tNso,Nsb))
    allocate(ColOffset(tNso,Nsb))
    RowOffset=0
    ColOffset=0
    !
    !
    if(allocated(SBleft_states))deallocate(SBleft_states)
    if(allocated(SBright_states))deallocate(SBright_states)
    if(allocated(A))deallocate(A)
    if(allocated(B))deallocate(B)
    if(allocated(Hleft))deallocate(Hleft)
    if(allocated(Hright))deallocate(Hright)
    if(allocated(isb2jsb))deallocate(isb2jsb)
    if(allocated(IsHconjg))deallocate(IsHconjg)
    if(allocated(Sl_n))deallocate(Sl_n)
    if(allocated(Sl_p))deallocate(Sl_p)
    if(allocated(Sr_n))deallocate(Sr_n)
    if(allocated(Sr_p))deallocate(Sr_p)
    !    
    allocate(SBleft_states(Nsb),SBright_states(Nsb))
    if(.not.direct_H_lazy)then
       allocate(A(tNso,Nsb),B(tNso,Nsb))
       allocate(Hleft(Nsb),Hright(Nsb))
       allocate(isb2jsb(tNso,Nsb));isb2jsb=0
       allocate(IsHconjg(tNso,Nsb));IsHconjg=0
    endif
    allocate(Sl_n(Nspin),Sr_n(Nspin))
    allocate(Sl_p(Nspin),Sr_p(Nspin))
    !
    if(MpiMaster)t0=t_start()
    !Main computation:
    do isb=1,Nsb 
       qn             = sb_sector%qn(index=isb)
       SBleft_states(isb)%states = sb2block_states(qn,'left')
       SBright_states(isb)%states = sb2block_states(qn,'right')
    enddo
    call setup_sector_filter_maps()
    if(MpiMaster)print*,"Get Filtered States:",t_stop()
    !

    !
    ! ROOT get basic operators from L/R blocks and bcast them
    if(MpiMaster)t0=t_start()
    if(MpiMaster)then
       do ispin=1,Nspin
          Sl_n(ispin) = left%operators%op("S"//left%okey(0,ispin,ilink='n'))
          Sr_n(ispin) = right%operators%op("S"//right%okey(0,ispin,ilink='n'))
          if(PBCdmrg)then
             Sl_p(ispin) = left%operators%op("S"//left%okey(0,ispin,ilink='p'))
             Sr_p(ispin) = right%operators%op("S"//right%okey(0,ispin,ilink='p'))
          endif
       enddo
       Hl = left%operators%op("H")
       Hr = right%operators%op("H")
    endif
#ifdef _MPI
    if(MpiStatus)then
       do ispin=1,Nspin
          call Sl_n(ispin)%bcast()
          call Sr_n(ispin)%bcast()
          if(PBCdmrg)then
             call Sl_p(ispin)%bcast()
             call Sr_p(ispin)%bcast()
          endif
       enddo
       call Hl%bcast()
       call Hr%bcast()
    endif
#endif
    if(MpiMaster)print*,"Build Operators:",t_stop()
    !
    !
    if(direct_H_lazy)then
       call setup_lazy_spin_operators(Sl_n,Sl_p,Sr_n,Sr_p,Hl,Hr)
    else
    !It is possible to MPI split the construction of the H_L,H_R,A,B ops
    if(MpiMaster)t0=t_start()
    isb2jsb=0
    do isb=1+MpiRank,Nsb,MpiSize
       !
       sizeA=size(SBleft_states(isb)%states)
       sizeB=size(SBright_states(isb)%states)
       if(MpiMaster.AND.sizeA>10)&
            write(LOGfile,*)"isb:"//str(isb)//"/"//str(Nsb)//&
               " N(isb):"//str(sizeA)//","//str(sizeB)
       !
       qn = sb_sector%qn(index=isb)
       !
       Hleft(isb) = sp_filter(Hl,SBleft_states(isb)%states)
       Hright(isb)= sp_filter(Hr,SBright_states(isb)%states)
       !
       !get  it=1: A = Jp*S_lz .x. B = S_rz + Row/Col Offsets (DIAGONAL)
       it = tMap(1,1,1)
       A(it,isb) = Hij(1,1)*sp_filter(Sl_n(1),SBleft_states(isb)%states,SBleft_states(isb)%states)
       B(it,isb) = sp_filter(Sr_n(1),SBright_states(isb)%states,SBright_states(isb)%states)
       qm  = qn
       jsb = isb
       Isb2Jsb(it,isb) =jsb
       IsHconjg(it,isb)=0!.false.
       RowOffset(it,isb)=Offset(isb)           
       ColOffset(it,isb)=Offset(isb)
       !
       !> get it=2: A = Jxy*S_l- .x. B = [S_r-]^+=S_r+ + Row/Col Offsets
       !> get it=3: A = Jxy*S_l+ .x. B = [S_r+]^+=S_r- + Row/Col Offsets 
       dq = [1d0]
       qm = qn - dq
       if(.not.sb_sector%has_qn(qm))cycle
       jsb = sb_sector%index(qn=qm)
       !
       !it=2
       it=tMap(2,1,1)
       A(it,isb) = Hij(2,2)*sp_filter(Sl_n(2),SBleft_states(isb)%states,SBleft_states(jsb)%states)
       B(it,isb) = sp_filter(hconjg(Sr_n(2)),SBright_states(isb)%states,SBright_states(jsb)%states)
       Isb2Jsb(it,isb)  =jsb
       IsHconjg(it,isb) =0!.false.
       RowOffset(it,isb)=Offset(isb)
       ColOffset(it,isb)=Offset(jsb)
       !
       !it=3
       it=tMap(3,1,1)
       A(it,isb) = hconjg(A(tMap(2,1,1),isb))
       B(it,isb) = hconjg(B(tMap(2,1,1),isb))
       Isb2Jsb(it,isb)  =jsb
       IsHconjg(it,isb) =1!.true.  !exchange jsb and isb
       RowOffset(it,isb)=Offset(jsb)
       ColOffset(it,isb)=Offset(isb)
       !
       if(PBCdmrg)then
          !get it=4: A = Jp*S_lz .x. B = S_rz + Row/Col Offsets (DIAGONAL)
          it = tMap(4,1,1)
          A(it,isb) = Hij(1,1)*sp_filter(Sl_p(1),SBleft_states(isb)%states,SBleft_states(isb)%states)
          B(it,isb) = sp_filter(Sr_p(1),SBright_states(isb)%states,SBright_states(isb)%states)
          qm  = qn
          jsb = isb
          Isb2Jsb(it,isb) =jsb
          IsHconjg(it,isb)=0!.false.
          RowOffset(it,isb)=Offset(isb)           
          ColOffset(it,isb)=Offset(isb)
          !
          !> get it=2: A = Jxy*S_l- .x. B = [S_r-]^+=S_r+ + Row/Col Offsets
          !> get it=3: A = Jxy*S_l+ .x. B = [S_r+]^+=S_r- + Row/Col Offsets 
          dq = [1d0]
          qm = qn - dq
          if(.not.sb_sector%has_qn(qm))cycle
          jsb = sb_sector%index(qn=qm)
          !
          it=tMap(5,1,1)
          A(it,isb) = Hij(2,2)*sp_filter(Sl_p(2),SBleft_states(isb)%states,SBleft_states(jsb)%states)
          B(it,isb) = sp_filter(hconjg(Sr_p(2)),SBright_states(isb)%states,SBright_states(jsb)%states)
          Isb2Jsb(it,isb)  =jsb
          IsHconjg(it,isb) =0!.false.
          RowOffset(it,isb)=Offset(isb)
          ColOffset(it,isb)=Offset(jsb)
          !
          it=tMap(6,1,1)
          A(it,isb) = hconjg(A(tMap(5,1,1),isb))
          B(it,isb) = hconjg(B(tMap(5,1,1),isb))
          Isb2Jsb(it,isb)  =jsb
          IsHconjg(it,isb) =1!.true.  !exchange jsb and isb
          RowOffset(it,isb)=Offset(jsb)
          ColOffset(it,isb)=Offset(isb)
       endif
    enddo
    if(MpiMaster)print*,"Get Op Blocks:",t_stop()
    !
#ifdef _MPI
    !This can possibly be improved workin on AllGather_MPI so to use
    !pack/unpack of the sparse matrix arrays.
    if(MpiStatus)then
       if(MpiMaster)t0=t_start()
       call AllGather_MPI(MpiComm,Hleft)
       call AllGather_MPI(MpiComm,Hright)
       call AllGather_MPI(MpiComm,A)
       call AllGather_MPI(MpiComm,B)
       call MPI_ALLREDUCE(Mpi_In_Place,Isb2Jsb,size(Isb2Jsb),&
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       call MPI_ALLREDUCE(Mpi_In_Place,RowOffset,size(RowOffset),&
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       call MPI_ALLREDUCE(Mpi_In_Place,ColOffset,size(ColOffset),&
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       call MPI_ALLREDUCE(Mpi_In_Place,IsHconjg,size(IsHconjg),&
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       !Root accumulate the exchanged data:
       do isb=1,Nsb
          kb_sb_setup_bcast = kb_sb_setup_bcast + Hleft(isb)%bytes() + Hright(isb)%bytes()
          do it=1,tNso
             kb_sb_setup_bcast = kb_sb_setup_bcast + A(it,isb)%bytes()  + B(it,isb)%bytes()
          enddo
       enddo
       if(MpiMaster)print*,"MpiComm Op Blocks:",t_stop()
    endif
#endif
    endif
    !
    !
    do ispin=1,Nspin
       call Sl_n(ispin)%free()
       call Sl_p(ispin)%free()
       call Sr_n(ispin)%free()
       call Sr_p(ispin)%free()
    enddo
    call Hl%free()
    call Hr%free()
    deallocate(Sl_n,Sr_n,Sl_p,Sr_p)
    !
    if(MpiMaster)call stop_timer("Setup SB Direct")
  end subroutine Setup_SuperBlock_Spin_Direct







  !##################################################################
  !                        FERMION CASE
  !##################################################################
  subroutine Setup_SuperBlock_Fermion_Direct()
    integer                                      :: Nso,Nsb
    integer                                      :: it,isb,jsb,ierr,ipr,fbc,sizeA,sizeB
    real(8),dimension(:),allocatable             :: qn,qm
    type(tstates),dimension(:),allocatable       :: Aj,Bj
    real(8),dimension(:),allocatable             :: dq
    real(8),dimension(:),allocatable             :: qnup,qndw
    integer,dimension(:,:,:),allocatable         :: tMap
    type(sparse_matrix)                          :: Hl,Hr
    type(sparse_matrix)                          :: P_n,P_p
    type(sparse_matrix),allocatable,dimension(:) :: Cl_n,Cr_n
    type(sparse_matrix),allocatable,dimension(:) :: Cl_p,Cr_p
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable        :: Hij
#else
    real(8),dimension(:,:),allocatable           :: Hij
#endif
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: Setup SB Direct - fermion"
#endif
    !
    if(MpiMaster)call start_timer("Setup SB Direct, Nsb: "//str(size(sb_sector)))
    !
    if(.not.left%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing left.H operator in the list"
    if(.not.right%operators%has_key("H"))&
         stop "Setup_SuperBlock_Direct ERROR: Missing right.H operator in the list"
    !
    !
    !> GET THE USER DEFINED MODEL HAMILTONIAN PARAMETERS:
    ! Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)
    !
    !
    fbc  = 2
    if(PBCdmrg)fbc=4
    Nso  = Nspin*Norb
    tNso = fbc*count(Hij/=zero)
    Nsb  = size(sb_sector)
    !
    !
    !Massive allocation
    if(allocated(tMap))deallocate(tMap)
    allocate(tMap(fbc,Nso,Nso))
    !Creating the sequence of operators A*_q, B*_q
    ! which decompose the term H^LR of the
    ! super-block Hamiltonian.
    it = 0
    do i=1,fbc
       do io=1,Nso
          do jo=1,Nso
             if(Hij(io,jo)==zero)cycle
             it = it+1
             tMap(i,io,jo)=it
          enddo
       enddo
    enddo
    !
    !
    allocate(RowOffset(tNso,Nsb))
    allocate(ColOffset(tNso,Nsb))
    RowOffset=0
    ColOffset=0
    !
    !
    if(allocated(SBleft_states))deallocate(SBleft_states)
    if(allocated(SBright_states))deallocate(SBright_states)
    if(allocated(A))deallocate(A)
    if(allocated(B))deallocate(B)
    if(allocated(Hleft))deallocate(Hleft)
    if(allocated(Hright))deallocate(Hright)
    if(allocated(isb2jsb))deallocate(isb2jsb)
    if(allocated(IsHconjg))deallocate(IsHconjg)
    if(allocated(Cl_n))deallocate(Cl_n)
    if(allocated(Cl_p))deallocate(Cl_p)
    if(allocated(Cr_n))deallocate(Cr_n)
    if(allocated(Cr_p))deallocate(Cr_p)
    !
    allocate(SBleft_states(Nsb),SBright_states(Nsb))
    if(.not.direct_H_lazy)then
       allocate(A(tNso,Nsb),B(tNso,Nsb))
       allocate(Hleft(Nsb),Hright(Nsb))
       allocate(isb2jsb(tNso,Nsb));isb2jsb=0
       allocate(IsHconjg(tNso,Nsb));IsHconjg=0
    endif
    allocate(Cl_n(Nso),Cr_n(Nso))
    allocate(Cl_p(Nso),Cr_p(Nso))
    !
    !
    if(MpiMaster)t0=t_start()
    !All nodes filter QN states:
    do isb=1,Nsb
       qn             = sb_sector%qn(index=isb)
       SBleft_states(isb)%states = sb2block_states(qn,'left')
       SBright_states(isb)%states = sb2block_states(qn,'right')
    enddo
    call setup_sector_filter_maps()
    if(MpiMaster)write(LOGfile,*)"Get Filtered States:",t_stop()
    !
    !
    ! ROOT get basic operators from L/R blocks and bcast them
    if(MpiMaster)t0=t_start()
    if(MpiMaster)then
       P_n = left%operators%op("P"//left%okey(0,0,ilink='n'))
       if(PBCdmrg) P_p = left%operators%op("P"//left%okey(0,0,ilink='p'))
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb+(ispin-1)*Norb
             Cl_n(io) = left%operators%op("C"//left%okey(iorb,ispin,ilink='n'))
             Cr_n(io) = right%operators%op("C"//right%okey(iorb,ispin,ilink='n'))
             if(PBCdmrg)then
                Cl_p(io) = left%operators%op("C"//left%okey(iorb,ispin,ilink='p'))
                Cr_p(io) = right%operators%op("C"//right%okey(iorb,ispin,ilink='p'))
             endif
          enddo
       enddo
       Hl = left%operators%op("H")
       Hr = right%operators%op("H")
    endif
#ifdef _MPI
    if(MpiStatus)then
       call P_n%bcast()
       if(PBCdmrg)call P_p%bcast()
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb+(ispin-1)*Norb
             call Cl_n(io)%bcast()
             call Cr_n(io)%bcast()
             if(PBCdmrg)then
                call Cl_p(io)%bcast()
                call Cr_p(io)%bcast()
             endif
          enddo
       enddo
       call Hl%bcast()
       call Hr%bcast()
    endif
#endif
    if(MpiMaster)write(LOGfile,*)"Build Operators:",t_stop()
    !
    !Setup dQ for each mode: dQ is Q after application of B=C(b,s)
    !the Hermitian cobjugate is obtained by symmetry
    allocate(qnup, mold=current_target_qn)
    allocate(qndw, mold=current_target_qn)
    select case(dmrg_mode)
    case default
       qnup = [1d0,0d0]         !Nup->Nup-1 (destroy electron spin-UP)
       qndw = [0d0,1d0]         !Ndw->Ndw-1 (destroy electron spin-DW)
    case("superc")
       qnup = [ 1d0]            !Sz->Sz-1 (destroy electron spin-UP)
       qndw = [-1d0]            !Sz->Sz+1 (destroy electron spin-DW)
    case("nonsu2")
       qnup = [ 1d0]            !N->N-1 (destroy electron whatever spin)
       qndw = [ 1d0]            !N->N-1 (destroy electron whatever spin)
    end select
    !
    if(direct_H_lazy)then
       call setup_lazy_fermion_operators(Cl_n,Cl_p,Cr_n,Cr_p,P_n,P_p,Hl,Hr)
    else
    !It is possible to MPI split the construction of the H_L,H_R,A,B ops
    if(MpiMaster)t0=t_start()
    isb2jsb=0
    do isb=1+MpiRank,Nsb,MpiSize
       !
       sizeA=size(SBleft_states(isb)%states)
       sizeB=size(SBright_states(isb)%states)
       if(MpiMaster.AND.sizeA>10)&
            write(LOGfile,*)"isb:"//str(isb)//"/"//str(Nsb)//&
            " N(isb):"//str(sizeA)//","//str(sizeB)
       !
       qn = sb_sector%qn(index=isb)
       !
       Hleft(isb) = sp_filter(Hl,SBleft_states(isb)%states)
       Hright(isb)= sp_filter(Hr,SBright_states(isb)%states)
       !
       !> get it=1,io,jo: A = H(a,b)*[Cl(a,s)^+@P] .x. B = Cr(b,s)  + Row/Col Offsets
       !> get it=2,io,jo: A = H(a,b)*[Cl(a,s)@P] .x. B = Cr(b,s)^+  + Row/Col Offsets
       do ispin=1,Nspin
          dq = qnup ; if(ispin==2)dq=qndw
          qm = qn - dq
          if(.not.sb_sector%has_qn(qm))cycle
          jsb = sb_sector%index(qn=qm)
          !
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb+(ispin-1)*Norb
                jo = jorb+(ispin-1)*Norb
                if(Hij(io,jo)==0d0)cycle
                !
                !it=1,a,b
                it=tMap(1,io,jo)
                A(it,isb) = Hij(io,jo)*sp_filter(matmul(Cl_n(io)%dgr(),P_n),SBleft_states(isb)%states,SBleft_states(jsb)%states)
                B(it,isb) = sp_filter(Cr_n(jo),SBright_states(isb)%states,SBright_states(jsb)%states)
                Isb2Jsb(it,isb)  =jsb
                IsHconjg(it,isb) =0!.false.
                RowOffset(it,isb)=Offset(isb)           
                ColOffset(it,isb)=Offset(jsb)
                !
                !it=2,a,b
                it=tMap(2,io,jo)
                A(it,isb) = hconjg(A(tMap(1,io,jo),isb))
                B(it,isb) = hconjg(B(tMap(1,io,jo),isb))
                Isb2Jsb(it,isb)  =jsb
                IsHconjg(it,isb) =1!.true.
                RowOffset(it,isb)=Offset(jsb)
                ColOffset(it,isb)=Offset(isb)
                !
                if(PBCdmrg)then
                   !it=3,a,b
                   it=tMap(3,io,jo)
                   A(it,isb) = Hij(io,jo)*sp_filter(matmul(Cl_p(io)%dgr(),P_p),SBleft_states(isb)%states,SBleft_states(jsb)%states)
                   B(it,isb) = sp_filter(Cr_p(jo),SBright_states(isb)%states,SBright_states(jsb)%states)
                   Isb2Jsb(it,isb)  =jsb
                   IsHconjg(it,isb) =0!.false.
                   RowOffset(it,isb)=Offset(isb)           
                   ColOffset(it,isb)=Offset(jsb)
                   !
                   !it=2,a,b
                   it=tMap(4,io,jo)
                   A(it,isb) = hconjg(A(tMap(3,io,jo),isb))
                   B(it,isb) = hconjg(B(tMap(3,io,jo),isb))
                   Isb2Jsb(it,isb)  =jsb
                   IsHconjg(it,isb) =1!.true.
                   RowOffset(it,isb)=Offset(jsb)
                   ColOffset(it,isb)=Offset(isb)
                end if
                !
             enddo
          enddo
          !
       enddo
    enddo
    if(MpiMaster)write(LOGfile,*)"Get Op Blocks:",t_stop()
    !
#ifdef _MPI
    if(MpiStatus)then
       if(MpiMaster)t0=t_start()
       call AllGather_MPI(MpiComm,Hleft)
       call AllGather_MPI(MpiComm,Hright)
       call AllGather_MPI(MpiComm,A)
       call AllGather_MPI(MpiComm,B)
       call MPI_ALLREDUCE(Mpi_In_Place,Isb2Jsb,size(Isb2Jsb),&
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       call MPI_ALLREDUCE(Mpi_In_Place,RowOffset,size(RowOffset),&
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       call MPI_ALLREDUCE(Mpi_In_Place,ColOffset,size(ColOffset),&
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       call MPI_ALLREDUCE(Mpi_In_Place,IsHconjg,size(IsHconjg),&
            MPI_INTEGER, MPI_SUM, MpiComm, ierr)
       !Root accumulate the exchanged data:
       do isb=1,Nsb
          kb_sb_setup_bcast = kb_sb_setup_bcast + Hleft(isb)%bytes() + Hright(isb)%bytes()
          do it=1,tNso
             kb_sb_setup_bcast = kb_sb_setup_bcast + A(it,isb)%bytes()  + B(it,isb)%bytes()
          enddo
       enddo
       if(MpiMaster)print*,"MpiComm Op Blocks:",t_stop()
    endif
#endif
    endif
    !
    call P_n%free()
    call P_p%free()
    call Hl%free()
    call Hr%free()
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb+(ispin-1)*Norb
          call Cl_n(io)%free()
          call Cr_n(io)%free()
          call Cl_p(io)%free()
          call Cr_p(io)%free()
       enddo
    enddo
    call Hl%free()
    call Hr%free()
    deallocate(Cl_n,Cr_n)
    deallocate(Cl_p,Cr_p)
    !
    if(MpiMaster)call stop_timer("Setup SB Direct")
  end subroutine Setup_SuperBlock_Fermion_Direct









  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_sparse_main(Nloc,v,Hv)
    integer                    :: Nloc
    integer                    :: i,j,jcol
#ifdef _CMPLX
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    complex(8)                 :: val
#else
    real(8),dimension(Nloc)    :: v
    real(8),dimension(Nloc)    :: Hv
    real(8)                    :: val
#endif
    t0=t_start()
    Hv=zero
    do i=1,Nloc
       matmul: do jcol=1, spHsb%row(i)%Size
          val = spHsb%row(i)%vals(jcol)
          j   = spHsb%row(i)%cols(jcol)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
    t_hxv_sparse=t_hxv_sparse+t_stop()
  end subroutine spMatVec_sparse_main


  subroutine setup_lazy_spin_operators(Sl_n,Sl_p,Sr_n,Sr_p,Hl,Hr)
    type(sparse_matrix),dimension(:),intent(in) :: Sl_n,Sl_p,Sr_n,Sr_p
    type(sparse_matrix),intent(in)              :: Hl,Hr
    integer                                     :: is
    !
    call Lazy_Hl%free()
    call Lazy_Hr%free()
    Lazy_Hl = Hl
    Lazy_Hr = Hr
    if(allocated(Lazy_Sl_n))then
       do is=1,size(Lazy_Sl_n)
          call Lazy_Sl_n(is)%free()
          call Lazy_Sr_n(is)%free()
          call Lazy_Sl_p(is)%free()
          call Lazy_Sr_p(is)%free()
       enddo
       deallocate(Lazy_Sl_n,Lazy_Sr_n,Lazy_Sl_p,Lazy_Sr_p)
    endif
    allocate(Lazy_Sl_n(size(Sl_n)),Lazy_Sr_n(size(Sr_n)))
    allocate(Lazy_Sl_p(size(Sl_p)),Lazy_Sr_p(size(Sr_p)))
    do is=1,size(Sl_n)
       Lazy_Sl_n(is) = Sl_n(is)
       Lazy_Sr_n(is) = Sr_n(is)
       Lazy_Sl_p(is) = Sl_p(is)
       Lazy_Sr_p(is) = Sr_p(is)
    enddo
  end subroutine setup_lazy_spin_operators


  subroutine setup_sector_filter_maps()
    integer :: isb,istate
    !
    if(allocated(SBleft_maps))deallocate(SBleft_maps)
    if(allocated(SBright_maps))deallocate(SBright_maps)
    allocate(SBleft_maps(size(sb_sector)),SBright_maps(size(sb_sector)))
    do isb=1,size(sb_sector)
       allocate(SBleft_maps(isb)%states(left%Dim))
       allocate(SBright_maps(isb)%states(right%Dim))
       SBleft_maps(isb)%states=0
       SBright_maps(isb)%states=0
       do istate=1,size(SBleft_states(isb)%states)
          SBleft_maps(isb)%states(SBleft_states(isb)%states(istate)) = istate
       enddo
       do istate=1,size(SBright_states(isb)%states)
          SBright_maps(isb)%states(SBright_states(isb)%states(istate)) = istate
       enddo
    enddo
  end subroutine setup_sector_filter_maps


  subroutine setup_lazy_fermion_operators(Cl_n,Cl_p,Cr_n,Cr_p,P_n,P_p,Hl,Hr)
    type(sparse_matrix),dimension(:),intent(in) :: Cl_n,Cl_p,Cr_n,Cr_p
    type(sparse_matrix),intent(in)              :: P_n,P_p,Hl,Hr
    integer                                     :: io
    !
    call Lazy_Hl%free()
    call Lazy_Hr%free()
    call Lazy_Pn%free()
    call Lazy_Pp%free()
    Lazy_Hl = Hl
    Lazy_Hr = Hr
    Lazy_Pn = P_n
    Lazy_Pp = P_p
    if(allocated(Lazy_Cl_n))then
       do io=1,size(Lazy_Cl_n)
          call Lazy_Cl_n(io)%free()
          call Lazy_Cr_n(io)%free()
          call Lazy_Cl_p(io)%free()
          call Lazy_Cr_p(io)%free()
          call Lazy_CdgP_n(io)%free()
          call Lazy_CdgP_p(io)%free()
       enddo
       deallocate(Lazy_Cl_n,Lazy_Cr_n,Lazy_Cl_p,Lazy_Cr_p)
       deallocate(Lazy_CdgP_n,Lazy_CdgP_p)
    endif
    allocate(Lazy_Cl_n(size(Cl_n)),Lazy_Cr_n(size(Cr_n)))
    allocate(Lazy_Cl_p(size(Cl_p)),Lazy_Cr_p(size(Cr_p)))
    allocate(Lazy_CdgP_n(size(Cl_n)),Lazy_CdgP_p(size(Cl_p)))
    do io=1,size(Cl_n)
       Lazy_Cl_n(io) = Cl_n(io)
       Lazy_Cr_n(io) = Cr_n(io)
       Lazy_Cl_p(io) = Cl_p(io)
       Lazy_Cr_p(io) = Cr_p(io)
       Lazy_CdgP_n(io) = matmul(Cl_n(io)%dgr(),P_n)
       if(PBCdmrg)Lazy_CdgP_p(io) = matmul(Cl_p(io)%dgr(),P_p)
    enddo
  end subroutine setup_lazy_fermion_operators


  subroutine apply_AxB_direct(Aop,Bop,row_offset,col_offset,v,Hv)
    type(sparse_matrix),intent(in) :: Aop,Bop
    integer,intent(in)             :: row_offset,col_offset
#ifdef _CMPLX
    complex(8),dimension(:)        :: v,Hv
    complex(8),dimension(:,:),allocatable :: C
    complex(8)                     :: val
#else
    real(8),dimension(:)           :: v,Hv
    real(8),dimension(:,:),allocatable :: C
    real(8)                        :: val
#endif
    integer                        :: ai,aj,bi,bj,ja,jb,j,ic,i,jc
    !
    if(.not.Aop%status.OR..not.Bop%status)return
    allocate(C(Bop%Nrow,Aop%Ncol));C=zero
    do aj=1,Aop%Ncol
       do bi=1,Bop%Nrow
          if(Bop%row(bi)%Size==0)cycle
          do jb=1,Bop%row(bi)%Size
             bj   = Bop%row(bi)%cols(jb)
             val  = Bop%row(bi)%vals(jb)
             jc   = bj + (aj-1)*Bop%Ncol
             j    = jc + col_offset
             C(bi,aj) = C(bi,aj) + val*v(j)
          enddo
       enddo
    enddo
    do bi=1,Bop%Nrow
       do ai=1,Aop%Nrow
          if(Aop%row(ai)%Size==0)cycle
          ic = bi + (ai-1)*Bop%Nrow
          i  = ic + row_offset
          do ja=1,Aop%row(ai)%Size
             aj  = Aop%row(ai)%cols(ja)
             val = Aop%row(ai)%vals(ja)
             Hv(i) = Hv(i) + val*C(bi,aj)
          enddo
       enddo
    enddo
    deallocate(C)
  end subroutine apply_AxB_direct


#ifdef _MPI
  subroutine apply_AxB_MPI_direct(Aop,Bop,k,q,is_hconjg,v,Hv)
    type(sparse_matrix),intent(in) :: Aop,Bop
    integer,intent(in)             :: k,q,is_hconjg
#ifdef _CMPLX
    complex(8),dimension(:)        :: v,Hv
    complex(8),dimension(:),allocatable   :: vt,Hvt
    complex(8),dimension(:,:),allocatable :: C,Ct
    complex(8)                     :: val
#else
    real(8),dimension(:)           :: v,Hv
    real(8),dimension(:),allocatable      :: vt,Hvt
    real(8),dimension(:,:),allocatable    :: C,Ct
    real(8)                        :: val
#endif
    integer                        :: ai,aj,bi,bj,ja,jb,jc,j,i
    integer                        :: mpiArow,mpiAcol,mpiBrow
    integer                        :: shift,abcomm,i_start,i_end
    !
    if(.not.Aop%status.OR..not.Bop%status)return
    mpiAcol = mpiDls(q)
    if(is_hconjg==1)mpiAcol=mpiDls(k)
    mpiArow = mpiDls(k)
    if(is_hconjg==1)mpiArow=mpiDls(q)
    mpiBrow = mpiDrs(k)
    if(is_hconjg==1)mpiBrow=mpiDrs(q)
    shift = mpiOffset(q)
    if(is_hconjg==1)shift = mpiOffset(k)
    !
    allocate(C(Bop%Nrow,mpiAcol));C=zero
    do aj=1,mpiAcol
       do bi=1,Bop%Nrow
          if(Bop%row(bi)%Size==0)cycle
          do jb=1,Bop%row(bi)%Size
             bj   = Bop%row(bi)%cols(jb)
             val  = Bop%row(bi)%vals(jb)
             jc   = bj + (aj-1)*Bop%Ncol
             j    = jc + shift
             C(bi,aj) = C(bi,aj) + val*v(j)
          enddo
       enddo
    enddo
    if(mpiNactive(q)>mpiNactive(k))then
       abcomm = mpiSBCOMM(q)
    else
       abcomm = mpiSBCOMM(k)
    endif
    allocate(Ct(Aop%Ncol,mpiBrow));Ct=zero
    call vector_transpose_MPI(Bop%Nrow,mpiAcol,C,Aop%Ncol,mpiBrow,Ct,abcomm)
    allocate(vt(mpiArow*Bop%Nrow));vt=zero
    allocate(Hvt(Aop%Nrow*mpiBrow));Hvt=zero
    do bi=1,mpiBrow
       do ai=1,Aop%Nrow
          if(Aop%row(ai)%Size==0)cycle
          i = ai + (bi-1)*Aop%Nrow
          do ja=1,Aop%row(ai)%Size
             aj  = Aop%row(ai)%cols(ja)
             val = Aop%row(ai)%vals(ja)
             Hvt(i) = Hvt(i) + val*Ct(aj,bi)
          enddo
       enddo
    enddo
    abcomm = mpiSBCOMM(k)
    if(is_hconjg==1)abcomm = mpiSBCOMM(q)
    call vector_transpose_MPI(Aop%Nrow,mpiBrow,Hvt,Bop%Nrow,mpiArow,vt,abcomm)
    i_start = 1 + mpiOffset(k)
    if(is_hconjg==1)i_start = 1 + mpiOffset(q)
    i_end = Bop%Nrow*mpiArow + mpiOffset(k)
    if(is_hconjg==1)i_end = Bop%Nrow*mpiArow + mpiOffset(q)
    Hv(i_start:i_end) = Hv(i_start:i_end) + vt
    deallocate(C,Ct,Hvt,vt)
  end subroutine apply_AxB_MPI_direct
#endif





  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_direct_main(Nloc,v,Hv)
    integer                               :: Nloc
#ifdef _CMPLX
    complex(8),dimension(Nloc)            :: v
    complex(8),dimension(Nloc)            :: Hv
    complex(8)                            :: val
    complex(8)                            :: aval,bval
    complex(8),dimension(:,:),allocatable :: C,Ct
#else
    real(8),dimension(Nloc)               :: v
    real(8),dimension(Nloc)               :: Hv
    real(8)                               :: val
    real(8)                               :: aval,bval
    real(8),dimension(:,:),allocatable    :: C,Ct
#endif
    integer                               :: i,j,k,q,n
    integer                               :: ir,il,jr,jl,it
    integer                               :: ai,aj,bi,bj,jcol
    integer                               :: ia,ib,ic,ja,jb,jc
    !
    Hv=zero
    t0=t_start()
    !> loop over all the SB sectors:
    sector: do k=1,size(sb_sector)
       !
       !> apply the 1^L x H^r
       t0 = t_start()
       do il=1,Dls(k)           !Fix the column il: v_il 
          !
          do ir=1,Drs(k)        !H^r.v_il
             i = ir + (il-1)*Drs(k) + offset(k)
             do jcol=1,Hright(k)%row(ir)%Size
                val = Hright(k)%row(ir)%vals(jcol)
                jr  = Hright(k)%row(ir)%cols(jcol)
                j   = jr + (il-1)*Drs(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             end do
          enddo
          !
       enddo
       t_hxv_1LxHR=t_hxv_1LxHR + t_stop()
       !
       !> apply the H^L x 1^r
       t0 = t_start()
       do ir=1,Drs(k)           !Fix the row ir: v_ir
          !
          do il=1,Dls(k)        !H^l.v_ir
             i = ir + (il-1)*Drs(k) + offset(k)
             do jcol=1,Hleft(k)%row(il)%Size
                val = Hleft(k)%row(il)%vals(jcol)
                jl  = Hleft(k)%row(il)%cols(jcol)
                j   = ir + (jl-1)*Drs(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             end do
          enddo
          !
       enddo
       t_hxv_HLx1R=t_hxv_HLx1R + t_stop()
       ! !
       !> apply the term sum_k sum_it A_it(k).x.B_it(k)
       !Hv = (A.x.B)vec(V) --> (A.x.B).V  -> vec(B.V.A^T)
       !
       !  B.V.A^T : [B.Nrow,B.Ncol].[B.Ncol,A.Ncol].[A.Ncol,A.Nrow]
       !            [Dr(k),Dr(k')].[Dr(k'),Dl(k')].[Dl(k'),Dl(k)]
       !   C.A^T  : [B.Nrow,A.Ncol].[A.Ncol,A.Nrow]
       !(A.C^T)^T : [ [A.Nrow,A.Ncol].[A.Ncol,B.Nrow] ]^T
       !              [B.Nrow,A.Nrow] = vec(Hv)
       t0 = t_start()
       do it=1,tNso
          q = isb2jsb(it,k)
          if(.not.A(it,k)%status.OR..not.B(it,k)%status)cycle
          !
          allocate(C(B(it,k)%Nrow,A(it,k)%Ncol));C=zero
          !
          !1. evaluate MMP: C = B.vec(V)
          !   \sum_bcol B(bi,bj)V_q(bj,aj)=C(bi,aj)
          !   j = bj+(aj-1)B.Ncol + ColOffset_q
          !   \sum_bcol B(bi,bj)v_q(j)=C(bi,aj)
          t0=t_start()         
          do aj=1,A(it,k)%Ncol             !
             do bi=1,B(it,k)%Nrow
                if(B(it,k)%row(bi)%Size==0)cycle
                !
                do jb=1,B(it,k)%row(bi)%Size
                   bj   = B(it,k)%row(bi)%cols(jb)
                   val  = B(it,k)%row(bi)%vals(jb)
                   jc   = bj + (aj-1)*B(it,k)%Ncol
                   j    = jc + ColOffset(it,k)
                   C(bi,aj) = C(bi,aj) + val*v(j)
                enddo
                !
             enddo
          enddo
          t_hxv_B=t_hxv_B+t_stop()
          !
          !2. evaluate MMP: C.A^t
          !   \sum_aj C(bi,aj)A^t(aj,ai)
          !  =\sum_aj [A(ai,aj)C^t(aj,bi)]^T
          t0=t_start()
          do bi=1,B(it,k)%Nrow
             !
             do ai=1,A(it,k)%Nrow
                if(A(it,k)%row(ai)%Size==0)cycle
                ic =  bi + (ai-1)*B(it,k)%Nrow
                i  =  ic + RowOffset(it,k)
                !
                do ja=1,A(it,k)%row(ai)%Size
                   aj  = A(it,k)%row(ai)%cols(ja)
                   val = A(it,k)%row(ai)%vals(ja)
                   Hv(i) = Hv(i) + val*C(bi,aj)
                enddo
                !
             enddo
          enddo
          t_hxv_B=t_hxv_B+t_stop()
          !
          deallocate(C)
       enddo
       !
       t_hxv_AxB=t_hxv_AxB+t_stop()
       !
    enddo sector
    !
    t_hxv_direct=t_hxv_direct + t_stop()
    !
  end subroutine spMatVec_direct_main


  subroutine spMatVec_direct_lazy_main(Nloc,v,Hv)
    integer                               :: Nloc
#ifdef _CMPLX
    complex(8),dimension(Nloc)            :: v
    complex(8),dimension(Nloc)            :: Hv
    complex(8)                            :: val
#else
    real(8),dimension(Nloc)               :: v
    real(8),dimension(Nloc)               :: Hv
    real(8)                               :: val
#endif
    type(sparse_matrix)                   :: Hlk,Hrk,Aop,Bop,tmpOp
    real(8),dimension(:),allocatable      :: qn,qm,dq,qnup,qndw
    character(len=:),allocatable          :: type
    integer                               :: k,q,ir,il,jr,jl,jcol,i,j
    integer                               :: ispin,iorb,jorb,io,jo
    !
    Hv=zero
    t0=t_start()
    type=str(left%type())
    sector: do k=1,size(sb_sector)
       Hrk = sp_filter(Lazy_Hr,SBright_states(k)%states,SBright_maps(k)%states,size(SBright_states(k)%states))
       do il=1,Dls(k)
          do ir=1,Drs(k)
             i = ir + (il-1)*Drs(k) + offset(k)
             do jcol=1,Hrk%row(ir)%Size
                val = Hrk%row(ir)%vals(jcol)
                jr  = Hrk%row(ir)%cols(jcol)
                j   = jr + (il-1)*Drs(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             enddo
          enddo
       enddo
       call Hrk%free()
       Hlk = sp_filter(Lazy_Hl,SBleft_states(k)%states,SBleft_maps(k)%states,size(SBleft_states(k)%states))
       do ir=1,Drs(k)
          do il=1,Dls(k)
             i = ir + (il-1)*Drs(k) + offset(k)
             do jcol=1,Hlk%row(il)%Size
                val = Hlk%row(il)%vals(jcol)
                jl  = Hlk%row(il)%cols(jcol)
                j   = ir + (jl-1)*Drs(k) + offset(k)
                Hv(i) = Hv(i) + val*v(j)
             enddo
          enddo
       enddo
       call Hlk%free()
       qn = sb_sector%qn(index=k)
       select case(to_lower(type(1:1)))
       case("s")
          Aop = HopH(1,1)*sp_filter(Lazy_Sl_n(1),SBleft_states(k)%states,SBleft_maps(k)%states,size(SBleft_states(k)%states))
          Bop = sp_filter(Lazy_Sr_n(1),SBright_states(k)%states,SBright_maps(k)%states,size(SBright_states(k)%states))
          call apply_AxB_direct(Aop,Bop,Offset(k),Offset(k),v,Hv)
          call Aop%free();call Bop%free()
          dq = [1d0]
          qm = qn - dq
          if(sb_sector%has_qn(qm))then
             q = sb_sector%index(qn=qm)
             Aop = HopH(2,2)*sp_filter(Lazy_Sl_n(2),SBleft_states(k)%states,SBleft_maps(q)%states,size(SBleft_states(q)%states))
             Bop = sp_filter(hconjg(Lazy_Sr_n(2)),SBright_states(k)%states,SBright_maps(q)%states,size(SBright_states(q)%states))
             call apply_AxB_direct(Aop,Bop,Offset(k),Offset(q),v,Hv)
             tmpOp = hconjg(Aop); call Aop%free(); Aop = tmpOp; call tmpOp%free()
             tmpOp = hconjg(Bop); call Bop%free(); Bop = tmpOp; call tmpOp%free()
             call apply_AxB_direct(Aop,Bop,Offset(q),Offset(k),v,Hv)
             call Aop%free();call Bop%free()
          endif
          if(PBCdmrg)then
             Aop = HopH(1,1)*sp_filter(Lazy_Sl_p(1),SBleft_states(k)%states,SBleft_maps(k)%states,size(SBleft_states(k)%states))
             Bop = sp_filter(Lazy_Sr_p(1),SBright_states(k)%states,SBright_maps(k)%states,size(SBright_states(k)%states))
             call apply_AxB_direct(Aop,Bop,Offset(k),Offset(k),v,Hv)
             call Aop%free();call Bop%free()
             qm = qn - [1d0]
             if(sb_sector%has_qn(qm))then
                q = sb_sector%index(qn=qm)
                Aop = HopH(2,2)*sp_filter(Lazy_Sl_p(2),SBleft_states(k)%states,SBleft_maps(q)%states,size(SBleft_states(q)%states))
                Bop = sp_filter(hconjg(Lazy_Sr_p(2)),SBright_states(k)%states,SBright_maps(q)%states,size(SBright_states(q)%states))
                call apply_AxB_direct(Aop,Bop,Offset(k),Offset(q),v,Hv)
                tmpOp = hconjg(Aop); call Aop%free(); Aop = tmpOp; call tmpOp%free()
                tmpOp = hconjg(Bop); call Bop%free(); Bop = tmpOp; call tmpOp%free()
                call apply_AxB_direct(Aop,Bop,Offset(q),Offset(k),v,Hv)
                call Aop%free();call Bop%free()
             endif
          endif
       case("f","e")
          allocate(qnup, mold=current_target_qn)
          allocate(qndw, mold=current_target_qn)
          select case(dmrg_mode)
          case default
             qnup = [1d0,0d0]
             qndw = [0d0,1d0]
          case("superc")
             qnup = [ 1d0]
             qndw = [-1d0]
          case("nonsu2")
             qnup = [1d0]
             qndw = [1d0]
          end select
          do ispin=1,Nspin
             dq = qnup ; if(ispin==2)dq=qndw
             qm = qn - dq
             if(.not.sb_sector%has_qn(qm))cycle
             q = sb_sector%index(qn=qm)
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb+(ispin-1)*Norb
                   jo = jorb+(ispin-1)*Norb
                   if(HopH(io,jo)==zero)cycle
                   Aop = HopH(io,jo)*sp_filter(Lazy_CdgP_n(io),SBleft_states(k)%states,SBleft_maps(q)%states,size(SBleft_states(q)%states))
                   Bop = sp_filter(Lazy_Cr_n(jo),SBright_states(k)%states,SBright_maps(q)%states,size(SBright_states(q)%states))
                   call apply_AxB_direct(Aop,Bop,Offset(k),Offset(q),v,Hv)
                   tmpOp = hconjg(Aop); call Aop%free(); Aop = tmpOp; call tmpOp%free()
                   tmpOp = hconjg(Bop); call Bop%free(); Bop = tmpOp; call tmpOp%free()
                   call apply_AxB_direct(Aop,Bop,Offset(q),Offset(k),v,Hv)
                   call Aop%free();call Bop%free()
                   if(PBCdmrg)then
                      Aop = HopH(io,jo)*sp_filter(Lazy_CdgP_p(io),SBleft_states(k)%states,SBleft_maps(q)%states,size(SBleft_states(q)%states))
                      Bop = sp_filter(Lazy_Cr_p(jo),SBright_states(k)%states,SBright_maps(q)%states,size(SBright_states(q)%states))
                      call apply_AxB_direct(Aop,Bop,Offset(k),Offset(q),v,Hv)
                      tmpOp = hconjg(Aop); call Aop%free(); Aop = tmpOp; call tmpOp%free()
                      tmpOp = hconjg(Bop); call Bop%free(); Bop = tmpOp; call tmpOp%free()
                      call apply_AxB_direct(Aop,Bop,Offset(q),Offset(k),v,Hv)
                      call Aop%free();call Bop%free()
                   endif
                enddo
             enddo
          enddo
          deallocate(qnup,qndw)
       end select
    enddo sector
    t_hxv_direct=t_hxv_direct + t_stop()
  end subroutine spMatVec_direct_lazy_main


#ifdef _MPI
  subroutine spMatVec_MPI_direct_lazy_main(Nloc,v,Hv)
    integer                               :: Nloc
#ifdef _CMPLX
    complex(8),dimension(Nloc)            :: v
    complex(8),dimension(Nloc)            :: Hv
    complex(8),dimension(:),allocatable   :: vt,Hvt
    complex(8)                            :: val
#else
    real(8),dimension(Nloc)               :: v
    real(8),dimension(Nloc)               :: Hv
    real(8),dimension(:),allocatable      :: vt,Hvt
    real(8)                               :: val
#endif
    type(sparse_matrix)                   :: Hlk,Hrk,Aop,Bop,tmpOp
    real(8),dimension(:),allocatable      :: qn,qm,dq,qnup,qndw
    character(len=:),allocatable          :: type
    integer                               :: k,q,ir,il,jr,jl,jcol,i,j
    integer                               :: ispin,iorb,jorb,io,jo
    integer                               :: i_start,i_end
    !
    if(.not.MpiStatus)stop "spMatVec_MPI_direct_lazy_main ERROR: MpiStatus = F"
    Hv=zero
    t0=t_start()
    type=str(left%type())
    sector: do k=1,size(sb_sector)
       Hrk = sp_filter(Lazy_Hr,SBright_states(k)%states,SBright_maps(k)%states,size(SBright_states(k)%states))
       do il=1,mpiDls(k)
          do ir=1,Drs(k)
             i = ir + (il-1)*Drs(k) + mpiOffset(k)
             do jcol=1,Hrk%row(ir)%Size
                val = Hrk%row(ir)%vals(jcol)
                jr  = Hrk%row(ir)%cols(jcol)
                j   = jr + (il-1)*Drs(k) + mpiOffset(k)
                Hv(i) = Hv(i) + val*v(j)
             enddo
          enddo
       enddo
       call Hrk%free()
       Hlk = sp_filter(Lazy_Hl,SBleft_states(k)%states,SBleft_maps(k)%states,size(SBleft_states(k)%states))
       allocate(vt(mpiDrs(k)*Dls(k)));vt=zero
       allocate(Hvt(mpiDrs(k)*Dls(k)));Hvt=zero
       i_start = 1 + mpiOffset(k)
       i_end   = Drs(k)*mpiDls(k)+mpiOffset(k)
       call vector_transpose_MPI(Drs(k),mpiDls(k),v(i_start:i_end),Dls(k),mpiDrs(k),vt,mpiSBCOMM(k))
       do il=1,mpiDrs(k)
          do ir=1,Dls(k)
             i = ir + (il-1)*Dls(k)
             do jcol=1,Hlk%row(ir)%Size
                val = Hlk%row(ir)%vals(jcol)
                jr  = Hlk%row(ir)%cols(jcol)
                j   = jr + (il-1)*Dls(k)
                Hvt(i) = Hvt(i) + val*vt(j)
             enddo
          enddo
       enddo
       deallocate(vt);allocate(vt(Drs(k)*mpiDls(k)));vt=zero
       call vector_transpose_MPI(Dls(k),mpiDrs(k),Hvt,Drs(k),mpiDls(k),vt,mpiSBCOMM(k))
       Hv(i_start:i_end) = Hv(i_start:i_end) + vt
       deallocate(vt,Hvt)
       call Hlk%free()
       qn = sb_sector%qn(index=k)
       select case(to_lower(type(1:1)))
       case("s")
          Aop = HopH(1,1)*sp_filter(Lazy_Sl_n(1),SBleft_states(k)%states,SBleft_maps(k)%states,size(SBleft_states(k)%states))
          Bop = sp_filter(Lazy_Sr_n(1),SBright_states(k)%states,SBright_maps(k)%states,size(SBright_states(k)%states))
          call apply_AxB_MPI_direct(Aop,Bop,k,k,0,v,Hv)
          call Aop%free();call Bop%free()
          dq = [1d0]
          qm = qn - dq
          if(sb_sector%has_qn(qm))then
             q = sb_sector%index(qn=qm)
             Aop = HopH(2,2)*sp_filter(Lazy_Sl_n(2),SBleft_states(k)%states,SBleft_maps(q)%states,size(SBleft_states(q)%states))
             Bop = sp_filter(hconjg(Lazy_Sr_n(2)),SBright_states(k)%states,SBright_maps(q)%states,size(SBright_states(q)%states))
             call apply_AxB_MPI_direct(Aop,Bop,k,q,0,v,Hv)
             tmpOp = hconjg(Aop); call Aop%free(); Aop = tmpOp; call tmpOp%free()
             tmpOp = hconjg(Bop); call Bop%free(); Bop = tmpOp; call tmpOp%free()
             call apply_AxB_MPI_direct(Aop,Bop,k,q,1,v,Hv)
             call Aop%free();call Bop%free()
          endif
          if(PBCdmrg)then
             Aop = HopH(1,1)*sp_filter(Lazy_Sl_p(1),SBleft_states(k)%states,SBleft_maps(k)%states,size(SBleft_states(k)%states))
             Bop = sp_filter(Lazy_Sr_p(1),SBright_states(k)%states,SBright_maps(k)%states,size(SBright_states(k)%states))
             call apply_AxB_MPI_direct(Aop,Bop,k,k,0,v,Hv)
             call Aop%free();call Bop%free()
             qm = qn - [1d0]
             if(sb_sector%has_qn(qm))then
                q = sb_sector%index(qn=qm)
                Aop = HopH(2,2)*sp_filter(Lazy_Sl_p(2),SBleft_states(k)%states,SBleft_maps(q)%states,size(SBleft_states(q)%states))
                Bop = sp_filter(hconjg(Lazy_Sr_p(2)),SBright_states(k)%states,SBright_maps(q)%states,size(SBright_states(q)%states))
                call apply_AxB_MPI_direct(Aop,Bop,k,q,0,v,Hv)
                tmpOp = hconjg(Aop); call Aop%free(); Aop = tmpOp; call tmpOp%free()
                tmpOp = hconjg(Bop); call Bop%free(); Bop = tmpOp; call tmpOp%free()
                call apply_AxB_MPI_direct(Aop,Bop,k,q,1,v,Hv)
                call Aop%free();call Bop%free()
             endif
          endif
       case("f","e")
          allocate(qnup, mold=current_target_qn)
          allocate(qndw, mold=current_target_qn)
          select case(dmrg_mode)
          case default
             qnup = [1d0,0d0]
             qndw = [0d0,1d0]
          case("superc")
             qnup = [ 1d0]
             qndw = [-1d0]
          case("nonsu2")
             qnup = [1d0]
             qndw = [1d0]
          end select
          do ispin=1,Nspin
             dq = qnup ; if(ispin==2)dq=qndw
             qm = qn - dq
             if(.not.sb_sector%has_qn(qm))cycle
             q = sb_sector%index(qn=qm)
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb+(ispin-1)*Norb
                   jo = jorb+(ispin-1)*Norb
                   if(HopH(io,jo)==zero)cycle
                   Aop = HopH(io,jo)*sp_filter(Lazy_CdgP_n(io),SBleft_states(k)%states,SBleft_maps(q)%states,size(SBleft_states(q)%states))
                   Bop = sp_filter(Lazy_Cr_n(jo),SBright_states(k)%states,SBright_maps(q)%states,size(SBright_states(q)%states))
                   call apply_AxB_MPI_direct(Aop,Bop,k,q,0,v,Hv)
                   tmpOp = hconjg(Aop); call Aop%free(); Aop = tmpOp; call tmpOp%free()
                   tmpOp = hconjg(Bop); call Bop%free(); Bop = tmpOp; call tmpOp%free()
                   call apply_AxB_MPI_direct(Aop,Bop,k,q,1,v,Hv)
                   call Aop%free();call Bop%free()
                   if(PBCdmrg)then
                      Aop = HopH(io,jo)*sp_filter(Lazy_CdgP_p(io),SBleft_states(k)%states,SBleft_maps(q)%states,size(SBleft_states(q)%states))
                      Bop = sp_filter(Lazy_Cr_p(jo),SBright_states(k)%states,SBright_maps(q)%states,size(SBright_states(q)%states))
                      call apply_AxB_MPI_direct(Aop,Bop,k,q,0,v,Hv)
                      tmpOp = hconjg(Aop); call Aop%free(); Aop = tmpOp; call tmpOp%free()
                      tmpOp = hconjg(Bop); call Bop%free(); Bop = tmpOp; call tmpOp%free()
                      call apply_AxB_MPI_direct(Aop,Bop,k,q,1,v,Hv)
                      call Aop%free();call Bop%free()
                   endif
                enddo
             enddo
          enddo
          deallocate(qnup,qndw)
       end select
    enddo sector
    t_hxv_direct=t_hxv_direct + t_stop()
  end subroutine spMatVec_MPI_direct_lazy_main
#endif


  







  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_MPI_direct_main(Nloc,v,Hv)
    integer                               :: Nloc
#ifdef _CMPLX
    complex(8),dimension(Nloc)            :: v
    complex(8),dimension(Nloc)            :: Hv
    complex(8),dimension(:),allocatable   :: vt,Hvt
    complex(8),dimension(:),allocatable   :: vin
    complex(8)                            :: val
    complex(8)                            :: aval,bval
    complex(8),dimension(:,:),allocatable :: C,Ct
#else
    real(8),dimension(Nloc)               :: v
    real(8),dimension(Nloc)               :: Hv
    real(8),dimension(:),allocatable      :: vt,Hvt
    real(8),dimension(:),allocatable      :: vin
    real(8)                               :: val
    real(8)                               :: aval,bval
    real(8),dimension(:,:),allocatable    :: C,Ct
#endif
    integer                               :: i,j,k,q,n,shift
    integer                               :: ir,il,jr,jl,it
    integer                               :: ai,aj,bi,bj,jcol
    integer                               :: ia,ib,ic,ja,jb,jc
    integer                               :: mpiArow,mpiAcol,mpiBrow,mpiBcol
    integer                               :: i_start,i_end, abcomm
    !
    if(.not.MpiStatus)stop "spMatVec_mpi_normal_main ERROR: MpiStatus = F"
    !
    !      
    Hv=zero
    t0=t_start()
    !> loop over all the SB sectors: k
    sector: do k=1,size(sb_sector)
      !
      ! if(MpiMaster)write(LOGfile,*)"SB sector:",k," qn:",sb_sector%qn(k)
      !
      !> apply the 1^L x H^r: share L columns
      ! if(MpiMaster)write(LOGfile,*)"Apply 1^L x H^R: share L columns"
      t0=t_start()
      do il=1,mpiDls(k)   !Fix the column il(q): v_il(q) for each thread
        !
        do ir=1,Drs(k)   !H^r.v_il
          i = ir + (il-1)*Drs(k) + mpiOffset(k)
          do jcol=1,Hright(k)%row(ir)%Size
            val = Hright(k)%row(ir)%vals(jcol)
            jr  = Hright(k)%row(ir)%cols(jcol)
            j   = jr + (il-1)*Drs(k) + mpiOffset(k)
            Hv(i) = Hv(i) + val*v(j)
          end do
        enddo
        !
      enddo
      t_hxv_1LxHR=t_hxv_1LxHR+t_stop()
      !       
      !> apply the H^L x 1^r
      !L part: non-contiguous in memory -> MPI transposition
      ! if(MpiMaster)write(LOGfile,*)"Apply H^L x 1^R: share R rows, MPI transpose L part"
      t0=t_start()
      allocate(vt(mpiDrs(k)*Dls(k))) ;vt=zero
      allocate(Hvt(mpiDrs(k)*Dls(k)));Hvt=zero
      i_start = 1 + mpiOffset(k)
      i_end   = Drs(k)*mpiDls(k)+mpiOffSet(k)
      call vector_transpose_MPI(Drs(k),mpiDls(k),v(i_start:i_end),Dls(k),mpiDrs(k),vt,mpiSBCOMM(k))
      do il=1,mpiDrs(k)  !Fix the *column* ir: v_ir(q). Transposed order
        do ir=1,Dls(k)  !go row-by-row H^l.v_ir: Transposed order
          i = ir + (il-1)*Dls(k)
          do jcol=1,Hleft(k)%row(ir)%Size
            val = Hleft(k)%row(ir)%vals(jcol)
            jr  = Hleft(k)%row(ir)%cols(jcol)
            j   = jr + (il-1)*Dls(k)
            Hvt(i) = Hvt(i) + val*vt(j)
          end do
        enddo
      enddo
      deallocate(vt) ; allocate(vt(Drs(k)*mpiDls(k))) ; vt=zero
      call vector_transpose_MPI(Dls(k),mpiDrs(k),Hvt,Drs(k),mpiDls(k),vt,mpiSBCOMM(k))
      Hv(i_start:i_end) = Hv(i_start:i_end) + Vt
      deallocate(vt,Hvt)
      t_hxv_HLx1R=t_hxv_HLx1R+t_stop()
      !
      !> apply the term sum_k sum_it A_it(k).x.B_it(k)
      !Hv = (A.x.B)vec(V) --> (A.x.B).V  -> vec(B.V.A^T)
      !
      !  B.V.A^T : [B.Nrow,B.Ncol].[B.Ncol,A.Ncol].[A.Ncol,A.Nrow]
      !            [Dr(k),Dr(k')].[Dr(k'),mpiDl(k')].[mpiDl(k'),Dl(k)]
      !   C.A^T  : [B.Nrow,A.Ncol].[A.Ncol,A.Nrow]
      !(A.C^T)^T : [ [A.Nrow,A.Ncol].[A.Ncol,B.Nrow] ]^T
      !              [B.Nrow,A.Nrow] = vec(Hv)
      !
      t0=t_start()
      do it=1,tNso
        if(.not.A(it,k)%status.OR..not.B(it,k)%status)cycle
      !   if(MpiMaster)write(LOGfile,*)"Apply A.x.B term it:",it," q:",isb2jsb(it,k)
        q = isb2jsb(it,k)
        !
        !1. evaluate MMP: C = B.vec(V)
        !   \sum_bcol B(bi,bj)V_q(bj,aj)=C(bi,aj)
        !   j = bj+(aj-1)B.Ncol + ColOffset_q
        !   \sum_bcol B(bi,bj)v_q(j)=C(bi,aj)
        !
        mpiAcol = mpiDls(q)
        if(isHconjg(it,k)==1)mpiAcol=mpiDls(k)
        !
        mpiArow = mpiDls(k)
        if(isHconjg(it,k)==1)mpiArow=mpiDls(q)
        !
        mpiBrow = mpiDrs(k)
        if(isHconjg(it,k)==1)mpiBrow=mpiDrs(q)
        !
        shift = mpiOffset(q)
        if(isHconjg(it,k)==1)shift = mpiOffset(k)
        !           
      !   if(MpiMaster)write(LOGfile,*)"MPI A.x.B: mpiArow,mpiAcol,mpiBrow:",mpiArow,mpiAcol,mpiBrow
        allocate(C(B(it,k)%Nrow,mpiAcol));C=zero
        t0=t_start()
        do aj=1,mpiAcol
          do bi=1,B(it,k)%Nrow
            if(B(it,k)%row(bi)%Size==0)cycle
            do jb=1,B(it,k)%row(bi)%Size
              bj   = B(it,k)%row(bi)%cols(jb)
              val  = B(it,k)%row(bi)%vals(jb)
              jc   = bj + (aj-1)*B(it,k)%Ncol
              j    = jc + shift
              C(bi,aj) = C(bi,aj) + val*v(j)
            enddo
          enddo
        enddo
        t_hxv_B=t_hxv_B+t_stop()
        !
        !Up to here we built "few", thread-related, columns of C(b,j*)
        !In the next step we will need to get [A.C^T]^T
        ! [C(b,j*)]^T = C(j*,b)
        ! [sum_j A_ij.[C^t]_jb]^T
        !
        !2. evaluate MMP: C.A^t
        !   \sum_aj C(bi,aj)A^t(aj,ai)
        !  =\sum_aj [A(ai,aj)C^t(aj,bi)]^T
        ! = [Hvt[A.Nrow,mpiBrow]]^T
        ! => Hv[B.Nrow,mpiArow]
        !use mpiSBCOMM(q) if mpiNactive(q)>mpiNactive(k) and mpiSBCOMM(k) otherwise
      !   if(MpiMaster)write(LOGfile,*)"MPI A.x.B: MPI transpose C"
        if(mpiNactive(q)>mpiNactive(k))then
           abcomm = mpiSBCOMM(q)
        else
           abcomm = mpiSBCOMM(k)
        endif
        allocate(Ct(A(it,k)%Ncol,mpiBrow));Ct=zero
        call vector_transpose_MPI(B(it,k)%Nrow,mpiAcol,C,A(it,k)%Ncol,mpiBrow,Ct,abcomm)
        !
        allocate(vt(mpiArow*B(it,k)%Nrow)) ; vt=zero
        allocate(Hvt(A(it,k)%Nrow*mpiBrow));Hvt=zero
        !
        t0=t_start()
        do bi=1,mpiBrow
          !
          do ai=1,A(it,k)%Nrow
            if(A(it,k)%row(ai)%Size==0)cycle
            i  =  ai + (bi-1)*A(it,k)%Nrow
            !
            do ja=1,A(it,k)%row(ai)%Size
              aj  = A(it,k)%row(ai)%cols(ja)
              val = A(it,k)%row(ai)%vals(ja)
              Hvt(i) = Hvt(i) + val*Ct(aj,bi)
            enddo
            !
          enddo
        enddo
        t_hxv_A=t_hxv_A+t_stop()
        !
      !   if(MpiMaster)write(LOGfile,*)"MPI A.x.B: MPI transpose A.C^T"
        abcomm = mpiSBCOMM(k)
        if(isHconjg(it,k)==1)abcomm = mpiSBCOMM(q) 
      !   if(MpiMaster)write(LOGfile,*)abcomm,MPI_COMM_NULL,MpiComm
        call vector_transpose_MPI(A(it,k)%Nrow,mpiBrow,Hvt,B(it,k)%Nrow,mpiArow,vt,abcomm)
        i_start = 1 + mpiOffset(k)
        if(isHconjg(it,k)==1)i_start = 1 + mpiOffset(q)
        i_end   = B(it,k)%Nrow*mpiArow+mpiOffSet(k)          
        if(isHconjg(it,k)==1)i_end   = B(it,k)%Nrow*mpiArow+mpiOffSet(q)
        !
        Hv(i_start:i_end) = Hv(i_start:i_end) + Vt
        !
        deallocate(C,Ct,Hvt,Vt)
      enddo
      t_hxv_AxB=t_hxv_AxB+t_stop()
      !
    enddo sector
    t_hxv_direct=t_hxv_direct + t_stop()
  end subroutine spMatVec_MPI_direct_main





  !##################################################################
  !              RETURN LEFT.states or RIGHT.states 
  !              contributing to the SUPERBLOCK with
  !              a given QN
  !. sb_sector: inherited from VARS_GLOBAL
  !##################################################################
  function sb2block_states(q,label) result(states)
    real(8),dimension(:)             :: q
    character(len=*)                 :: label
    integer,dimension(:),allocatable :: tmp,states,sb_map
    integer                          :: i,istate,l,r,isb,m,Rdim
    !
    if(.not.associated(sb_sector%root))&
         stop "sb2block_states error: sb_sector is not allocated"
    !
    if(allocated(states))deallocate(states)
    !
    !> get the map from the QN of the sector:
    sb_map = sb_sector%map(qn=q)
    !
    !> left,right, sb_sector and sb_states have to be known at this time:
    ! add a check
#ifdef _MPI
    if(MpiStatus)then
       if(MpiMaster)rDim=right%Dim
       call Bcast_MPI(MpiComm,rDim)
    else
       rDim=right%Dim
    endif
#else
    rDim=right%Dim
#endif
    allocate(tmp(size(sb_map)))
    select case(to_lower(str(label)))
    case("left","l","sys","s")
       do i=1,size(sb_map)
          istate = sb_states(sb_map(i))
          l = (istate-1)/rDim+1
          tmp(i) = l
       enddo
    case("right","r","env","e")
       do i=1,size(sb_map)
          istate = sb_states(sb_map(i))
          r = mod(istate,rDim);if(r==0)r=rDim
          tmp(i)=r             
       enddo
    end select
    allocate(states, source=uniq(tmp))
    deallocate(tmp)
  end function sb2block_states




END MODULE DMRG_SUPERBLOCK_SETUP













!   subroutine spMatVec_direct_main_(Nloc,v,Hv)
!     integer                    :: Nloc
! #ifdef _CMPLX
!     complex(8),dimension(Nloc) :: v
!     complex(8),dimension(Nloc) :: Hv
!     complex(8)                 :: val
!     complex(8)                 :: aval,bval
! #else
!     real(8),dimension(Nloc)    :: v
!     real(8),dimension(Nloc)    :: Hv
!     real(8)                    :: val
!     real(8)                    :: aval,bval
! #endif
!     integer                    :: i,j,k,n
!     integer                    :: ir,il,jr,jl,it
!     integer                    :: arow,brow,acol,bcol,jcol
!     integer                    :: ia,ib,ic,ja,jb,jc,sum
!     !
!     Hv=zero
!     !> loop over all the SB sectors:
!     sector: do k=1,size(sb_sector)
!        !
!        !> apply the 1^L x H^r
!        do il=1,Dls(k)           !Fix the column il: v_il 
!           !
!           do ir=1,Drs(k)        !H^r.v_il
!              i = ir + (il-1)*Drs(k) + offset(k)
!              do jcol=1,Hright(k)%row(ir)%Size
!                 val = Hright(k)%row(ir)%vals(jcol)
!                 jr  = Hright(k)%row(ir)%cols(jcol)
!                 j   = jr + (il-1)*Drs(k) + offset(k)
!                 Hv(i) = Hv(i) + val*v(j)
!              end do
!           enddo
!           !
!        enddo
!        !
!        !> apply the H^L x 1^r
!        do ir=1,Drs(k)           !Fix the row ir: v_ir
!           !
!           do il=1,Dls(k)        !H^l.v_ir
!              i = ir + (il-1)*Drs(k) + offset(k)
!              do jcol=1,Hleft(k)%row(il)%Size
!                 val = Hleft(k)%row(il)%vals(jcol)
!                 jl  = Hleft(k)%row(il)%cols(jcol)
!                 j   = ir + (jl-1)*Drs(k) + offset(k)
!                 Hv(i) = Hv(i) + val*v(j)
!              end do
!           enddo
!           !
!        enddo
!        !
!        !> apply the term sum_k sum_it A_it(k).x.B_it(k)
!        do it=1,tNso          
!           if(.not.A(it,k)%status.OR..not.B(it,k)%status)cycle
!           do ic=1,A(it,k)%Nrow*B(it,k)%Nrow
!              arow = (ic-1)/B(it,k)%Nrow+1
!              brow = mod(ic,B(it,k)%Nrow);if(brow==0)brow=B(it,k)%Nrow
!              if(A(it,k)%row(arow)%Size==0.OR.B(it,k)%row(brow)%Size==0)cycle
!              i = ic + RowOffset(it,k)
!              do ja=1,A(it,k)%row(arow)%Size
!                 acol = A(it,k)%row(arow)%cols(ja)
!                 aval = A(it,k)%row(arow)%vals(ja)
!                 do jb=1,B(it,k)%row(brow)%Size
!                    bcol = B(it,k)%row(brow)%cols(jb)
!                    bval = B(it,k)%row(brow)%vals(jb)
!                    j = bcol+(acol-1)*B(it,k)%Ncol + ColOffset(it,k)
!                    Hv(i) = Hv(i) + aval*bval*v(j)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo sector
!   end subroutine spMatVec_direct_main_
