MODULE DMRG_SUPERBLOCK_SETUP_DIRECT_SPIN
  USE DMRG_GLOBAL
  USE DMRG_SUPERBLOCK_SETUP_GLOBAL
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private


  public :: Setup_SuperBlock_Spin_Direct




contains


  !##################################################################
  !         SETUP THE SUPERBLOCK HAMILTONIAN for SPIN PROBLEMS
  !  the filtered operator blocks are stored in the memory as Lazy_*
  !  These are used either:
!In lazy mode the unfiltered operators remain in the Lazy_* pool.
!In cache mode they are only temporary sources for A/B/Hleft/Hright.
  !  - lazy mode: the unfiltered operators remain in memory and blocks 
  !    are extracted  on-the-fly at each H*v call.
  !  - cache mode: the filtered operators are tmp sources to A/B/Hleft/Hright 
  !    used in H*v operations
  !##################################################################
  subroutine Setup_SuperBlock_Spin_Direct()
    integer                                      :: Nso,Nsb
    integer                                      :: i,is,isb,it,jsb,ierr,sizeA,sizeB,qDim
    real(8),dimension(:),allocatable             :: qn,qm,dq
    integer,dimension(:,:,:),allocatable         :: tMap
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable        :: Hij
#else
    real(8),dimension(:,:),allocatable           :: Hij
#endif
    character(len=:),allocatable                 :: lkey,rkey
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
    !
    Nsb  = size(sb_sector)
    !
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
    !    
    allocate(SBleft_states(Nsb),SBright_states(Nsb))
    if(.not.direct_H_lazy)then
       allocate(A(tNso,Nsb),B(tNso,Nsb))
       allocate(Hleft(Nsb),Hright(Nsb))
       allocate(isb2jsb(tNso,Nsb));isb2jsb=0
       allocate(IsHconjg(tNso,Nsb));IsHconjg=0
    endif
    !
    if(MpiMaster)t0=t_start()
    !Main computation:
    !>SETUP STATES/MAPS FILTER.  
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
    !
    !free the memory first:
    call free_lazy_operators()
    qDim=size(current_target_QN)
    !
    allocate(Lazy_Sl_n(Nspin),Lazy_Sr_n(Nspin))
    allocate(Lazy_Sl_p(Nspin),Lazy_Sr_p(Nspin))
    allocate(Lazy_dql_n(qDim,Nspin),Lazy_dqr_n(qDim,Nspin))
    allocate(Lazy_dql_p(qDim,Nspin),Lazy_dqr_p(qDim,Nspin))
    if(MpiMaster)then
       !is=1 => Sz
       !is=2 => S+
       do is=1,Nspin
          lkey            = "S"//left%okey(0,is,ilink='n')
          rkey            = "S"//right%okey(0,is,ilink='n')
          Lazy_Sl_n(is)   = left%operators%op(key=lkey)
          Lazy_Sr_n(is)   = right%operators%op(key=rkey)
          Lazy_dql_n(:,is)= left%operators%dq(key=lkey)           
          Lazy_dqr_n(:,is)= right%operators%dq(key=rkey)
          if(PBCdmrg)then
             lkey             = "S"//left%okey(0,is,ilink='p')
             rkey             = "S"//right%okey(0,is,ilink='p')
             Lazy_Sl_p(is)   = left%operators%op(key=lkey)
             Lazy_Sr_p(is)   = right%operators%op(key=rkey)
             Lazy_dql_p(:,is)= left%operators%dq(key=lkey)             
             Lazy_dqr_p(:,is)= right%operators%dq(key=rkey)
          endif
       enddo
       Lazy_Hl = left%operators%op("H")
       Lazy_Hr = right%operators%op("H")
    endif
#ifdef _MPI
    if(MpiStatus)then
       do is=1,Nspin
          call Lazy_Sl_n(is)%bcast()
          call Lazy_Sr_n(is)%bcast()
          if(PBCdmrg)then
             call Lazy_Sl_p(is)%bcast()
             call Lazy_Sr_p(is)%bcast()
          endif
       enddo
       call Lazy_Hl%bcast()
       call Lazy_Hr%bcast()
       call Bcast_MPI(MpiComm,Lazy_dql_n)
       call Bcast_MPI(MpiComm,Lazy_dqr_n)
       if(PBCdmrg)then
          call Bcast_MPI(MpiComm,Lazy_dql_p)
          call Bcast_MPI(MpiComm,Lazy_dqr_p) 
       endif 
    endif
#endif
    if(MpiMaster)write(LOGfile,*),"Load Setup Operators:",t_stop()
    !
    !In lazy mode the unfiltered operators remain in the Lazy_* pool.
    !In cache mode they are only temporary sources for A/B/Hleft/Hright.
    if(.not.direct_H_lazy)then
      if(MpiMaster)t0=t_start()
      isb2jsb=0
      do isb=1+MpiRank,size(sb_sector),MpiSize
        sizeA=size(SBleft_states(isb)%states)
        sizeB=size(SBright_states(isb)%states)
        if(MpiMaster.AND.sizeA>10)&
              write(LOGfile,*)"isb:"//str(isb)//"/"//str(size(sb_sector))//&
              " N(isb):"//str(sizeA)//","//str(sizeB)
        !
        qn = sb_sector%qn(index=isb)
        Hleft(isb) = filter_left_operator(Lazy_Hl,isb,isb)
        Hright(isb)= filter_right_operator(Lazy_Hr,isb,isb)
        !
        it = tMap(1,1,1)
        A(it,isb) = HopH(1,1)*filter_left_operator(Lazy_Sl_n(1),isb,isb)
        B(it,isb) = filter_right_operator(Lazy_Sr_n(1),isb,isb)
        jsb = isb
        Isb2Jsb(it,isb) = jsb
        IsHconjg(it,isb)= 0
        RowOffset(it,isb)=Offset(isb)
        ColOffset(it,isb)=Offset(isb)
        !
        ! formerly hard coded: dq = [1d0]
        dq = Lazy_dql_n(:,2) !spin=1: S_z, spin=2: S_+
        qm = qn - dq
        if(sb_sector%has_qn(qm))then
          jsb = sb_sector%index(qn=qm)
          !
          it=tMap(2,1,1)
          A(it,isb) = HopH(2,2)*filter_left_operator(Lazy_Sl_n(2),isb,jsb)
          B(it,isb) = filter_right_operator(hconjg(Lazy_Sr_n(2)),isb,jsb)
          Isb2Jsb(it,isb)  = jsb
          IsHconjg(it,isb) = 0
          RowOffset(it,isb)= Offset(isb)
          ColOffset(it,isb)= Offset(jsb)
          !
          it=tMap(3,1,1)
          A(it,isb) = hconjg(A(tMap(2,1,1),isb))
          B(it,isb) = hconjg(B(tMap(2,1,1),isb))
          Isb2Jsb(it,isb)  = jsb
          IsHconjg(it,isb) = 1
          RowOffset(it,isb)= Offset(jsb)
          ColOffset(it,isb)= Offset(isb)
        endif
        !
        if(PBCdmrg)then
            it = tMap(4,1,1)
            A(it,isb) = HopH(1,1)*filter_left_operator(Lazy_Sl_p(1),isb,isb)
            B(it,isb) = filter_right_operator(Lazy_Sr_p(1),isb,isb)
            jsb = isb
            Isb2Jsb(it,isb) = jsb
            IsHconjg(it,isb)= 0
            RowOffset(it,isb)= Offset(isb)
            ColOffset(it,isb)= Offset(isb)
            !
            dq = Lazy_dql_p(:,2)
            qm = qn - dq
            if(sb_sector%has_qn(qm))then
              jsb = sb_sector%index(qn=qm)
              !
              it=tMap(5,1,1)
              A(it,isb) = HopH(2,2)*filter_left_operator(Lazy_Sl_p(2),isb,jsb)
              B(it,isb) = filter_right_operator(hconjg(Lazy_Sr_p(2)),isb,jsb)
              Isb2Jsb(it,isb)  = jsb
              IsHconjg(it,isb) = 0
              RowOffset(it,isb)= Offset(isb)
              ColOffset(it,isb)= Offset(jsb)
              !
              it=tMap(6,1,1)
              A(it,isb) = hconjg(A(tMap(5,1,1),isb))
              B(it,isb) = hconjg(B(tMap(5,1,1),isb))
              Isb2Jsb(it,isb)  = jsb
              IsHconjg(it,isb) = 1
              RowOffset(it,isb)= Offset(jsb)
              ColOffset(it,isb)= Offset(isb)
            endif
        endif
      enddo
      if(MpiMaster)write(LOGfile,*)"Get Cache Op.blocks:",t_stop()
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
        do isb=1,size(sb_sector)
            kb_sb_setup_bcast = kb_sb_setup_bcast + Hleft(isb)%bytes() + Hright(isb)%bytes()
            do it=1,tNso
              kb_sb_setup_bcast = kb_sb_setup_bcast + A(it,isb)%bytes()  + B(it,isb)%bytes()
            enddo
        enddo
        if(MpiMaster)write(LOGfile,*)"MpiComm Cache Op.blocks:",t_stop()
      endif
#endif
      !Free Lazy operators in cache mode
      call free_lazy_operators()
    endif
    !
    if(MpiMaster)call stop_timer("Setup SB Direct")
    !
  end subroutine Setup_SuperBlock_Spin_Direct


END MODULE DMRG_SUPERBLOCK_SETUP_DIRECT_SPIN
