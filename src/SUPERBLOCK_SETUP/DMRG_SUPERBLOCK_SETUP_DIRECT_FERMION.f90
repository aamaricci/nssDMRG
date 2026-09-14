MODULE DMRG_SUPERBLOCK_SETUP_DIRECT_FERMION
  USE DMRG_GLOBAL
  USE DMRG_SUPERBLOCK_SETUP_GLOBAL
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private


  public :: Setup_SuperBlock_Fermion_Direct




  contains





  !##################################################################
  !         SETUP THE SUPERBLOCK HAMILTONIAN for FERMION PROBLEMS
  !  the filtered operator blocks are stored in the memory as Lazy_*
  !  These are used either:
  !In lazy mode the unfiltered operators remain in the Lazy_* pool.
  !In cache mode they are only temporary sources for A/B/Hleft/Hright.
  !  - lazy mode: the unfiltered operators remain in memory and blocks 
  !    are extracted  on-the-fly at each H*v call.
  !  - cache mode: the filtered operators are tmp sources to A/B/Hleft/Hright 
  !    used in H*v operations
  !##################################################################
  subroutine Setup_SuperBlock_Fermion_Direct()
    integer                               :: Nso,Nsb
    integer                               :: it,jsb,ierr,ipr,fbc,sizeA,sizeB,qDim
    real(8),dimension(:),allocatable      :: qn,qm,dq
    integer,dimension(:,:,:),allocatable  :: tMap
    integer                               :: i,io,jo,is,isb,iorb_,ispin_,io_,jo_
    type(sparse_matrix)                   :: Ctmp,Pn,Pp
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hij
#else
    real(8),dimension(:,:),allocatable    :: Hij
#endif
    character(len=:),allocatable          :: lkey,rkey,pkey
    real(8),allocatable                   :: dqC(:),dqPn(:),dqPp(:)
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
    tNso = fbc*count(Hij/=zero) !Nspin=2*# of non-zero terms
    !
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
    !
    allocate(SBleft_states(Nsb),SBright_states(Nsb))
    if(.not.direct_H_lazy)then
       allocate(A(tNso,Nsb),B(tNso,Nsb))
       allocate(Hleft(Nsb),Hright(Nsb))
       allocate(isb2jsb(tNso,Nsb));isb2jsb=0
       allocate(IsHconjg(tNso,Nsb));IsHconjg=0
    endif
    !
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
    !free the memory first:s
    call free_lazy_operators()
    qDim = size(current_target_qn)
    !
    allocate(Lazy_CdgP_n(Nspin*Norb),Lazy_Cr_n(Nspin*Norb))
    allocate(Lazy_CdgP_p(Nspin*Norb),Lazy_Cr_p(Nspin*Norb))
    allocate(Lazy_dql_n(qDim,Nspin*Norb),Lazy_dqr_n(qDim,Nspin*Norb))
    allocate(Lazy_dql_p(qDim,Nspin*Norb),Lazy_dqr_p(qDim,Nspin*Norb))    
    if(MpiMaster)then
       lkey = "P"//left%okey(0,0,ilink='n')
       Pn   = left%operators%op(key=lkey)
       dqPn = left%operators%dq(key=lkey)
       if(PBCdmrg)then  
          lkey =  "P"//left%okey(0,0,ilink='p')
          Pp   = left%operators%op(key=lkey)
          dqPp = left%operators%dq(key=lkey)
       endif
       do ispin_=1,Nspin
          do iorb_=1,Norb
             io_  = iorb_+(ispin_-1)*Norb
             lkey = "C"//left%okey(iorb_,ispin_,ilink='n')
             rkey = "C"//right%okey(iorb_,ispin_,ilink='n')
             Ctmp = left%operators%op(key=lkey)
             dqC  = left%operators%dq(key=lkey)
             Lazy_CdgP_n(io_)  = matmul(Ctmp%dgr(),Pn)
             Lazy_Cr_n(io_)    = right%operators%op(key=rkey)
             Lazy_dql_n(:,io_) = -dqC + dqPn
             Lazy_dqr_n(:,io_) = right%operators%dq(key=rkey)
             call Ctmp%free()
             if(PBCdmrg)then
                lkey = "C"//left%okey(iorb_,ispin_,ilink='p')
                rkey = "C"//right%okey(iorb_,ispin_,ilink='p')              
                Ctmp = left%operators%op(key=lkey)
                dqC  = left%operators%dq(key=lkey)
                Lazy_CdgP_p(io_)  = matmul(Ctmp%dgr(),Pp)
                Lazy_Cr_p(io_)    = right%operators%op(key=rkey)
                Lazy_dql_p(:,io_) = -dqC + dqPp
                Lazy_dqr_p(:,io_) = right%operators%dq(key=rkey)
                call Ctmp%free()
             endif
          enddo
       enddo
       Lazy_Hl = left%operators%op("H")
       Lazy_Hr = right%operators%op("H")
       !
       call Pn%free()
       if(PBCdmrg)call Pp%free()
       !
    endif
#ifdef _MPI
    if(MpiStatus)then
       do is=1,Nspin*Norb
          call Lazy_CdgP_n(is)%bcast()
          call Lazy_Cr_n(is)%bcast()
          if(PBCdmrg)then
             call Lazy_CdgP_p(is)%bcast()
             call Lazy_Cr_p(is)%bcast()
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
        do io_=1,Nso
          dq = Lazy_dql_n(:,io_)
          qm = qn - dq
          if(sb_sector%has_qn(qm))then
            jsb = sb_sector%index(qn=qm)
            do jo_=1,Nso
              if(HopH(io_,jo_)==zero)cycle
              !
              it=tMap(1,io_,jo_)
              A(it,isb) = HopH(io_,jo_)*filter_left_operator(Lazy_CdgP_n(io_),isb,jsb)
              B(it,isb) = filter_right_operator(Lazy_Cr_n(jo_),isb,jsb)
              Isb2Jsb(it,isb)  = jsb
              IsHconjg(it,isb) = 0
              RowOffset(it,isb)= Offset(isb)
              ColOffset(it,isb)= Offset(jsb)
              !
              it=tMap(2,io_,jo_)
              A(it,isb) = hconjg(A(tMap(1,io_,jo_),isb))
              B(it,isb) = hconjg(B(tMap(1,io_,jo_),isb))
              Isb2Jsb(it,isb)  = jsb
              IsHconjg(it,isb) = 1
              RowOffset(it,isb)= Offset(jsb)
              ColOffset(it,isb)= Offset(isb)
            enddo
          endif
        enddo
        !
        !PBC:
        if(PBCdmrg)then
          do io_=1,Nso
            dq = Lazy_dql_p(:,io_)
            qm = qn - dq
            if(sb_sector%has_qn(qm))then
              jsb = sb_sector%index(qn=qm)
              do jo_=1,Nso
                if(HopH(io_,jo_)==zero)cycle
                it=tMap(3,io_,jo_)
                A(it,isb) = HopH(io_,jo_)*filter_left_operator(Lazy_CdgP_p(io_),isb,jsb)
                B(it,isb) = filter_right_operator(Lazy_Cr_p(jo_),isb,jsb)
                Isb2Jsb(it,isb)  = jsb
                IsHconjg(it,isb) = 0
                RowOffset(it,isb)= Offset(isb)
                ColOffset(it,isb)= Offset(jsb)
                !
                it=tMap(4,io_,jo_)
                A(it,isb) = hconjg(A(tMap(3,io_,jo_),isb))
                B(it,isb) = hconjg(B(tMap(3,io_,jo_),isb))
                Isb2Jsb(it,isb)  = jsb
                IsHconjg(it,isb) = 1
                RowOffset(it,isb)= Offset(jsb)
                ColOffset(it,isb)= Offset(isb)
              enddo
            endif
          enddo
        endif
        !
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
      !
    endif
    !
    if(MpiMaster)call stop_timer("Setup SB Direct")
    !
  end subroutine Setup_SuperBlock_Fermion_Direct

  




END MODULE DMRG_SUPERBLOCK_SETUP_DIRECT_FERMION
