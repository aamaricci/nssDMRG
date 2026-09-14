MODULE DMRG_SUPERBLOCK_HXV
  USE DMRG_GLOBAL
  USE DMRG_SUPERBLOCK_SETUP_GLOBAL
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private

  public :: spMatVec_sparse_main
  public :: spMatVec_cache_main
  public :: spMatVec_lazy_main
#ifdef _MPI
  public :: spMatVec_MPI_cache_main
  public :: spMatVec_MPI_lazy_main
#endif


contains



!##################################################################
!              SuperBlock MATRIX-VECTOR PRODUCTS
!                          SPARSE OP 
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










  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !                     CACHE OP (serial + MPI)
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_cache_main(Nloc,v,Hv)
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
  end subroutine spMatVec_cache_main




#ifdef _MPI
  subroutine spMatVec_MPI_cache_main(Nloc,v,Hv)
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
  end subroutine spMatVec_MPI_cache_main
#endif

















  !##################################################################
  !              SuperBlock MATRIX-VECTOR PRODUCTS
  !                     LAZY OP (serial + MPI)
  !              using shared quantities in GLOBAL
  !##################################################################
  subroutine spMatVec_lazy_main(Nloc,v,Hv)
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
    type(sparse_matrix)                   :: Hlk,Hrk,Aop,Bop
    real(8),dimension(:),allocatable      :: qn,qm,dq
    character(len=:),allocatable          :: type
    integer                               :: k,q,ir,il,jr,jl,jcol,i,j,Nso
    integer                               :: ispin,iorb,jorb,io,jo
    !
    Hv=zero
    t0=t_start()
    Nso=Nspin*Norb
    !
    !Detect the type of operators
    type=str(left%type())
    !
    !Loop over the sectors:
    sector: do k=1,size(sb_sector)
       !apply the 1^L x H^r
       Hrk = filter_right_operator(Lazy_Hr,k,k)
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
       !
       !apply the H^L x 1^r
       Hlk = filter_left_operator(Lazy_Hl,k,k)
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
       !
       !apply sum_k A^L_k x B^R_k
       qn = sb_sector%qn(index=k)
       select case(to_lower(type(1:1)))
       case("s")                !SPIN
          !This part is Sz.Sz
          Aop = HopH(1,1)*filter_left_operator(Lazy_Sl_n(1),k,k)
          Bop = filter_right_operator(Lazy_Sr_n(1),k,k)
          call apply_AxB_direct(Aop,Bop,Offset(k),Offset(k),v,Hv)
          call Aop%free()
          call Bop%free()
          !
          !This is Sp.S- (+ H.c.)
          dq = Lazy_dql_n(:,2) !spin=1: S_z, spin=2: S_+
          qm = qn - dq
          if(sb_sector%has_qn(qm))then
             q = sb_sector%index(qn=qm)
             Aop = HopH(2,2)*filter_left_operator(Lazy_Sl_n(2),k,q)
             Bop = filter_right_operator(hconjg(Lazy_Sr_n(2)),k,q)
             call apply_AxB_direct(Aop,Bop,Offset(k),Offset(q),v,Hv)
             call apply_AxB_direct(hconjg(Aop),hconjg(Bop),Offset(q),Offset(k),v,Hv)
             call Aop%free()
             call Bop%free()
          endif
          if(PBCdmrg)then
            !This part is Sz.Sz
             Aop = HopH(1,1)*filter_left_operator(Lazy_Sl_p(1),k,k)
             Bop = filter_right_operator(Lazy_Sr_p(1),k,k)
             call apply_AxB_direct(Aop,Bop,Offset(k),Offset(k),v,Hv)
             call Aop%free()
             call Bop%free()
             !
             !This is Sp.S- (+ H.c.)
             qm = qn - Lazy_dql_p(:,2) ![1d0]
             if(sb_sector%has_qn(qm))then
                q = sb_sector%index(qn=qm)
                Aop = HopH(2,2)*filter_left_operator(Lazy_Sl_p(2),k,q)
                Bop = filter_right_operator(hconjg(Lazy_Sr_p(2)),k,q)
                call apply_AxB_direct(Aop,Bop,Offset(k),Offset(q),v,Hv)
                call apply_AxB_direct(hconjg(Aop),hconjg(Bop),Offset(q),Offset(k),v,Hv)
                call Aop%free()
                call Bop%free()
             endif
          endif
          !
       case("f","e")            !FERMIONS
         !
         do io=1,Nso
            dq = Lazy_dql_n(:,io)
            qm = qn - dq
            if(sb_sector%has_qn(qm))then
              q = sb_sector%index(qn=qm)
              do jo=1,Nso
                if(HopH(io,jo)==zero)cycle
                Aop = HopH(io,jo)*filter_left_operator(Lazy_CdgP_n(io),k,q)
                Bop = filter_right_operator(Lazy_Cr_n(jo),k,q)
                call apply_AxB_direct(Aop,Bop,Offset(k),Offset(q),v,Hv)
                call apply_AxB_direct(hconjg(Aop),hconjg(Bop),Offset(q),Offset(k),v,Hv)
                call Aop%free()
                call Bop%free()
              enddo
            endif
         enddo
         !PBC:
         if(PBCdmrg)then
          do io=1,Nso
            dq = Lazy_dql_p(:,io)
            qm = qn - dq
            if(sb_sector%has_qn(qm))then
              q = sb_sector%index(qn=qm)
              do jo=1,Nso
                if(HopH(io,jo)==zero)cycle
                Aop = HopH(io,jo)*filter_left_operator(Lazy_CdgP_p(io),k,q)
                Bop = filter_right_operator(Lazy_Cr_p(jo),k,q)
                call apply_AxB_direct(Aop,Bop,Offset(k),Offset(q),v,Hv)
                call apply_AxB_direct(hconjg(Aop),hconjg(Bop),Offset(q),Offset(k),v,Hv)
                call Aop%free();call Bop%free()
              enddo
            endif
          enddo
        endif
        !
       end select
    enddo sector
    t_hxv_direct=t_hxv_direct + t_stop()
  end subroutine spMatVec_lazy_main



#ifdef _MPI
  subroutine spMatVec_MPI_lazy_main(Nloc,v,Hv)
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
    type(sparse_matrix)                   :: Hlk,Hrk,Aop,Bop
    real(8),dimension(:),allocatable      :: qn,qm,dq
    character(len=:),allocatable          :: type
    integer                               :: k,q,ir,il,jr,jl,jcol,i,j,Nso
    integer                               :: ispin,iorb,jorb,io,jo
    integer                               :: i_start,i_end
    !
    if(.not.MpiStatus)stop "spMatVec_MPI_direct_lazy_main ERROR: MpiStatus = F"
    Hv=zero
    t0=t_start()
    Nso=Nspin*Norb
    !
    !Detect the type of operators
    type=str(left%type())
    !
    !Loop over the sectors:
    sector: do k=1,size(sb_sector)
       !> apply the 1^L x H^r: share L columns
       Hrk = filter_right_operator(Lazy_Hr,k,k)
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
       !
       !> apply the H^L x 1^r: share R rows, requires MPI transpose
       Hlk = filter_left_operator(Lazy_Hl,k,k)
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
       !
       !apply the sum_k A_L^k x B_R^k
       !See the corresponding part in the cache procedure to get a better idea
       !of what is going on here. apply_AxB_MPI_direct does all the work 
       qn = sb_sector%qn(index=k)
       select case(to_lower(type(1:1)))
       case("s")                !SPIN
          !This part is Sz.Sz
          Aop = HopH(1,1)*filter_left_operator(Lazy_Sl_n(1),k,k)
          Bop = filter_right_operator(Lazy_Sr_n(1),k,k)
          call apply_AxB_MPI_direct(Aop,Bop,k,k,0,v,Hv)
          call Aop%free()
          call Bop%free()
          !
          !This is Sp.S- (+ H.c.)
          dq = Lazy_dql_n(:,2) !spin=1: S_z, spin=2: S_+
          qm = qn - dq
          if(sb_sector%has_qn(qm))then
             q = sb_sector%index(qn=qm)
             Aop = HopH(2,2)*filter_left_operator(Lazy_Sl_n(2),k,q)
             Bop = filter_right_operator(hconjg(Lazy_Sr_n(2)),k,q)
             call apply_AxB_MPI_direct(Aop,Bop,k,q,0,v,Hv)
             call apply_AxB_MPI_direct(hconjg(Aop),hconjg(Bop),k,q,1,v,Hv)
             call Aop%free()
             call Bop%free()
          endif
          if(PBCdmrg)then
             !This part is Sz.Sz
             Aop = HopH(1,1)*filter_left_operator(Lazy_Sl_p(1),k,k)
             Bop = filter_right_operator(Lazy_Sr_p(1),k,k)
             call apply_AxB_MPI_direct(Aop,Bop,k,k,0,v,Hv)
             call Aop%free()
             call Bop%free()
             !
             !This is Sp.S- (+ H.c.)
             qm = qn - Lazy_dql_p(:,2)
             if(sb_sector%has_qn(qm))then
                q = sb_sector%index(qn=qm)
                Aop = HopH(2,2)*filter_left_operator(Lazy_Sl_p(2),k,q)
                Bop = filter_right_operator(hconjg(Lazy_Sr_p(2)),k,q)
                call apply_AxB_MPI_direct(Aop,Bop,k,q,0,v,Hv)
                call apply_AxB_MPI_direct(hconjg(Aop),hconjg(Bop),k,q,1,v,Hv)
                call Aop%free()
                call Bop%free()
             endif
          endif
          !
       case("f","e")            !FERMIONS
          do io=1,Nso
          dq = Lazy_dql_n(:,io)
          qm = qn - dq
          if(sb_sector%has_qn(qm))then
            q = sb_sector%index(qn=qm)
            do jo=1,Nso
              if(HopH(io,jo)==zero)cycle
              Aop = HopH(io,jo)*filter_left_operator(Lazy_CdgP_n(io),k,q)
              Bop = filter_right_operator(Lazy_Cr_n(jo),k,q)
              call apply_AxB_MPI_direct(Aop,Bop,k,q,0,v,Hv)
              call apply_AxB_MPI_direct(hconjg(Aop),hconjg(Bop),k,q,1,v,Hv)
              call Aop%free()
              call Bop%free()
            enddo
          endif
        enddo
        !PBC:
        if(PBCdmrg)then
          do io=1,Nso
            dq = Lazy_dql_p(:,io)
            qm = qn - dq
            if(sb_sector%has_qn(qm))then
              q = sb_sector%index(qn=qm)
              do jo=1,Nso
                if(HopH(io,jo)==zero)cycle
                Aop = HopH(io,jo)*filter_left_operator(Lazy_CdgP_p(io),k,q)
                Bop = filter_right_operator(Lazy_Cr_p(jo),k,q)
                call apply_AxB_MPI_direct(Aop,Bop,k,q,0,v,Hv)
                call apply_AxB_MPI_direct(hconjg(Aop),hconjg(Bop),k,q,1,v,Hv)
                call Aop%free()
                call Bop%free()
              enddo
            endif
          enddo
        endif
        !
       end select
    enddo sector
    !
    t_hxv_direct=t_hxv_direct + t_stop()
    !
  end subroutine spMatVec_MPI_lazy_main
#endif


  




  !##################################################################
  !                       APPLY AxB 
  !                  DIRECT (serial + MPI)
  !##################################################################  
  subroutine apply_AxB_direct(Aop,Bop,row_offset,col_offset,v,Hv)
    type(sparse_matrix),intent(in)        :: Aop,Bop
    integer,intent(in)                    :: row_offset,col_offset
#ifdef _CMPLX
    complex(8),dimension(:)               :: v,Hv
    complex(8),dimension(:,:),allocatable :: C
    complex(8)                            :: val
#else
    real(8),dimension(:)                  :: v,Hv
    real(8),dimension(:,:),allocatable    :: C
    real(8)                               :: val
#endif
    integer                               :: ai,aj,bi,bj,ja,jb,j,ic,i,jc
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
    type(sparse_matrix),intent(in)        :: Aop,Bop
    integer,intent(in)                    :: k,q,is_hconjg
#ifdef _CMPLX
    complex(8),dimension(:)               :: v,Hv
    complex(8),dimension(:),allocatable   :: vt,Hvt
    complex(8),dimension(:,:),allocatable :: C,Ct
    complex(8)                            :: val
#else
    real(8),dimension(:)                  :: v,Hv
    real(8),dimension(:),allocatable      :: vt,Hvt
    real(8),dimension(:,:),allocatable    :: C,Ct
    real(8)                               :: val
#endif
    integer                               :: ai,aj,bi,bj,ja,jb,jc,j,i
    integer                               :: mpiArow,mpiAcol,mpiBrow
    integer                               :: shift,abcomm,i_start,i_end
    !
    if(.not.Aop%status.OR..not.Bop%status)return
    !
    mpiAcol = mpiDls(q)
    if(is_hconjg==1)mpiAcol=mpiDls(k)
    !
    mpiArow = mpiDls(k)
    if(is_hconjg==1)mpiArow=mpiDls(q)
    !
    mpiBrow = mpiDrs(k)
    if(is_hconjg==1)mpiBrow=mpiDrs(q)
    !
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
    !
    if(mpiNactive(q)>mpiNactive(k))then
       abcomm = mpiSBCOMM(q)
    else
       abcomm = mpiSBCOMM(k)
    endif
    !
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
    !
    i_start = 1 + mpiOffset(k)
    if(is_hconjg==1)i_start = 1 + mpiOffset(q)
    !
    i_end = Bop%Nrow*mpiArow + mpiOffset(k)
    if(is_hconjg==1)i_end = Bop%Nrow*mpiArow + mpiOffset(q)
    !
    Hv(i_start:i_end) = Hv(i_start:i_end) + vt
    !
    deallocate(C,Ct,Hvt,vt)
  end subroutine apply_AxB_MPI_direct
#endif


END MODULE DMRG_SUPERBLOCK_HXV
