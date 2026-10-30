! Original vector_transpose_MPI from commit 502d934.
! Preserved for review before the packed single-collective replacement.
! Every line below is commented out; this file is not compiled.
!   subroutine vector_transpose_MPI(nrow,qcol,a,ncol,qrow,b,sub_comm)
!     !
!     integer                             :: nrow !Global number of rows
!     integer                             :: ncol !Global number of columns
!     integer                             :: qrow !Local number of rows on each thread
!     integer                             :: qcol !Local number of columns on each thread
!     integer                             :: sub_comm !Sub-communicator for active ranks
! #ifdef _CMPLX
!     complex(8)                          :: a(nrow,qcol) ! Input vector to be transposed
!     complex(8)                          :: b(ncol,qrow) ! Output vector :math:`b = v^T`
!     complex(8),dimension(:),allocatable :: Vtmp
! #else
!     real(8)                             :: a(nrow,qcol) ! Input vector to be transposed
!     real(8)                             :: b(ncol,qrow) ! Output vector :math:`b = v^T`
!     real(8),dimension(:),allocatable    :: Vtmp
! #endif
!     integer,allocatable,dimension(:,:)  :: send_counts,send_offset
!     integer,allocatable,dimension(:,:)  :: recv_counts,recv_offset
!     integer                             :: counts,Ntot
!     integer                             :: i,j,irank,ierr
!     !
!     if (sub_comm == MPI_COMM_NULL) return
! 	  !
! 		!
! 		!Round robin distribution of columns: hardcoded for now, can be easily changed to block distribution if needed
!     counts = Nrow/MpiSize
!     Ntot   = Ncol/MpiSize
!     if(mod(Ncol,MpiSize)/=0)Ntot=Ntot+1
!     !
!     allocate(send_counts(0:MpiSize-1,Ntot));send_counts=0
!     allocate(send_offset(0:MpiSize-1,Ntot));send_offset=0
!     allocate(recv_counts(0:MpiSize-1,Ntot));recv_counts=0
!     allocate(recv_offset(0:MpiSize-1,Ntot));recv_offset=0
!     !
! 	 !Round robin distribution of rows: hardcoded for now, can be easily changed to block distribution if needed
!
!     do i=1,qcol
!        do irank=0,MpiSize-1
!           if(irank < mod(Nrow,MpiSize))then
!              send_counts(irank,i) = counts+1
!           else
!              send_counts(irank,i) = counts
!           endif
!        enddo
!     enddo
!     !
!     do i=1,Ntot
!        call MPI_AllToAll(&
!             send_counts(0:,i),1,MPI_INTEGER,&
!             recv_counts(0:,i),1,MPI_INTEGER,&
!             sub_comm,ierr) !not sure if this is the right communicator to use here, but it should be fine since only active ranks have non-zero counts
!     enddo
!     !
!     do i=1,Ntot
!        do irank=1,MpiSize-1
!           send_offset(irank,i) = send_counts(irank-1,i) + send_offset(irank-1,i)
!        enddo
!     enddo
!     !
!     !Get the irank=0 elements, i.e. first entries:
!     recv_offset(0,1) = 0
!     do i=2,Ntot
!        recv_offset(0,i) = sum(recv_counts(0,:i-1))
!     enddo
!     !the rest of the entries:
!     do i=1,Ntot
!        do irank=1,MpiSize-1
!           recv_offset(irank,i) = recv_offset(irank-1,i) + sum(recv_counts(irank-1,:))
!        enddo
!     enddo
!     !
!     !
!     t0=t_start()
!     do j=1,Ntot
!        !Fix issue with empty columns arising from having few MPI nodes
!        if(j<=size(A,2))then
!           Vtmp = A(:,j)            !automatic allocation
!        else
!           allocate(Vtmp(0))
!        endif
!
! #ifdef _CMPLX
!        call MPI_AllToAllV(& ! A(:,j),send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
!             Vtmp,send_counts(:,j),send_offset(:,j),MPI_DOUBLE_COMPLEX,&
!             B(:,:),recv_counts(:,j),recv_offset(:,j),MPI_DOUBLE_COMPLEX,&
!             sub_comm,ierr)
! #else
!        call MPI_AllToAllV(& ! A(:,j),send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
!             Vtmp,send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
!             B(:,:),recv_counts(:,j),recv_offset(:,j),MPI_DOUBLE_PRECISION,&
!             sub_comm,ierr)
! #endif
!        deallocate(Vtmp)
!     enddo
!     t_a2av=t_a2av+t_stop()
!     !
!     call local_transpose(b,ncol,qrow)
!     !
!     !
!     !Accumulate onto global variable (initialized or reset to 0 somewhere else)
!     kb_sent_a2av = kb_sent_a2av + sum(send_counts)*DATA_kb_size
!     kb_recv_a2av = kb_recv_a2av + sum(recv_counts)*DATA_kb_size
!     !
!     return
!   end subroutine vector_transpose_MPI
