! Compile with the matching real or complex DMRG library and -D_MPI.
program test_vector_transpose_mpi
  use MPI
  use DMRG_GLOBAL, only: vector_transpose_MPI, MpiRank, MpiSize
  implicit none
  integer :: ierr,rank,nproc,icase,nrow,ncol,nactive,sub_comm,color
  integer :: qcol,qrow,col_first,row_first,i,j,fail,all_fail
  integer,parameter :: cases(2,7)=reshape([7,11, 2,3, 1,1, 5,2, 3,2, 8,8, 2,9],[2,7])
#ifdef _CMPLX
  complex(8),allocatable :: a(:,:),b(:,:)
  complex(8) :: expected
#else
  real(8),allocatable :: a(:,:),b(:,:)
  real(8) :: expected
#endif

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  MpiRank=rank
  MpiSize=nproc
  fail=0
  do icase=1,size(cases,2)
     nrow=cases(1,icase)
     ncol=cases(2,icase)
     nactive=max(min(nrow,nproc),min(ncol,nproc))
     color=MPI_UNDEFINED
     if(rank<nactive)color=1
     call MPI_Comm_split(MPI_COMM_WORLD,color,rank,sub_comm,ierr)
     if(sub_comm/=MPI_COMM_NULL)then
        qcol=ncol/nproc+merge(1,0,rank<mod(ncol,nproc))
        qrow=nrow/nproc+merge(1,0,rank<mod(nrow,nproc))
        col_first=rank*(ncol/nproc)+min(rank,mod(ncol,nproc))+1
        row_first=rank*(nrow/nproc)+min(rank,mod(nrow,nproc))+1
        allocate(a(nrow,qcol),b(ncol,qrow))
        do j=1,qcol
           do i=1,nrow
#ifdef _CMPLX
              a(i,j)=cmplx(1000*(col_first+j-1)+i,-i,kind=8)
#else
              a(i,j)=dble(1000*(col_first+j-1)+i)
#endif
           enddo
        enddo
        call vector_transpose_MPI(nrow,qcol,a,ncol,qrow,b,sub_comm)
        do i=1,qrow
           do j=1,ncol
#ifdef _CMPLX
              expected=cmplx(1000*j+row_first+i-1,-(row_first+i-1),kind=8)
#else
              expected=dble(1000*j+row_first+i-1)
#endif
              if(b(j,i)/=expected)fail=fail+1
           enddo
        enddo
        deallocate(a,b)
        call MPI_Comm_free(sub_comm,ierr)
     endif
  enddo
  call MPI_Allreduce(fail,all_fail,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(rank==0)then
     if(all_fail==0)then
        write(*,*)'vector_transpose_MPI: all cases passed'
     else
        write(*,*)'vector_transpose_MPI: mismatches=',all_fail
     endif
  endif
  call MPI_Finalize(ierr)
  if(all_fail/=0)error stop 1
end program test_vector_transpose_mpi
