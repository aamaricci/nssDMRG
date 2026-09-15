module REGRESSION_UTILS
  implicit none
  private
  public :: assert_table,remove_file

contains

  subroutine assert_table(reference,computed,ncol,atol,rtol)
    character(len=*),intent(in) :: reference,computed
    integer,intent(in)          :: ncol
    real(8),intent(in)          :: atol,rtol
    real(8),allocatable         :: ref(:),val(:)
    integer                     :: uref,uval,ios_ref,ios_val,irow
    logical                     :: found_ref,found_val
    !
    open(newunit=uref,file=reference,status='old',action='read',iostat=ios_ref)
    if(ios_ref/=0)error stop "assert_table ERROR: missing reference file"
    open(newunit=uval,file=computed,status='old',action='read',iostat=ios_val)
    if(ios_val/=0)error stop "assert_table ERROR: missing computed file"
    allocate(ref(ncol),val(ncol));irow=0
    do
       call next_row(uref,ncol,ref,found_ref)
       call next_row(uval,ncol,val,found_val)
       if(.not.found_ref.OR..not.found_val)exit
       irow=irow+1
       if(any(abs(val-ref)>atol+rtol*abs(ref)))then
          write(*,*)"Regression failure: ",trim(computed)
          write(*,*)"Row              : ",irow
          write(*,*)"Reference        : ",ref
          write(*,*)"Computed         : ",val
          write(*,*)"Absolute error   : ",abs(val-ref)
          error stop 2
       endif
    enddo
    if(found_ref.neqv.found_val)error stop "assert_table ERROR: different row counts"
    if(irow==0)error stop "assert_table ERROR: empty tables"
    close(uref);close(uval)
    write(*,*)trim(computed),": regression test PASSED"
  end subroutine assert_table


  subroutine remove_file(file)
    character(len=*),intent(in) :: file
    integer                     :: unit,ios
    logical                     :: exists
    !Remove output left by a previous CTest invocation.  DMRG history
    !files are append-only, while every regression run must start clean.
    inquire(file=file,exist=exists)
    if(.not.exists)return
    open(newunit=unit,file=file,status='old',iostat=ios)
    if(ios/=0)error stop "remove_file ERROR: can not open output file"
    close(unit,status='delete',iostat=ios)
    if(ios/=0)error stop "remove_file ERROR: can not delete output file"
  end subroutine remove_file


  subroutine next_row(unit,ncol,row,found)
    integer,intent(in)  :: unit,ncol
    real(8),intent(out) :: row(ncol)
    logical,intent(out) :: found
    character(len=2048) :: line
    integer             :: ios
    !
    found=.false.
    do
       read(unit,'(A)',iostat=ios)line
       if(ios<0)return
       if(ios>0)error stop "assert_table ERROR: failed reading table"
       line=adjustl(line)
       if(len_trim(line)==0)cycle
       if(line(1:1)=='#'.OR.line(1:1)=='!')cycle
       read(line,*,iostat=ios)row
       if(ios/=0)error stop "assert_table ERROR: malformed numeric row"
       found=.true.
       return
    enddo
  end subroutine next_row

end module REGRESSION_UTILS
