program wf_prediction_indices
  implicit none
  integer :: ml,mr,d,il,ir,sl,sr,link
  integer :: xl,xr,il2,ir2,sl2,sr2
  integer :: pl,pr,psl,psr,phase

  ml=3; mr=4; d=2
  do link=1,2
     do il=1,ml
        do ir=1,mr
           do sl=1,d
              do sr=1,d
                 call enlarged_indices(link,il,ir,sl,sr,ml,mr,d,xl,xr)
                 call decode_indices(link,xl,xr,ml,mr,d,il2,ir2,sl2,sr2)
                 if(any([il2,ir2,sl2,sr2]/=[il,ir,sl,sr])) error stop 'index round trip'
                 do pl=0,1
                    do pr=0,1
                       do psl=0,1
                          do psr=0,1
                             phase=graded_phase(link,pl,pr,psl,psr)
                             if(abs(phase)/=1) error stop 'invalid fermion phase'
                             if(mod(psl+psr,2)==0 .and. mod(pl+pr,2)==0) then
                                if(phase/=1) error stop 'even pair phase'
                             endif
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  print *, 'wf prediction index and parity checks passed'
contains
  subroutine enlarged_indices(link,il,ir,sl,sr,ml,mr,d,xl,xr)
    integer,intent(in) :: link,il,ir,sl,sr,ml,mr,d
    integer,intent(out) :: xl,xr
    select case(link)
    case(1) ! n: L.o.sL and sR.o.R; second factor is fastest
       xl=sl+(il-1)*d
       xr=ir+(sr-1)*mr
    case(2) ! p: sL.o.L and R.o.sR
       xl=il+(sl-1)*ml
       xr=sr+(ir-1)*d
    case default
       error stop 'invalid link'
    end select
  end subroutine

  subroutine decode_indices(link,xl,xr,ml,mr,d,il,ir,sl,sr)
    integer,intent(in) :: link,xl,xr,ml,mr,d
    integer,intent(out) :: il,ir,sl,sr
    select case(link)
    case(1)
       il=(xl-1)/d+1; sl=mod(xl-1,d)+1
       sr=(xr-1)/mr+1; ir=mod(xr-1,mr)+1
    case(2)
       sl=(xl-1)/ml+1; il=mod(xl-1,ml)+1
       ir=(xr-1)/d+1; sr=mod(xr-1,d)+1
    case default
       error stop 'invalid link'
    end select
  end subroutine

  integer function graded_phase(link,pl,pr,psl,psr) result(phase)
    integer,intent(in) :: link,pl,pr,psl,psr
    ! Starting order is L,R,sL,sR. Move the site pair to its
    ! position in the enlarged superblock, counting odd crossings.
    select case(link)
    case(1) ! L,sL,sR,R
       phase=1-2*mod(pr*(psl+psr),2)
    case(2) ! sL,L,R,sR
       phase=1-2*mod(psl*(pl+pr),2)
    case default
       error stop 'invalid link'
    end select
  end function
end program
