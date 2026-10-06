!Small-chain driver used for ED comparisons and cross-version restart checks.
program qn_shift_test
  use SCIFOR
  use DMRG
  use DMRG_GLOBAL, only: gs_energy,current_target_qn,current_L
#ifdef _MPI
  use MPI
#endif
  implicit none
  character(len=256) :: finput
  integer :: model,i,unit,n,comm
  logical :: master=.true.,do_run
  type(site) :: sites(1)
  type(sparse_matrix) :: op(2),c
#ifdef _CMPLX
  complex(8),allocatable :: hopping(:,:)
#else
  real(8),allocatable :: hopping(:,:)
#endif
  real(8),allocatable :: values(:,:)
  real(8) :: ebond,eloc,energy
#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  master=get_Master_MPI(comm)
#endif
  call parse_cmd_variable(finput,"FINPUT",default="DMRG.conf")
  call parse_input_variable(model,"TEST_MODEL",finput,default=0)
  call parse_input_variable(do_run,"IRUN",finput,default=.true.)
  call read_input(finput)
  if(model/=1)then
     if(model==0)then
        sites(1)=spin_site(sun=2)
     else
        sites(1)=spin_site(sun=3)
     endif
     hopping=diag([1d0,0.5d0])
     op(1)=sites(1)%operators%op("S"//sites(1)%okey(0,1,ilink="n"))
  else
     sites(1)=electron_site()
     hopping=diag([-1d0,-1d0])
     do i=1,2
        c=sites(1)%operators%op("C"//sites(1)%okey(1,i,ilink="n"))
        op(i)=matmul(c%dgr(),c)
     enddo
  endif
  call init_dmrg(hopping,ModelDot=sites)
  if(do_run)call run_DMRG()
  n=2*Ldmrg
  if(n>4.and.to_lower(DMRGtype)=="i")then
     if(model/=1)then
        call Measure_DMRG(op(:1),pos=arange(1,n),avOp=values)
     else
        call Measure_DMRG(op,pos=arange(1,n),avOp=values)
     endif
     call Measure_Energy_DMRG(hopping,ebond,eloc,energy)
  else
     !The existing measurement initializer requires rotation history beyond L=4.
     !At L=4 no block truncation has occurred: use the diagonalized energy/sector.
     energy=gs_energy(1)
  endif
  if(master)then
     open(newunit=unit,file="qn_result.out",status="replace")
     write(unit,*)energy
     if(n>4.and.to_lower(DMRGtype)=="i")then
        write(unit,*)sum(values,dim=2)
     else
        write(unit,*)current_target_qn
     endif
     write(unit,*)current_L
     close(unit)
  endif
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif
end program qn_shift_test
