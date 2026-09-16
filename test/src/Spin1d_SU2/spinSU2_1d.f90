program dmrg_spin_1d
  USE SCIFOR
  USE DMRG
  USE REGRESSION_UTILS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=64)                     :: finput
  character(len=:),allocatable          :: run_label
  integer                               :: i,j,unit,Nsites,comm
  type(site),allocatable                :: MyDot(:)
  type(sparse_matrix)                   :: Sz,Sz2,Jij,Hi
#ifdef _CMPLX
  complex(8),allocatable                :: Hlr(:,:)
  complex(8)                            :: corr
#else
  real(8),allocatable                   :: Hlr(:,:)
  real(8)                               :: corr
#endif
  real(8),allocatable                   :: avSz(:),avSz2(:)
  real(8)                               :: Espin,Eloc,Etotal
  real(8),parameter                     :: atol=1d-8,rtol=1d-7
  real(8),parameter                     :: observable_atol=1d-6
  logical                               :: master=.true.

#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  master=get_Master_MPI(comm)
#endif

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call read_input(finput)
  if(to_lower(DMRGtype)/='i')error stop "Spin regression test requires iDMRG"
  run_label=label_DMRG(DMRGtype)
  if(master)then
     call remove_file("energyVSblock.length"//run_label)
     call remove_file("SentropyVSblock.length"//run_label)
     call remove_file("spin_local.out")
     call remove_file("spin_nn.out")
     call remove_file("spin_1j.out")
  endif
#ifdef _MPI
  call MPI_BARRIER(comm,i)
#endif
  allocate(MyDot(1));MyDot=spin_site(sun=2)
  Hlr=diag([Jp,Jx/2d0])
  call init_dmrg(Hlr,ModelDot=MyDot)
  call run_DMRG()

  !The final iDMRG superblock contains two blocks of length Ldmrg.
  Nsites=2*Ldmrg
  Sz =MyDot(1)%operators%op(key="S"//MyDot(1)%okey(0,1,ilink="n"))
  Sz2=matmul(Sz,Sz)
  call Measure_DMRG(Sz ,pos=arange(1,Nsites),avOp=avSz)
  call Measure_DMRG(Sz2,pos=arange(1,Nsites),avOp=avSz2)

  if(master)then
     open(newunit=unit,file="spin_local.out",status="replace")
     do i=1,Nsites
        write(unit,*)i,avSz(i),avSz2(i)
     enddo
     close(unit)
     open(newunit=unit,file="spin_nn.out",status="replace")
  endif
  do i=1,Nsites-1
     corr=Measure_SpinSpin_DMRG(i,i+1)
     if(master)write(unit,*)i,i+1,real(corr,8)
#ifdef _CMPLX
     if(abs(aimag(corr))>atol)error stop "spin_nn ERROR: non-real correlation"
#endif
  enddo
  if(master)then
     close(unit)
     open(newunit=unit,file="spin_1j.out",status="replace")
  endif
  do j=1,Nsites
     corr=Measure_SpinSpin_DMRG(1,j)
     if(master)write(unit,*)1,j,real(corr,8)
#ifdef _CMPLX
     if(abs(aimag(corr))>atol)error stop "spin_1j ERROR: non-real correlation"
#endif
  enddo
  if(master)close(unit)
  call Measure_Energy_DMRG(Hlr,Espin,Eloc,Etotal,Jij,Hi)
  if(master)write(*,*)"Measured energies [Espin,Eloc,Etotal]:",Espin,Eloc,Etotal
  call End_Measure_DMRG()

  if(master)then
     call assert_table("energy.check","energyVSblock.length"//run_label,3,atol,rtol)
     call assert_table("entropy.check","SentropyVSblock.length"//run_label,4,atol,rtol)
     call assert_table("spin_local.check","spin_local.out",3,observable_atol,rtol)
     call assert_table("spin_nn.check","spin_nn.out",3,observable_atol,rtol)
     call assert_table("spin_1j.check","spin_1j.out",3,observable_atol,rtol)
  endif

  call Sz%free();call Sz2%free()
  call Jij%free();call Hi%free()
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif
end program dmrg_spin_1d
