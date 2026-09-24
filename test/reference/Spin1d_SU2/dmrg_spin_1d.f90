program dmrg_spin_1d
  USE SCIFOR
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=64)                  :: finput
  integer                            :: i,j,SUN,Unit,pos,Nsites
  real(8)                            :: Hvec,Noise,Sij,Espin,Eloc,Etotal
  type(site)                         :: MyDot
  type(sparse_matrix)                :: Sz,Sz2,Hi
  real(8),dimension(:,:),allocatable :: Hlr
  real(8),dimension(:),allocatable   :: avSz,avSz2
  integer                            :: irank,comm,rank,ierr
  logical                            :: master=.true.,irun,imeasure
  character(len=:),allocatable       :: key_Sz,run_label
  
#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(imeasure,"imeasure",finput,default=.true.,&
       comment="Bool to perform measurements. T for post-processing.")
  call parse_input_variable(irun,"irun",finput,default=.true.,&
       comment="Bool to run DMRG. F for post-processing")       
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(Noise,"NOISE",finput,default=0d0,&
       comment="Magnetic field noise amplitude")
  call parse_input_variable(Hvec,"Hvec",finput,default=0d0,&
       comment="Magnetic field direction")

  call read_input(finput)

  if(Imeasure)then
     save_block=.true.
     save_umat=.true.
  endif

  Nsites=2*Ldmrg

  MyDot = spin_site(sun=SUN,Hz=Hvec)
  Hlr   = diag([Jp,Jx/2d0])

  !Init DMRG
  call init_dmrg(Hlr,ModelDot=[MyDot])


  !Run DMRG algorithm
  if(Irun)call run_DMRG()



  if(imeasure)then
     !Post-processing and measure quantities:
     !Measure <Sz(i)>, <Sz(i).Sz(i)>
     key_Sz="S"//mydot%okey(0,1,ilink="n")
     Sz =MyDot%operators%op(key_Sz)
     Sz2=matmul(Sz,Sz)
     !
     call Measure_DMRG(Sz ,pos=arange(1,Nsites),avOp=avSz)
     call Measure_DMRG(Sz2,pos=arange(1,Nsites),avOp=avSz2)
     if(master)then
        unit=fopen("spin_local.check",append=.false.)
        do i=1,Nsites
           write(unit,*)i,avSz(i),avSz2(i)
        enddo
        close(unit)
     endif

     !Nearest-neighbour reference correlations: i, i+1, <S_i.S_(i+1)>.
     if(master)unit=fopen("spin_nn.check",append=.false.)
     do i=1,Nsites-1
        Sij=Measure_SpinSpin_DMRG(i,i+1)
        if(master)write(unit,*)i,i+1,Sij
     enddo
     if(master)close(unit)

     !Long-range reference correlations: 1, j, <S_1.S_j>.
     if(master)unit=fopen("spin_1j.check",append=.false.)
     do j=1,Nsites
        Sij=Measure_SpinSpin_DMRG(1,j)
        if(master)write(unit,*)1,j,Sij
     enddo
     if(master)close(unit)

     !Exchange, local and total energies from the final effective
     !Hamiltonian.  Hi stores the site-resolved local contribution.
     call Measure_Energy_DMRG(Hlr,Espin,Eloc,Etotal,Hi)
     if(master)then
        unit=fopen("Ecomponents"//str(label_DMRG('u')),append=.true.)
        write(unit,*)Espin,Eloc,Etotal
        close(unit)
     endif
     call End_Measure_DMRG()
     call Sz%free()
     call Sz2%free()
     call Hi%free()
  endif

  if(master)then
     run_label=label_DMRG(DMRGtype)
     call copy_table("energyVSblock.length"//run_label,"energy.check")
     call copy_table("SentropyVSblock.length"//run_label,"entropy.check")
  endif

  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif

contains

  subroutine copy_table(source,destination)
    character(len=*),intent(in) :: source,destination
    character(len=4096)         :: line
    integer                     :: source_unit,destination_unit,ios
    open(newunit=source_unit,file=source,status="old",action="read")
    open(newunit=destination_unit,file=destination,status="replace",action="write")
    do
       read(source_unit,'(A)',iostat=ios)line
       if(ios/=0)exit
       write(destination_unit,'(A)')trim(line)
    enddo
    close(source_unit);close(destination_unit)
  end subroutine copy_table

end program dmrg_spin_1d
