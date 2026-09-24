program hubbard_1d
  USE SCIFOR
  USE DMRG
  USE REGRESSION_UTILS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=64)                              :: finput
  character(len=:),allocatable                   :: key,run_label
  integer                                        :: i,j,iorb,ispin
  integer                                        :: unit,Nsites,Nso,comm
  real(8)                                        :: ts(2),Mh(2)
  type(site),allocatable                         :: MyDot(:)
#ifdef _CMPLX
  complex(8),allocatable                         :: Hlr(:,:),corr(:,:)
#else
  real(8),allocatable                            :: Hlr(:,:),corr(:,:)
#endif
  type(sparse_matrix),allocatable                :: Cop(:,:),Nop(:,:)
  type(sparse_matrix)                            :: Docc,H0loc,Hint,Hshift
  type(sparse_matrix)                            :: Cdag
  real(8),allocatable                            :: dqC(:),dqs(:,:)
#ifdef _CMPLX
  complex(8)                                     :: product
#else
  real(8)                                        :: product
#endif
  real(8),allocatable                            :: avLocal(:,:),values(:)
  real(8)                                        :: Ekin,Eloc,Etotal,E0loc,Eint,Eshift
  real(8),parameter                              :: atol=1d-8,rtol=1d-7
  real(8),parameter                              :: observable_atol=1d-6
  logical                                        :: master=.true.

#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  master=get_Master_MPI(comm)
#endif

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -1d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(Mh,"MH",finput,default=(/(0d0,i=1,2 )/),&
       comment="Crystal field splittings")
  call read_input(finput)

  if(to_lower(DMRGtype)/='i')error stop "Hubbard regression test requires iDMRG"
  if(Norb/=1.OR.Nspin/=2)error stop "Hubbard regression test requires Norb=1 and Nspin=2"
  run_label=label_DMRG(DMRGtype)
  if(master)then
     call remove_file("energyVSblock.length"//run_label)
     call remove_file("SentropyVSblock.length"//run_label)
     call remove_file("hubbard_local.out")
     call remove_file("density_nn.out")
     call remove_file("density_1j.out")
  endif
#ifdef _MPI
  call MPI_BARRIER(comm,i)
#endif
  Nso=Nspin*Norb
  allocate(MyDot(1))
  MyDot=electron_site(H0loc=H0loc,Hint=Hint,Hshift=Hshift)
  !
  allocate(Hlr(Nso,Nso));Hlr=diag([ts(1:Norb),ts(1:Norb)])
  call init_dmrg(Hlr,ModelDot=MyDot)
  call run_DMRG()

  Nsites=2*Ldmrg
  allocate(Cop(Norb,Nspin),Nop(Norb,Nspin))
  do ispin=1,Nspin
     do iorb=1,Norb
        key="C"//MyDot(1)%okey(iorb,ispin,ilink="n")
        Cop(iorb,ispin)=MyDot(1)%operators%op(key)
        Nop(iorb,ispin)=matmul(Cop(iorb,ispin)%dgr(),Cop(iorb,ispin))
     enddo
  enddo
  Docc=matmul(Nop(1,1),Nop(1,2))
  call Measure_DMRG([Nop(1,1),Nop(1,2),Docc],&
       pos=arange(1,Nsites),avOp=avLocal)

  if(master)then
     open(newunit=unit,file="hubbard_local.out",status="replace")
     do i=1,Nsites
        write(unit,*)i,avLocal(1,i),avLocal(2,i),avLocal(3,i)
     enddo
     close(unit)
     open(newunit=unit,file="density_nn.out",status="replace")
  endif
  allocate(values(Nso*Nso))
  do i=1,Nsites-1
     corr=Measure_DensityDensity_DMRG(i,i+1)
     call correlation_values(corr,values)
     if(master)write(unit,*)i,i+1,values
  enddo
  if(master)then
     close(unit)
     open(newunit=unit,file="density_1j.out",status="replace")
  endif
  do j=1,Nsites
     corr=Measure_DensityDensity_DMRG(1,j)
     call correlation_values(corr,values)
     if(master)write(unit,*)1,j,values
  enddo
  if(master)close(unit)
  !
  !Four odd factors form n_up(1)n_up(Nsites), including the L/R cut.
  key="C"//MyDot(1)%okey(1,1,ilink="n")
  dqC=MyDot(1)%operators%dq(key)
  allocate(dqs(size(dqC),4))
  dqs(:,1)=-dqC;dqs(:,2)=dqC
  dqs(:,3)=-dqC;dqs(:,4)=dqC
  Cdag=Cop(1,1)%dgr()
  !
  !check <C^+C> with two methods
  product=Measure_Product_DMRG([Cdag,Cop(1,1)],dqs(:,1:2),&
       ["fermionic","fermionic"],[1,Nsites])
  if(abs(product-Measure_Corr_DMRG(Cdag,-dqC,Cop(1,1),dqC,&
       1,Nsites,"fermionic","fermionic"))>observable_atol)&
       error stop "Hubbard product correlation ERROR: fermionic L/R string"
  !
  !check <C.C^+> = 1-<C^+C>=1-n
  product=Measure_Product_DMRG([Cop(1,1),Cdag],-dqs(:,1:2),&
       ["fermionic","fermionic"],[1,1])
  if(abs(product-(1d0-avLocal(1,1)))>observable_atol)&
       error stop "Hubbard product correlation ERROR: factor order"
  !
  !check <C^+C.C^+C> = <n(1)n(Nsites)>
  product=Measure_Product_DMRG([Cdag,Cop(1,1),Cdag,Cop(1,1)],&
       dqs,["fermionic","fermionic","fermionic","fermionic"],&
       [1,1,Nsites,Nsites])
  corr=Measure_DensityDensity_DMRG(1,Nsites)
  if(abs(product-corr(1,1))>observable_atol)&
       error stop "Hubbard product correlation ERROR: density product"
  call Cdag%free()
  !
  call Measure_Energy_DMRG(Hlr,Ekin,Eloc,Etotal,&
       E0loc=E0loc,Eint=Eint,Eshift=Eshift,&
       H0loc=H0loc,Hint=Hint,Hshift=Hshift)
  if(abs(Eloc-E0loc-Eint-Eshift)>observable_atol*max(1d0,abs(Eloc)))&
       error stop "Hubbard local-energy decomposition ERROR"
  if(master)write(*,*)"Measured energies [Ekin,Eloc,Etotal]:",Ekin,Eloc,Etotal
  if(master)write(*,*)"Local components [E0loc,Eint,Eshift]:",E0loc,Eint,Eshift
  call End_Measure_DMRG()

  if(master)then
     call assert_table("energy.check","energyVSblock.length"//run_label,2,atol,rtol)
     call assert_table("entropy.check","SentropyVSblock.length"//run_label,4,atol,rtol)
     call assert_table("hubbard_local.check","hubbard_local.out",4,observable_atol,rtol)
     call assert_table("density_nn.check","density_nn.out",2+Nso*Nso,observable_atol,rtol)
     call assert_table("density_1j.check","density_1j.out",2+Nso*Nso,observable_atol,rtol)
  endif

  do ispin=1,Nspin
     do iorb=1,Norb
        call Cop(iorb,ispin)%free()
        call Nop(iorb,ispin)%free()
     enddo
  enddo
  call Docc%free()
  call H0loc%free()
  call Hint%free()
  call Hshift%free()
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif

contains

  subroutine correlation_values(Cij,vals)
#ifdef _CMPLX
    complex(8),intent(in) :: Cij(:,:)
#else
    real(8),intent(in)    :: Cij(:,:)
#endif
    real(8),intent(out)   :: vals(:)
    integer               :: ia,ib,k
    !Use io as the outer index and jo as the inner index in reference files.
    if(size(vals)/=size(Cij))error stop "correlation_values ERROR: wrong size"
#ifdef _CMPLX
    if(maxval(abs(aimag(Cij)))>atol)error stop "density correlation has a finite imaginary part"
#endif
    k=0
    do ia=1,size(Cij,1)
       do ib=1,size(Cij,2)
          k=k+1;vals(k)=real(Cij(ia,ib),8)
       enddo
    enddo
  end subroutine correlation_values

end program hubbard_1d
