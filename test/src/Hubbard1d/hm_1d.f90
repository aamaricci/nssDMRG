program hubbard_1d
  USE SCIFOR
  USE DMRG
  USE ASSERTING
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,unit,iorb,ispin,L
  real(8)                                        :: ts(2),Mh(2)
  type(site),dimension(:),allocatable            :: MyDot
  real(8),dimension(:,:),allocatable             :: Hloc
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable          :: Hlr
#else
  real(8),dimension(:,:),allocatable             :: Hlr
#endif
  type(sparse_matrix),dimension(:,:),allocatable :: Nop,Cop
  type(sparse_matrix),dimension(:),allocatable   :: dens,docc,s2z
  type(sparse_matrix)                            :: C,Cdg
  real(8),dimension(:),allocatable               :: dqC
  real(8),dimension(:,:),allocatable             :: avO
  real(8),dimension(:),allocatable               :: x,data,data_
  real(8),dimension(:),allocatable               :: n,d,m,e,s
#ifdef _CMPLX
  complex(8)                                     :: Gii,Gll,Gll_,Glr,Glr_,Alr
#else
  real(8)                                        :: Gii,Gll,Gll_,Glr,Glr_,Alr
#endif
  real(8)                                        :: nii
  integer                             :: irank,comm,rank,ierr
  logical                             :: master=.true.

#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=(/( -1d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(Mh,"MH",finput,default=(/(0d0,i=1,2 )/),&
       comment="Crystal field splittings")
  call read_input(finput)


  Nso = Nspin*Norb

  ! allocate(Hloc(Nso,Nso))
  Hloc = diag([Mh(1:Norb),Mh(1:Norb)])


  allocate(MyDot(1))
  MyDot = electron_site()

  if(allocated(Hlr))deallocate(Hlr)
  allocate(Hlr(Nso,Nso))
  Hlr = diag([ts(1:Norb),ts(1:Norb)])

  call init_dmrg(Hlr,ModelDot=MyDot)

  !Run DMRG algorithm
  call run_DMRG()


  !Post-processing and measure quantities:
  !Measure <Sz(i)>
  allocate(Cop(Norb,Nspin),Nop(Norb,Nspin))
  do ispin=1,Nspin
     do iorb=1,Norb
        Cop(iorb,ispin) = myDot(1)%operators%op(key="C"//myDot(1)%okey(iorb,ispin))
        Nop(iorb,ispin) = matmul(Cop(iorb,ispin)%dgr(),Cop(iorb,ispin))
     enddo
  enddo
  allocate(dens(Norb),docc(Norb),s2z(Norb))
  do iorb=1,Norb
     dens(iorb) = Nop(iorb,1)+Nop(iorb,2)
     docc(iorb) = matmul(Nop(iorb,1),Nop(iorb,2))
     s2z(iorb)  = matmul((Nop(iorb,1)-Nop(iorb,2)),(Nop(iorb,1)-Nop(iorb,2)))
  enddo


  call Measure_DMRG([dens,docc,s2z],pos=arange(1,Ldmrg),avOp=avO)


  !Odd-fermion correlations. Each two-point function is parity even,
  !but its representation contains the Jordan-Wigner string between
  !the two endpoints. Test equal-site composition, hermiticity and the
  !anticommutation sign, both within one block and across the L/R cut.
  C   = myDot(1)%operators%op(key="C"//myDot(1)%okey(1,1))
  Cdg = C%dgr()
  dqC = myDot(1)%operators%dq(key="C"//myDot(1)%okey(1,1))
  Gii = Measure_Corr_DMRG(Cdg,-dqC,C,dqC,1,1,"fermionic","fermionic")
  nii = Measure_Op_DMRG(Nop(1,1),1)
  Gll = Measure_Corr_DMRG(Cdg,-dqC,C,dqC,1,2,"fermionic","fermionic")
  Gll_= Measure_Corr_DMRG(Cdg,-dqC,C,dqC,2,1,"fermionic","fermionic")
  Glr = Measure_Corr_DMRG(Cdg,-dqC,C,dqC,Ldmrg,Ldmrg+1,"fermionic","fermionic")
  Glr_= Measure_Corr_DMRG(Cdg,-dqC,C,dqC,Ldmrg+1,Ldmrg,"fermionic","fermionic")
  Alr = Measure_Corr_DMRG(C,dqC,Cdg,-dqC,Ldmrg+1,Ldmrg,"fermionic","fermionic")
  if(master)then
#ifdef _CMPLX
     call assert(Gii,cmplx(nii,0d0,8),"equal-site <Cdg.C>",tol=1d-10)
     call assert(Gll,conjg(Gll_),"same-block fermion hermiticity",tol=1d-10)
     call assert(Glr,conjg(Glr_),"left/right fermion hermiticity",tol=1d-10)
#else
     call assert(Gii,nii,"equal-site <Cdg.C>",tol=1d-10)
     call assert(Gll,Gll_,"same-block fermion hermiticity",tol=1d-10)
     call assert(Glr,Glr_,"left/right fermion hermiticity",tol=1d-10)
#endif
     call assert(Alr,-Glr,"left/right fermion anticommutation",tol=1d-10)
  endif
  call C%free()
  call Cdg%free()
  call End_Measure_DMRG()


  if(master)then
     !
     call save_array("n.out",avO(1,:))
     call save_array("d.out",avO(2,:))
     call save_array("s2z.out",avO(3,:))
     !
     !Check energy:
     L = file_length("energyVSleft.length_L40_M40_iDMRG.dmrg")
     allocate(x(L),data(L))
     call sread("energyVSleft.length_L40_M40_iDMRG.dmrg",x,data)
     L = file_length("energy.check")
     allocate(data_(L))
     call read_array("energy.check",data_)
     if(size(data)/=size(data_))stop "Energy files have different sizes"
     call assert(data,data_,"E",1d-6)
     deallocate(x,data,data_)
     !
     i = 1
     L = file_length("n.check")
     allocate(data_(L))
     call read_array("n.check",data_)
     if(size(avO,2)/=size(data_))stop "N files have different sizes"
     call assert(avO(i,:),data_,"N",1d-6)
     deallocate(data_)
     !
     i = 2
     L = file_length("d.check")
     allocate(data_(L))
     call read_array("d.check",data_)
     if(size(avO,2)/=size(data_))stop "N files have different sizes"
     call assert(avO(i,:),data_,"D",1d-6)
     deallocate(data_)
     !
     i = 3
     L = file_length("s2z.check")
     allocate(data_(L))
     call read_array("s2z.check",data_)
     if(size(avO,2)/=size(data_))stop "S2z files have different sizes"
     call assert(avO(i,:),data_,"S2z",1d-6)
     deallocate(data_)
     !
     ! L = file_length("SentropyVSleft.length_L40_M40_iDMRG.dmrg")
     ! allocate(x(L),data(L))
     ! call sread("SentropyVSleft.length_L40_M40_iDMRG.dmrg",x,data)
     ! L = file_length("Sentropy.check")
     ! allocate(data_(L))
     ! call read_array("Sentropy.check",data_)
     ! if(size(data)/=size(data_))stop "Sentropy files have different sizes"
     ! call assert(data,data_,"S",1d-6)
     ! deallocate(x,data,data_)
     !
  endif



  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif


end program hubbard_1d




