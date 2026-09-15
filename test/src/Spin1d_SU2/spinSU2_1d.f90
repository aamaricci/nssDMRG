program dmrg_spin_1d
  USE SCIFOR
  USE DMRG
  USE ASSERTING
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=64)                   :: finput
  integer                             :: i,Unit,L
  type(site),dimension(:),allocatable :: myDot
  type(sparse_matrix)                 :: bSz,bSp,SiSj
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable :: Hlr
#else
  real(8),dimension(:,:),allocatable    :: Hlr
#endif
  real(8),dimension(:),allocatable    :: avSz,x,data,data_
  real(8),dimension(:),allocatable    :: e,s
#ifdef _CMPLX
  complex(8)                          :: Sii,Sll,Sll_,Slr,Slr_
#else
  real(8)                             :: Sii,Sll,Sll_,Slr,Slr_
#endif
  integer                             :: irank,comm,rank,ierr
  logical                             :: master

#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call read_input(finput)


  !Init DMRG
  allocate(MyDot(1))              !Homogeneous system
  MyDot = spin_site(sun=2)        !spin_site handles array!
  Hlr = diag([Jp,Jx/2d0])       !implicit fortran allocation
  call init_dmrg(Hlr,ModelDot=MyDot)

  call run_DMRG()

  !Post-processing and measure quantities:
  call Measure_DMRG(myDot(1)%operators%op(key="S"//myDot(1)%okey(0,1,ilink="n")),&
       pos=arange(1,Ldmrg),avOp=avSz)

  !Generic static correlations: same site, same block and across the
  !left/right superblock cut. The reversed correlators test both the
  !growth-order construction and the sector-changing L/R contraction.
  Sii  = Measure_SpinSpin_DMRG(1,1)
  Sll  = Measure_SpinSpin_DMRG(1,2)
  Sll_ = Measure_SpinSpin_DMRG(2,1)
  Slr  = Measure_SpinSpin_DMRG(Ldmrg,Ldmrg+1)
  Slr_ = Measure_SpinSpin_DMRG(Ldmrg+1,Ldmrg)
  if(master)then
#ifdef _CMPLX
     call assert(Sii,cmplx(0.75d0,0d0,8),"<S_i.S_i>",tol=1d-10)
     call assert(Sll,conjg(Sll_),"same-block <S_i.S_j> symmetry",tol=1d-10)
     call assert(Slr,conjg(Slr_),"left/right <S_i.S_j> symmetry",tol=1d-10)
#else
     call assert(Sii,0.75d0,"<S_i.S_i>",tol=1d-10)
     call assert(Sll,Sll_,"same-block <S_i.S_j> symmetry",tol=1d-10)
     call assert(Slr,Slr_,"left/right <S_i.S_j> symmetry",tol=1d-10)
#endif
  endif
  call End_Measure_DMRG()

  if(master)then
     !Check energy:
     L = file_length("energyVSleft.length_L40_M20_iDMRG.dmrg")
     allocate(x(L),data(L))
     call sread("energyVSleft.length_L40_M20_iDMRG.dmrg",x,data)
     L = file_length("energy.check")
     deallocate(x)
     allocate(x(L),data_(L))
     call sread("energy.check",x,data_)
     if(size(data)/=size(data_))stop "Energy files have different sizes"
     call assert(data,data_,"E")
     deallocate(x,data,data_)
     !
     L = file_length("sz.check")
     allocate(data_(L))
     call read_array("sz.check",data_)
     if(size(avSz)/=size(data_))stop "Sz files have different sizes"
     call assert(avSz,data_,"Sz",tol=1d-8)
     deallocate(avSz,data_)
     !
     ! L = file_length("SentropyVSleft.length_L40_M20_iDMRG.dmrg")
     ! allocate(x(L),data(L))
     ! call sread("SentropyVSleft.length_L40_M20_iDMRG.dmrg",x,data)
     ! L = file_length("Sentropy.check")
     ! allocate(data_(L))
     ! call read_array("Sentropy.check",data_)
     ! if(size(data)/=size(data_))stop "Sentropy files have different sizes"
     ! call assert(data,data_,"S")
     ! deallocate(x,data,data_)
     !
  endif


  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif


end program dmrg_spin_1d
