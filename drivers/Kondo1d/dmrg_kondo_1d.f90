program dmrg_kondo_1d
  USE SCIFOR, only: parse_cmd_variable,parse_input_variable,init_MPI,&
       StartMsg_MPI,get_Rank_MPI,get_Master_MPI,finalize_MPI,zero,one
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none

  character(len=64)                       :: finput
  integer                                 :: comm,rank,unit,unit_zz
  integer                                 :: i,ispin,Nsites,impurity_site
  logical                                 :: master,irun,imeasure
  real(8)                                 :: ts,Jz,Jxy,Vk
  integer,dimension(:),allocatable        :: growth_to_physical
  type(site)                              :: bath_dot,kondo_dot
  type(site),dimension(:),allocatable     :: ModelDot
  type(sparse_matrix)                     :: Cup,Cdown,Nup,Ndown
  type(sparse_matrix)                     :: Se_z,Se_p,Se_m
  type(sparse_matrix)                     :: Simp_z,Simp_p,Simp_m
  type(sparse_matrix)                     :: Kondo_Se_z,Kondo_Se_p,Kondo_Se_m
  type(sparse_matrix)                     :: Nimp_op
  real(8),dimension(2)                    :: dq0,dqplus,dqminus
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable   :: Hlr
  complex(8)                              :: correlation,correlation_zz
#else
  real(8),dimension(:,:),allocatable      :: Hlr
  real(8)                                 :: correlation,correlation_zz
#endif
  !
  ! Initialize MPI and parse model-specific input.
  master=.true.
#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank=get_Rank_MPI(comm)
  master=get_Master_MPI(comm)
#endif
  call parse_cmd_variable(finput,'FINPUT',default='DMRG.conf')
  call parse_input_variable(irun,'IRUN',finput,default=.true.,&
       comment='Run the DMRG calculation')
  call parse_input_variable(imeasure,'IMEASURE',finput,default=.true.,&
       comment='Measure the impurity-bath spin correlation')
  call parse_input_variable(ts,'TS',finput,default=-1d0,&
       comment='Uniform nearest-neighbour bath hopping')
  call parse_input_variable(Jz,'JZ',finput,default=1d0,&
       comment='Longitudinal Kondo exchange')
  call parse_input_variable(Jxy,'JXY',finput,default=1d0,&
       comment='Transverse Kondo exchange')
  call parse_input_variable(Vk,'VK',finput,default=0d0,&
       comment='Local conduction-electron potential at the Kondo site')
  call parse_input_variable(impurity_site,'KONDO_SITE',finput,default=0,&
       comment='Physical impurity position; zero selects the left PBC seed')
  call read_input(finput)
  !
  ! Validate the vanilla fixed-impurity contract.
  if(to_lower(DMRGtype)/='i')&
       error stop "dmrg_kondo_1d: DMRGTYPE must be 'i'"
  if(.not.PBCdmrg)error stop 'dmrg_kondo_1d: PBCDMRG must be true'
  if(Nspin/=2)error stop 'dmrg_kondo_1d: Nspin must be 2'
  if(Norb/=1)error stop 'dmrg_kondo_1d: NORB must be 1'
  if(to_lower(trim(dmrg_mode))/='normal')&
       error stop "dmrg_kondo_1d: DMRG_MODE must be 'normal'"
  if(QNdim/=2)error stop 'dmrg_kondo_1d: QNDIM must be 2'
  if(any(abs(QN_density-[0.5d0,0.5d0])>1d-12))&
       error stop 'dmrg_kondo_1d: QN_DENSITY must be 0.5,0.5'
  if(any(abs(QN_offset-[1d0,0d0])>1d-12))&
       error stop 'dmrg_kondo_1d: QN_OFFSET must be 1,0'
  !
  ! Build ModelDot in physical chain order. The default impurity position is
  ! the first left seed selected by the PBC growth geometry.
  Nsites=fNlat
  growth_to_physical=Build_DMRG_Growth_Map(Nsites,PBCdmrg)
  if(impurity_site==0)impurity_site=growth_to_physical(1)
  if(impurity_site<1.or.impurity_site>Nsites)&
       error stop 'dmrg_kondo_1d: KONDO_SITE outside the physical chain'
  if(impurity_site/=growth_to_physical(1))&
       error stop 'dmrg_kondo_1d: fixed impurity must occupy the left PBC seed'
  bath_dot=electron_site()
  kondo_dot=kondo_fixed_site(Jz,Jxy,Vk,Simp_z=Simp_z,Simp_p=Simp_p,&
       Sele_z=Kondo_Se_z,Sele_p=Kondo_Se_p,Nimp_op=Nimp_op)
  allocate(ModelDot(Nsites))
  do i=1,Nsites
     ModelDot(i)=bath_dot
  enddo
  ModelDot(impurity_site)=kondo_dot
  !
  ! Construct the uniform bath hopping. Bond-dependent hopping is deliberately
  ! postponed to the Wilson-chain extension.
  allocate(Hlr(Nspin*Norb,Nspin*Norb))
  Hlr=zero
  do ispin=1,Nspin
     Hlr(ispin,ispin)=ts*one
  enddo
  if(imeasure)then
     save_block=.true.
     save_umat=.true.
  endif
  call init_dmrg(Hlr,ModelDot=ModelDot)
  !
  ! Run iDMRG and measure the Kondo cloud.
  if(irun)call run_DMRG()
  if(imeasure)then
     Cup=bath_dot%operators%op(&
          'C'//bath_dot%okey(1,1,ilink='n'))
     Cdown=bath_dot%operators%op(&
          'C'//bath_dot%okey(1,2,ilink='n'))
     Nup=matmul(Cup%dgr(),Cup)
     Ndown=matmul(Cdown%dgr(),Cdown)
     Se_z=0.5d0*(Nup-Ndown)
     Se_p=matmul(Cup%dgr(),Cdown)
     Se_m=Se_p%dgr()
     Simp_m=Simp_p%dgr()
     Kondo_Se_m=Kondo_Se_p%dgr()
     dq0=[0d0,0d0]
     dqplus=[1d0,-1d0]
     dqminus=[-1d0,1d0]
     call Init_Measure_DMRG('Kondo cloud')
     if(master)then
        unit=fopen('SimpDotSeVSsite.dat',append=.false.)
        unit_zz=fopen('SimpSzSeSzVSsite.dat',append=.false.)
     endif
     do i=1,Nsites
        if(i==impurity_site)then
           correlation_zz=Measure_Corr_DMRG(Simp_z,dq0,Kondo_Se_z,dq0,&
                impurity_site,i)
           correlation=correlation_zz
           correlation=correlation+0.5d0*&
                Measure_Corr_DMRG(Simp_p,dqplus,Kondo_Se_m,dqminus,&
                impurity_site,i)
           correlation=correlation+0.5d0*&
                Measure_Corr_DMRG(Simp_m,dqminus,Kondo_Se_p,dqplus,&
                impurity_site,i)
        else
           correlation_zz=Measure_Corr_DMRG(Simp_z,dq0,Se_z,dq0,&
                impurity_site,i)
           correlation=correlation_zz
           correlation=correlation+0.5d0*&
                Measure_Corr_DMRG(Simp_p,dqplus,Se_m,dqminus,&
                impurity_site,i)
           correlation=correlation+0.5d0*&
                Measure_Corr_DMRG(Simp_m,dqminus,Se_p,dqplus,&
                impurity_site,i)
        endif
#ifdef _CMPLX
        if(master)write(unit,*)i,real(correlation,8),aimag(correlation)
        if(master)write(unit_zz,*)i,real(correlation_zz,8),aimag(correlation_zz)
#else
        if(master)write(unit,*)i,correlation
        if(master)write(unit_zz,*)i,correlation_zz
#endif
     enddo
     if(master)then
        close(unit)
        close(unit_zz)
     endif
     call End_Measure_DMRG()
  endif
  !
  ! Release DMRG and local driver data.
  call finalize_dmrg()
  do i=1,size(ModelDot)
     call ModelDot(i)%free()
  enddo
  call bath_dot%free()
  call kondo_dot%free()
  call Cup%free()
  call Cdown%free()
  call Nup%free()
  call Ndown%free()
  call Se_z%free()
  call Se_p%free()
  call Se_m%free()
  call Simp_z%free()
  call Simp_p%free()
  call Simp_m%free()
  call Kondo_Se_z%free()
  call Kondo_Se_p%free()
  call Kondo_Se_m%free()
  call Nimp_op%free()
#ifdef _MPI
  call finalize_MPI()
#endif


end program dmrg_kondo_1d
