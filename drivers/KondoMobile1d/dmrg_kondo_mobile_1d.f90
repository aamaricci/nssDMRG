program dmrg_kondo_mobile_1d
  USE SCIFOR, only: parse_cmd_variable,parse_input_variable,init_MPI,&
       StartMsg_MPI,get_Rank_MPI,get_Master_MPI,finalize_MPI,zero,one
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none

  character(len=64)                       :: finput
  integer                                 :: comm,rank
  integer                                 :: i,j,ispin,length,Nimp_length,Ntotal
  integer                                 :: Nimp,unit_nimp,unit_ie,unit_ii
  integer                                 :: unit_ie_zz,unit_ii_zz
  logical                                 :: master,irun,imeasure,icorrelations
  real(8)                                 :: ts,timp,Jz,Jxy,Vk,Eimp,pinning
  real(8)                                 :: Eimp_site,nimp_value
  integer,dimension(:),allocatable        :: pinning_sites
  integer,dimension(:),allocatable        :: growth_to_physical
  real(8),dimension(:,:),allocatable      :: TargetQNs
  real(8),dimension(3)                    :: dq0,dqplus,dqminus
  type(site)                              :: mobile_dot
  type(site),dimension(:),allocatable     :: ModelDot
  type(sparse_matrix)                     :: Simp_z,Simp_p,Simp_m
  type(sparse_matrix)                     :: Sele_z,Sele_p,Sele_m,Nimp_op
#ifdef _CMPLX
  complex(8),dimension(:,:),allocatable   :: Hlr
  complex(8)                              :: corr_ie,corr_ii
  complex(8)                              :: corr_ie_zz,corr_ii_zz
#else
  real(8),dimension(:,:),allocatable      :: Hlr
  real(8)                                 :: corr_ie,corr_ii
  real(8)                                 :: corr_ie_zz,corr_ii_zz
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
       comment='Measure the local impurity occupation')
  call parse_input_variable(icorrelations,'ICORRELATIONS',finput,&
       default=.true.,comment='Measure the full spatial spin correlations')
  call parse_input_variable(ts,'TS',finput,default=-1d0,&
       comment='Uniform conduction-electron hopping')
  call parse_input_variable(timp,'TIMP',finput,default=-1d0,&
       comment='Uniform projected impurity hopping')
  call parse_input_variable(Jz,'JZ',finput,default=1d0,&
       comment='Longitudinal Kondo exchange')
  call parse_input_variable(Jxy,'JXY',finput,default=1d0,&
       comment='Transverse Kondo exchange')
  call parse_input_variable(Vk,'VK',finput,default=0d0,&
       comment='Local conduction-electron potential')
  call parse_input_variable(Eimp,'EIMP',finput,default=0d0,&
       comment='Local impurity energy')
  call parse_input_variable(Nimp,'NIMP',finput,default=1,&
       comment='Global number of projected impurity electrons')
  call parse_input_variable(pinning,'PINNING',finput,default=0d0,&
       comment='Energy lowering on the selected impurity sites')
  if(Nimp>0)then
     allocate(pinning_sites(Nimp))
     call parse_input_variable(pinning_sites,'PINNING_SITES',finput,&
          default=[(i,i=1,Nimp)],comment='Physical sites lowered by PINNING')
  endif
  call read_input(finput)
  !
  ! Validate the mobile Kondo-lattice contract.
  if(to_lower(DMRGtype)/='i')&
       error stop "dmrg_kondo_mobile_1d: DMRGTYPE must be 'i'"
  if(.not.PBCdmrg)&
       error stop 'dmrg_kondo_mobile_1d: PBCDMRG must be true'
  if(Nspin/=2)error stop 'dmrg_kondo_mobile_1d: Nspin must be 2'
  if(Norb/=2)error stop 'dmrg_kondo_mobile_1d: NORB must be 2'
  if(to_lower(trim(dmrg_mode))/='normal')&
       error stop "dmrg_kondo_mobile_1d: DMRG_MODE must be 'normal'"
  if(QNdim/=3)error stop 'dmrg_kondo_mobile_1d: QNDIM must be 3'
  if(Nimp<0.or.Nimp>fNlat)&
       error stop 'dmrg_kondo_mobile_1d: NIMP must lie in [0,fNlat]'
  if(Nimp==0.and.pinning/=0d0)&
       error stop 'dmrg_kondo_mobile_1d: PINNING requires NIMP > 0'
  if(Nimp>0.and.pinning/=0d0)then
     if(any(pinning_sites<1).or.any(pinning_sites>fNlat))&
          error stop 'dmrg_kondo_mobile_1d: PINNING_SITES outside the chain'
     do i=1,Nimp
        if(count(pinning_sites==pinning_sites(i))/=1)&
             error stop 'dmrg_kondo_mobile_1d: PINNING_SITES must be distinct'
     enddo
  endif
  if(master.and.Nimp>0.and.timp==0d0.and.pinning==0d0)then
     write(LOGfile,*)&
          'dmrg_kondo_mobile_1d WARNING: TIMP=0 without pinning leaves ',&
          'disconnected impurity configurations.'
  endif
  !
  ! Build a consistent warm-up schedule.  The conduction band is currently
  ! fixed at half filling, Nc(L)=L.  A later extension may instead derive its
  ! spin-resolved filling from QN_DENSITY(1:2).  Use CEILING so that every
  ! impurity species present in the final system already occurs in the first
  ! compatible warm-up superblock and can spread over the renormalized blocks.
  allocate(TargetQNs(QNdim,fNlat))
  growth_to_physical=Build_DMRG_Growth_Map(fNlat,PBCdmrg)
  do length=1,fNlat
     if(Nimp>0.and.timp==0d0.and.pinning/=0d0)then
        !Static pinned impurities enter with their physical sites.
        Nimp_length=0
        do i=1,length
           if(any(pinning_sites==growth_to_physical(i)))&
                Nimp_length=Nimp_length+1
        enddo
     else
        Nimp_length=(length*Nimp+fNlat-1)/fNlat
     endif
     Ntotal=length+Nimp_length
     TargetQNs(1,length)=dble((Ntotal+1)/2)
     TargetQNs(2,length)=dble(Ntotal/2)
     TargetQNs(3,length)=dble(Nimp_length)
  enddo
  !
  ! Construct the physical site array.  Positive PINNING lowers the local
  ! impurity energy on the selected sites by Eimp(i)=Eimp-PINNING.
  mobile_dot=kondo_mobile_site(Jz,Jxy,Vk,Eimp,Simp_z=Simp_z,&
       Simp_p=Simp_p,Sele_z=Sele_z,Sele_p=Sele_p,Nimp_op=Nimp_op)
  allocate(ModelDot(fNlat))
  do i=1,fNlat
     Eimp_site=Eimp
     if(Nimp>0.and.pinning/=0d0)then
        if(any(pinning_sites==i))Eimp_site=Eimp-pinning
     endif
     ModelDot(i)=kondo_mobile_site(Jz,Jxy,Vk,Eimp_site)
  enddo
  !
  ! Construct diagonal hopping for the conduction and projected impurity
  ! orbitals.
  allocate(Hlr(Nspin*Norb,Nspin*Norb))
  Hlr=zero
  do ispin=1,Nspin
     Hlr(1+(ispin-1)*Norb,1+(ispin-1)*Norb)=ts*one
     Hlr(2+(ispin-1)*Norb,2+(ispin-1)*Norb)=timp*one
  enddo
  if(imeasure)then
     save_block=.true.
     save_umat=.true.
  endif
  call init_dmrg(Hlr,ModelDot=ModelDot,TargetQNs=TargetQNs)
  !
  ! Run iDMRG and measure the impurity distribution and spin correlations.
  if(irun)call run_DMRG()
  if(imeasure)then
     Simp_m=Simp_p%dgr()
     Sele_m=Sele_p%dgr()
     dq0=0d0
     dqplus=[1d0,-1d0,0d0]
     dqminus=-dqplus
     call Init_Measure_DMRG('Mobile Kondo observables')
     if(master)then
        unit_nimp=fopen('NimpVSsite.dat',append=.false.)
        if(icorrelations)then
           unit_ie=fopen('SimpDotSeleVSsites.dat',append=.false.)
           unit_ie_zz=fopen('SimpSzSeleSzVSsites.dat',append=.false.)
           unit_ii=fopen('SimpDotSimpVSsites.dat',append=.false.)
           unit_ii_zz=fopen('SimpSzSimpSzVSsites.dat',append=.false.)
        endif
     endif
     do i=1,fNlat
        nimp_value=Measure_Op_DMRG(Nimp_op,i)
        if(master)write(unit_nimp,*)i,nimp_value
        if(.not.icorrelations)cycle
        do j=1,fNlat
           corr_ie_zz=Measure_Corr_DMRG(Simp_z,dq0,Sele_z,dq0,i,j)
           corr_ie=corr_ie_zz
           corr_ie=corr_ie+0.5d0*&
                Measure_Corr_DMRG(Simp_p,dqplus,Sele_m,dqminus,i,j)
           corr_ie=corr_ie+0.5d0*&
                Measure_Corr_DMRG(Simp_m,dqminus,Sele_p,dqplus,i,j)
           corr_ii_zz=Measure_Corr_DMRG(Simp_z,dq0,Simp_z,dq0,i,j)
           corr_ii=corr_ii_zz
           corr_ii=corr_ii+0.5d0*&
                Measure_Corr_DMRG(Simp_p,dqplus,Simp_m,dqminus,i,j)
           corr_ii=corr_ii+0.5d0*&
                Measure_Corr_DMRG(Simp_m,dqminus,Simp_p,dqplus,i,j)
#ifdef _CMPLX
           if(master)then
              write(unit_ie,*)i,j,real(corr_ie,8),aimag(corr_ie)
              write(unit_ie_zz,*)i,j,real(corr_ie_zz,8),aimag(corr_ie_zz)
              write(unit_ii,*)i,j,real(corr_ii,8),aimag(corr_ii)
              write(unit_ii_zz,*)i,j,real(corr_ii_zz,8),aimag(corr_ii_zz)
           endif
#else
           if(master)then
              write(unit_ie,*)i,j,corr_ie
              write(unit_ie_zz,*)i,j,corr_ie_zz
              write(unit_ii,*)i,j,corr_ii
              write(unit_ii_zz,*)i,j,corr_ii_zz
           endif
#endif
        enddo
     enddo
     if(master)then
        close(unit_nimp)
        if(icorrelations)then
           close(unit_ie)
           close(unit_ie_zz)
           close(unit_ii)
           close(unit_ii_zz)
        endif
     endif
     call End_Measure_DMRG()
  endif
  !
  ! Release DMRG and local driver data.
  call finalize_dmrg()
  do i=1,size(ModelDot)
     call ModelDot(i)%free()
  enddo
  call mobile_dot%free()
  call Simp_z%free()
  call Simp_p%free()
  call Simp_m%free()
  call Sele_z%free()
  call Sele_p%free()
  call Sele_m%free()
  call Nimp_op%free()
#ifdef _MPI
  call finalize_MPI()
#endif


end program dmrg_kondo_mobile_1d
