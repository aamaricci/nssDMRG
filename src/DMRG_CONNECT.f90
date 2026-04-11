MODULE DMRG_CONNECT
  USE DMRG_GLOBAL
  implicit none
  private


  public :: enlarge_block
  public :: connect_fermion_blocks
  public :: connect_spin_blocks

contains


  !-----------------------------------------------------------------!
  ! Purpose: enlarge a given BLOCK "SELF" growing it left/right with
  ! a SITE "DOT" (specified in the init)
  !-----------------------------------------------------------------!
  subroutine enlarge_block(self,dot,label,link)
    type(block)                  :: self
    type(site)                   :: dot
    character(len=*),optional    :: label
    character(len=*),optional    :: link
    character(len=1)             :: label_
    character(len=1)             :: ilink,olink
    character(len=:),allocatable :: key,dtype,otype,error_otype
    type(tbasis)                 :: self_basis,dot_basis,enl_basis
    type(sparse_matrix)          :: Hb,Hd,H2,eO
    integer                      :: i,l
    !
    label_='l'; if(present(label))label_=to_lower(str(label(1:1)))
    iLink ="n";if(present(Link))iLink=to_lower(str(Link(1:1)))

    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: ENLARGE block "//to_upper(label_)
#endif
    !
    if(MpiMaster)call start_timer("Enlarge blocks "//to_upper(label_))
    t0=t_start()
    !
    dtype=dot%type()

    select case(label_)
    case default;stop "Enlarge_Block ERROR: label_ not in ['l','r']"
    case ("l","r");continue
    end select
    !
    select case(iLink)
    case default;stop "Enlarge_Block ERROR: iLink not in ['n','p']"
    case ("n","p");continue
    end select
    !
    if(.not.self%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing self.H operator in the list"
    if(.not.dot%operators%has_key("H"))&
         stop "Enlarge_Block ERROR: Missing dot.H operator in the list"
    if(dtype/=self%type())&
         stop "Enlarge_Block ERROR: Dot.Type != Self.Type"
    !
    dtype=to_lower(dtype(1:1))
    select case(dtype)
    case default;stop "Enlarge_Block ERROR: wrong dot.Type"
    case ("s","f","e","b");continue
    end select
    !> Update Hamiltonian:
    if(MpiMaster)then
#ifdef _DEBUG
       write(LOGfile,*)"DEBUG: ENLARGE block: update H"
#endif
       t0=t_start()
       select case(iLink)          
       case ("n")  !OBC/PBC_next L:[o-o...o-]-*; R:*-[-o...o-o]
          select case(label_)
          case ("l")
             Hb = self%operators%op("H").x.id(dot%dim)
             Hd = id(self%dim).x.dot%operators%op("H")
             select case(dtype)
             case ("s")    ;H2 = connect_spin_blocks(self,as_block(dot),link=iLink)
             case ("f","e");H2 = connect_fermion_blocks(self,as_block(dot),link=iLink)
             end select
          case ("r")
             Hb = id(dot%dim).x.self%operators%op("H")
             Hd = dot%operators%op("H").x.id(self%dim)
             select case(dtype)
             case ("s")     ;H2 = connect_spin_blocks(as_block(dot),self,link=iLink)
             case ("f","e") ;H2 = connect_fermion_blocks(as_block(dot),self,link=iLink)
             end select
          end select
       case ("p")  !PBC_prev L:@-[-o...o-o]; R:[o-o...o-]-@
          select case(label_)
          case ("l")
             Hb = id(dot%dim).x.self%operators%op("H")
             Hd = dot%operators%op("H").x.id(self%dim)
             select case(dtype)
             case ("s")     ;H2 = connect_spin_blocks(as_block(dot),self,link=iLink)
             case ("f","e") ;H2 = connect_fermion_blocks(as_block(dot),self,link=iLink)
             end select
          case ("r")
             Hb = self%operators%op("H").x.id(dot%dim)
             Hd = id(self%dim).x.dot%operators%op("H")
             select case(dtype)
             case ("s")    ;H2 = connect_spin_blocks(self,as_block(dot),link=iLink)
             case ("f","e");H2 = connect_fermion_blocks(self,as_block(dot),link=iLink)
             end select
          end select
       end select
       call self%put_op("H", sp_add3(Hb,Hd,H2), type="bosonic")
       write(LOGfile,*)"Build&Put H*",t_stop()
       !
       !
       !
       !> Update all the other operators in the list involved in the enlargement procedure:
#ifdef _DEBUG
       write(LOGfile,*)"DEBUG: ENLARGE block: update Op list"
#endif
       t0=t_start()
       do i=1,size(self%operators)
          key   = str(self%operators%key(index=i)) !respect the case
          if(key=="H")cycle
          otype = to_lower(str(self%operators%type(index=i)))
          otype = otype(1:1)
          olink = to_lower(key(len(key):len(key)))
          !
          error_otype="Enlarge_BLock ERROR: wrong self.operators.type L: !\in['Bosonic','Fermionic','PSign']"
          select case(iLink)
          case ("n")                    !OBC/PBC_next
             select case(label_)
             case ("l")
                select case(olink)
                case ("n")
                   !       [o-o...o-]-@ 
                   !Bosons   : O_L -> i_L.x.O_d
                   !Fermions : O_L -> p_L.x.O_d
                   !FermiSign: O_L -> p_L.x.P_d
                   select case(otype)
                   case default;stop error_otype
                   case('b');eO = Id(self%dim)
                   case('f');eO = self%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   case('p');eO = self%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   end select
                   call self%put_op(key, eO.x.dot%operators%op(key), type=otype)
                case ("p")
                   !        [@-o...o-]-o
                   !Bosons   : O_L -> O_L.x.i_d
                   !Fermions : O_L -> O_L.x.i_d
                   !FermiSign: O_L -> P_L.x.p_d
                   select case(otype)
                   case default;stop error_otype
                   case('b');eO = Id(dot%dim)
                   case('f');eO = Id(dot%dim)
                   case('p');eO = dot%operators%op(key="P"//dot%okey(0,0,ilink=olink))
                   end select
                   call self%put_op(key, self%operators%op(key).x.eO, type=otype )
                end select
             case ("r")
                select case(olink)
                case ("n")
                   !         @-[-o...o-o]
                   !Bosons   : O_R -> O_d.x.i_R
                   !Fermions : O_R -> O_d.x.i_R
                   !FermiSign: O_R -> P_d.x.p_R
                   select case(otype)
                   case default;stop error_otype
                   case('b');eO = Id(self%dim)
                   case('f');eO = Id(self%dim)
                   case('p');eO = self%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   end select
                   call self%put_op(key, dot%operators%op(key).x.eO, type=otype )
                case ("p")
                   !        o-[-o...o-@]
                   !Bosons   : O_R -> i_d.x.O_R
                   !Fermions : O_R -> p_d.x.O_R
                   !FermiSign: O_R -> p_d.x.P_R
                   select case(otype)
                   case default;stop error_otype
                   case('b');eO = Id(dot%dim)
                   case('f');eO = dot%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   case('p');eO = dot%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   end select
                   call self%put_op(key, eO.x.self%operators%op(key), type=otype)
                end select
             end select
          case ("p")                    !PBC_prev
             select case(label_)
             case ("l")
                select case(olink)
                case ("n")
                   !      o-[-o...o-@]
                   !Bosons   : O_L -> i_d.x.O_L
                   !Fermions : O_L -> p_d.x.O_L
                   !FermiSign: P_L -> p_d.x.P_L
                   select case(otype)
                   case default;stop error_otype
                   case('b');eO = Id(dot%dim)
                   case('f');eO = dot%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   case('p');eO = dot%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   end select
                   call self%put_op(key, eO.x.self%operators%op(key), type=otype)
                case ("p")
                   !      @-[-o...o-o]
                   !Bosons   : O_L -> O_d.x.i_L
                   !Fermions : O_L -> O_d.x.i_L
                   !FermiSign: P_L -> P_d.x.p_L
                   select case(otype)
                   case default;stop error_otype
                   case('b');eO = Id(self%dim)
                   case('f');eO = Id(self%dim)
                   case('p');eO = self%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   end select
                   call self%put_op(key, dot%operators%op(key).x.eO, type=otype )
                end select
             case ("r")
                select case(olink)
                case ("n")
                   !     [@-o...o-]-o
                   !Bosons   : O_R -> O_R.x.i_d
                   !Fermions : O_R -> O_R.x.i_d
                   !FermiSign: P_R -> P_R.x.p_d
                   select case(otype)
                   case default;stop error_otype
                   case('b');eO = Id(dot%dim)
                   case('f');eO = Id(dot%dim)
                   case('p');eO = dot%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   end select
                   call self%put_op(key, self%operators%op(key).x.eO, type=otype )
                case ("p")
                   !     [o-o...o-]-@
                   !Bosons   : O_R -> i_R.x.O_d
                   !Fermions : O_R -> p_R.x.O_d
                   !FermiSign: P_R -> p_R.x.P_d
                   select case(otype)
                   case default;stop error_otype
                   case('b');eO = Id(self%dim)
                   case('f');eO = self%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   case('p');eO = self%operators%op(key="P"//self%okey(0,0,ilink=olink))
                   end select
                   call self%put_op(key, eO.x.dot%operators%op(key), type=otype)
                end select
             end select
          end select
       enddo
       write(LOGfile,*)"Build&Put O*",t_stop()
    endif
    !
    !> Enlarge dimensions
    self%length = self%length + 1
    self%Dim    = self%Dim*dot%Dim
    !
    !> Enlarge the basis states
    if(MpiMaster)t0=t_start()
    call self%get_basis(self_basis)
    call dot%get_basis(dot_basis)
    select case(iLink)
    case ("n")  !OBC/PBC_next L:[o-o...o-]-*; R:*-[-o...o-o]
       select case(label_)
       case ("l");enl_basis = (self_basis.o.dot_basis)
       case ("r");enl_basis = (dot_basis.o.self_basis)
       end select
    case ("p")  !PBC_prev L:@-[-o...o-o]; R:[o-o...o-]-@
       select case(label_)
       case ("l");enl_basis = (dot_basis.o.self_basis)
       case ("r");enl_basis = (self_basis.o.dot_basis)
       end select
    end select
    call self%set_basis( basis=enl_basis )
    if(MpiMaster)write(LOGfile,*)"Build basis",t_stop()
    !
    if(MpiMaster)then 
       if(.not.self%is_valid())then
          write(LOGfile,*)"ENLARGE_BLOCK error: not valid enlarged_block "//to_upper(label_)
          call self%show(file="corrupted_enl_block_"//label_//".dat")
          stop
       endif
    endif
    !
    !Free the memory:
    call Hb%free()
    call Hd%free()
    call H2%free()
    call self_basis%free()
    call dot_basis%free()
    call enl_basis%free()
    !
#ifdef _MPI
    call Barrier_MPI(MpiComm)
#endif
    if(MpiMaster)call stop_timer("Enlarge blocks")
    t_enlarge_blocks=t_enlarge_blocks+t_stop()
    t_connect_blocks=0d0  !reset: this is counted in here
    !
  end subroutine enlarge_block



  !H_lr = \sum_{a}h_aa*(C^+_{left,a}@P_left) x C_{right,a}] + H.c.
  function connect_fermion_blocks(left,right,states,link) result(H2)
    type(block)                               :: left
    type(block)                               :: right
    integer,dimension(:),optional             :: states
    character(len=*),optional                 :: link
    character(len=1)                          :: ilink
    type(sparse_matrix),dimension(Nspin*Norb) :: Cl,Cr
    type(sparse_matrix)                       :: Pl,Pr
    type(sparse_matrix)                       :: H2,Hh
    integer                                   :: ispin,iorb,jorb,io,jo
    integer,dimension(2)                      :: Hdims,dleft,dright
    character(len=:),allocatable              :: key
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable     :: Hij
    complex(8)                                :: Tr,Tl
#else
    real(8),dimension(:,:),allocatable        :: Hij
    real(8)                                   :: Tr,Tl
#endif
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: connect Fermion blocks"
#endif
    !
    iLink ="n";if(present(Link))iLink=to_lower(str(Link(1:1)))
    !
    t0=t_start()
    !
    !Hij is shared:
    !Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)
    !
    !> Get H2 dimensions:
    dleft = shape(left%operators)
    dright= shape(right%operators)
    Hdims = dleft*dright
    if(present(states))Hdims = size(states)
    call H2%init(Hdims(1),Hdims(2))
    !
    !
    !FERMION SPECIFIC:
    !>Retrieve operators:
    do ispin=1,Nspin
       do iorb=1,Norb
          key = "C"//left%okey(iorb,ispin,ilink=iLink)
          io  = iorb + (ispin-1)*Norb
          Cl(io) = left%operators%op(key)
          Cr(io) = right%operators%op(key)
       enddo
    enddo
    !
    !
    !>Build H2:
    ! A = (C_{La}x1_R)^+ 
    ! B = (P_LxC_{Rb})
    ! H2 = \sum_{ab} t_{ab}A^+ @B + t*_{ab} B^+ @A
    ! (AxB)@(CxD)=(A@C)x(B@D) + (AxB)^+=(A^+xB^+)
    !    = \sum_{ab} t_{ab} (C^+_{La}@P_L)x(1_R@C_{Rb}) 
    !              + t*_{ab}(P^+_L@C_{La})x(C^+_{Rb}@1_R) 
    !
    !P^+ = P ?
    key = "P"//left%okey(0,0,ilink=iLink)
    Pl  = left%operators%op(key)
    do io=1,Nspin*Norb
       do jo=1,Nspin*Norb
          if(Hij(io,jo)==zero)cycle
          Tr = Hij(io,jo)
#ifdef _CMPLX
          Tl = conjg(Tr) !== H(jo,io) by Hermiticity
#else
          Tl = Tr
#endif
          if(present(states))then
             H2 = H2 + Tr*sp_kron(matmul(Cl(io)%dgr(),Pl),Cr(jo),states)
             H2 = H2 + Tl*sp_kron(matmul(Pl,Cl(io)),Cr(jo)%dgr(),states)
          else
             H2 = H2 + Tr*(matmul(Cl(io)%dgr(),Pl).x.Cr(jo))
             H2 = H2 + Tl*(matmul(Pl,Cl(io)).x.Cr(jo)%dgr())
          endif
       enddo
    enddo
    !
    !> free memory
    call Pl%free
    call Pr%free
    do io=1,Nspin*Norb
       call Cl(io)%free
       call Cr(io)%free
    enddo
    !
    t_connect_blocks=t_connect_blocks+t_stop()
    !
  end function connect_fermion_blocks




  function connect_spin_blocks(left,right,states,link) result(H2)
    type(block)                           :: left
    type(block)                           :: right
    integer,dimension(:),optional         :: states
    character(len=*),optional             :: link
    character(len=1)                      :: ilink
    type(sparse_matrix)                   :: Sl(Nspin)![Sz,Sp]
    type(sparse_matrix)                   :: Sr(Nspin)![Sz,Sp]
    type(sparse_matrix)                   :: H2
    integer,dimension(2)                  :: Hdims
    integer                               :: ispin
#ifdef _CMPLX
    complex(8),dimension(:,:),allocatable :: Hij
#else
    real(8),dimension(:,:),allocatable    :: Hij
#endif
    !
#ifdef _DEBUG
    if(MpiMaster)write(LOGfile,*)"DEBUG: connect Spin blocks"
#endif
    !
    iLink ="n";if(present(Link))iLink=to_lower(str(Link(1:1)))  
    !
    t0=t_start()
    !
    !Hij is shared:
    !Hij = Hmodel(left,right)
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij, source=HopH)    
    !
    !> Get H2 dimensions:
    Hdims = shape(left%operators)*shape(right%operators)
    if(present(states))Hdims = [size(states),size(states)]
    call H2%init(Hdims(1),Hdims(2))
    !
    !>Retrieve operators:
    do ispin=1,Nspin
       Sl(ispin) = left%operators%op("S"//left%okey(iorb=0,ispin=ispin,ilink=iLink))
       Sr(ispin) = right%operators%op("S"//right%okey(iorb=0,ispin=ispin,ilink=iLink))
    enddo
    !
    !>Build H2:
    if(present(states))then
       H2 = H2 + Hij(1,1)*sp_kron(Sl(1),Sr(1),states) + &
            Hij(2,2)*sp_kron(Sl(2),Sr(2)%dgr(),states)+ &
            Hij(2,2)*sp_kron(Sl(2)%dgr(),Sr(2),states)
    else
       H2 = H2 + Hij(1,1)*(Sl(1).x.Sr(1)) + &
            Hij(2,2)*(Sl(2).x.Sr(2)%dgr())+ &
            Hij(2,2)*(Sl(2)%dgr().x.Sr(2))
    endif
    !
    !> Free memory
    do ispin=1,Nspin
       call Sl(ispin)%free
       call Sr(ispin)%free
    enddo
    !
    t_connect_blocks=t_connect_blocks+t_stop()
    !
  end function connect_spin_blocks



END MODULE DMRG_CONNECT
