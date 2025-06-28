! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! SET UP THE MESH AT BEGINNING OF THE RUN
      module Mesh
      use globals
      use physical_constants
      use arrays
      use RandomNumbers
      use general_functions
      use atomic_physics
      implicit none

      real(8) :: TotMass,TotNi56
      real(8) :: deltav(3)
      integer , allocatable , save :: ind(:,:,:)
      real(8) , allocatable , save :: fracNi56(:,:,:)
      real(8) , allocatable , save :: v1(:),v2(:),v3(:)  !!! grid velocity
      real(8) , allocatable , save :: tau(:)
      integer , allocatable , save :: pellets(:,:,:)


      
      contains

      subroutine init_mesh
      real(8) :: rho,vol,vmin(3)
      integer :: i,j,k,n,ieff,jeff,keff,dilute,z,a
      real(8) :: tmpmass,tmptotmass,tmpiso(50)
      real(8) :: Mni0,Mfe0,Mime0,ekin,dmix
      real(8) :: Mej,Vej,Rej,Vejmax,fac,vejfacmax,ni_fac
      integer :: nmix,i1,i2
      integer :: ejecta_type,mesh_type
      logical :: ismix
      character(50) :: filename
      real(8) :: isomass(max_isotops)
!     1 - constant density
!     2 - exponential
!     3 - read from file
!     10 - read w7
      integer :: ino

      namelist /mesh/ dim,nc1,nc2,nc3,ejecta_type,Mni0,Mfe0,Mime0,&
                      ekin,vej,vejmax,mesh_type,ismix,dmix,nmix,vejfacmax,ni_fac

!!    first define all isotops
      call add_full_isotop_list

      ejecta_type=1
      mesh_type=1
      mej=1.38d0
      ekin=1.1780d0
      Vej=1.d9
      mfe0=0.1d0
      mni0=0.7d0
      mime0=0.3d0
      vejmax=-1.0d0
      vejfacmax=10.0d0
      ni_fac=1d0

      ismix=.false.
      dmix=0.02d0
      nmix=3
 
      open(unit=5,file=data_file)
      read(5,nml=mesh,iostat=ino)
      close(5)

      ekin=ekin*1d51
      mej=mej*solar_mass
      if (ejecta_type.eq.2) then
        if (ekin.gt.0.0d0) vej=sqrt(ekin/mej*2.0d0*exp_dist_F2(vejfacmax)/exp_dist_F4(vejfacmax))
        vejmax=vejfacmax*vej
      endif

      dmix=dmix*solar_mass

      write(fout,nml=mesh)

!     dilute=4
      deltav=0.0d0
!     nc1=320/dilute
!     nc2=640/dilute
      nctot=nc1*nc2*nc3
      nc1p=nc1+1
      nc2p=nc2+1
      nc3p=nc3+1

      n=1
      allocate(ind(nc1,nc2,nc3))
      do k=1,nc3
        do j=1,nc2
          do i=1,nc1
            ind(i,j,k)=n
            n=n+1
          enddo
        enddo
      enddo
           
      allocate(atoms(niso,nctot))
      allocate(mass(nctot))
      allocate(rhov(nctot))

      allocate(fracNi56(nc1,nc2,nc3))
      allocate(pellets(nc1,nc2,nc3))
      allocate(v1(nc1p),v2(nc2p),v3(nc3p))
      allocate(tau(nc1p))
      mass=0.0d0
      atoms=0.0d0
      rhov=0.0d0
      fracNi56=0.0d0
      pellets=0
      v1=0.0d0
      v2=0.0d0
      v3=0.0d0

      TotMass=0.0d0
      TotNi56=0.0d0

      if (dim.eq.1) then

!     set density distribution
      if (ejecta_type.eq.1) then   !! rho=const
        rho=mej/(4.0d0*pi/3.0d0*Vej**3.0d0) !rho defined actually as M/(V/t^3), so to get real density must devide by t^3.
        totmass=0.0d0
        do i=1,nc1
          if (mesh_type.eq.1) then         ! constant dv per shell
            v1(i+1)=Vej/dble(nc1)*dble(i)
          elseif (mesh_type.eq.2) then     ! constant dm per shell
            mass(ind(i,1,1))=mej/dble(nc1)
            totmass=totmass+mass(ind(i,1,1))
            v1(i+1)=(3.0d0/4.0d0/pi*totmass/rho)**(1.0d0/3.0d0)
          endif
          vol=4.0d0*pi/3.0d0*(v1(i+1)**3.0d0-v1(i)**3.0d0)  !vol defined as V/t^3, so must multiply by t^3 to get real volume
          mass(ind(i,1,1))=rho*vol
          rhov(ind(i,1,1))=1.0d0/vol
        enddo
      elseif (ejecta_type.eq.2) then !! rho=exp(-v/vej)
        do i=1,nc1
          if (mesh_type.eq.1) then 
! rho is defined as ro_0*t_0^3, which is equal to M(all)/(4*pi*v_e^3*F2(zmax)), and then m(z)=4*pi*rho*ve^3*(F2(z(i+1))-F2(z(i)))
            rho=mej/(4.0d0*pi*vej**3.0d0*(exp_dist_F2(vejmax/vej)-exp_dist_F2(0.0d0)))
            v1(i+1)=Vejmax/dble(nc1)*dble(i)
            mass(ind(i,1,1))=4.0d0*pi*rho*vej**3.0d0*(exp_dist_F2(v1(i+1)/vej)-exp_dist_F2(v1(i)/vej))
          elseif (mesh_type.eq.2) then
! here rho defined differently such that rho*F2(z(i)) gives the mass up to the end of the i shell.
            rho=mej/exp_dist_F2(vejmax/vej)
            mass(ind(i,1,1))=mej/dble(nc1)
            v1(i+1)=vej*solve_exp_dist_F2(mej*dble(i)/dble(nc1)/rho)
          endif
          vol=4.0d0*pi/3.0d0*(v1(i+1)**3.0d0-v1(i)**3.0d0)
          rhov(ind(i,1,1))=1.0d0/vol
        enddo
      elseif (ejecta_type.eq.3) then !! read from file
        open(unit=10,file='ejecta')
        do i=1,nc1p
          if (i.le.nc1) then
            read(10,*) v1(i),mass(ind(i,1,1)),fracni56(i,1,1)
          else
            read(10,*) v1(i)
          endif
        enddo
        do i=1,nc1
          vol=4.0d0*pi/3.0d0*(v1(i+1)**3.0d0-v1(i)**3.0d0)
          rhov(ind(i,1,1))=1.0d0/vol
        enddo
      endif

      if (ejecta_type.lt.10) then
!     set isotop distribution

      totmass=0.0d0
      do n=1,nc1
        i=ind(n,1,1)
        totmass=totmass+mass(i)
        if (ejecta_type.le.2) then
          fac=totmass/solar_mass
! Lucy profile
!         if (fac.le.0.5d0) then
!           atoms(ind_ni56,i)=1.0d0
!         elseif (totmass/solar_mass.le.0.75d0) then
!           atoms(ind_ni56,i)=min(1.0d0,(max(0.0d0,(0.75d0-fac)/0.25d0)))
!         else
!           atoms(ind_ni56,i)=0.0d0
!         endif
!         atoms(ind_si28,i)=1.0d0-atoms(ind_ni56,i)
! Woosley & Kasen
          if (fac.le.mfe0) then
            atoms(indiso(26,54),i)=1.0d0     ! iron
          elseif(fac.gt.mfe0 .and. fac.le.(mfe0+mni0)) then
            atoms(indiso(28,56),i)=1.0d0*ni_fac     ! nickel
            atoms(indiso(14,28),i)=0.535d0*(1d0-ni_fac)   ! Si
            atoms(indiso(16,32),i)=0.32d0*(1d0-ni_fac)    ! S
            atoms(indiso(18,36),i)=0.062d0*(1d0-ni_fac)   ! Ar
            atoms(indiso(20,40),i)=0.083d0*(1d0-ni_fac)   ! Ca
          elseif(fac.gt.(mfe0+mni0) .and. fac.le.(mfe0+mni0+mime0)) then
            atoms(indiso(14,28),i)=0.535d0   ! Si
            atoms(indiso(16,32),i)=0.32d0    ! S
            atoms(indiso(18,36),i)=0.062d0   ! Ar
            atoms(indiso(20,40),i)=0.083d0   ! Ca
          else
            atoms(indiso(6,12),i)=0.5d0      ! C
            atoms(indiso(8,16),i)=0.5d0      ! O
          endif
        endif
        atoms(:,i)=mass(i)*atoms(:,i)/iso(1:niso)%A
      enddo

      if (ismix) call mix1d(dmix,nmix)

      endif


      if (ejecta_type.eq.10) then !! read profile from file.
      filename='XXX'
! write the rest of this part according to the file format.
      endif

      elseif (dim.eq.2) then

!     AA(1)=m_h1
!     ZZ(1)=1.0d0
!     AA(2)=m_he3
!     ZZ(2)=2.0d0
!     AA(3)=m_he4
!     ZZ(3)=2.0d0
!     AA(4)=m_c12
!     ZZ(4)=6.0d0
!     AA(5)=m_n14
!     ZZ(5)=7.0d0
!     AA(6)=m_o16
!     ZZ(6)=8.0d0
!     AA(7)=m_ne20
!     ZZ(7)=10.0d0
!     AA(8)=m_mg24
!     ZZ(8)=12.0d0
!     AA(9)=m_si28
!     ZZ(9)=14.0d0
!     AA(10)=m_s32
!     ZZ(10)=16.0d0
!     AA(11)=m_ar36
!     ZZ(11)=18.0d0
!     AA(12)=m_ca40
!     ZZ(12)=20.0d0
!     AA(13)=m_ti44
!     ZZ(13)=22.0d0
!     AA(14)=m_cr48
!     ZZ(14)=24.0d0
!     AA(15)=m_fe52
!     ZZ(15)=26.0d0
!     AA(16)=m_fe54
!     ZZ(16)=26.0d0
!     AA(17)=m_ni56
!     ZZ(17)=28.0d0
      ind_ni56=17

      deltav(1)=1.226334d7
      deltav(2)=1.226334d7
      vmin(1)=0.0d0
      vmin(2)=-3.918139d9-deltav(2)/2.0d0

      do i=1,nc1p
        v1(i)=vmin(1)+deltav(1)*dilute*dble(i-1)
      enddo
      do i=1,nc2p
        v2(i)=vmin(2)+deltav(2)*dilute*dble(i-1)
      enddo


      
      open(unit=10,file='ejecta.dat')
      do n=1,320*640
        read(10,*) i,j,tmpmass,tmpiso(1:niso)
        ieff=(i-1)/dilute+1
        jeff=(j-1)/dilute+1
        tmptotmass=sum(tmpiso(1:niso))
        if (abs(tmpmass-tmptotmass)/tmpmass.gt.1.d-6) print*,'prob',abs(tmpmass-tmptotmass)/tmpmass

        tmpmass=sum(tmpiso(1:niso))
!       mass(ieff,jeff,1)=mass(ieff,jeff,1)+tmpmass
!       atoms(1:niso,ieff,jeff,1)=atoms(1:niso,ieff,jeff,1)+tmpiso(1:niso)/AA(1:niso)

        totmass=totmass+tmpmass

      enddo
      close(10)

      do j=1,nc2
        do i=1,nc1
        vol=pi*(v1(i+1)**2.0d0-v1(i)**2.0d0)*(v2(j+1)-v2(j))
        if (vol.lt.0.0d0) print*,'neg vol !!!'
        rhov(ind(i,j,1))=1.0d0/vol
!       fracNi56(i,j,1)=atoms(ind_ni56,i,j,1)*AA(ind_ni56)/mass(i,j,1)
!       TotNi56=TotNi56+mass(i,j,1)*fracNi56(i,j,1)
        enddo
      enddo

      endif


      isomass=0.0d0
      do i=1,niso
        isomass(i)=sum(atoms(i,:)*iso(i)%a)/solar_mass
      enddo

      totmass=sum(isomass(:))
      TotNi56=sum(atoms(indiso(28,56),:))*iso(indiso(28,56))%A
      write(fout,1001) totmass

      call reduce_isotop_list(isomass(:))

      do i=1,niso
        z=iso(i)%z
        a=nint(iso(i)%a*avogadro)
        write(fout,1002) indiso(z,a),z,a,sum(atoms(i,:)*iso(i)%a)/solar_mass
      enddo
1001  format('Total ejecta mass = ',1pe14.6,' solar mass')
1002  format('Index=',I4,' Total isotop (Z,A)=',2I4,' mass = ',&
              1pe14.6,' solar mass')
      ind_fe56=indiso(26,56)
      ind_co56=indiso(27,56)
      ind_ni56=indiso(28,56)
      write(fout,1005) ind_fe56,ind_co56,ind_ni56
1005  format('radioactive isotop indexs: Fe56=',I4,' Co56=',I4,' Ni56=',I4)


      if (dim.eq.1) then
        open(unit=50,file='ejecta_profile',POSITION='APPEND')
        write(50,1003) v1(:)
        write(50,1003) 0.0d0,mass(1:nc1)
        do n=1,niso
          write(50,1004) iso(n)%sym,atoms(n,1:nc1)*iso(n)%A/mass(1:nc1)
        enddo
        close(50)
      endif

1003  format(1000(1pe14.6))
1004  format(I5,1000(1pe14.6))

      return
      end subroutine init_mesh

      subroutine findijk(r,t,i,j,k)
      real(8) , intent(in) :: r(3),t
      integer , intent(out) :: i,j,k
      real(8) :: v(3) , vnorm
      integer :: ierr

      i=1
      j=1
      k=1

      v(:)=vofr(r,t)

      if (dim.eq.1) then
        vnorm=norm(v)
        i=find_index(vnorm,v1,ierr)
        if (ierr.ne.0) print*,'ierr prob'
      endif

      return
      end subroutine findijk


      function random_location(i,j,k,t)
      integer , intent(in) :: i,j,k
      real(8) , intent(in) :: t
      real(8) :: random_location(3)
      real(8) :: z(3),v

      call random_number(z)

      if (dim.eq.1) then

        v=(v1(i)**3.0d0*z(1)+v1(i+1)**3.0d0*(1.0d0-z(1)))**(1.0d0/3.0d0)
        random_location=v*random_unit_vec1(3)

      elseif (dim.eq.2) then

        v=(v1(i)**2.0d0*z(1)+v1(i+1)**2.0d0*(1.0d0-z(1)))**(1.0d0/2.0d0)
        random_location(1:2)=v*random_unit_vec1(2)
        random_location(3)=v2(j)*z(2)+v2(j+1)*(1.0d0-z(2))

      elseif (dim.eq.3) then

        random_location(1)=v1(i)*z(1)+v1(i+1)*(1.0d0-z(1))
        random_location(2)=v2(j)*z(2)+v2(j+1)*(1.0d0-z(2))
        random_location(3)=v3(k)*z(3)+v3(k+1)*(1.0d0-z(3))

      endif

      random_location=rofv(random_location,t)

      return
      end function random_location

      function rofv(v,t)
      real(8) , intent(in) :: v(3) , t
      real(8) :: rofv(3)

      rofv=v(:)*t

      return
      end function rofv

      function vofr(r,t)
      real(8) , intent(in) :: r(3) , t
      real(8) :: vofr(3)

      vofr=r(:)/t

      return
      end function vofr

      function rhooft(rho,t)
      real(8) , intent(in) :: rho , t
      real(8) :: rhooft

      rhooft=rho/t**3.0d0

      return
      end function rhooft


      function exp_dist_f2(z)
!     function returns integral from 0 to z of
!     z^2*exp(-z)
      real(8) , intent(in) :: z
      real(8) :: exp_dist_f2

      exp_dist_f2=2.0d0-(z**2.0d0+2.0d0*z+2.0d0)*exp(-z)

      return
      end function exp_dist_f2

      function solve_exp_dist_f2(f2)
!     function returns integral from 0 to z of
!     z^2*exp(-z)
      real(8) , intent(in) :: f2
      real(8) :: solve_exp_dist_f2
      real(8) :: conv,z,zold
      integer :: niter

      conv=1.0d0
      niter=0
      z=1.0d0
      do while (conv.gt.1.d-12 .and. niter.lt.1000)
      zold=z

      z=-log(-(f2-2.0d0)/(z**2.0d0+2.0d0*z+2.0d0))
      conv=abs(z-zold)/zold
      niter=niter+1

      enddo
      solve_exp_dist_f2=z

      return
      end function solve_exp_dist_f2

      function exp_dist_f4(z)
!     function returns integral from 0 to z of
!     z^4*exp(-z)
      real(8) , intent(in) :: z
      real(8) :: exp_dist_f4

      exp_dist_f4=24.0d0-(z**4.0d0+4.0d0*(z**3.0d0+3.0d0*(z**2.0d0+2.0d0*(z+1.0d0))))*exp(-z)

      return
      end function exp_dist_f4

      subroutine mix1d(dmix,nmix)
      real(8) , intent(in) :: dmix
      integer , intent(in) :: nmix
      real(8) , allocatable :: mixfrac(:),summass(:)
      real(8) :: mix1,mix2,tota
      integer :: nm,n,i,j,i1,i2

      allocate(mixfrac(nctot))
      allocate(summass(nc1p))

      summass=0.0d0

      do i=2,nctot
        summass(i)=summass(i-1)+mass(i-1)
      enddo

   
      do nm=1,nmix

      do n=1,nc1
!       mix1=dble(n-1)*dmix
!       mix2=dble(n)*dmix
        mix1=summass(n)
        mix2=summass(n)+dmix
        if (mix2.gt.summass(nctot)) cycle

        mixfrac=0.0d0
        i1=2
        do while (summass(i1).lt.mix1)
          i1=i1+1
        enddo
        i1=i1-1

        i2=2
        do while (summass(i2).lt.mix2 .and. i2.lt.nc1p)
          i2=i2+1
        enddo
        i2=i2-1
   
        if (i1.eq.i2) cycle

        mixfrac(i1)=(summass(i1+1)-mix1)/mass(i1)
        mixfrac(i2)=(mix2-summass(i2))/mass(i2)
        if (i2-i1.gt.1) mixfrac(i1+1:i2-1)=1.0d0


        do j=1,niso
          tota=sum(atoms(j,:)*iso(j)%A*mixfrac(:))
          atoms(j,:)=(atoms(j,:)*iso(j)%A-atoms(j,:)*iso(j)%A*mixfrac(:)+&
            tota*mixfrac(:)*mass(:)/sum(mixfrac(:)*mass(:)))/iso(j)%A
        enddo

      enddo

      enddo
        
      deallocate(mixfrac)
      deallocate(summass)

      return
      end subroutine mix1d

      subroutine reduce_isotop_list(totmass)
      real(8) :: totmass(max_isotops)
      real(8) , allocatable :: tempatoms(:,:)
      type (isotop_data) , allocatable :: tempiso(:)
      integer :: n,nisonew,i1,i2,a,z

      nisonew=0

!!    make sure fe56,co56,ni56 always present
      totmass(indiso(26,56))=1.d-10
      totmass(indiso(27,56))=1.d-10
      totmass(indiso(28,56))=1.d-10

      do n=1,niso
        if (totmass(n).gt.0.0d0) nisonew=nisonew+1
      enddo

      allocate(tempatoms(nisonew,nctot))
      allocate(tempiso(nisonew))
      tempatoms=0.0d0

      nisonew=0
      do n=1,niso
        z=iso(n)%z
        a=nint(iso(n)%a*avogadro)
        if (totmass(n).eq.0.0d0) then
          indiso(z,a)=-1
        else
          nisonew=nisonew+1
          tempatoms(nisonew,:)=atoms(n,:)
          tempiso(nisonew)=iso(n)
          indiso(z,a)=nisonew
        endif
      enddo
      niso=nisonew
      deallocate(atoms)
      deallocate(iso)
      allocate(atoms(niso,nctot))
      allocate(iso(niso))
      atoms(:,:)=tempatoms(:,:)
      iso(:)=tempiso(:)
      deallocate(tempatoms)
      deallocate(tempiso)

      return
      end subroutine reduce_isotop_list

      end module Mesh
