! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! SUBROUTINES FOR READING THE ATOMIC DATA (LINES LISTS) AND FOR
! CALCULATING IONIZATION WITH SAHA EQUATION
      Module Atomic_Physics
      use globals , only : niso , fout , data_file , max_isotops
      use arrays , only : iso , indiso
      use physical_constants
      use general_functions
      implicit none

      integer , parameter :: max_ion_levels=8
      integer , parameter :: max_atoms=99

      type level_data
        real(8) :: g=0.0d0
        real(8) :: energy=1.d90
      end type level_data

      type line_data
        integer :: Z=0
        integer :: ion=0
        real(8) :: lambda=0.0d0
        integer :: l1,l2
        real(8) :: fij
      end type line_data

      type ion_data
        integer :: nlevels=0
        integer :: nlines=0
        real(8) , allocatable :: partition(:)
        real(8) , allocatable :: ptemp(:)
        type (level_data) , allocatable :: level(:)
      end type ion_data

      type atom_data
        logical :: active=.false.
        integer :: Z=0
        real(8) :: ionization_energy(0:max_ion_levels)=0.0d0
        type (ion_data) :: ion(0:max_ion_levels)
      end type atom_data

      type (atom_data) , save :: atom(max_atoms)

      type (line_data) , allocatable , save :: line(:)


      character(20) , save :: linelist='kuruczCD23'
      character(20) , save :: pfdata='kuruczCD23'
      integer , save :: nlines
      real(8) , save :: etherm=1.0d0
      real(8) , save :: gfcut=-20d0

      contains

      subroutine init_atomic_data
      integer :: i,j,k,n,ino,i1,i2
      integer :: eof,ZZ,LL,NN,nlevels
      real(8) :: gg,en
      namelist /atomic_physics/ etherm,linelist,pfdata,gfcut

      open(unit=5,file=data_file)
      read(5,nml=atomic_physics,iostat=ino)
      close(5)

      do i=1,max_atoms
        atom(i)%z=i
      enddo

      do i=1,niso
        atom(iso(i)%z)%active=.true.
      enddo

      write(fout,nml=atomic_physics)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                       !!
!!    Atomic data is taken from the following databases:                 !!
!!                                                                       !!
!!    Ionization data: from 'atom_data' of the publicly available        !!
!!    TARDIS code written by WE Kerzendorf and SA Sim , arXiv:1401.5469  !!
!!    and SA Sim , arXiv:1401.5469                                       !!
!!                                                                       !!
!!    Line list and level data from Kurucz gfall.dat file                !!
!!                                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      open (unit=10,file='XXX/YYY/ionization_data')
      eof=0
      do while (eof.eq.0)
        read(10,*,IOSTAT=eof) ZZ,LL,EN
        if (eof.eq.0) then
          if (zz.le.max_atoms .and. ll.le.max_ion_levels) then
            atom(zz)%ionization_energy(ll)=en*electron_volt
          endif
        endif
      enddo
      close (10)

      if (pfdata.eq.'kuruczCD23') then
        call prepare_partition_functions(1)
      elseif (pfdata.eq.'kurucz_table') then
        call prepare_partition_functions(1)
        call prepare_partition_functions(2)
      else
        print*,'no partition function method defined'
        stop
      endif
      

      call read_kurucz

      return
      end subroutine init_atomic_data

      subroutine SahaIonization(ni,zi,temp,ne,nij,partition,zavg)
      real(8) , intent(in) :: ni(:),temp
      integer , intent(in) :: zi(:)
      real(8) , intent(out) :: ne,zavg
      real(8) , dimension (0:max_ion_levels,size(ni)) , intent(out) :: nij,partition
      real(8) :: A(1:max_ion_levels,size(ni))
      real(8) :: phi(1:max_ion_levels,size(ni))
      real(8) :: tempmin,eps,g,en,kbt,parte,conv,nenew,convfac,neold,tmax
      integer :: niso,z,i,j,k,max1i,niter,nlevels

!!    subroutine calculates the ionization equilibrium statr for a mixture of isotops with
!!    atomic density ni, atomic number zi and temperature temp and returns the electron density ne and average
!!    ionization level <z>=ne/sum(ni)
!!
!!    Solution method is based on the paper "Ionization Equilibrium Equation of State. II. Mixtures",
!!    by Rouse, C.A. Astrophysical Journal, vol. 135, p.599.
!!
!!    No attempt is made at this time to optimize solution for execution time

      eps=1.0d-12
      tempmin=1000.0d0
      niso=size(ni)


      tmax=max(temp,tempmin)
      kbt=kboltz*tmax
      parte=2.0d0*(2.0d0*pi*electron_mass*kbt)**1.5d0/planck**3.0d0
      partition=0.0d0
      phi=0.0d0

      do i=1,niso
        if (ni(i).lt.eps) cycle
        do j=0,max_ion_levels
!     step 1: calculate partition functions for all isotops ionization levels
          partition(j,i)=calc_partition(tmax,zi(i),j,2)
!     step 2: calcualte phi(i,j)=n(i,j+1)/n(i,j)*ne
          if (j.gt.0) then
            en=atom(zi(i))%ionization_energy(j)
            phi(j,i)=parte*partition(j,i)/(partition(j-1,i)+eps)*exp(-en/kbt)
          endif
        enddo
      enddo

!     step 3: iterate to find nij

!     first guess - assume all ions are fully ionized
      ne=sum(ni) 

      max1i=100       !! maximum number of iterations
      convfac=0.5d0   !! mixing factor between old and new value of ne

      niter=0
      conv=1.0d0
      do while (niter.lt.max1i .and. conv.gt.1.d-6)
      neold=ne
      nenew=0.0d0
      nij=0.0d0
        do i=1,niso
          if (ni(i).lt.eps) cycle
!     step 3.1 - calculate netural ion n(i,j=0) density from 
!                conservation equation ni(i)=sum(nij(i,j=0...max_ion_levels))
          A(1,i)=phi(1,i)/ne
          do j=2,max_ion_levels
            A(j,i)=A(j-1,i)*phi(j,i)/ne
          enddo
          nij(0,i)=ni(i)/(1.0d0+sum(A(:,i)))
!     step 3.2 - calculate netural ion n(i,j=0) density from 
!                conservation equation ni(i)=sum(nij(i,j=0...max_ion_levels))
          nij(1:max_ion_levels,i)=A(1:max_ion_levels,i)*nij(0,i)
!     step 3.3 - calculate contribution to electron density
          do j=1,max_ion_levels
            nenew=nenew+dble(j)*nij(j,i)
          enddo
        enddo

      ne=convfac*nenew+(1.0d0-convfac)*neold
      conv=abs(ne-neold)/(neold+eps)
      niter=niter+1
      enddo

      if (niter.ge.max1i) then
        print*,'saha problem - iterations not converging'
        print*,'conv,temp=',conv,max(temp,tempmin)
      endif

      ne=max(ne,eps)
      zavg=ne/sum(ni(:))

      return
      end subroutine SahaIonization

      subroutine prepare_partition_functions(method)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Subroutine prepares partition functions tables as     !!
!!    a function of temperatures for all relevant ions      !!
!!    using one of several methos:                          !!
!!                                                          !!
!!    method = 1 - analyticaly from kurucz CD 23 line list  !!
!!    method = 2 - ...                                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer , intent(in) :: method
      integer :: i,j,i1,i2,n,totlines,eof
      character*50 :: filename
      real(8) :: lambda,fij,elo,ehi,glo,ghi
      real(8) :: temp,part(0:5)
      integer :: z,ion
      type (level_data) :: lvl
      type (level_data) , allocatable :: tmplvl(:)
      character*10 :: label1,label2
      integer :: nloc,npart(max_atoms),nlevels
      logical :: add,zlist(max_atoms)


      if (method.eq.1) then

      filename='XXX/YYY/gfall.dat'
      zlist=.false.

!!!   temporary allocation for big level array per isotop
      do z=1,max_atoms
        if (.not.atom(z)%active) cycle
        zlist(z)=.true.
        do ion=0,max_ion_levels
          allocate(atom(z)%ion(ion)%partition(100))
          allocate(atom(z)%ion(ion)%ptemp(100))
          allocate(atom(z)%ion(ion)%level(5000))
        enddo
      enddo

!!    read energy levels and sort
      open (unit=10,file=filename)
      nloc=1
      nlines=0
      eof=0
      do while (eof.eq.0)
        call get_kurucz_23_next_line(10,nloc,zlist(:),z,ion,lambda,fij,elo,glo,ehi,ghi,eof)
        if (eof.eq.0) then
          nlines=nlines+1
          lvl%energy=elo
          lvl%g=glo
          call add_level(lvl,atom(z)%ion(ion)%level(:),atom(z)%ion(ion)%nlevels)
          lvl%energy=ehi
          lvl%g=ghi
          call add_level(lvl,atom(z)%ion(ion)%level(:),atom(z)%ion(ion)%nlevels)
        endif
      enddo

      close(10)

!!    calculate partition functions
      do z=1,max_atoms
        if (.not.atom(z)%active) cycle
        do ion=0,max_ion_levels
          do n=1,100
            atom(z)%ion(ion)%ptemp(n)=1000.0d0*dble(n)
            atom(z)%ion(ion)%partition(n)=calc_partition(atom(z)%ion(ion)%ptemp(n),z,ion,1)
          enddo
          deallocate(atom(z)%ion(ion)%level)
        enddo
      enddo

      elseif (method.eq.2) then

      filename='XXX/YYY/partfnall'

!!!   temporary allocation for big level array per isotop
      do z=1,max_atoms
        if (.not.atom(z)%active) cycle
        do ion=0,5
          if (allocated(atom(z)%ion(ion)%partition)) then
            deallocate(atom(z)%ion(ion)%ptemp)
            deallocate(atom(z)%ion(ion)%partition)
          endif
          allocate(atom(z)%ion(ion)%ptemp(201))
          allocate(atom(z)%ion(ion)%partition(201))
        enddo
      enddo

!!    read energy levels and sort
      open (unit=10,file=filename)
      eof=0
      npart=0
      do while (eof.eq.0)
        read(10,*,iostat=eof) z,temp,part(0:5)
        if (eof.eq.0 .and. atom(z)%active) then
          npart(z)=npart(z)+1
          do ion=0,5
            atom(z)%ion(ion)%ptemp(npart(z))=temp
            atom(z)%ion(ion)%partition(npart(z))=part(ion)
          enddo
        endif
      enddo

      close(10)

      endif

      return
      end subroutine prepare_partition_functions

      subroutine read_kurucz
      integer :: i,j,i1,i2,n,totlines,eof
      character*50 :: filename1,filename2
      real(8) :: lambda(2),fij(2),elo(2),ehi(2),glo(2),ghi(2)
      integer :: zz,z(2),ion(2),nloc(2),lentape(2)
      type (level_data) :: lvl
      type (level_data) , allocatable :: tmplvl(:)
      integer :: nlevels
      integer :: nlist1,nlist2
      logical , dimension(max_atoms) :: zlist1,zlist2
      integer , dimension(max_atoms) :: wlist2
      real(8) :: facread
      logical :: readall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Subroutine reads kurucz CD23 linelist and builds      !!
!!    sorted energy level list per ion , partition function !!
!!    per ion for predetermined temperature grid and        !!
!!    combind linelist with references to the energy level  !!
!!    indices.                                              !!
!!                                                          !!
!!                                                          !!
!!    opt = 1 - read and build the complete lists           !!
!!    opt = 2 - build only partition functions              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      filename1='XXX/YYY/gfall.dat'
      filename2='XXX/YYYlowhighlines'
      lentape(1)=564986
      lentape(2)=31627892+10268732

!!!   temporary allocation for big level array per isotop
      do zz=1,max_atoms
        if (.not.atom(zz)%active) cycle
        do i=0,max_ion_levels
          allocate(atom(zz)%ion(i)%level(5000))
        enddo
      enddo

      nlist1=0
      nlist2=0
      zlist1(:)=.false.
      zlist2(:)=.false.
      wlist2(:)=0
      do zz=1,max_atoms
        if (.not.atom(zz)%active) cycle
        select case (linelist)
          case ('kuruczCD23')
            nlist1=nlist1+1
            zlist1(zz)=.true.
          case ('kuruczCD1')
            nlist2=nlist2+1
            zlist2(zz)=.true.
            wlist2(zz)=1
          case ('kuruczCD1+23')
            nlist1=nlist1+1
            zlist1(zz)=.true.
            nlist2=nlist2+1
            zlist2(zz)=.true.
            wlist2(zz)=-1
        end select 
      enddo

!!    read energy levels and sort
      open (unit=11,file=filename1)
      open (unit=12,file=filename2)
      nloc=1
      if (nlist1.gt.0) then
      call get_kurucz_23_next_line(11,nloc(1),zlist1(:),z(1),ion(1),&
                    lambda(1),fij(1),elo(1),glo(1),ehi(1),ghi(1),eof)
      else
        lambda(1)=1.d99
        nloc(1)=lentape(1)+1
      endif
      if (nlist2.gt.0) then
      call get_kurucz_1_next_line(12,nloc(2),zlist2(:),&
                         wlist2(:),z(2),ion(2),&
                       lambda(2),fij(2),elo(2),glo(2),ehi(2),ghi(2),eof)
      else
        lambda(2)=1.d99
        nloc(2)=lentape(2)+1
      endif

      readall=.false.
      nlines=0
      facread=0.0d0
      do while (.not.readall)
        if (lambda(1).lt.lambda(2)) then
          i=1
        else
          i=2
        endif
        nlines=nlines+1
        atom(z(i))%ion(ion(i))%nlines=atom(z(i))%ion(ion(i))%nlines+1
        lvl%energy=elo(i)
        lvl%g=glo(i)
        call add_level(lvl,atom(z(i))%ion(ion(i))%level(:),atom(z(i))%ion(ion(i))%nlevels)
!       lvl%energy=ehi
!       lvl%g=ghi
!       call add_level(lvl,atom(z)%ion(ion)%level(:),atom(z)%ion(ion)%nlevels)
        if (i.eq.1) call get_kurucz_23_next_line(11,nloc(1),zlist1(:),z(1),ion(1),&
                    lambda(1),fij(1),elo(1),glo(1),ehi(1),ghi(1),eof)
        if (i.eq.2) call get_kurucz_1_next_line(12,nloc(2),zlist2(:),&
                         wlist2(:),z(2),ion(2),&
                       lambda(2),fij(2),elo(2),glo(2),ehi(2),ghi(2),eof)
        if (eof.ne.0) lambda(i)=1.d99

        if (minval(lambda(:)).gt.1.d90) readall=.true.

        if (dble(sum(nloc(:)-1))/dble(sum(lentape)).gt.facread) then
          write(fout,1000) 100.0d0*facread,nlines
          print*,'%=',facread
          facread=facread+0.01d0
        endif
      enddo
1000  format ('1st read, completed ',F6.2,'%, nlines=',I8)

!!    resize level arrays
      allocate(tmplvl(5000))
      do zz=1,max_atoms
        if (.not.atom(zz)%active) cycle
        write(fout,*) 'z=',zz
        do i=0,max_ion_levels
          nlevels=atom(zz)%ion(i)%nlevels
          tmplvl(1:nlevels)=atom(zz)%ion(i)%level(1:nlevels)
          deallocate(atom(zz)%ion(i)%level)
          allocate(atom(zz)%ion(i)%level(nlevels))
          atom(zz)%ion(i)%level(1:nlevels)=tmplvl(1:nlevels)
          write(fout,1001) i,nlevels,atom(zz)%ion(i)%nlines
        enddo
      enddo
1001  format('ion=',I2,' nlevels=',I6,' nlines=',I8)
      deallocate(tmplvl)

      write(fout,*) 'total lines=',nlines
      
      allocate(line(nlines))

!!    read the line list
      rewind(11)
      rewind(12)
      nloc=1
      if (nlist1.gt.0) then
      call get_kurucz_23_next_line(11,nloc(1),zlist1(:),z(1),ion(1),&
                    lambda(1),fij(1),elo(1),glo(1),ehi(1),ghi(1),eof)
      else
        lambda(1)=1.d99
      endif
      if (nlist2.gt.0) then
      call get_kurucz_1_next_line(12,nloc(2),zlist2(:),&
                         wlist2(:),z(2),ion(2),&
                       lambda(2),fij(2),elo(2),glo(2),ehi(2),ghi(2),eof)
      else
        lambda(2)=1.d99
      endif


      eof=0
      j=0
      facread=0.0d0
      readall=.false.
      do while (.not.readall)
        if (lambda(1).lt.lambda(2)) then
          i=1
        else
          i=2
        endif
        j=j+1
        line(j)%z=z(i)
        line(j)%ion=ion(i)
        line(j)%lambda=lambda(i)

        lvl%energy=elo(i)
        lvl%g=glo(i)
        i1=find_level(lvl,atom(z(i))%ion(ion(i))%level(:),atom(z(i))%ion(ion(i))%nlevels)
        line(j)%l1=i1

!       lvl%energy=ehi(i)
!       lvl%g=ghi(i)
!       i2=find_level(lvl,atom(z(i))%ion(ion(i))%level(:),atom(z(i))%ion(ion(i))%nlevels)
!       line(j)%l2=i2

        line(j)%fij=fij(i)
        if (i.eq.1) call get_kurucz_23_next_line(11,nloc(1),zlist1(:),z(1),ion(1),&
                    lambda(1),fij(1),elo(1),glo(1),ehi(1),ghi(1),eof)
        if (i.eq.2) call get_kurucz_1_next_line(12,nloc(2),zlist2(:),&
                         wlist2(:),z(2),ion(2),&
                       lambda(2),fij(2),elo(2),glo(2),ehi(2),ghi(2),eof)
        if (eof.ne.0) lambda(i)=1.d99

        if (minval(lambda(:)).gt.1.d90) readall=.true.
!         write(fout,1001) i,line(i)%lambda/angstrom,zz,ion,line(i)%fij,line(i)%l1,line(i)%l2
        if (dble(sum(nloc(:)-1))/dble(sum(lentape)).gt.facread) then
          write(fout,1100) 100.0d0*facread,j,nlines
          facread=facread+0.01d0
        endif
1100  format ('2nd read, completed ',F6.2,'%, nlines=',I8,' out of ',I8)

      enddo

      close(10)

      return
      end subroutine read_kurucz

      subroutine get_kurucz_23_next_line(nfile,nline,zlist,z,ion,lambda,&
                                         fij,elo,glo,ehi,ghi,eof)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Subroutine reads kurucz CD23 line list and returns    !!
!!    next line according to the rules:                     !!
!!    1. line for isotop with given z in zlist(:)           !!
!!    2. ...                                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer , intent(in) :: nfile
      logical , intent(in) :: zlist(max_atoms)
      integer , intent(inout) :: nline
      integer , intent(out) :: z,ion,eof
      real(8) , intent(out) :: lambda,fij,elo,ehi,glo,ghi
      integer :: i1,i2,n
      real(8) :: lam,gf,zzii,en(2),j(2)
      character*10 :: label1,label2
      logical :: found

      eof=0
      found=.false.
      do while (eof.eq.0 .and. .not.found)
        read(nfile,1000,IOSTAT=eof) lam,gf,zzii,en(1),j(1),label1,en(2),j(2),label2

        if (eof.ne.0) cycle
        nline=nline+1

        z=nint(zzii)
        ion=nint(100.0d0*(zzii-z))
        if (.not.zlist(z)) cycle
        if (ion.gt.max_ion_levels) cycle
        if (gf.lt.gfcut) cycle
        lambda=lam*1d-7 !! translate from nm
        en(:)=planck*clight*abs(en(:)) !! translate from cm-1 to erg
        i1=1
        i2=2
        if (en(2).lt.en(1)) then
          i1=2
          i2=1
        endif
        elo=en(i1)
        ehi=en(i2)
        glo=2.0d0*j(i1)+1.0d0
        ghi=2.0d0*j(i2)+1.0d0
        fij=10.0d0**gf/glo
        found=.true.
      enddo
1000  format(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10)

      return
      end subroutine get_kurucz_23_next_line

      subroutine get_kurucz_1_next_line(nfile,nline,zlist,wlist,z,ion,lambda,&
                                        fij,elo,glo,ehi,ghi,eof)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Subroutine reads kurucz CD1 big line list and returns !!
!!    next line according to the rules:                     !!
!!    1. line for isotop with given z in zlist(:)           !!
!!    2. ...                                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer , intent(in) :: nfile,wlist(max_atoms)
      logical , intent(in) :: zlist(max_atoms)
      integer , intent(inout) :: nline
      integer , intent(out) :: z,ion,eof
      real(8) , intent(out) :: lambda,fij,elo,ehi,glo,ghi
      integer(4) :: iwl
      integer(2) :: ielion,ielo,igflog,igr,igs,igw
      real(8) :: ratiolog
      real(8) :: wlvac,code,facread,gf
      real(8) , save :: elem(840)
      integer :: n,nelion,nread(2),lentape
      logical :: found
      data elem/&
        1.00, 1.01,&
        2.00, 2.01, 2.02,&
        3.00, 3.01, 3.02, 3.03,&
        4.00, 4.01, 4.02, 4.03, 4.04,&
        5.00, 5.01, 5.02, 5.03, 5.04, 5.05,&
        6.00, 6.01, 6.02, 6.03, 6.04, 6.05, 6.06,&
        7.00, 7.01, 7.02, 7.03, 7.04, 7.05, 7.06, 7.07,&
        8.00, 8.01, 8.02, 8.03, 8.04, 8.05, 8.06, 8.07, 8.08,&
        9.00, 9.01, 9.02, 9.03, 9.04, 9.05, 9.06, 9.07, 9.08, 9.09,&
       10.00,10.01,10.02,10.03,10.04,10.05,10.06,10.07,10.08,10.09,&
       10.10,&
       11.00,11.01,11.02,11.03,11.04,11.05,11.06,11.07,11.08,11.09,&
       11.10,11.11,&
       12.00,12.01,12.02,12.03,12.04,12.05,12.06,12.07,12.08,12.09,&
       12.10,12.11,12.12,&
       13.00,13.01,13.02,13.03,13.04,13.05,13.06,13.07,13.08,13.09,&
       13.10,13.11,13.12,13.13,&
       14.00,14.01,14.02,14.03,14.04,14.05,14.06,14.07,14.08,14.09,&
       14.10,14.11,14.12,14.13,14.14,&
       15.00,15.01,15.02,15.03,15.04,15.05,15.06,15.07,15.08,15.09,&
       15.10,15.11,15.12,15.13,15.14,15.15,&
       16.00,16.01,16.02,16.03,16.04,16.05,16.06,16.07,16.08,16.09,&
       16.10,16.11,16.12,16.13,16.14,16.15,16.16,&
       17.00,17.01,17.02,17.03,17.04,17.05,17.06,17.07,17.08,17.09,&
       17.10,17.11,17.12,17.13,17.14,17.15,17.16,17.17,&
       18.00,18.01,18.02,18.03,18.04,18.05,18.06,18.07,18.08,18.09,&
       18.10,18.11,18.12,18.13,18.14,18.15,18.16,18.17,18.18,&
       19.00,19.01,19.02,19.03,19.04,19.05,19.06,19.07,19.08,19.09,&
       19.10,19.11,19.12,19.13,19.14,19.15,19.16,19.17,19.18,19.19,&
       20.00,20.01,20.02,20.03,20.04,20.05,20.06,20.07,20.08,20.09,&
       20.10,20.11,20.12,20.13,20.14,20.15,20.16,20.17,20.18,20.19,&
       20.20,&
       21.00,21.01,21.02,21.03,21.04,21.05,21.06,21.07,21.08,21.09,&
       21.10,21.11,21.12,21.13,21.14,21.15,21.16,21.17,21.18,21.19,&
       21.20,21.21,&
       22.00,22.01,22.02,22.03,22.04,22.05,22.06,22.07,22.08,22.09,&
       22.10,22.11,22.12,22.13,22.14,22.15,22.16,22.17,22.18,22.19,&
       22.20,22.21,22.22,&
       23.00,23.01,23.02,23.03,23.04,23.05,23.06,23.07,23.08,23.09,&
       23.10,23.11,23.12,23.13,23.14,23.15,23.16,23.17,23.18,23.19,&
       23.20,23.21,23.22,23.23,&
       24.00,24.01,24.02,24.03,24.04,24.05,24.06,24.07,24.08,24.09,&
       24.10,24.11,24.12,24.13,24.14,24.15,24.16,24.17,24.18,24.19,&
       24.20,24.21,24.22,24.23,24.24,&
       25.00,25.01,25.02,25.03,25.04,25.05,25.06,25.07,25.08,25.09,&
       25.10,25.11,25.12,25.13,25.14,25.15,25.16,25.17,25.18,25.19,&
       25.20,25.21,25.22,25.23,25.24,25.25,&
       26.00,26.01,26.02,26.03,26.04,26.05,26.06,26.07,26.08,26.09,&
       26.10,26.11,26.12,26.13,26.14,26.15,26.16,26.17,26.18,26.19,&
       26.20,26.21,26.22,26.23,26.24,26.25,26.26,&
       27.00,27.01,27.02,27.03,27.04,27.05,27.06,27.07,27.08,27.09,&
       27.10,27.11,27.12,27.13,27.14,27.15,27.16,27.17,27.18,27.19,&
       27.20,27.21,27.22,27.23,27.24,27.25,27.26,27.27,&
       28.00,28.01,28.02,28.03,28.04,28.05,28.06,28.07,28.08,28.09,&
       28.10,28.11,28.12,28.13,28.14,28.15,28.16,28.17,28.18,28.19,&
       28.20,28.21,28.22,28.23,28.24,28.25,28.26,28.27,28.28,&
       29.00,29.01,29.02,29.03,29.04,29.05,29.06,29.07,29.08,29.09,&
       29.10,29.11,29.12,29.13,29.14,29.15,29.16,29.17,29.18,29.19,&
       29.20,29.21,29.22,29.23,29.24,29.25,29.26,29.27,29.28,29.29,&
       30.00,30.01,30.02,30.03,30.04,30.05,30.06,30.07,30.08,30.09,&
       30.10,30.11,30.12,30.13,30.14,30.15,30.16,30.17,30.18,30.19,&
       30.20,30.21,30.22,30.23,30.24,30.25,30.26,30.27,30.28,30.29,&
       30.30,&
       31.00,31.01,31.02,31.03,31.04,&
       32.00,32.01,32.02,32.03,32.04,&
       33.00,33.01,33.02,33.03,33.04,&
       34.00,34.01,34.02,34.03,34.04,&
       35.00,35.01,35.02,35.03,35.04,&
       36.00,36.01,36.02,36.03,36.04,&
       37.00,37.01,37.02,37.03,37.04,&
       38.00,38.01,38.02,38.03,38.04,&
       39.00,39.01,39.02,39.03,39.04,&
       40.00,40.01,40.02,40.03,40.04,&
       41.00,41.01,41.02,41.03,41.04,&
       42.00,42.01,42.02,42.03,42.04,&
       43.00,43.01,43.02,43.03,43.04,&
       44.00,44.01,44.02,44.03,44.04,&
       45.00,45.01,45.02,45.03,45.04,&
       46.00,46.01,46.02,46.03,46.04,&
       47.00,47.01,47.02,47.03,47.04,&
       48.00,48.01,48.02,48.03,48.04,&
       49.00,49.01,49.02,49.03,49.04,&
       50.00,50.01,50.02,50.03,50.04,&
       51.00,51.01,51.02,51.03,51.04,&
       52.00,52.01,52.02,52.03,52.04,&
       53.00,53.01,53.02,53.03,53.04,&
       54.00,54.01,54.02,54.03,54.04,&
       55.00,55.01,55.02,55.03,55.04,&
       56.00,56.01,56.02,56.03,56.04,&
       57.00,57.01,57.02,57.03,57.04,&
       58.00,58.01,58.02,58.03,58.04,&
       59.00,59.01,59.02,59.03,59.04,&
       60.00,60.01,60.02,60.03,60.04,&
       61.00,61.01,61.02,61.03,61.04,&
       62.00,62.01,62.02,62.03,62.04,&
       63.00,63.01,63.02,63.03,63.04,&
       64.00,64.01,64.02,64.03,64.04,&
       65.00,65.01,65.02,65.03,65.04,&
       66.00,66.01,66.02,66.03,66.04,&
       67.00,67.01,67.02,67.03,67.04,&
       68.00,68.01,68.02,68.03,68.04,&
       69.00,69.01,69.02,69.03,69.04,&
       70.00,70.01,70.02,70.03,70.04,&
       71.00,71.01,71.02,71.03,71.04,&
       72.00,72.01,72.02,72.03,72.04,&
       73.00,73.01,73.02,73.03,73.04,&
       74.00,74.01,74.02,74.03,74.04,&
       75.00,75.01,75.02,75.03,75.04,&
       76.00,76.01,76.02,76.03,76.04,&
       77.00,77.01,77.02,77.03,77.04,&
       78.00,78.01,78.02,78.03,78.04,&
       79.00,79.01,79.02,79.03,79.04,&
       80.00,80.01,80.02,80.03,80.04,&
       81.00,81.01,81.02,81.03,81.04,&
       82.00,82.01,82.02,82.03,82.04,&
       83.00,83.01,83.02,83.03,83.04,&
       84.00,84.01,84.02,84.03,84.04,&
       85.00,85.01,85.02,85.03,85.04,&
       86.00,86.01,86.02,86.03,86.04,&
       87.00,87.01,87.02,87.03,87.04,&
       88.00,88.01,88.02,88.03,88.04,&
       89.00,89.01,89.02,89.03,89.04,&
       90.00,90.01,90.02,90.03,90.04,&
       91.00,91.01,91.02,91.03,91.04,&
       92.00,92.01,92.02,92.03,92.04,&
       93.00,93.01,93.02,93.03,93.04,&
       94.00,94.01,94.02,94.03,94.04,&
       95.00,95.01,95.02,95.03,95.04,&
       96.00,96.01,96.02,96.03,96.04,&
       97.00,97.01,97.02,97.03,97.04,&
       98.00,98.01,98.02,98.03,98.04,&
       99.00,99.01,99.02,99.03,99.04/

      lentape=31627892+10268732

      ratiolog=log(1.d0+1.d0/2000000.d0)

      eof=0
      found=.false.
      do while (eof.eq.0 .and. .not.found)
        read(nfile,*,IOSTAT=eof) iwl,ielion,ielo,igflog,igr,igs,igw
        if (eof.ne.0) cycle
        nline=nline+1

        nelion=abs(ielion)/10
        code=elem(nelion)
        z=nint(code)
        ion=nint(100.0d0*(code-z))
        if (.not.zlist(z)) cycle
        if (ion.gt.max_ion_levels) cycle
        gf=tablog(igflog)
        if (log10(gf).lt.gfcut) cycle
        if (wlist(z).eq.-1 .and. ielion.gt.0) cycle
        lambda=exp(iwl*ratiolog)*1.d-7 !! translate from nm
        elo=planck*clight*tablog(ielo) !! translate from cm-1 to erg
        ehi=elo+planck*clight/lambda
!!  fiducial values as no degeneracy values are given for low/hi energy levels
        glo=1.0d0
        ghi=1.0d0
        fij=gf/glo
        found=.true.
      enddo

      return
      end subroutine get_kurucz_1_next_line

      subroutine add_level(lvl,level,length)
      type (level_data) , intent(in) :: lvl
      type (level_data) , intent(inout) :: level(:)
      integer , intent(inout) :: length
!     integer , intent(out) :: ierr
      integer :: i

!     ierr=0
      
      if (length.eq.0) then
        level(1)=lvl
        length=1
        return
      endif

      i=1
      do while (level(i)%energy.lt.lvl%energy)
        i=i+1
      enddo

      if (level(i)%energy.eq.lvl%energy) then
        do while (level(i)%energy.eq.lvl%energy .and. level(i)%g.lt.lvl%g)
          i=i+1
        enddo
        if (level(i)%energy.eq.lvl%energy .and. level(i)%g.eq.lvl%g) return
      endif

      level(i+1:length+1)=level(i:length)
      level(i)=lvl
      length=length+1

      return
      end subroutine add_level

      integer function find_level(lvl,level,length)
      type (level_data) , intent(in) :: lvl,level(:)
      integer , intent(in) :: length
      integer :: i
      logical :: found

      found=.false.
      i=1
      do while (.not.found .and. i.le.length)
        if (level(i)%energy.eq.lvl%energy .and. level(i)%g.eq.lvl%g) then 
          found=.true.
        else 
          i=i+1
        endif
      enddo
      if (.not.found) print*,'find prob'

      find_level=i

      return
      end function find_level

      real(8) function calc_partition(temp,z,ion,method)
      real(8) , intent(in) :: temp
      integer , intent(in) :: z,ion,method
      real(8) :: part,kbt,en,g
      integer :: i,k,ierr

      kbt=kboltz*temp

      part=0.0d0

      if (method.eq.1) then

        do k=1,atom(z)%ion(ion)%nlevels
          g=atom(z)%ion(ion)%level(k)%g
          en=atom(z)%ion(ion)%level(k)%energy
          part=part+g*exp(-en/kbt)
        enddo

      elseif (method.eq.2) then

        i=find_index(temp,atom(z)%ion(ion)%ptemp(:),ierr)
        part=interp1(atom(z)%ion(ion)%ptemp(i:i+1),atom(z)%ion(ion)%partition(i:i+1),temp)

      endif

      calc_partition=part

      return
      end function calc_partition

      real(4) function tablog(ii)
      integer(2) , intent(in) :: ii

      tablog=10.**((ii-16384)*.001)
      
      return
      end function tablog

      subroutine add_full_isotop_list
      integer :: ni(2,99)
      integer :: n,i,j

      data ni(:,1)/1,3/
      data ni(:,2)/3,6/
      data ni(:,3)/6,8/
      data ni(:,4)/7,11/
      data ni(:,5)/8,12/
      data ni(:,6)/10,15/
      data ni(:,7)/12,17/
      data ni(:,8)/14,20/
      data ni(:,9)/17,21/
      data ni(:,10)/18,25/
      data ni(:,11)/20,26/
      data ni(:,12)/21,28/
      data ni(:,13)/23,30/
      data ni(:,14)/25,33/
      data ni(:,15)/27,35/
      data ni(:,16)/29,38/
      data ni(:,17)/31,40/
      data ni(:,18)/33,44/
      data ni(:,19)/35,46/
      data ni(:,20)/37,49/
      data ni(:,21)/40,50/
      data ni(:,22)/42,52/
      data ni(:,23)/44,54/
      data ni(:,24)/46,56/
      data ni(:,25)/48,58/
      data ni(:,26)/50,62/
      data ni(:,27)/52,63/
      data ni(:,28)/54,67/
      data ni(:,29)/57,69/
      data ni(:,30)/59,72/
      data ni(:,31)/61,74/
      data ni(:,32)/63,78/
                 
      allocate(indiso(max_atoms,100))
      allocate(iso(max_isotops))
      indiso=-1
      n=0
      do i=1,32
        do j=ni(1,i),ni(2,i)
          n=n+1
          iso(n)%z=i
          iso(n)%a=j
          iso(n)%sym=100*iso(n)%z+nint(iso(n)%a)
          indiso(iso(n)%z,nint(iso(n)%a))=n
          iso(n)%a=iso(n)%a/avogadro
        enddo
      enddo
      niso=n
                 
      return
      end subroutine add_full_isotop_list

      end module Atomic_Physics
