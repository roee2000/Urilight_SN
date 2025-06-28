! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! ACUTAL GAMMA-RAY TRANSFER WITH MONTE-CARLO

      Module GammaTransfer
      use globals
      use arrays
      use general_functions
      use transport_general_functions
      use Radioactive_Decay
      use RandomNumbers
      use Mesh
      use diagnostics
      use Gamma_Physics
      implicit none

      logical , save :: ispe=.true.
      logical , save :: iscompton=.true.
      logical , save :: ispp=.false.

      real(8) , save :: EpelletGamma

      contains

      subroutine init_gamma
      integer :: i,j,k,totpel,ino
      real(8) :: Ni56Edep,Co56Edep
      real(8) :: xmin,xmax,deltax,mni
      integer :: nbins
      namelist /gamma/ N_GammaPellets,spect_type_gamma,xmin,xmax,deltax,nbins,ispe,iscompton,ispp

!!    default
      spect_type_gamma=1 !! bins by energy (in mev)
      xmin=0.0d0
      xmax=4.0d0
      nbins=100
     
      open(unit=5,file=data_file)
      read(5,nml=gamma,iostat=ino)
      close(5)
      write(fout,nml=gamma)

      allocate (Edep_gamma(ntimes,nctot),edep_pos(ntimes,nctot))
      allocate (ecr_gamma(ntimes))
      allocate (spect_gamma(nbins,ntimes))
      allocate (spect_bins_gamma(nbins+1))
      Edep_gamma=0.0d0
      edep_pos=0.0d0
      ecr_gamma=0.0d0
      spect_gamma=0.0d0
      spect_bins_gamma=0.0d0
      
      do i=1,nbins+1
        spect_bins_gamma(i)=xmin+(xmax-xmin)/dble(nbins)*dble(i-1)
      enddo

      totpel=0
      do k=1,nc3
        do j=1,nc2
          do i=1,nc1
            mni=atoms(ind_ni56,ind(i,j,k))*iso(ind_ni56)%A
            pellets(i,j,k)=nint(dble(N_GammaPellets)*mni/TotNi56)
            totpel=totpel+pellets(i,j,k)
          enddo
        enddo
      enddo

      call RadioactiveEnergyDepo(sum(atoms(ind_ni56,:)),0.0d0,tfinal,Ni56Edep,Co56Edep)
      EpelletGamma=(Ni56Edep+Co56Edep)*mev/totpel
       
      write(fout,1001) (Ni56Edep+Co56Edep)*mev
      write(fout,1002) Co56Edep*PosEfrac*mev
      write(fout,1003) totpel
1001  format('Gamma energy source (not including positrons) is ',1pe14.6,' erg')
1002  format('Positron energy deposited for uvoir transfer is ',1pe14.6,' erg')
1003  format('Total number of gamma pellets is ',I10)

      return
      end subroutine init_gamma

      subroutine gamma_transport
      integer :: i,j,k,nt,np,ip,jp,kp,NiCo,ierr
      logical :: inmesh,intime,isabs,isprob
      type (epacket) :: p
      integer :: nprob,ndep,nout,ndirect,nscat
      real(8) :: Ni56Edep,Co56Edep,Etotm,vm(3)

!     call calc_rodr
      if (onlyrodr) stop

      nprob=0
      nout=0
      ndep=0
      ndirect=0
      nscat=0
     
      do k=1,nc3
        do j=1,nc2
          do i=1,nc1
            print*,'i,j,k,np=',i,j,k,pellets(i,j,k)
            do np=1,pellets(i,j,k)

            p=new_gamma(i,j,k,EpelletGamma,NiCo)

            nt=max(find_index(p%t,times,ierr),1) !! all energy deposited before tinit is given to first time step

            ecr_gamma(nt)=ecr_gamma(nt)+EpelletGamma

            if (NiCo.eq.2) then  !! deposit Co56 positron locally
!             Edep_pos(nt,ind(i,j,k))=Edep_pos(nt,ind(i,j,k))+p%Etot*posEfrac
              edep_pos(nt,ind(i,j,k))=edep_pos(nt,ind(i,j,k))+EpelletGamma*posefrac*1.0d0
            endif


            ip=i
            jp=j
            kp=k
            call track_gamma(ip,jp,kp,p,inmesh,intime,isabs,isprob)

            if (isprob) then
              nprob=nprob+1
            endif

            if (.not.inmesh) then
              call write_to_spectrum(p,spect_gamma,spect_bins_gamma,spect_type_gamma)
              call integrate_bolometric(p,gout,2)
              nout=nout+1
              if (p%direct) then
                ndirect=ndirect+1
              else
                nscat=nscat+1
              endif
            elseif (isabs) then
              nt=max(find_index(p%t,times,ierr),1) !! all energy deposited before tinit is given to first time step
!!  deposit in comoving frame !
              vm=vofr(p%r,teff(nt))
              Etotm=comoving2lab_transform_E(p%Etot,p%n,vm,2)
              Edep_gamma(nt,ind(ip,jp,kp))=Edep_gamma(nt,ind(ip,jp,kp))+Etotm
              ndep=ndep+1
            endif

            enddo
          enddo
        enddo
      enddo
      print*,'nprob=',nprob
      print*,'nout=',nout
      print*,'ndirect,nscat,tot=',ndirect,nscat,ndirect+nscat
      print*,'ndep=',ndep

      return
      end subroutine gamma_transport

      function new_gamma(i,j,k,Etot,NiCo)
      integer , intent (in) :: i,j,k
      real(8) , intent (in) :: Etot
      integer , intent (out) :: NiCo
      type (epacket) :: new_gamma
      real(8) :: v(3),fac
      integer :: nt,ierr

10    call get_new_pellet(new_gamma%t,new_gamma%hnu,NiCo)
      if (new_gamma%t.gt.tfinal) goto 10
        
      fac=1.0d0
      if (new_gamma%t.lt.tinit) then
        fac=new_gamma%t/tinit
        new_gamma%t=tinit
      endif

      nt=find_index(new_gamma%t,times,ierr)
      
      new_gamma%r=random_location(i,j,k,teff(nt))
      new_gamma%n=random_unit_vec1(3)
      v=vofr(new_gamma%r,teff(nt))

      new_gamma%hnu=comoving2lab_transform_E(new_gamma%hnu,new_gamma%n,v,1)
      new_gamma%Etot=comoving2lab_transform_E(Etot,new_gamma%n,v,1)*fac
      new_gamma%n=comoving2lab_transform_n(new_gamma%n,v,1)

      new_gamma%direct=.true.

      return
      end function new_gamma

      subroutine track_gamma(i,j,k,p,inmesh,intime,isabs,isprob)
      integer , intent (inout) :: i,j,k
      type (epacket) , intent (inout) :: p
      logical , intent (out) :: inmesh,intime,isabs,isprob
      integer :: n,nt,i1,j1,k1,ierr
      real(8) :: v(3),vm(3),nm(3),rhom,Etotm,rcut(3),tcut,hnum,z,dt,rad
      real(8) :: q,f,cost,phi
      real(8) :: ds,ds_time,ds_edge,ds_event
      real(8) :: sig_s,sig_pe,sig_pp,sig_tot


      inmesh=.true.
      intime=.true.
      isprob=.false.
      isabs=.false.

!     print*,'start track'

      n=0

      nt=find_index(p%t,times,ierr)

      do while (inmesh .and. intime .and. .not.isabs)

        ds_time=(times(nt+1)-p%t)*clight

!       print*,'in cell i,j,k,nt=',i,j,k,nt
        call cut_nearest_surface(p%r,p%n,p%t,nt,i,j,k,ds_edge,rcut,tcut,i1,j1,k1,isprob)

        if (isprob) exit
        rhom=rhooft(rhov(ind(i,j,k)),teff(nt))
        vm=vofr(p%r,teff(nt))

        hnum=comoving2lab_transform_E(p%hnu,p%n,vm,2)
        Etotm=comoving2lab_transform_E(p%Etot,p%n,vm,2)
        sig_s=0.0d0
        sig_pe=0.0d0
        sig_pp=0.0d0
        if (iscompton) then
          call total_compton_scattering(hnum,atoms(1:Niso,ind(i,j,k)),iso(1:Niso)%z,sig_s)
        endif
        if (ispe) then
          call total_photoelectric_absorption(hnum,atoms(1:Niso,ind(i,j,k)),iso(1:Niso)%z,sig_pe)
        endif
        if (ispp) then
          call total_pair_production(hnum,atoms(1:Niso,ind(i,j,k)),iso(1:Niso)%z,sig_pe)
        endif

        sig_tot=(hnum/p%hnu)*(sig_s+sig_pe+sig_pp)

        call random_number(z)

        ds_event=-log(z)/(rhom*sig_tot)

        ds=min(ds_time,ds_edge,ds_event)

        p%r=p%r+p%n*ds
        p%t=p%t+ds/clight

        if (ds.eq.ds_time) then
          nt=nt+1
          if (nt.gt.ntimes) then
            intime=.false.
          else
!           call findijk(p%r,teff(nt),i,j,k)
            rad=sqrt(sum(p%r(:)**2.0d0))
            do i1=1,nc1
              if (rad.ge.v1(i1)*teff(nt) .and. rad.le.v1(i1+1)*teff(nt)) i=i1
            enddo
          endif
        elseif (ds.eq.ds_edge) then
          i=i1
          j=j1
          k=k1
          if (i.gt.nc1 .or. j.gt.nc2 .or. k.gt.nc3) inmesh=.false.
          if (i.lt.1   .or. j.lt.1   .or. k.lt.1)   inmesh=.false.
        elseif (ds_event.eq.ds) then
          n=n+1
          p%direct=.false.

          vm=vofr(p%r,teff(nt))
          hnum=comoving2lab_transform_E(p%hnu,p%n,vm,2)
          Etotm=comoving2lab_transform_E(p%Etot,p%n,vm,2)

          call random_number(z)
          if (z.lt.sig_s/(sig_s+sig_pe+sig_pp)) then
            
            call sample_compton_scattering(hnum,f,cost,phi)
            call random_number(z)
            if (z.lt.f) then  !! scatter
              nm(:)=comoving2lab_transform_n(p%n(:),vm,2)
              nm=new_direction(nm,cost,phi)
              hnum=f*hnum
              Etotm=Etotm
              p%hnu=comoving2lab_transform_E(hnum,nm,vm,1)
              p%Etot=comoving2lab_transform_E(Etotm,nm,vm,1)
              p%n=comoving2lab_transform_n(nm,vm,1)
            else              !! absorb
              isabs=.true.
            endif
          else
            call random_number(z)
            if (z.lt.sig_pe/(sig_pe+sig_pp)) then
              isabs=.true.
            else
              call random_number(z)
              q=mec2mev/hnum
              if (z.lt.(0.5d0+q)) then
                nm=random_unit_vec1(3)
                hnum=mec2mev
                Etotm=Etotm
                p%hnu=comoving2lab_transform_E(hnum,nm,vm,1)
                p%Etot=comoving2lab_transform_E(Etotm,nm,vm,1)
                p%n=comoving2lab_transform_n(nm,vm,1)
              else
                isabs=.true.
              endif
            endif
          endif
        endif
            
      enddo

      return
      end subroutine track_gamma

      subroutine calc_rodr
      integer :: i,j,k,np,ip,jp,kp,ip1,jp1,kp1,n
      logical :: inmesh,isabs,isprob
      type (epacket) :: p
      integer :: nprob,ndep,nout,ndirect,nscat
      real(8) :: rcut(3),tcut,ds,vnorm,rhom
      real(8) :: rodr,rodr1

      rodr=0.0d0

      nprob=0
      nout=0
      ndep=0
      ndirect=0
      nscat=0

!     slow mesh
      vnorm=v1(nc1p)
      v1=v1/vnorm
     
      n=0
      do k=1,nc3
        do j=1,nc2
          do i=1,nc1
!           print*,'i,j,k,np=',i,j,k,pellets(i,j,k)
            do np=1,pellets(i,j,k)/10

            p%t=1.0d0
            p%r=random_location(i,j,k,p%t)
            p%n=random_unit_vec1(3)

!           p%r=p%r/sqrt(sum(p%r**2.0d0))*v1(51)*1.0001d0
!           p%n=p%r/sqrt(sum(p%r**2.0d0))

            ip=i
            jp=j
            kp=k


            inmesh=.true.
            isprob=.false.
            isabs=.false.
          
            rodr1=0.0d0
            do while (inmesh .and. .not.isabs)

!             call cut_nearest_surface(p%r,p%n,p%t,ip,jp,kp,ds,rcut,tcut,ip1,jp1,kp1,isprob)

              if (isprob) exit

              p%r=rcut
              rhom=rhooft(rhov(ind(ip,jp,kp)),p%t)

              rodr=rodr+ds*vnorm*mass(ind(ip,jp,kp))*rhom
              rodr1=rodr1+ds*vnorm*mass(ind(ip,jp,kp))*rhom

              p%t=tcut
              ip=ip1
              jp=jp1
              kp=kp1
              if (ip.gt.nc1 .or. jp.gt.nc2 .or. kp.gt.nc3) inmesh=.false.
              if (ip.lt.1   .or. jp.lt.1   .or. kp.lt.1)   inmesh=.false.

            enddo

            if (isprob) then
              nprob=nprob+1
            endif

            if (.not.isprob) n=n+1

            enddo

          enddo
        enddo
      enddo
      rodr=rodr/dble(n)

      print*,'rodr,tday=',rodr,sqrt(rodr*0.025d0)/day

      v1=v1*vnorm

      return
      end subroutine calc_rodr

      end Module GammaTransfer
