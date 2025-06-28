! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! CALCULATING UVOIR INTERACTIONS
      Module Uvoir_Physics
      use physical_constants
      use atomic_physics
      implicit none

      real(8) , parameter :: bigexp=200.0d0

      real(8) , save :: mintemp=100.0d0

      contains

      subroutine calc_planck_int(bp,fnorm,reslow,reshigh,wbins,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!A!!!!!!!!!!!!!!!!!!!!
!!    This subourinte calculates Integral on planck function !!
!!    for all wavelength bins.                               !!
!!    At the moment a very simple interation is implemented. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) , intent(in) :: wbins(:),temp
      real(8) , intent(out) :: bp(:),fnorm,reslow,reshigh
      integer :: i,j,nb
      real(8) :: lam1,lam2,x1,x2,x,dx,fx,kbt,x11,x22

      kbt=kboltz*temp
      fnorm=2.0d0*planck/clight**2.0d0*(kbt/planck)**4.0d0
      fx=planck*clight/kbt

      nb=100

      do i=1,size(bp)
        lam1=wbins(i)
        lam2=wbins(i+1)

        x1=fx/lam2
        x2=fx/lam1

        dx=(x2-x1)

        bp(i)=0.0d0
        do j=1,nb
          x11=x1+dx/dble(nb)*dble(j-1)
          x22=x1+dx/dble(nb)*dble(j)
          x=(x11+x22)/2.0d0
          x=min(x,bigexp)  !! exp(bigexp)~1e99 already
          bp(i)=bp(i)+x**3.0d0/(exp(x)-1.0d0)
        enddo
        bp(i)=bp(i)*dx/dble(nb)
      enddo
      
      bp=fnorm*bp

      x1=fx/wbins(1)
      reslow=fnorm*(bradley_clark_pi(1000.0d0)-bradley_clark_pi(x1))
      x1=fx/wbins(size(bp)+1)
      reshigh=fnorm*bradley_clark_pi(x1)

      return
      end subroutine calc_planck_int

      subroutine calc_freefree_abs(alpha,temp,ne,nij,spect_bins)
      real(8) , intent(in) :: temp,ne,nij(0:,:),spect_bins(:)
      real(8) , intent(out) :: alpha(:)
      real(8) :: coeff,z2neni,hokt,nu,lam
      integer :: i,j,niso,nions,nwave
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    This subroutine calculates the free-free absorption coefficient [1/cm]
!!    according to the Bremsstrahlung formula from "Radiative Processes in Astrophysics"
!!    by Rybicki & Lightman, equation 5.18a.
!!    guant factor assums to be 1 for the moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      niso=size(nij,2)
      nions=size(nij,1)-1
      nwave=size(alpha)

      hokt=planck/kboltz/temp

      coeff=4.0d0*electron_charge**6.0d0/(3.0d0*electron_mass*planck*clight)
      coeff=coeff*(2.0d0*pi/(3.0d0*kboltz*electron_mass))**0.5d0
      coeff=coeff/temp**0.5d0

      z2neni=0.0d0
      do i=1,niso
        do j=1,nions
          z2neni=z2neni+dble(j)**2.0d0*nij(j,i)
        enddo
      enddo
      z2neni=z2neni*ne

      do i=1,nwave
        lam=(spect_bins(i)+spect_bins(i+1))/2.0d0
        nu=clight/lam
        alpha(i)=(1.0d0-exp(-hokt*nu))/nu**3.0d0
      enddo

      alpha(:)=coeff*z2neni*alpha(:)

      return
      end subroutine calc_freefree_abs

      subroutine expansion_opacity_LTE(alpha_abs,alpha_scat,time,temp,nij,partition,spect_bins)
      real(8) , intent(in) :: time,temp,nij(0:,:),partition(0:,:),spect_bins(:)
      real(8) , intent(out) :: alpha_abs(:),alpha_scat(:)
      integer :: n,i,j,z,k1,k2,niso,nions,nwave,nbins,lambin
      real(8) :: nijk,nu,g,fij,kbt,en1,en2,eline,stimfac,sigtot,alpha,onexp,pa,ps
      real(8) :: dlam,lambda,ct,tausob,psob
      real(8) :: nop(0:max_ion_levels,max_atoms)
      logical :: activez(max_atoms)

      alpha_abs(:)=0.0d0
      alpha_scat(:)=0.0d0

      niso=size(nij,2)
      nions=size(nij,1)-1
      nwave=size(alpha_abs)
      nbins=size(spect_bins)

      nop=0.0d0
      activez=.false.

!! note: There can be several isotops in a given cell with the same Z,
!!       for example Fe54 and Fe56. Both have the same lines but differenet
!!       densities and both contribute to line opacity. The energy level density
!!       should be therefore first summed up as the opacity is exp(-tau) and so
!!       exp(-(n1+n2)*...) is the correct expression and not exp(-n1*...)+exp(-n2*...)
!!       if differnet isotops of the same element would have been treated seperatly !
      do i=1,niso
        do j=0,max_ion_levels

        if (nij(j,i).gt.0.0d0) then
          nop(j,iso(i)%z)=nop(j,iso(i)%z)+nij(j,i)/partition(j,i)
          activez(iso(i)%z)=.true.
        endif
         
        enddo
      enddo

      kbt=kboltz*temp
      sigtot=pi*electron_charge**2.0d0/(electron_mass*clight)
      ct=clight*time

      lambin=1
      do n=1,nlines

        lambda=line(n)%lambda

        if (lambda.lt.spect_bins(1) .or. lambda.gt.spect_bins(nbins)) cycle

        do while (lambda.gt.spect_bins(lambin+1))
          lambin=lambin+1
        enddo
        if (lambda.lt.spect_bins(lambin) .or. lambda.gt.spect_bins(lambin+1)) print*,'bin prob sob'

        z=line(n)%z
        j=line(n)%ion
        if (.not.activez(z)) cycle

        k1=line(n)%l1
        k2=line(n)%l2
        fij=line(n)%fij
        nu=clight/lambda

        g=atom(z)%ion(j)%level(k1)%g
        en1=atom(z)%ion(j)%level(k1)%energy
!         en2=atom(z)%ion(j)%level(k2)%energy
!         eline=en2-en1
        eline=planck*nu

!       nijk=nij(j,i)*g*exp(-en1/kbt)/partition(j,i)
        nijk=nop(j,z)*g*exp(-en1/kbt)

        stimfac=1.0d0-exp(-eline/kbt)

        alpha=nijk*sigtot*fij*stimfac

        tausob=alpha*ct/nu

        if (tausob.lt.1d-20) cycle

        onexp=(1.0d0-exp(-tausob))
        psob=onexp/tausob
       
        pa=etherm/(etherm+psob*(1.0d0-etherm))

        if (z.eq.20) pa=0.0d0 !! Kasen 2006 ?!
        ps=1.0d0-pa

        alpha_abs(lambin)=alpha_abs(lambin)+pa*lambda*onexp
        alpha_scat(lambin)=alpha_scat(lambin)+ps*lambda*onexp

      enddo

      do lambin=1,nwave
        dlam=spect_bins(lambin+1)-spect_bins(lambin)
        alpha_abs(lambin)=alpha_abs(lambin)/dlam/ct
        alpha_scat(lambin)=alpha_scat(lambin)/dlam/ct
      enddo

      return
      end subroutine expansion_opacity_LTE

      real(8) function radiation_temperature(Jnu,model)
      real(8) , intent(in) :: jnu
      integer , intent(in) :: model

      radiation_temperature=(pi*jnu/sigma_sb)**0.25d0

      return
      end function radiation_temperature

      real(8) function color_temperature(jnu,nujnu)
      real(8) , intent(in) :: jnu , nujnu
      real(8) , save :: avgx
      real(8) :: avgnu
      data avgx/3.8322d0/

      avgnu=nujnu/(jnu+1.d-32)

      color_temperature=planck*avgnu/avgx/kboltz

      return
      end function color_temperature

      real(8) function plasma_temperature(Edep,alpha_abs,wbins,tguess)
      real(8) , intent(in) :: Edep,alpha_abs(:),wbins(:),tguess
      real(8) :: bpx(size(alpha_abs))
      integer :: nb,niter
      real(8) :: fnorm,f1,f2,told,tnew,conv,reslow,reshigh,int_alpha_bp

      nb=size(alpha_abs)

      tnew=max(tguess,mintemp)
      conv=1.0d0
      niter=0
      do while (conv.gt.1.d-4 .and. niter.lt.100)
        told=tnew

        call calc_planck_int(bpx(:),fnorm,reslow,reshigh,wbins(:),told)

        f1=edep
        int_alpha_bp=sum(alpha_abs(:)*bpx(:))+&
                     alpha_abs(1)*reslow+alpha_abs(nb)*reshigh
        int_alpha_bp=int_alpha_bp/fnorm
        f2=2.0d0*planck/clight**2.0d0*int_alpha_bp

        tnew=planck/kboltz*(f1/f2)**(1.0d0/4.0d0)

        tnew=0.5d0*tnew+0.5d0*told
        tnew=max(tnew,mintemp)

        conv=abs(tnew-told)/(told+50.0d0)

        niter=niter+1
        if (niter.gt.80) then
          write(66,1212) niter,edep,tguess,told,tnew,conv
1212  format('not going to converge...=',I4,10(1pe14.6))
        endif
      enddo

      if (niter.eq.100) write(66,*) &
          'plasma temperature iter problem , conv=',conv

      plasma_temperature=tnew

      return
      end function plasma_temperature

      real(8) function bradley_clark_pin(x)
      real(8) , intent(in) :: x
      integer :: i
      real(8) :: res
      real(8) :: coef(13)
      data coef(1:13)/0d0,0d0,3d0,8d0,60d0,&
      0d0,5040d0,0d0,272160d0,0d0,13305600d0,&
      0d0,622702080d0/

      res=0.0d0

      do i=3,5
        res=res+x**dble(i)/coef(i)
      enddo
      do i=7,13,2
        res=res+x**dble(i)/coef(i)
      enddo
     
      bradley_clark_pin=res

      return
      end function bradley_clark_pin

      real(8) function bradley_clark_pil(x)
      real(8) , intent(in) :: x
      integer :: i,order,l
      real(8) :: res,expl(4),sig(4)

      res=0.0d0

      order=5

      sig=0.0d0
      do l=1,order
        expl(1)=exp(-l*x)/dble(l)
        expl(2)=expl(1)/dble(l)
        expl(3)=expl(2)/dble(l)
        expl(4)=expl(3)/dble(l)
        sig(1:4)=sig(1:4)+expl(1:4)
      enddo
        
      res=-x**3.0d0*sig(1)-3.0d0*x**2.0d0*sig(2)-6.0d0*x*sig(3)-6.0d0*sig(4)

      bradley_clark_pil=pi**4.0d0/15.0d0+res

      return
      end function bradley_clark_pil

      real(8) function bradley_clark_pi(x)
      real(8) , intent(in) :: x
      real(8) :: pil,pin

      pil=bradley_clark_pil(x)
      pin=bradley_clark_pin(x)
      
      if (pin.gt.pil) then
        bradley_clark_pi=pil!-pi**4.0d0/15.0d0
      else
        bradley_clark_pi=pin
      endif

      return
      end function bradley_clark_pi

      end module Uvoir_Physics
