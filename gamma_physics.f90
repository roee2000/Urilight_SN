! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! CALCULATING GAMMA-RAYS INTERACTIONS

       module Gamma_Physics
       use physical_constants
       implicit none







       contains


       subroutine total_compton_scattering(Ein,xNi,xZi,sigtot)
       real(8) , intent (in) :: Ein,xNi(:)
       integer , intent (in) :: xZi(:)
       real(8) , intent (out) :: sigtot
       real(8) :: x,x2,xp1,x2p1,logx2p1
       real(8) :: Kx,A1,A2,B1,C1

       x=Ein/mec2mev
       x2=2*x
       xp1=1.0d0+x
       x2p1=1.0d0+x2
       logx2p1=log(x2p1)

       A1=xp1/x**3.0d0
       A2=x2*xp1/x2p1-logx2p1
       B1=1.0d0/x2*logx2p1
       C1=-(1.0d0+3.0d0*x)/x2p1**2.0d0

       Kx=3.0d0/4.0d0*(A1*A2+B1+C1)

       sigtot=sigma_thomson*Kx*sum(xNi(:)*dble(xZi(:)))

       return
       end subroutine total_compton_scattering
       
       subroutine sample_compton_scattering(Ein,x,cost,phi)
       real(8) , intent(in) :: Ein
       real(8) , intent(out) :: x,cost,phi
       real(8) :: E0,x0,alpha1,alpha2,z(3),zp,f,sin2t,ge
       logical :: accept

       E0=Ein/mec2mev
       x0=1.0d0/(1.0d0+2.0d0*E0)

       accept=.false.
       do while (.not. accept) 

         call random_number(z(:))
         alpha1=log(1./x0)
         alpha2=(1.0d0-x0**2.0d0)/2.0d0
         if (z(1).lt.alpha1/(alpha1+alpha2)) then
           x=x0**z(2)
         else
           x=sqrt(x0**2.0d0+(1.0d0-x0**2.0d0)*z(2))
         endif
         f=(1.0d0-x)/x/E0
         sin2t=f*(2.0d0-f)
         ge=1.0d0-x/(1+x**2.0d0)*sin2t
         if (ge.gt.z(3)) accept=.true.

       enddo

       cost=1.0d0-f
      
       call random_number(zp)

       phi=2.0d0*pi*zp

       return
       end subroutine sample_compton_scattering

       subroutine total_photoelectric_absorption(Ein,xNi,xZi,sigtot)
       real(8) , intent (in) :: Ein,xNi(:)
       integer , intent (in) :: xZi(:)
       real(8) , intent (out) :: sigtot
       real(8) :: x,x2,xp1,x2p1,logx2p1
       real(8) :: Kx,A1,A2,B1,C1

       x=Ein/mec2mev

       Kx=fine_structure**4.0d0*8.0d0*sqrt(2.0d0)*x**(-3.5d0)

       sigtot=sigma_thomson*Kx*sum(xNi(:)*dble(xZi(:))**5.0d0)

       return
       end subroutine total_photoelectric_absorption

       subroutine total_pair_production(Ein,xNi,xZi,sigtot)
       real(8) , intent (in) :: Ein,xNi(:)
       integer , intent (in) :: xZi(:)
       real(8) , intent (out) :: sigtot

       sigtot=0.0d0
       
       if (Ein.lt.2.0d0*mec2mev) return

       if (Ein.gt.1.5d0) then
         sigtot=0.0481d0+0.301d0*(Ein-1.5d0)
       else
         sigtot=1.0063d0*(Ein-2.0d0*mec2mev)
       endif

       sigtot=sigtot*1.d-27*sum(xNi(:)*dble(xZi(:))**2.0d0)

       return
       end subroutine total_pair_production

       end module Gamma_Physics
