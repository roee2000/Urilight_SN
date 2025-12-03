! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! NUCLEAR DATA FOR CALCULATING RADIOACTIVE DECAY

      module radioactive_decay
      use physical_constants
      use RandomNumbers
      implicit none

      type gamma_line
        real(8) :: energy=0.0d0
        real(8) :: frac=0.0d0
        real(8) :: energy_frac=0.0d0
      end type gamma_line

      integer :: Ni56lines,Co56lines
      real(8) :: Ni56t12,Co56t12
      real(8) :: Ni56tau,Co56tau
      real(8) :: Ni56Etot,Co56Etot
      real(8) :: Ni56EtotFrac
      real(8) :: PositronEnergy,PosEfrac
      type (gamma_line) , save :: Ni56(30) , Co56(30)


      contains

      subroutine init_nuclear_data
      integer :: i
      
!     default data from Ambwani&Sutherland (Apj,325,820-827,1988)
!     half lifes from Milne 2004

      Ni56t12=6.075d0*day
      Ni56tau=Ni56t12/log(2.0d0)
      Ni56lines=6
      Ni56(1)%energy=0.158d0
      Ni56(2)%energy=0.270d0
      Ni56(3)%energy=0.480d0
      Ni56(4)%energy=0.750d0
      Ni56(5)%energy=0.812d0
      Ni56(6)%energy=1.562d0

      Ni56(1)%frac=1.0d0
      Ni56(2)%frac=0.36d0
      Ni56(3)%frac=0.36d0
      Ni56(4)%frac=0.5d0
      Ni56(5)%frac=0.87d0
      Ni56(6)%frac=0.14d0

      Co56t12=77.233d0*day
      Co56tau=Co56t12/log(2.0d0)
      Co56lines=23
      Co56(1)%energy=0.511d0
      Co56(2)%energy=0.734d0
      Co56(3)%energy=0.788d0
      Co56(4)%energy=0.847d0
      Co56(5)%energy=0.978d0
      Co56(6)%energy=1.038d0
      Co56(7)%energy=1.140d0
      Co56(8)%energy=1.175d0
      Co56(9)%energy=1.238d0
      Co56(10)%energy=1.360d0
      Co56(11)%energy=1.443d0
      Co56(12)%energy=1.772d0
      Co56(13)%energy=1.811d0
      Co56(14)%energy=1.964d0
      Co56(15)%energy=2.015d0
      Co56(16)%energy=2.035d0
      Co56(17)%energy=2.213d0
      Co56(18)%energy=2.598d0
      Co56(19)%energy=3.010d0
      Co56(20)%energy=3.202d0
      Co56(21)%energy=3.254d0
      Co56(22)%energy=3.273d0
      Co56(23)%energy=3.452d0

      Co56(1)%frac=0.3800d0
      Co56(2)%frac=0.0021d0
      Co56(3)%frac=0.0030d0
      Co56(4)%frac=0.9998d0
      Co56(5)%frac=0.0144d0
      Co56(6)%frac=0.1408d0
      Co56(7)%frac=0.0015d0
      Co56(8)%frac=0.0224d0
      Co56(9)%frac=0.6758d0
      Co56(10)%frac=0.0428d0
      Co56(11)%frac=0.0020d0
      Co56(12)%frac=0.1600d0
      Co56(13)%frac=0.0048d0
      Co56(14)%frac=0.0072d0
      Co56(15)%frac=0.0309d0
      Co56(16)%frac=0.0795d0
      Co56(17)%frac=0.0063d0
      Co56(18)%frac=0.1672d0
      Co56(19)%frac=0.0100d0
      Co56(20)%frac=0.0303d0
      Co56(21)%frac=0.0743d0
      Co56(22)%frac=0.0176d0
      Co56(23)%frac=0.0086d0

      PositronEnergy=0.632d0*0.19d0
   


      Ni56Etot=sum(Ni56(:)%energy*Ni56(:)%frac)
      Ni56(:)%energy_frac=Ni56(:)%energy*Ni56(:)%frac/Ni56Etot
      Co56Etot=sum(Co56(:)%energy*Co56(:)%frac)
      Co56(:)%energy_frac=Co56(:)%energy*Co56(:)%frac/Co56Etot
      Ni56EtotFrac=Ni56Etot/(Ni56Etot+Co56Etot)

      print*,'Ni56Etot=',Ni56Etot
      print*,'Co56Etot=',Co56Etot

      PosEfrac=PositronEnergy/Co56Etot


      return
      end subroutine init_nuclear_data

      subroutine get_new_pellet(ti,energy,NiCo)
      real(8) , intent(out) :: ti,energy
      integer , intent(out) :: NiCo
      integer :: decay_type,i
      real(8) :: frac,zrand(2)

      
      decay_type=1 !! set default decay to Ni56

      call random_number(zrand(1))
      if (zrand(1).gt.Ni56EtotFrac) decay_type=2 !! decay is Co56

      if (decay_type.eq.1) then 
         call random_number(zrand(1))
         ti=-Ni56tau*log(zrand(1))
         i=choose_from_probability_distribution(Ni56(1:Ni56lines)%energy_frac)
         energy=Ni56(i)%energy
         NiCo=1
      else
         call random_number(zrand)
         ti=-Ni56tau*log(zrand(1))-Co56tau*log(zrand(2))
         i=choose_from_probability_distribution(Co56(1:Co56lines)%energy_frac)
         energy=Co56(i)%energy
         NiCo=2
      endif

      return
      end subroutine get_new_pellet

      subroutine RadioactiveEnergyDepo(nnuc,t1,t2,Ni56Edep,Co56Edep)
      implicit none
      real(8) , intent (in) :: nnuc,t1,t2
      real(8) , intent (out) :: Ni56Edep,Co56Edep
      real(8) :: totni56,totco56,nni1,nco1,nfe1,nni2,nco2,nfe2

      totni56=-Ni56tau*(exp(-t2/Ni56tau)-exp(-t1/Ni56tau))
      totco56=-Co56tau*(exp(-t2/Co56tau)-exp(-t1/Co56tau))

      Ni56edep=nnuc*Ni56Etot/Ni56tau*totni56
      Co56edep=nnuc*Co56Etot/(Ni56tau-Co56tau)*(totni56-totco56)

!     call Ni56DecayChain(nnuc,t1,nni1,nco1,nfe1)
!     call Ni56DecayChain(nnuc,t2,nni2,nco2,nfe2)

      return
      end subroutine RadioactiveEnergyDepo

      subroutine Ni56DecayChain(nnuc,t,nni,nco,nfe)
      implicit none
      real(8) , intent (in) :: nnuc,t
      real(8) , intent (out) :: nni,nco,nfe

      nni=nnuc*exp(-t/Ni56tau)
      nco=nnuc*Co56tau/(Ni56tau-Co56tau)*(exp(-t/Ni56tau)-exp(-t/Co56tau))
      nfe=nnuc-nni-nco

      return
      end subroutine Ni56DecayChain


      end module radioactive_Decay
