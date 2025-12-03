! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
      Module physical_constants
      implicit none

      real(8) , parameter :: pi=3.14159265358979323846d0
      real(8) , parameter :: angstrom=1.d-8
      real(8) , parameter :: avogadro=6.02214129d23
      real(8) , parameter :: clight=29979245800.0d0     ! [cm/s]
      real(8) , parameter :: planck=6.62606957d-27
      real(8) , parameter :: kboltz=1.3806488d-16
      real(8) , parameter :: electron_mass=9.10938215d-28
      real(8) , parameter :: electron_charge=4.80320425d-10
      real(8) , parameter :: electron_volt=1.602176565d-12
      real(8) , parameter :: mev=1.d6*electron_volt
      real(8) , parameter :: mec2=electron_mass*clight**2.0d0
      real(8) , parameter :: mec2mev=mec2/mev
      real(8) , parameter :: electron_radius=electron_charge**2.0d0/mec2
      real(8) , parameter :: sigma_thomson=8.0d0*pi/3.0d0*electron_radius**2.0d0
      real(8) , parameter :: sigma_sb=2.0d0/15.0d0*pi**5.0d0*kboltz**4.0d0/planck**3.0d0/clight**2.0d0
      real(8) , parameter :: arad=4.0d0*sigma_sb/clight
      real(8) , parameter :: fine_structure=electron_charge**2.0d0/(clight*planck/2.0d0/pi)

      real(8) , parameter :: solar_mass=1.98855d33 ! [g]
      real(8) , parameter :: solar_radius=6.955d10 ! [cm]
      real(8) , parameter :: parsec=3.08567758d18 ! [cm]

      real(8) , parameter :: day=86400.0d0

      end module physical_constants
