! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
      Module Globals
      implicit none


!     simulation data
      integer :: fout
      character(30) :: data_file
      integer :: ntimes
      real(8) :: tinit , tfinal


!     gamma 
      logical :: onlyrodr=.false.
      integer :: N_GammaPellets=500000
      integer :: spect_type_gamma=1
!     uvoir
      logical :: isuvoir=.true.
      integer :: N_UvoirPellets=200000
      integer :: spect_type_uvoir=2
      integer :: nwavelengths=0

!     geometry
      integer , save :: dim=1
!     Mesh data
      integer , save :: nc1=1,nc2=1,nc3=1,nctot=1
      integer , save :: nc1p,nc2p,nc3p
!     ejecta data
      real(8) , save :: Minit,Vej,Rej,t0ej,Vejmax
      integer , save :: ejecta_type=1
!     1 - constant density
!     2 - exponential
!     3 - read from file

      integer , save :: Nitermax=10

      real(8) , parameter :: eps=1.d-32
      real(8) , parameter :: epstemp=100.0d0

      integer , save :: niso
      integer , parameter :: max_isotops=350
      integer , save :: ind_ni56,ind_co56,ind_fe56,ind_fe54,&
                        ind_si28,ind_s32,ind_ar36,ind_ca40,&
                        ind_c12,ind_o16
      type isotop_data
        integer :: Z=0.0d0
        real(8) :: A=0.0d0
        integer :: sym=0
      end type isotop_data


      end module globals
