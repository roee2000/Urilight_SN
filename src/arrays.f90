! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! ARRAYS ALLOCATION
      Module arrays
      use globals , only : isotop_data , max_isotops
      implicit none

      type epacket
        real(8) :: r(3)=0.0d0
        real(8) :: n(3)=0.0d0
        real(8) :: Etot=0.0d0
        real(8) :: hnu=0.0d0
        real(8) :: lam=0.0d0
        real(8) :: t=1.0d99
        logical :: direct=.true.
        logical :: active=.true.
      end type epacket

      real(8) , allocatable :: atoms(:,:)
      real(8) , allocatable :: mass(:)
      real(8) , allocatable :: rhov(:)
      real(8) , allocatable :: nelec(:),zavg(:)

      real(8) , allocatable :: Times(:) , Teff(:) , dt(:)

      real(8) , allocatable :: Edep_gamma(:,:) , Edep_pos(:,:) , Ecr_gamma(:)
      real(8) , allocatable :: spect_gamma(:,:) , spect_bins_gamma(:)

      real(8) , allocatable :: spect_uvoir(:,:) , spect_bins_uvoir(:) ,&
                               dspect_bins_uvoir(:)
      real(8) , allocatable :: jnudnu(:),nujnudnu(:),edep(:),esca(:)
      real(8) , allocatable :: kappa_abs(:),kappa_scat(:)
      real(8) , allocatable :: temp(:),temp_old(:),trad(:),tcolor(:),tplasma(:)

      real(8) , allocatable :: bp(:,:)
      real(8) , allocatable :: emissivity(:,:)
      real(8) , allocatable :: alpha_ff(:,:),alpha_abs_exp(:,:),alpha_scat_exp(:,:)
      real(8) , allocatable :: alpha_scat(:)


      integer , allocatable :: ntracks(:)

      type(epacket) , allocatable :: photon(:)

      type (isotop_data) , allocatable :: iso(:)
      integer , allocatable :: indiso(:,:)

      end module arrays
