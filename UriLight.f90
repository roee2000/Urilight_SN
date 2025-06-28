! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! THIS IS THE PROGRAM, REFER TO THE VARIOUS MODULES FOR FURTHER
! DESCRIPTIONS

      program UriLight
      use RandomNumbers
      use Atomic_Physics
      use Radioactive_Decay
      use Mesh
      use GammaTransfer
      use UvoirTransfer
      use diagnostics
      implicit none

      call init_simulation
      call init_random_numbers
      call init_nuclear_data
      call init_mesh
      call init_atomic_data
      call init_gamma
      call init_diagnostics

      call gamma_transport

      if (isuvoir) then
        call init_uvoir
        call uvoir_transport
      endif

      call diag_write_totals

      contains

      subroutine init_simulation
      use globals
      use arrays , only : times , indiso
      implicit none
      integer :: i,ino
      namelist /simulation/tinit,tfinal,ntimes,isuvoir

      fout=66
      data_file='data_file'
      open (unit=fout,file='out',position='append')

      write(fout,*) 'Initializing simualtion'
      write(fout,*) '-----------------------'

      tinit=2.0d0
      tfinal=50.0d0
      ntimes=100

      open(unit=5,file=data_file)
      read(5,nml=simulation,iostat=ino)
      close(5)
      write(fout,nml=simulation)

      tinit=tinit*day
      tfinal=tfinal*day

      allocate(times(ntimes+1),teff(ntimes),dt(ntimes))
      times=0.0d0
      
      do i=1,ntimes+1
        times(i)=log(tinit)+(log(tfinal)-log(tinit))/dble(ntimes)*dble(i-1)
        times(i)=exp(times(i))
      enddo

      do i=1,ntimes
        teff(i)=sqrt(times(i)*times(i+1))
        dt(i)=times(i+1)-times(i)
      enddo

      return
      end subroutine init_simulation

      end program UriLight
