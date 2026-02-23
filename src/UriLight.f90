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
      use Logger
      use globals, only: fout
      implicit none
      integer :: count_start, count_end, count_rate, count_max
      real(8) :: elapsed_sec

      call system_clock(count_start, count_rate, count_max)

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

      call system_clock(count_end, count_rate, count_max)
      elapsed_sec = merge(dble(count_end - count_start) / dble(count_rate), 0.0d0, count_rate > 0)
      call diag_write_totals(elapsed_sec=elapsed_sec)

      contains

      subroutine init_simulation
      use globals
      use arrays , only : times , indiso
      implicit none
      integer :: i,ino
      integer :: ios
      namelist /simulation/tinit,tfinal,ntimes,isuvoir,freeze_composition,output_dir

      fout=66
      data_file='data_file'
      output_dir='output'
      tinit=2.0d0
      tfinal=50.0d0
      ntimes=100
      open(unit=5,file=data_file)
      read(5,nml=simulation,iostat=ino)
      close(5)

      call execute_command_line('mkdir -p '//trim(output_dir), wait=.true., exitstat=ios)
      if (ios /= 0) then
        print *, 'Error creating output directory'
        stop
      endif
      open (unit=fout,file=trim(output_dir)//'/out',position='append')

      call log_info('Simulation started')
      call log_info('Initializing simulation')
      call log_info('-----------------------')

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
