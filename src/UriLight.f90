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
      character(8) :: date
      character(10) :: time

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

      call diag_write_totals

      call system_clock(count_end, count_rate, count_max)
      if (count_rate > 0) then
        elapsed_sec = dble(count_end - count_start) / dble(count_rate)
      else
        elapsed_sec = 0.0d0
      endif
      call date_and_time(date, time)
      write(fout,'(A)') ''
      write(fout,'(A)') '=============================================================================='
      write(fout,'(A)') ' Run complete'
      write(fout,'(A)') '=============================================================================='
      write(fout,'(A,I0.2,A,I0.2,A,I0.2)') ' Total runtime = ', int(elapsed_sec)/3600, ':', mod(int(elapsed_sec),3600)/60, ':', mod(int(elapsed_sec),60)
      write(fout,'(A,F12.2,A)') '             (', elapsed_sec, ' s)'
      write(fout,'(A,A4,A,A2,A,A2,A,A2,A,A2,A,A2)') ' Ended: ', date(1:4), '-', date(5:6), '-', date(7:8), ' ', time(1:2), ':', time(3:4), ':', time(5:6)
      write(fout,'(A)') '=============================================================================='

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
