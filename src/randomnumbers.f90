! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
       Module RandomNumbers
       use physical_constants , only : pi
       use omp_lib
       implicit none

       integer(8) , save :: rng_state = 1_8
       !$omp threadprivate(rng_state)

       interface rand_number
         module procedure rand_number_0 , rand_number_1
       end interface

       contains

       subroutine init_random_numbers
       integer :: tid
       integer(8) :: seed0 , t

       call system_clock(t)
       if (t.le.0_8) t = 123456789_8
       seed0 = t

       !$omp parallel private(tid)
       tid = omp_get_thread_num()
       rng_state = seed0 + int(tid + 1 , 8) * 978641_8
       if (rng_state.eq.0_8) rng_state = 1_8
       !$omp end parallel

       end subroutine init_random_numbers

       subroutine rand_number_0(x)
       real(8) , intent(out) :: x
       integer(8) :: s
       integer(8) , parameter :: mask53 = 9007199254740991_8

       s = rng_state
       s = ieor(s , ishft(s , 13))
       s = ieor(s , ishft(s , -7))
       s = ieor(s , ishft(s , 17))
       rng_state = s
       x = dble(iand(s , mask53)) * (1.0d0 / 9007199254740992.0d0)
       if (x.le.0.0d0) x = 1.0d-16
       if (x.ge.1.0d0) x = 1.0d0 - 1.0d-16

       end subroutine rand_number_0

       subroutine rand_number_1(x)
       real(8) , intent(out) :: x(:)
       integer :: i
       do i = 1 , size(x)
         call rand_number_0(x(i))
       enddo
       end subroutine rand_number_1

       integer function choose_from_probability_distribution(pdf)
       real(8) , intent(in) , dimension (:) :: pdf
       real(8) :: zrand,cumprob
       integer n,i

       n=size(pdf)
       i=0
       cumprob=0.0d0

       call rand_number(zrand)

       do while (zrand>cumprob)
         i=i+1
         if (i.gt.n) then
           print*,'err'
           stop
         endif
         cumprob=cumprob+pdf(i)
       enddo

       choose_from_probability_distribution=i

       return
       end function choose_from_probability_distribution

       function random_unit_vec1(n)
       integer , intent(in) :: n
       real(8) :: random_unit_vec1(n)
       real(8) :: z(2) , vec(3) , vnorm , phi , teta

       call rand_number(z)

       phi=2.0d0*pi*z(1)
       z(2)=1.0d0-2.0d0*z(2)

       if (n.eq.2) z(2)=0.0d0

       teta=acos(z(2))

       vec(1)=sin(teta)*cos(phi)
       vec(2)=sin(teta)*sin(phi)
       vec(3)=cos(teta)

       vnorm=sqrt(sum(vec(:)**2.0d0))
       random_unit_vec1=vec(1:n)/vnorm

       return
       end function random_unit_vec1


       end Module RandomNumbers
