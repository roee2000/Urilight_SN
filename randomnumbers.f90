! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
       Module RandomNumbers
       use physical_constants , only : pi
       implicit none



       contains

       subroutine init_random_numbers
       real(8) :: zrand
       real(8) :: zrandvec(1:3)
      
!      call random_seed

       call random_number(zrand)

       end subroutine init_random_numbers



       integer function choose_from_probability_distribution(pdf)
       real(8) , intent(in) , dimension (:) :: pdf
       real(8) :: zrand,cumprob
       integer n,i

       n=size(pdf)
       i=0
       cumprob=0.0d0
       

       call random_number(zrand)

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

       call random_number(z)

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
