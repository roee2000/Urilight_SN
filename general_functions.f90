! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
       module general_functions
       implicit none




       contains
 
       integer function find_index(x,vec,ierr)
       implicit none
       real(8) , intent(in) :: x,vec(:)
       integer , intent(out) :: ierr
       integer :: i,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     given x and a vector vec, this function returns index i for which
!!     vec(i)<x<=vec(i+1).
!!     if x<vec(1)          function returns i=1         and ierr=1
!!     if x>vec(size(vec))  function returns i=size(vec) and ierr=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ierr=0

       n=size(vec)

       if (x.lt.vec(1)) then
         ierr=1
         find_index=1
       elseif (x.gt.vec(n)) then
         ierr=2
         find_index=n-1
       endif

       if (ierr.ne.0) return

       i=1
       do while (x.gt.vec(i+1))
         i=i+1
       enddo

       find_index=i

       return
       end function find_index

       recursive function find_index1(x,vec,ierr) result(indx)
       implicit none
       real(8) , intent(in) :: x,vec(:)
       integer , intent(out) :: ierr
       integer :: indx
       integer :: i,n,imid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     given x and a vector vec, this function returns index i for which
!!     vec(i)<x<=vec(i+1).
!!     if x<vec(1)          function returns i=1         and ierr=1
!!     if x>vec(size(vec))  function returns i=size(vec) and ierr=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ierr=0

       n=size(vec)

       if (x.lt.vec(1)) then
         ierr=1
         indx=1
       elseif (x.gt.vec(n)) then
         ierr=2
         indx=n-1
       endif

       if (ierr.ne.0) return

       if (n.eq.2) then
         indx=1
         return
       endif

       imid=n/2+1

       if (x.lt.vec(imid)) then
         indx=find_index1(x,vec(1:imid),ierr)
       else
         indx=(imid-1)+find_index1(x,vec(imid:n),ierr)
       endif
         
       return
       end function find_index1

       real(8) function norm(x)
       implicit none
       real(8) :: x(:)

       norm=sqrt(sum(x(:)**2.0d0))

       return
       end function norm

      integer function count_lines(num)
      integer , intent(in) :: num
      integer :: eof
      character*1 :: nothing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    simple io function that counts the number of lines in a file
!!    with unit=num
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rewind(num)
      
      eof=0
      count_lines=0
      do while (eof.eq.0)
        read(num,*,IOSTAT=eof) nothing
        if (eof.eq.0) then
        count_lines=count_lines+1
        endif
      enddo

      rewind(num)

      return
      end function count_lines

      integer function find(x,vec,length,ierr)
      real(8) , intent(in) :: x,vec(:)
      integer , intent(in) :: length
      integer , intent(out) :: ierr
      integer :: i
      logical :: found

      ierr=0

      i=1
      do while (i.le.length .and. .not.found)
        if (x.eq.vec(i)) found=.true.
        i=i+1
      enddo

      if (found) then
        find=i
      else
        find=-1
        ierr=1
      endif

      return
      end function find

      real(8) function interp1(x,y,x1)
      real(8) , intent(in) :: x(2),y(2),x1

      interp1=y(1)+(y(2)-y(1))/(x(2)-x(1))*(x1-x(1))

      return
      end function interp1

      end module general_functions
    

       
