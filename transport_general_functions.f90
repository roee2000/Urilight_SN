! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
      Module transport_general_functions
      use globals , only : dim 
      use arrays , only : teff
      use physical_constants , only : clight
      use general_functions
      use mesh , only : v1,v2,v3
      implicit none

      contains

      subroutine cut_nearest_surface(rp,np,t,nt,i,j,k,ds,rcut,tcut,i1,j1,k1,isprob)
      real(8) , intent (in) :: rp(3),np(3),t
      integer , intent (in) :: nt,i,j,k
      real(8) , intent (out) :: rcut(3),tcut,ds
      integer , intent (out) :: i1,j1,k1
      logical , intent (out) :: isprob
      real(8) :: dsi(6),vi(6),vin
      integer :: n,ind
      integer , save :: addind(3,6)
      data addind(1,:)/-1,1,0,0,0,0/
      data addind(2,:)/0,0,-1,1,0,0/
      data addind(3,:)/0,0,0,0,-1,1/

      isprob=.false.
      dsi=1.d99
      ds=1.d99
      vi=0.0d0
      i1=i
      j1=j
      k1=k

      if (dim.eq.1) then

      if (i.gt.1) then
        vi(1)=v1(i)*teff(nt)
        dsi(1)=cut_line_sphere(rp(1:3),np(1:3),vi(1))
      endif
        vi(2)=v1(i+1)*teff(nt)
        dsi(2)=cut_line_sphere(rp(1:3),np(1:3),vi(2))


      elseif (dim.eq.2) then
      if (i.gt.1) then
        vi(1)=v1(i)*teff(nt)
        dsi(1)=cut_line_sphere(rp(1:2),np(1:2),v1(i))
      endif
        vi(2)=v1(i+1)*teff(nt)
        dsi(2)=cut_line_sphere(rp(1:2),np(1:2),v1(i+1))

!     vi(3)=v2(j)*times
!     dsi(3)=(rp(3)-v2(j)*t)/(v2(j)/clight-np(3))

!     vi(4)=v2(j+1)*times
!     dsi(4)=(rp(3)-v2(j+1)*t)/(v2(j+1)/clight-np(3))

      elseif (dim.eq.3) then

      endif

      ind=0
      do n=1,2*dim
        if (dsi(n).lt.ds .and. dsi(n).gt.1.d-8*abs(vi(n))) then
          ds=dsi(n)
          ind=n
        endif
      enddo
!     print*,'in r1,r2,rp=',vi(1),vi(2),norm(rp)
      if (ind.eq.0) then
!       print*,'prob'
!       print*,'i,j,k=',i,j,k
!       print*,'r1,r2,rp=',vi(1),vi(2),norm(rp)
!       print*,'rp(:),np(:)=',rp(:),np(:)
!       print*,'dsi=',dsi(1:4)
        isprob=.true.
!       stop
        return
      endif
      i1=i1+addind(1,ind)
      j1=j1+addind(2,ind)
      k1=k1+addind(3,ind)
      

      rcut=rp(:)+np(:)*ds
      tcut=t+ds/clight
!     print*,'out rcut=',norm(rcut)

      return
      end subroutine cut_nearest_surface

      real(8) function cut_line_sphere(rp,np,rin,ierr)
      real(8) , intent (in) :: rp(:),np(:),rin
      real(8) :: rnn,rnrn,nn,AA,BB,CC,det
      real(8) :: ds1,ds2,eps
      integer , optional , intent(in) :: ierr
      
      cut_line_sphere=-1.0d0
      eps=1.d-12*rin

      nn=sum(np(:)*np(:))
      rnn=sum(rp(:)*np(:))
      rnrn=sum(rp(:)*rp(:))
      
      AA=nn
      BB=2.0d0*rnn
      CC=rnrn-rin**2.0d0
      det=BB**2.0d0-4.0d0*AA*CC

      if (det.lt.0d0) then
!       print*,'nocut'
      else
        ds1=(-BB-sqrt(det))/2.0d0/AA
        ds2=(-BB+sqrt(det))/2.0d0/AA
!       print*,'ds1,ds2,epsvdt=',ds1,ds2
        if (present(ierr)) then
!         print*,'ds1,ds2,epsvdt=',ds1,ds2,1.d-14*vi*t
        endif
       
        if (min(ds1,ds2).gt.eps) then
          cut_line_sphere=min(ds1,ds2)
        else
          cut_line_sphere=max(ds1,ds2)
        endif

      endif

      return
      end function cut_line_sphere

      function comoving2lab_transform_E(E,n,v,cmf)
      real(8) , intent (in) :: E
      real(8) , intent (in) :: n(3) , v(3)
      integer , intent (in) :: cmf
      real(8) :: nvdot,Etag,gamma
      real(8) :: comoving2lab_transform_E

      nvdot=sum(n(:)*v(:))
      gamma=1.0d0/sqrt(1.0d0-sum(v(:)*v(:))/clight**2.0d0)

      if (cmf.eq.1) then  !!! comoving to lab transform

!       Etag=E*(1.0d0+nvdot/clight)

        Etag=gamma*E*(1.0d0+nvdot/clight)

      elseif (cmf.eq.2) then !!! lab to comoving transform
  
!       Etag=E*(1.0d0-nvdot/clight)

        Etag=gamma*E*(1.0d0-nvdot/clight)
   
      endif

      comoving2lab_transform_E=Etag

      return
      end function comoving2lab_transform_E

      function comoving2lab_transform_n(n,v,cmf)
      real(8) , intent (in) :: n(3) , v(3)
      integer , intent (in) :: cmf
      real(8) :: nvdot,ntag(3),gamma
      real(8) :: comoving2lab_transform_n(3)

      nvdot=sum(n(:)*v(:))
      gamma=1.0d0/sqrt(1.0d0-sum(v(:)*v(:))/clight**2.0d0)

      if (cmf.eq.1) then  !!! comoving to lab transform

!       ntag=(n+v/clight)/(1.0d0+nvdot/clight)

        ntag=n+gamma*v/clight*(1.0d0+gamma/(1.0d0+gamma)*nvdot/clight)
        ntag=ntag/(gamma*(1.0d0+nvdot/clight))

      elseif (cmf.eq.2) then !!! lab to comoving transform
  
!       ntag=(n-v/clight)/(1.0d0-nvdot/clight)

        ntag=n-gamma*v/clight*(1.0d0-gamma/(1.0d0+gamma)*nvdot/clight)
        ntag=ntag/(gamma*(1.0d0-nvdot/clight))
   
      endif

      comoving2lab_transform_n=ntag/sqrt(sum(ntag(:)**2.0d0))

      return
      end function comoving2lab_transform_n

      function new_direction(n,cost,phi)
      real(8) , intent(in) :: n(3),cost,phi
      real(8) :: new_direction(3),ni(3)
      real(8) :: sint,cosphi,sinphi
      real(8) :: eps,sgn,det

      eps=1.d0-1.d-6

      sgn=n(3)/abs(n(3))

      sint=sqrt(1.0d0-cost**2.0d0)
      cosphi=cos(phi)
      sinphi=sin(phi)

      if (abs(n(3)).gt.eps) then
        ni(1)=sint*cosphi
        ni(2)=sgn*sint*sinphi
        ni(3)=sgn*cost
      else
        det=sqrt(1.0d0-n(3)**2.0d0)
        ni(1)=n(1)*cost+(n(1)*n(3)*cosphi-n(2)*sinphi)*sint/det
        ni(2)=n(2)*cost+(n(2)*n(3)*cosphi+n(1)*sinphi)*sint/det
        ni(3)=n(3)*cost-cosphi*sint*det
      endif

      new_direction=ni
      return
      end function new_direction

      end Module transport_general_functions
