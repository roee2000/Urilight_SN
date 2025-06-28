! URILIGHT - MONTE CARLO RADIATIVE TRANSFER SIMULATION
! WRITTEN BY Y. ELBAZ
! CODE DESCRIPTION AND EXAMPLES IN WYGODA, ELBAZ AND KATZ 2018.
!
! DIAGNOSTICS MODULE
       Module diagnostics
       use globals
       use arrays , only : epacket , times , dt , Edep_gamma , &
                           edep_pos , temp , nelec , emissivity
       use general_functions
       use mesh
       implicit none

       integer , save :: nuvoir_bands=8
       character(3) , save :: uvoir_band_names(8)

       integer :: Nbins,Ntetas
       real(8) :: emin_lc,emax_lc    !! minimal and maximal LightCurve energies [Mev]
       real(8) , allocatable , save :: eobs(:) , ctetaobs(:)
       real(8) , allocatable , save :: Gout(:),Bolout(:)
       real(8) , allocatable , save :: uvoir_f(:,:)
       real(8) , allocatable , save :: directspectrum(:,:,:),indirectspectrum(:,:,:)
       integer , allocatable , save :: counts(:,:,:)

       contains

       subroutine init_diagnostics
       integer :: i
       real(8) :: logemin,logemax

       nuvoir_bands=8
       uvoir_band_names(1)='U'
       uvoir_band_names(2)='B'
       uvoir_band_names(3)='V'
       uvoir_band_names(4)='R'
       uvoir_band_names(5)='I'
       uvoir_band_names(6)='J'
       uvoir_band_names(7)='H'
       uvoir_band_names(8)='K'

       emin_lc=1.d-32
       emax_lc=4.0d0

       logemin=log10(emin_lc)
       logemax=log10(emax_lc)
    
       Nbins=200
       Ntetas=10
       allocate(eobs(Nbins+1),ctetaobs(Ntetas+1))
       allocate(directspectrum(Ntimes,Nbins,Ntetas))
       allocate(indirectspectrum(Ntimes,Nbins,Ntetas))
       allocate(counts(Ntimes,Nbins,Ntetas))
       allocate(Gout(Ntimes),Bolout(Ntimes))
       allocate(uvoir_f(Ntimes,nuvoir_bands))
       eobs=0.0d0
       ctetaobs=0.0d0
       directspectrum=0.0d0
       indirectspectrum=0.0d0
       counts=0.0d0
       Gout=0.0d0
       Bolout=0.0d0
       uvoir_f=0.0d0
 
       do i=1,Nbins+1
         eobs(i)=emin_lc+(emax_lc-emin_lc)/dble(Nbins)*dble(i-1)
       enddo
       do i=1,Ntetas+1
         ctetaobs(i)=-1.0d0+2.0d0/dble(Ntetas)*dble(i-1)
       enddo

       return
       end subroutine init_diagnostics

       subroutine write_to_spectrum(p,spect,bins,spect_type)
       type (epacket) , intent (in) :: p
       real(8) , intent(in) :: bins(:)
       integer , intent(in) :: spect_type
       real(8) , intent(inout) :: spect(:,:)
       real(8) :: tau,ctet,delta
       integer :: i1,i2,i3,ierr,nbins

       nbins=size(bins)-1

       tau=p%t-sum(p%r(:)*p%n(:))/clight
       ctet=p%n(3)/sqrt(sum(p%n(:)**2.0d0))

       if (spect_type.eq.1) then
         i1=find_index(p%hnu,bins(1:nbins+1),ierr)
       elseif (spect_type.eq.2) then
         i1=find_index(p%lam,bins(1:nbins+1),ierr)
       endif

       if (ierr.gt.0) return

       i2=find_index(tau,times(1:ntimes+1),ierr)

       if (ierr.gt.0) return
!      i3=1
!      do while (ctet.gt.ctetaobs(i3))
!        i3=i3+1
!      enddo
!      i3=i3-1
      
       delta=(bins(i1+1)-bins(i1))*(times(i2+1)-times(i2))

       spect(i1,i2)=spect(i1,i2)+p%Etot/delta

!      if (p%direct) then
!        directspectrum(i1,i2,i3)=directspectrum(i1,i2,i3)+p%Etot/p%hnu
!      else
!        indirectspectrum(i1,i2,i3)=indirectspectrum(i1,i2,i3)+p%Etot/p%hnu
!      endif
!      counts(i1,i2,1)=counts(i1,i2,1)+1

       return
       end subroutine write_to_spectrum

       subroutine integrate_bolometric(p,vec,method)
       type (epacket) , intent (in) :: p
       real(8) , intent (inout) :: vec(:)
       integer , intent (in) :: method
       real(8) :: tau,delta
       integer :: i1,i2,i3,ierr

       if (method.eq.1) then
         tau=p%t
       else
         tau=p%t-sum(p%r(:)*p%n(:))/clight
       endif

       i1=find_index(tau,times(1:ntimes+1),ierr)

       if (ierr.gt.0) return
      
       delta=(times(i1+1)-times(i1))
       vec(i1)=vec(i1)+p%Etot/delta

       return
       end subroutine integrate_bolometric

       subroutine diag_integrate_uvoir_bands(p,spect_bins,method)
       type (epacket) , intent (in) :: p
       real(8) , intent (in) :: spect_bins(:)
       integer , intent (in) :: method
       real(8) :: tau,delta,dlam,phi,phin
       integer :: nb,i1,ierr,ilam

       if (method.eq.1) then
         tau=p%t
       else
         tau=p%t-sum(p%r(:)*p%n(:))/clight
       endif

       i1=find_index1(tau,times(1:ntimes+1),ierr)
 
       if (ierr.gt.0) return
       delta=(times(i1+1)-times(i1))

       ilam=find_index1(p%lam,spect_bins(:),ierr)
       dlam=(spect_bins(ilam+1)-spect_bins(ilam))/angstrom

       do nb=1,size(uvoir_band_names)
         phi=spectroscopic_filter(p%lam,uvoir_band_names(nb),phin)
         uvoir_f(i1,nb)=uvoir_f(i1,nb)+p%Etot/delta*phi
       enddo
      
       return
       end subroutine diag_integrate_uvoir_bands

      subroutine diag_write_spectrum(namef,spect,bins)
      real(8) , intent(in) :: spect(:,:),bins(:)
      character*20 , intent(in) :: namef
      integer :: n,nbins

      nbins=size(bins)-1
      
      open (unit=101,file=namef,position='append')
      write(101,1000) -1.0d0,times(1:ntimes+1)
      do n=1,nbins
        write(101,1000) bins(n),bins(n+1),spect(n,1:ntimes)
      enddo
      close (101)

1000  format(1000(1pe14.6))
      return
      end subroutine diag_write_spectrum

      subroutine diag_write_totals
      integer :: i
      real(8) :: edep,epos,delta,gcr

      open (unit=101,file='totals',position='append')
      open (unit=102,file='luminocity',position='append')
      open (unit=103,file='magnitude',position='append')
      do i=1,Ntimes
        delta=(times(i+1)-times(i))
        edep=sum(edep_gamma(i,:))/delta
        epos=sum(edep_pos(i,:))/delta
        gcr=ecr_gamma(i)/delta
        write(101,1000) times(i),times(i+1),gcr,edep,epos,Gout(i),Bolout(i)
        write(102,1000) times(i),times(i+1),Bolout(i),uvoir_f(i,:)
        write(103,1000) times(i),times(i+1),bolometric_magnitude(Bolout(i)),&
            magnitude(uvoir_f(i,:),uvoir_band_names(:))
      enddo
      close(101)
      close(102)
      close(103)

1000  format(1000(1pe14.6))

      return
      end subroutine diag_write_totals

      subroutine diag_write_profiles(nt)
      integer , intent(in) :: nt
      integer :: nb
      character*20 :: namefile,chf

      IF (NT.LT.10) THEN
        WRITE(CHF,1101) NT
      ELSE IF (NT.LT.100) THEN
        WRITE(CHF,1102) NT
      ELSE
        WRITE(CHF,1103) NT
      ENDIF
        namefile='profile'//chf

      OPEN(UNIT=50,file=namefile,POSITION='APPEND')

      write(50,1000) (edep_gamma(nt,1:nc1)+edep_pos(nt,1:nc1))/dt(nt)
      write(50,1000) nelec(1:nc1)
      write(50,1000) zavg(1:nc1)
      write(50,1000) trad(1:nc1)
      write(50,1000) tcolor(1:nc1)
      write(50,1000) tplasma(1:nc1)
      write(50,1000) kappa_abs(1:nc1)
      write(50,1000) kappa_scat(1:nc1)
      do nb=1,nuvoir_bands
      write(50,1000) emissivity(nb,1:nc1)
      enddo

      close(50)
1000  format(1000(1pe14.6))
1101  FORMAT('-',I3.3)
1102  FORMAT('-',I3.3)
1103  FORMAT('-',I3.3)

      return
      end subroutine diag_write_profiles

      real(8) function spectroscopic_filter(lam,band,totint)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    A function calculating spectroscopic bands transmission     !!
!!    fractions.                                                  !!
!!    Sources for bands data:                                     !!
!!      UBVRI - Bessell, M. S. 1990, PASP, 102, 1181              !!
!!      JHK   - Persson, S. E. et al.  1998, AJ, 116, 2475        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) , intent(in) :: lam
      character(3) , intent(in) :: band
      real(8) , intent(out) :: totint
      real(8) :: la
      integer :: i,ierr

      real(8) , save :: I_filter(23,2) , R_filter(24,2) , V_filter(24,2), &
                        B_filter(21,2) , U_filter(25,5)
      real(8) , save :: J_filter(32,2) , H_filter(62,2) , K_filter(73,2)
      real(8) , save :: Uint,Bint,Vint,Rint,Iint,Jint,Hint,Kint

      data I_filter(1:23,1)/7000.d0,7100.d0,7200.d0,7300.d0,7400.d0,&
                            7500.d0,7600.d0,7700.d0,7800.d0,7900.d0,&
                            8000.d0,8100.d0,8200.d0,8300.d0,8400.d0,&
                            8500.d0,8600.d0,8700.d0,8800.d0,8900.d0,&
                            9000.d0,9100.d0,9200.d0/
      data I_filter(1:23,2)/0.000d0,0.024d0,0.232d0,0.555d0,0.785d0,&
                            0.910d0,0.965d0,0.985d0,0.990d0,0.995d0,&
                            1.000d0,1.000d0,0.990d0,0.980d0,0.950d0,&
                            0.910d0,0.860d0,0.750d0,0.560d0,0.330d0,&
                            0.150d0,0.030d0,0.000d0/
      data R_filter(1:24,1)/5500.d0,5600.d0,5700.d0,5800.d0,5900.d0,&
                            6000.d0,6100.d0,6200.d0,6300.d0,6400.d0,&
                            6500.d0,6600.d0,6700.d0,6800.d0,6900.d0,&
                            7000.d0,7100.d0,7200.d0,7300.d0,7400.d0,&
                            7500.d0,8000.d0,8500.d0,9000.d0/
      data R_filter(1:24,2)/0.00d0,0.23d0,0.74d0,0.91d0,0.98d0,&
                            1.00d0,0.98d0,0.96d0,0.93d0,0.90d0,&
                            0.86d0,0.81d0,0.78d0,0.72d0,0.67d0,&
                            0.61d0,0.56d0,0.51d0,0.46d0,0.40d0,&
                            0.35d0,0.14d0,0.03d0,0.00d0/
      data V_filter(1:24,1)/4700.d0,4800.d0,4900.d0,5000.d0,5100.d0,&
                            5200.d0,5300.d0,5400.d0,5500.d0,5600.d0,&
                            5700.d0,5800.d0,5900.d0,6000.d0,6100.d0,&
                            6200.d0,6300.d0,6400.d0,6500.d0,6600.d0,&
                            6700.d0,6800.d0,6900.d0,7000.d0/
      data V_filter(1:24,2)/0.000d0,0.030d0,0.163d0,0.458d0,0.780d0,&
                            0.967d0,1.000d0,0.973d0,0.898d0,0.792d0,&
                            0.684d0,0.574d0,0.461d0,0.359d0,0.270d0,&
                            0.197d0,0.135d0,0.081d0,0.045d0,0.025d0,&
                            0.017d0,0.013d0,0.009d0,0.000d0/
      data B_filter(1:21,1)/3600.d0,3700.d0,3800.d0,3900.d0,4000.d0,&
                            4100.d0,4200.d0,4300.d0,4400.d0,4500.d0,&
                            4600.d0,4700.d0,4800.d0,4900.d0,5000.d0,&
                            5100.d0,5200.d0,5300.d0,5400.d0,5500.d0,&
                            5600.d0/
      data B_filter(1:21,2)/0.000d0,0.030d0,0.134d0,0.567d0,0.920d0,&
                            0.978d0,1.000d0,0.978d0,0.935d0,0.853d0,&
                            0.740d0,0.640d0,0.536d0,0.424d0,0.325d0,&
                            0.235d0,0.150d0,0.095d0,0.043d0,0.009d0,&
                            0.000d0/
      data U_filter(1:25,1)/3000.d0,3050.d0,3100.d0,3150.d0,3200.d0,&
                            3250.d0,3300.d0,3350.d0,3400.d0,3450.d0,&
                            3500.d0,3550.d0,3600.d0,3650.d0,3700.d0,&
                            3750.d0,3800.d0,3850.d0,3900.d0,3950.d0,&
                            4000.d0,4050.d0,4100.d0,4150.d0,4200.d0/
      data U_filter(1:25,2)/0.000d0,0.016d0,0.068d0,0.167d0,0.287d0,&
                            0.423d0,0.560d0,0.673d0,0.772d0,0.841d0,&
                            0.905d0,0.943d0,0.981d0,0.993d0,1.000d0,&
                            0.989d0,0.916d0,0.804d0,0.625d0,0.423d0,&
                            0.238d0,0.114d0,0.051d0,0.019d0,0.000d0/
      data J_filter(1:32,1)/10800.d0,10900.d0,11000.d0,11100.d0,11200.d0,&
                            11300.d0,11400.d0,11500.d0,11600.d0,11700.d0,&
                            11800.d0,11900.d0,12000.d0,12100.d0,12200.d0,&
                            12300.d0,12400.d0,12500.d0,12600.d0,12700.d0,&
                            12800.d0,12900.d0,13000.d0,13100.d0,13200.d0,&
                            13300.d0,13400.d0,13500.d0,13600.d0,13700.d0,&
                            13800.d0,13900.d0/
      data J_filter(1:32,2)/0.000d0,0.025d0,0.065d0,0.135d0,0.320d0,&
                            0.555d0,0.720d0,0.785d0,0.850d0,0.860d0,&
                            0.860d0,0.865d0,0.875d0,0.883d0,0.892d0,&
                            0.900d0,0.905d0,0.905d0,0.900d0,0.902d0,&
                            0.906d0,0.910d0,0.915d0,0.918d0,0.927d0,&
                            0.925d0,0.875d0,0.495d0,0.140d0,0.045d0,&
                            0.015d0,0.000d0/
      data H_filter(1:62,1)/13100.d0,13200.d0,13300.d0,13400.d0,13500.d0,&
                            13600.d0,13700.d0,13800.d0,13900.d0,14000.d0,&
                            14100.d0,14200.d0,14300.d0,14400.d0,14500.d0,&
                            14600.d0,14700.d0,14800.d0,14900.d0,15000.d0,&
                            15100.d0,15200.d0,15300.d0,15400.d0,15500.d0,&
                            15600.d0,15700.d0,15800.d0,15900.d0,16000.d0,&
                            16100.d0,16200.d0,16300.d0,16400.d0,16500.d0,&
                            16600.d0,16700.d0,16800.d0,16900.d0,17000.d0,&
                            17100.d0,17200.d0,17300.d0,17400.d0,17500.d0,&
                            17600.d0,17700.d0,17800.d0,17900.d0,18000.d0,&
                            18100.d0,18200.d0,18300.d0,18400.d0,18500.d0,&
                            18600.d0,18700.d0,18800.d0,18900.d0,19000.d0,&
                            19100.d0,19200.d0/
      data H_filter(1:62,2)/0.000d0,0.001d0,0.002d0,0.003d0,0.004d0,&
                            0.004d0,0.005d0,0.006d0,0.007d0,0.008d0,&
                            0.012d0,0.016d0,0.024d0,0.032d0,0.063d0,&
                            0.087d0,0.131d0,0.178d0,0.261d0,0.396d0,&
                            0.594d0,0.653d0,0.693d0,0.721d0,0.725d0,&
                            0.729d0,0.737d0,0.741d0,0.746d0,0.745d0,&
                            0.739d0,0.737d0,0.737d0,0.737d0,0.737d0,&
                            0.721d0,0.709d0,0.713d0,0.729d0,0.745d0,&
                            0.748d0,0.733d0,0.661d0,0.622d0,0.594d0,&
                            0.598d0,0.634d0,0.665d0,0.673d0,0.665d0,&
                            0.436d0,0.218d0,0.099d0,0.020d0,0.014d0,&
                            0.012d0,0.009d0,0.008d0,0.006d0,0.004d0,&
                            0.002d0,0.000d0/
      data K_filter(1:73,1)/18400.d0,18500.d0,18600.d0,18700.d0,18800.d0,&
                            18900.d0,19000.d0,19100.d0,19200.d0,19300.d0,&
                            19400.d0,19500.d0,19600.d0,19700.d0,19800.d0,&
                            19900.d0,20000.d0,20100.d0,20200.d0,20300.d0,&
                            20400.d0,20500.d0,20600.d0,20700.d0,20800.d0,&
                            20900.d0,21000.d0,21100.d0,21200.d0,21300.d0,&
                            21400.d0,21500.d0,21600.d0,21700.d0,21800.d0,&
                            21900.d0,22000.d0,22100.d0,22200.d0,22300.d0,&
                            22400.d0,22500.d0,22600.d0,22700.d0,22800.d0,&
                            22900.d0,23000.d0,23100.d0,23200.d0,23300.d0,&
                            23400.d0,23500.d0,23600.d0,23700.d0,23800.d0,&
                            23900.d0,24000.d0,24100.d0,24200.d0,24300.d0,&
                            24400.d0,24500.d0,24600.d0,24700.d0,24800.d0,&
                            24900.d0,25000.d0,25100.d0,25200.d0,25300.d0,&
                            25400.d0,25500.d0,25600.d0/
      data K_filter(1:73,2)/0.000d0,0.003d0,0.004d0,0.006d0,0.008d0,&
                            0.012d0,0.016d0,0.020d0,0.024d0,0.038d0,&
                            0.044d0,0.059d0,0.071d0,0.095d0,0.115d0,&
                            0.139d0,0.178d0,0.277d0,0.356d0,0.491d0,&
                            0.566d0,0.610d0,0.653d0,0.657d0,0.653d0,&
                            0.645d0,0.661d0,0.693d0,0.713d0,0.748d0,&
                            0.772d0,0.784d0,0.792d0,0.800d0,0.800d0,&
                            0.802d0,0.802d0,0.802d0,0.802d0,0.806d0,&
                            0.810d0,0.806d0,0.814d0,0.812d0,0.810d0,&
                            0.804d0,0.796d0,0.786d0,0.772d0,0.744d0,&
                            0.729d0,0.705d0,0.709d0,0.713d0,0.721d0,&
                            0.744d0,0.760d0,0.768d0,0.693d0,0.554d0,&
                            0.475d0,0.345d0,0.198d0,0.107d0,0.079d0,&
                            0.048d0,0.020d0,0.020d0,0.012d0,0.010d0,&
                            0.008d0,0.004d0,0.000d0/
      data uint/640.4d-8/,bint/959.2d-8/,vint/893.1d-8/,rint/1591.0d-8/,&
           iint/1495.1d-8/,jint/2027.3d-8/,hint/2304.9d-8/,kint/3289.3d-8/

      spectroscopic_filter=0.0d0

      la=lam/angstrom

      select case (band)
      
        case('I')
          totint=iint
          i=find_index(la,I_filter(:,1),ierr)
          if (ierr.gt.0) return
          spectroscopic_filter=interp1(i_filter(i:i+1,1),i_filter(i:i+1,2),la)
        case('R')
          totint=rint
          i=find_index(la,R_filter(:,1),ierr)
          if (ierr.gt.0) return
          spectroscopic_filter=interp1(r_filter(i:i+1,1),r_filter(i:i+1,2),la)
        case('V')
          totint=vint
          i=find_index(la,V_filter(:,1),ierr)
          if (ierr.gt.0) return
          spectroscopic_filter=interp1(v_filter(i:i+1,1),v_filter(i:i+1,2),la)
        case('B')
          totint=bint
          i=find_index(la,B_filter(:,1),ierr)
          if (ierr.gt.0) return
          spectroscopic_filter=interp1(b_filter(i:i+1,1),b_filter(i:i+1,2),la)
        case('U')
          totint=uint
          i=find_index(la,U_filter(:,1),ierr)
          if (ierr.gt.0) return
          spectroscopic_filter=interp1(u_filter(i:i+1,1),u_filter(i:i+1,2),la)
        case('J')
          totint=jint
          i=find_index(la,J_filter(:,1),ierr)
          if (ierr.gt.0) return
          spectroscopic_filter=interp1(j_filter(i:i+1,1),j_filter(i:i+1,2),la)
        case('H')
          totint=hint
          i=find_index(la,H_filter(:,1),ierr)
          if (ierr.gt.0) return
          spectroscopic_filter=interp1(h_filter(i:i+1,1),h_filter(i:i+1,2),la)
        case('K')
          totint=kint
          i=find_index(la,k_filter(:,1),ierr)
          if (ierr.gt.0) return
!         if (ierr.gt.0 .or. la.gt.24200d0) return
          spectroscopic_filter=interp1(k_filter(i:i+1,1),k_filter(i:i+1,2),la)

      end select
 
      return
      end function spectroscopic_filter

      function magnitude(L,band) result(mag)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Given total luminocity L filtered in band (or total luminocity)  !!
!!    the function returns the absolute magnitude.                     !!
!!    zero point fluxes are taken from M.S. Bessell et. al.,           !!
!!    Astron. Astrophys. 333, 231.250 (1998).                          !!
!!    NOTE: there is a type in table A2 and zp(f_lam)<->zp(f_nu) !     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) , intent(in) :: L(:)
      real(8) :: mag(size(L))
      character(3) :: band(:)
      real(8) :: surf,lam,totphi,flux,phi,zpf
      integer :: i
      real(8) , save :: zpfu,zpfb,zpfv,zpfr,zpfi,zpfj,zpfh,zpfk
      real(8) , save :: vega0
      data vega0/21.1d0/
      data zpfu/-0.152d0/,zpfb/-0.602d0/,zpfv/0.0d0/,zpfr/0.555d0/,&
           zpfi/1.271d0/,zpfj/2.655d0/,zpfh/3.760d0/,zpfk/4.906d0/
      
      mag=0.0d0
      lam=1000.0d0*angstrom
      surf=4.0d0*pi*(10.0d0*parsec)**2.0d0

      do i=1,size(L)

      phi=spectroscopic_filter(lam,band(i),totphi)

      select case (band(i))
        case('I')
          zpf=zpfi
        case('R')
          zpf=zpfr
        case('V')
          zpf=zpfv
        case('B')
          zpf=zpfb
        case('U')
          zpf=zpfu
        case('J')
          zpf=zpfj
        case('H')
          zpf=zpfh
        case('K')
          zpf=zpfk
      end select
!!    translate to flux density in units erg/s/cm^2/Angstrom
      flux=L(i)/surf/(totphi/angstrom)

      if (flux.gt.1.d-90) then
        mag(i)=-2.5d0*log10(flux)-vega0-zpf
      else
        mag(i)=0.0d0
      endif

      enddo

      return
      end function magnitude

      real(8) function bolometric_magnitude(L)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Given total bolometric luminocity L                              !!
!!    the function returns the absolute magnitude.                     !!
!!    zero point fluxes are taken from M.S. Bessell et. al.,           !!
!!    Astron. Astrophys. 333, 231.250 (1998).                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) , intent(in) :: L
      real(8) , save :: Lsun,msun
      data lsun/3.855d33/,msun/4.74d0/

      if (l.gt.1.d-90) then
        bolometric_magnitude=msun-2.5d0*log10(L/Lsun)
      else
        bolometric_magnitude=100.0d0
      endif

      return
      end function bolometric_magnitude

      end module diagnostics
