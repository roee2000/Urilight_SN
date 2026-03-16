! NIST atomic DB reader (out_gg format only).
! Levels from out_gg/outglv_*; lines from out_gg/outggf_sorted_*; each line connected to closest level by energy.
      module read_nist_db_mod
      use atomic_data_types
      use physical_constants
      implicit none
      integer, parameter :: max_levels_per_ion = 50000
      logical, save :: ii_iii_log_initialized = .false.
      private
      public read_nist_load
      contains

      subroutine read_nist_load(path, atom, line, nlines_out, fout, gfcut, spect_bins)
      character(*), intent(in) :: path
      type(atom_data), intent(inout) :: atom(max_atoms)
      type(line_data), allocatable, intent(inout) :: line(:)
      integer, intent(out) :: nlines_out
      integer, intent(in) :: fout
      real(8), intent(in) :: gfcut
      real(8), intent(in) :: spect_bins(:)
      integer :: z, ion, ZZ, LL, eof, nlines_tot, j, iu
      real(16) :: EN
      character(256) :: ionfile, fgf, flv
      character(2) :: sym
      character(4) :: roman

      ii_iii_log_initialized = .false.

      ionfile = trim(path)//'/ionization_data'
      open(newunit=iu, file=ionfile, status='old', action='read', iostat=eof)
      if (eof.eq.0) then
        do
          read(iu, *, iostat=eof) ZZ, LL, EN
          if (eof.ne.0) exit
          if (ZZ.le.max_atoms .and. LL.le.max_ion_levels) &
            atom(ZZ)%ionization_energy(LL) = EN * electron_volt
        enddo
        close(iu)
      endif

      do z = 1, max_atoms
        if (.not. atom(z)%active) cycle
        do ion = 0, max_ion_levels
          if (allocated(atom(z)%ion(ion)%level)) deallocate(atom(z)%ion(ion)%level)
          atom(z)%ion(ion)%nlevels = 0
          atom(z)%ion(ion)%nlines = 0
        enddo
      enddo

      nlines_tot = 0
      do z = 1, max_atoms
        if (.not. atom(z)%active) cycle
        call z_to_sym(z, sym)
        do ion = 0, max_ion_levels
          call ion_to_roman(ion, roman)
          call load_one_ion(path, z, ion, sym, roman, atom, gfcut, nlines_tot, fout)
        enddo
      enddo

      nlines_out = nlines_tot
      allocate(line(nlines_out))
      j = 0
      do z = 1, max_atoms
        if (.not. atom(z)%active) cycle
        call z_to_sym(z, sym)
        do ion = 0, max_ion_levels
          if (atom(z)%ion(ion)%nlines.le.0) cycle
          call ion_to_roman(ion, roman)
          call fill_lines_one_ion(path, z, ion, sym, roman, atom, line, j, gfcut, spect_bins, fout)
        enddo
      enddo

      ! Combined list over different Z/ion is not globally sorted by lambda.
      ! Keep file-order line ingestion; downstream can bin per frequency directly.

      return
      end subroutine read_nist_load

      subroutine load_one_ion(path, z, ion, sym, roman, atom, gfcut, nlines_tot, fout)
      character(*), intent(in) :: path
      integer, intent(in) :: z, ion
      character(2), intent(in) :: sym
      character(4), intent(in) :: roman
      type(atom_data), intent(inout) :: atom(max_atoms)
      real(8), intent(in) :: gfcut
      integer, intent(inout) :: nlines_tot
      integer, intent(in) :: fout
      character(256) :: fgf, flv
      integer :: eof, iu, iu2, nlev, ilev
      real(16) :: wl, elo_1000, gf, e_1000, jval
      type(level_data) :: lvl
      type(level_data), allocatable :: buf(:)

      flv = trim(path)//'/out_gg/outglv_0_'//trim(sym)//'_'//trim(roman)
      fgf = trim(path)//'/out_gg/outggf_sorted_'//trim(sym)//'_'//trim(roman)

      open(newunit=iu2, file=flv, status='old', action='read', iostat=eof)
      if (eof.ne.0) return

      nlev = 0
      allocate(buf(max_levels_per_ion))
      do
        read(iu2, *, iostat=eof) e_1000, jval
        if (eof.ne.0) exit
        if (nlev.ge.max_levels_per_ion) exit
        lvl%energy = planck * clight * abs(e_1000) * 1000.0d0
        lvl%g = 2.0d0 * jval + 1.0d0
        call add_level_sorted(lvl, buf, nlev, max_levels_per_ion)
      enddo
      close(iu2)


      atom(z)%ion(ion)%nlevels = nlev
      if (nlev.gt.0) then
        allocate(atom(z)%ion(ion)%level(nlev))
        atom(z)%ion(ion)%level(1:nlev) = buf(1:nlev)
      endif
      deallocate(buf)
      ! write all level energies to output file (one level per line)
      if (nlev.gt.0) then
        write(fout, *) 'level energies for Z=', z, ' ion=', ion, ' nlevels=', nlev
        do ilev = 1, nlev
          write(fout,'(I8,1X,ES24.16E3)') ilev, atom(z)%ion(ion)%level(ilev)%energy
        enddo
      endif
      if (nlev.eq.0) return

      open(newunit=iu, file=fgf, status='old', action='read', iostat=eof)
      if (eof.ne.0) return
      do
        read(iu, *, iostat=eof) wl, elo_1000, gf
        if (eof.ne.0) exit
        if (wl.le.0.0d0) cycle
        if (wl.lt.100.0d0 .or. wl.gt.1.d7) cycle
        if (gfcut.gt.-99.0d0 .and. gf.lt.10.0d0**gfcut) cycle
        atom(z)%ion(ion)%nlines = atom(z)%ion(ion)%nlines + 1
        nlines_tot = nlines_tot + 1
      enddo
      close(iu)

      return
      end subroutine load_one_ion

      subroutine fill_lines_one_ion(path, z, ion, sym, roman, atom, line, jstart, gfcut, spect_bins, fout)
      character(*), intent(in) :: path
      integer, intent(in) :: z, ion
      character(2), intent(in) :: sym
      character(4), intent(in) :: roman
      type(atom_data), intent(in) :: atom(max_atoms)
      type(line_data), intent(inout) :: line(:)
      integer, intent(inout) :: jstart
      real(8), intent(in) :: gfcut
      real(8), intent(in) :: spect_bins(:)
      integer, intent(in) :: fout
      character(256) :: fgf
      integer :: eof, iu, i1, i2, nlevels, lambin, nbins, iu_log, row_idx, accepted_idx
      real(16) :: wl, elo_1000, gf, elo_erg, ehi_erg, lambda

      fgf = trim(path)//'/out_gg/outggf_sorted_'//trim(sym)//'_'//trim(roman)
      nlevels = atom(z)%ion(ion)%nlevels
      if (nlevels.le.0) return

      nbins = size(spect_bins)
      lambin = 1
      row_idx = 0
      accepted_idx = 0

      open(newunit=iu, file=fgf, status='old', action='read')
      do
        read(iu, *, iostat=eof) wl, elo_1000, gf
        if (eof.ne.0) exit
        row_idx = row_idx + 1
        if (wl.le.0.0d0) cycle
        if (wl.lt.100.0d0 .or. wl.gt.1.d7) cycle
        if (gfcut.gt.-99.0d0 .and. gf.lt.10.0d0**gfcut) cycle
        accepted_idx = accepted_idx + 1
        elo_erg = planck * clight * abs(elo_1000) * 1000.0d0
        ehi_erg = elo_erg + planck * clight / (wl * 1.d-8)
        i1 = find_closest_level(elo_erg, atom(z)%ion(ion)%level, nlevels)
        i2 = find_closest_level(ehi_erg, atom(z)%ion(ion)%level, nlevels)
        jstart = jstart + 1
        line(jstart)%z = z
        line(jstart)%ion = ion
        line(jstart)%lambda = wl * 1.d-8
        line(jstart)%elo = elo_erg
        line(jstart)%l1 = i1
        line(jstart)%l2 = i2
        line(jstart)%fij = gf / atom(z)%ion(ion)%level(i1)%g
        ! Bin assignment during read (lines per file are sorted by lambda)
        if (nbins.ge.2) then
          lambda = line(jstart)%lambda
          if (lambda.ge.spect_bins(1) .and. lambda.le.spect_bins(nbins)) then
            do while (lambin.lt.nbins-1 .and. lambda.gt.spect_bins(lambin+1))
              lambin = lambin + 1
            enddo
            line(jstart)%ibin = lambin
          else
            line(jstart)%ibin = 0
          endif
        else
          line(jstart)%ibin = 0
        endif

        if (z.eq.52 .and. (ion.eq.1 .or. ion.eq.2)) then
          if (.not.ii_iii_log_initialized) then
            open(newunit=iu_log, file='output/urilight_lines_II_III.log', status='replace', action='write')
            write(iu_log,'(A)') '# z ion roman src_file row_idx accepted_idx wl_A lambda_cm e_low_1000 gf e_low_raw_erg e_up_raw_erg i1 i2 level1_energy level2_energy g1 g2 delta_low delta_high fij ibin'
            close(iu_log)
            ii_iii_log_initialized = .true.
          endif
          open(newunit=iu_log, file='output/urilight_lines_II_III.log', status='old', position='append', action='write')
          write(iu_log,'(I3,1X,I3,1X,A,1X,A,1X,I8,1X,I8,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,I8,1X,I8,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,I8)') &
            z, ion, trim(roman), trim(fgf), row_idx, accepted_idx, wl, line(jstart)%lambda, elo_1000, gf, &
            elo_erg, ehi_erg, i1, i2, atom(z)%ion(ion)%level(i1)%energy, atom(z)%ion(ion)%level(i2)%energy, &
            atom(z)%ion(ion)%level(i1)%g, atom(z)%ion(ion)%level(i2)%g, &
            abs(atom(z)%ion(ion)%level(i1)%energy-elo_erg), abs(atom(z)%ion(ion)%level(i2)%energy-ehi_erg), &
            line(jstart)%fij, line(jstart)%ibin
          close(iu_log)
        endif
      enddo
      close(iu)

      write(fout, *) 'number of lines in ion ', ion, ' of atom ', z, ' = ', jstart
      write(fout, *) 'last line = ', line(jstart)%lambda
      write(fout, *) 'last line = ', line(jstart)%l1
      write(fout, *) 'last line = ', line(jstart)%l2
      write(fout, *) 'last line = ', line(jstart)%fij
      write(fout, *) 'last line = ', line(jstart)%ibin
      write(fout, *) 'last line = ', line(jstart)%z
      write(fout, *) 'last line = ', line(jstart)%ion

      return
      end subroutine fill_lines_one_ion

      subroutine z_to_sym(z, sym)
      integer, intent(in) :: z
      character(2), intent(out) :: sym
      character(2) :: s(70)
      data s / 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb' /
      if (z.ge.1 .and. z.le.70) then
        sym = s(z)
      else
        sym = '??'
      endif
      end subroutine z_to_sym

      subroutine ion_to_roman(ion, roman)
      integer, intent(in) :: ion
      character(4), intent(out) :: roman
      character(4) :: r(0:7)
      data r / 'I   ','II  ','III ','IV  ','V   ','VI  ','VII ','VIII' /
      if (ion.ge.0 .and. ion.le.7) then
        roman = r(ion)
      else
        roman = ' '
      endif
      end subroutine ion_to_roman

      subroutine add_level_sorted(lvl, level, length, maxlen)
      type(level_data), intent(in) :: lvl
      type(level_data), intent(inout) :: level(:)
      integer, intent(inout) :: length
      integer, intent(in) :: maxlen
      integer :: i
      if (length.ge.maxlen) return
      if (length.eq.0) then
        level(1) = lvl
        length = 1
        return
      endif
      i = 1
      do while (i.le.length .and. level(i)%energy .lt. lvl%energy)
        i = i + 1
      enddo
      if (i.le.length .and. level(i)%energy .eq. lvl%energy) then
        do while (i.le.length .and. level(i)%energy .eq. lvl%energy .and. level(i)%g .lt. lvl%g)
          i = i + 1
        enddo
        if (i.le.length .and. level(i)%energy .eq. lvl%energy .and. level(i)%g .eq. lvl%g) return
      endif
      if (length.ge.maxlen) return
      level(i+1:length+1) = level(i:length)
      level(i) = lvl
      length = length + 1
      end subroutine add_level_sorted

      integer function find_closest_level(e_erg, level, length)
      real(16), intent(in) :: e_erg
      type(level_data), intent(in) :: level(:)
      integer, intent(in) :: length
      integer :: lo, hi, mid
      if (length.le.0) then
        find_closest_level = 1
        return
      endif
      if (e_erg.le.level(1)%energy) then
        find_closest_level = 1
        return
      endif
      if (e_erg.ge.level(length)%energy) then
        find_closest_level = length
        return
      endif

      lo = 1
      hi = length
      do while (hi - lo .gt. 1)
        mid = (lo + hi) / 2
        if (level(mid)%energy.le.e_erg) then
          lo = mid
        else
          hi = mid
        endif
      enddo

      if (abs(level(lo)%energy - e_erg).le.abs(level(hi)%energy - e_erg)) then
        find_closest_level = lo
      else
        find_closest_level = hi
      endif
      end function find_closest_level

      end module read_nist_db_mod
