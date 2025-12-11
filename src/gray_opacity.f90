! Gray Opacity Module for UriLight
! Handles reading and interpolation of pre-computed Rosseland mean opacity tables

module GrayOpacity
  use globals
  use physical_constants
  implicit none

  integer, parameter :: max_gray_materials = 10
  integer, parameter :: max_table_size = 10000

  type gray_opacity_table_type
    integer :: nT, nrho, nt_time
    real(8), allocatable :: T_vals(:)
    real(8), allocatable :: rho_vals(:)
    real(8), allocatable :: time_vals(:)
    real(8), allocatable :: kappa_P_abs(:,:,:)
    real(8), allocatable :: kappa_R_tot(:,:,:)
    real(8), allocatable :: f_abs(:,:,:)
    logical :: loaded = .false.
  end type gray_opacity_table_type

  type(gray_opacity_table_type), save :: gray_tables(max_gray_materials)
  integer, save :: n_gray_materials_loaded = 0

  contains

  subroutine load_gray_opacity_table(filename, material_id)
    implicit none
    character(*), intent(in) :: filename
    integer, intent(in) :: material_id
    integer :: iunit, line_idx, j, n, iostat, nlines
    integer :: i_temp, i_rho, i_time
    character(200) :: line
    real(8) :: T_val, rho_val, t_days_val, kR_abs_dummy, kR_scat_dummy, kP_abs, kR_tot, f_abs_val
    real(8), allocatable :: T_grid(:), rho_grid(:), time_grid(:)
    integer :: nT, nrho, nt_time
    logical :: T_found, rho_found, t_time_found

    if (material_id < 1 .or. material_id > max_gray_materials) then
      write(fout,*) 'ERROR: Material ID out of range:', material_id
      stop
    endif

    write(fout,*) 'Loading gray opacity table for material', material_id, ': ', trim(filename)

    iunit = 10
    open(unit=iunit, file=filename, status='old', iostat=iostat)
    if (iostat /= 0) then
      write(fout,*) 'ERROR: Cannot open gray opacity table: ', trim(filename)
      stop
    endif

    ! Count lines and find unique grid values
    nlines = 0
    do
      read(iunit, '(A)', iostat=iostat) line
      if (iostat /= 0) exit
      if (line(1:1) == '#') cycle
      nlines = nlines + 1
    enddo
    rewind(iunit)

    ! Skip header lines
    do
      read(iunit, '(A)') line
      if (line(1:1) /= '#') then
        backspace(iunit)
        exit
      endif
    enddo

    ! Read all data to determine grid
    allocate(T_grid(max_table_size), rho_grid(max_table_size), time_grid(max_table_size))
    nT = 0
    nrho = 0
    nt_time = 0

    do line_idx = 1, nlines
      ! Read all columns but only use: T, rho, t, kP_abs, kR_tot, f_abs
      ! Skip kR_abs and kR_scat (columns 4 and 5) - will derive from kR_tot and f_abs
      read(iunit, *, iostat=iostat) T_val, rho_val, t_days_val, kR_abs_dummy, kR_scat_dummy, kP_abs, kR_tot, f_abs_val
      if (iostat /= 0) exit

      ! Check if T is new
      T_found = .false.
      do j = 1, nT
        if (abs(T_val - T_grid(j)) < 1e-10) then
          T_found = .true.
          exit
        endif
      enddo
      if (.not. T_found) then
        nT = nT + 1
        T_grid(nT) = T_val
      endif

      ! Check if rho is new
      rho_found = .false.
      do j = 1, nrho
        if (abs(rho_val - rho_grid(j)) < 1e-20) then
          rho_found = .true.
          exit
        endif
      enddo
      if (.not. rho_found) then
        nrho = nrho + 1
        rho_grid(nrho) = rho_val
      endif

      ! Check if t is new
      t_time_found = .false.
      do j = 1, nt_time
        if (abs(t_days_val - time_grid(j)) < 1e-10) then
          t_time_found = .true.
          exit
        endif
      enddo
      if (.not. t_time_found) then
        nt_time = nt_time + 1
        time_grid(nt_time) = t_days_val
      endif
    enddo

    ! Sort grid values
    call sort_array(T_grid(1:nT))
    call sort_array(rho_grid(1:nrho))
    call sort_array(time_grid(1:nt_time))

    gray_tables(material_id)%nT = nT
    gray_tables(material_id)%nrho = nrho
    gray_tables(material_id)%nt_time = nt_time

    allocate(gray_tables(material_id)%T_vals(nT))
    allocate(gray_tables(material_id)%rho_vals(nrho))
    allocate(gray_tables(material_id)%time_vals(nt_time))
    allocate(gray_tables(material_id)%kappa_P_abs(nT, nrho, nt_time))
    allocate(gray_tables(material_id)%kappa_R_tot(nT, nrho, nt_time))
    allocate(gray_tables(material_id)%f_abs(nT, nrho, nt_time))

    gray_tables(material_id)%T_vals = T_grid(1:nT)
    gray_tables(material_id)%rho_vals = rho_grid(1:nrho)
    gray_tables(material_id)%time_vals = time_grid(1:nt_time)

    ! Read data into arrays
    rewind(iunit)
    ! Skip header
    do
      read(iunit, '(A)') line
      if (line(1:1) /= '#') then
        backspace(iunit)
        exit
      endif
    enddo

    do n = 1, nlines
      ! Read all columns but only use: T, rho, t, kP_abs, kR_tot, f_abs
      ! Skip kR_abs and kR_scat (columns 4 and 5) - will derive from kR_tot and f_abs
      read(iunit, *, iostat=iostat) T_val, rho_val, t_days_val, kR_abs_dummy, kR_scat_dummy, kP_abs, kR_tot, f_abs_val
      if (iostat /= 0) exit

      ! Find grid indices - use exact match
      i_temp = find_index_in_array(T_val, gray_tables(material_id)%T_vals, gray_tables(material_id)%nT)
      i_rho = find_index_in_array(rho_val, gray_tables(material_id)%rho_vals, gray_tables(material_id)%nrho)
      i_time = find_index_in_array(t_days_val, gray_tables(material_id)%time_vals, gray_tables(material_id)%nt_time)

      if (i_temp > 0 .and. i_rho > 0 .and. i_time > 0) then
        gray_tables(material_id)%kappa_P_abs(i_temp, i_rho, i_time) = kP_abs
        gray_tables(material_id)%kappa_R_tot(i_temp, i_rho, i_time) = kR_tot
        gray_tables(material_id)%f_abs(i_temp, i_rho, i_time) = f_abs_val
      else
        write(fout,*) 'WARNING: Could not find grid indices for T=', T_val, ' rho=', rho_val, ' t=', t_days_val
        write(fout,*) '  i_temp=', i_temp, ' i_rho=', i_rho, ' i_time=', i_time
      endif
    enddo

    close(iunit)
    gray_tables(material_id)%loaded = .true.
    n_gray_materials_loaded = max(n_gray_materials_loaded, material_id)

    write(fout,*) 'Loaded gray opacity table: nT=', nT, ' nrho=', nrho, ' nt=', nt_time
    write(fout,*) 'T range: ', gray_tables(material_id)%T_vals(1), ' to ', &
                  gray_tables(material_id)%T_vals(nT), ' eV'
    write(fout,*) 'rho range: ', gray_tables(material_id)%rho_vals(1), ' to ', &
                  gray_tables(material_id)%rho_vals(nrho), ' g/cm³'
    write(fout,*) 't range: ', gray_tables(material_id)%time_vals(1), ' to ', &
                  gray_tables(material_id)%time_vals(nt_time), ' days'

    deallocate(T_grid, rho_grid, time_grid)
  end subroutine load_gray_opacity_table

  function get_gray_opacity_rosseland_abs(material_id, T_eV, rho, t_days) result(kappa)
    integer, intent(in) :: material_id
    real(8), intent(in) :: T_eV, rho, t_days
    real(8) :: kappa
    kappa = interpolate_table(material_id, T_eV, rho, t_days, 1)
  end function get_gray_opacity_rosseland_abs

  function get_gray_opacity_rosseland_scat(material_id, T_eV, rho, t_days) result(kappa)
    integer, intent(in) :: material_id
    real(8), intent(in) :: T_eV, rho, t_days
    real(8) :: kappa
    kappa = interpolate_table(material_id, T_eV, rho, t_days, 2)
  end function get_gray_opacity_rosseland_scat

  function get_gray_opacity_planck_abs(material_id, T_eV, rho, t_days) result(kappa)
    integer, intent(in) :: material_id
    real(8), intent(in) :: T_eV, rho, t_days
    real(8) :: kappa
    kappa = interpolate_table(material_id, T_eV, rho, t_days, 3)
  end function get_gray_opacity_planck_abs

  function get_gray_opacity_rosseland_tot(material_id, T_eV, rho, t_days) result(kappa)
    integer, intent(in) :: material_id
    real(8), intent(in) :: T_eV, rho, t_days
    real(8) :: kappa
    kappa = interpolate_table(material_id, T_eV, rho, t_days, 4)
  end function get_gray_opacity_rosseland_tot

  function get_gray_opacity_f_abs(material_id, T_eV, rho, t_days) result(f_abs)
    integer, intent(in) :: material_id
    real(8), intent(in) :: T_eV, rho, t_days
    real(8) :: f_abs
    f_abs = interpolate_table(material_id, T_eV, rho, t_days, 5)
  end function get_gray_opacity_f_abs

  function interpolate_table(material_id, T_eV, rho, t_days, opac_type) result(kappa)
    integer, intent(in) :: material_id, opac_type
    real(8), intent(in) :: T_eV, rho, t_days
    real(8) :: kappa
    integer :: iT, irho, it_idx
    real(8) :: wT, wrho, wt_val
    real(8) :: kappa_000, kappa_001, kappa_010, kappa_011
    real(8) :: kappa_100, kappa_101, kappa_110, kappa_111

    if (.not. gray_tables(material_id)%loaded) then
      write(fout,*) 'ERROR: Gray opacity table not loaded for material', material_id
      kappa = 0.0d0
      return
    endif

    ! Find grid indices (log space interpolation)
    iT = find_index_log(T_eV, gray_tables(material_id)%T_vals, gray_tables(material_id)%nT, wT)
    irho = find_index_log(rho, gray_tables(material_id)%rho_vals, gray_tables(material_id)%nrho, wrho)
    it_idx = find_index_log(t_days, gray_tables(material_id)%time_vals, gray_tables(material_id)%nt_time, wt_val)

    ! Trilinear interpolation
    ! For cases 1 and 2, derive from kappa_R_tot and f_abs
    select case(opac_type)
    case(1)  ! kappa_R_abs = f_abs * kappa_R_tot
      ! Interpolate f_abs and kappa_R_tot separately, then multiply
      kappa_000 = gray_tables(material_id)%f_abs(iT, irho, it_idx) * gray_tables(material_id)%kappa_R_tot(iT, irho, it_idx)
      kappa_001 = gray_tables(material_id)%f_abs(iT, irho, it_idx+1) * gray_tables(material_id)%kappa_R_tot(iT, irho, it_idx+1)
      kappa_010 = gray_tables(material_id)%f_abs(iT, irho+1, it_idx) * gray_tables(material_id)%kappa_R_tot(iT, irho+1, it_idx)
      kappa_011 = gray_tables(material_id)%f_abs(iT, irho+1, it_idx+1) * gray_tables(material_id)%kappa_R_tot(iT, irho+1, it_idx+1)
      kappa_100 = gray_tables(material_id)%f_abs(iT+1, irho, it_idx) * gray_tables(material_id)%kappa_R_tot(iT+1, irho, it_idx)
      kappa_101 = gray_tables(material_id)%f_abs(iT+1, irho, it_idx+1) * gray_tables(material_id)%kappa_R_tot(iT+1, irho, it_idx+1)
      kappa_110 = gray_tables(material_id)%f_abs(iT+1, irho+1, it_idx) * gray_tables(material_id)%kappa_R_tot(iT+1, irho+1, it_idx)
      kappa_111 = gray_tables(material_id)%f_abs(iT+1, irho+1, it_idx+1) * gray_tables(material_id)%kappa_R_tot(iT+1, irho+1, it_idx+1)
    case(2)  ! kappa_R_scat = (1 - f_abs) * kappa_R_tot
      kappa_000 = (1.0d0 - gray_tables(material_id)%f_abs(iT, irho, it_idx)) * gray_tables(material_id)%kappa_R_tot(iT, irho, it_idx)
      kappa_001 = (1.0d0 - gray_tables(material_id)%f_abs(iT, irho, it_idx+1)) * gray_tables(material_id)%kappa_R_tot(iT, irho, it_idx+1)
      kappa_010 = (1.0d0 - gray_tables(material_id)%f_abs(iT, irho+1, it_idx)) * gray_tables(material_id)%kappa_R_tot(iT, irho+1, it_idx)
      kappa_011 = (1.0d0 - gray_tables(material_id)%f_abs(iT, irho+1, it_idx+1)) * gray_tables(material_id)%kappa_R_tot(iT, irho+1, it_idx+1)
      kappa_100 = (1.0d0 - gray_tables(material_id)%f_abs(iT+1, irho, it_idx)) * gray_tables(material_id)%kappa_R_tot(iT+1, irho, it_idx)
      kappa_101 = (1.0d0 - gray_tables(material_id)%f_abs(iT+1, irho, it_idx+1)) * gray_tables(material_id)%kappa_R_tot(iT+1, irho, it_idx+1)
      kappa_110 = (1.0d0 - gray_tables(material_id)%f_abs(iT+1, irho+1, it_idx)) * gray_tables(material_id)%kappa_R_tot(iT+1, irho+1, it_idx)
      kappa_111 = (1.0d0 - gray_tables(material_id)%f_abs(iT+1, irho+1, it_idx+1)) * gray_tables(material_id)%kappa_R_tot(iT+1, irho+1, it_idx+1)
    case(3)  ! kappa_P_abs
      kappa_000 = gray_tables(material_id)%kappa_P_abs(iT, irho, it_idx)
      kappa_001 = gray_tables(material_id)%kappa_P_abs(iT, irho, it_idx+1)
      kappa_010 = gray_tables(material_id)%kappa_P_abs(iT, irho+1, it_idx)
      kappa_011 = gray_tables(material_id)%kappa_P_abs(iT, irho+1, it_idx+1)
      kappa_100 = gray_tables(material_id)%kappa_P_abs(iT+1, irho, it_idx)
      kappa_101 = gray_tables(material_id)%kappa_P_abs(iT+1, irho, it_idx+1)
      kappa_110 = gray_tables(material_id)%kappa_P_abs(iT+1, irho+1, it_idx)
      kappa_111 = gray_tables(material_id)%kappa_P_abs(iT+1, irho+1, it_idx+1)
    case(4)  ! kappa_R_tot
      kappa_000 = gray_tables(material_id)%kappa_R_tot(iT, irho, it_idx)
      kappa_001 = gray_tables(material_id)%kappa_R_tot(iT, irho, it_idx+1)
      kappa_010 = gray_tables(material_id)%kappa_R_tot(iT, irho+1, it_idx)
      kappa_011 = gray_tables(material_id)%kappa_R_tot(iT, irho+1, it_idx+1)
      kappa_100 = gray_tables(material_id)%kappa_R_tot(iT+1, irho, it_idx)
      kappa_101 = gray_tables(material_id)%kappa_R_tot(iT+1, irho, it_idx+1)
      kappa_110 = gray_tables(material_id)%kappa_R_tot(iT+1, irho+1, it_idx)
      kappa_111 = gray_tables(material_id)%kappa_R_tot(iT+1, irho+1, it_idx+1)
    case(5)  ! f_abs
      kappa_000 = gray_tables(material_id)%f_abs(iT, irho, it_idx)
      kappa_001 = gray_tables(material_id)%f_abs(iT, irho, it_idx+1)
      kappa_010 = gray_tables(material_id)%f_abs(iT, irho+1, it_idx)
      kappa_011 = gray_tables(material_id)%f_abs(iT, irho+1, it_idx+1)
      kappa_100 = gray_tables(material_id)%f_abs(iT+1, irho, it_idx)
      kappa_101 = gray_tables(material_id)%f_abs(iT+1, irho, it_idx+1)
      kappa_110 = gray_tables(material_id)%f_abs(iT+1, irho+1, it_idx)
      kappa_111 = gray_tables(material_id)%f_abs(iT+1, irho+1, it_idx+1)
    end select

    ! Trilinear interpolation
    kappa = (1.0d0-wT) * (1.0d0-wrho) * (1.0d0-wt_val) * kappa_000 + &
            (1.0d0-wT) * (1.0d0-wrho) * wt_val * kappa_001 + &
            (1.0d0-wT) * wrho * (1.0d0-wt_val) * kappa_010 + &
            (1.0d0-wT) * wrho * wt_val * kappa_011 + &
            wT * (1.0d0-wrho) * (1.0d0-wt_val) * kappa_100 + &
            wT * (1.0d0-wrho) * wt_val * kappa_101 + &
            wT * wrho * (1.0d0-wt_val) * kappa_110 + &
            wT * wrho * wt_val * kappa_111
  end function interpolate_table

  function find_index_log(x, arr, n, weight) result(idx)
    integer, intent(in) :: n
    real(8), intent(in) :: x, arr(n)
    real(8), intent(out) :: weight
    integer :: idx
    integer :: i
    real(8) :: logx, logarr1, logarr2

    if (x <= arr(1)) then
      idx = 1
      weight = 0.0d0
      return
    endif
    if (x >= arr(n)) then
      idx = n - 1
      weight = 1.0d0
      return
    endif

    logx = log(x)
    do i = 1, n-1
      if (x >= arr(i) .and. x <= arr(i+1)) then
        logarr1 = log(arr(i))
        logarr2 = log(arr(i+1))
        if (abs(logarr2 - logarr1) < 1e-10) then
          weight = 0.5d0
        else
          weight = (logx - logarr1) / (logarr2 - logarr1)
        endif
        idx = i
        return
      endif
    enddo
    idx = 1
    weight = 0.0d0
  end function find_index_log

  function find_index_in_array(x, arr, n) result(idx)
    integer, intent(in) :: n
    real(8), intent(in) :: x, arr(n)
    integer :: idx, i
    real(8) :: tol
    ! Use relative tolerance for better matching
    tol = max(abs(x) * 1e-8, 1e-20)
    do i = 1, n
      if (abs(x - arr(i)) < tol) then
        idx = i
        return
      endif
    enddo
    idx = 0
  end function find_index_in_array

  subroutine sort_array(arr)
    real(8), intent(inout) :: arr(:)
    integer :: i, j, n
    real(8) :: temp
    n = size(arr)
    do i = 1, n-1
      do j = i+1, n
        if (arr(i) > arr(j)) then
          temp = arr(i)
          arr(i) = arr(j)
          arr(j) = temp
        endif
      enddo
    enddo
  end subroutine sort_array

end module GrayOpacity
