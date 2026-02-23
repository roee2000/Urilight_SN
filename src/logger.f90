! Lightweight timestamped logging to the main "out" file.
module Logger
  use globals, only: fout
  implicit none

contains

  subroutine log_info(message)
    character(*), intent(in) :: message
    character(8) :: date
    character(10) :: time

    call date_and_time(date, time)
    write(fout,'("[",A2,"-",A2,"-",A2," ",A2,":",A2,":",A2,"] ",A)') &
         date(3:4), date(5:6), date(7:8), &
         time(1:2), time(3:4), time(5:6), trim(message)
  end subroutine log_info

end module Logger
