module timing 

real :: timetab(20), timetab_start(20), timetab_stop(20)
public :: timetab, timetab_start,  timetab_stop

contains 

subroutine start_time(item)
  integer item 
!  real, save ::  timetab_start(20)
!  real, save ::  timetab_stop(20)
!  real ::  timetab_start(20)
!  real ::  timetab_stop(20)
  real       :: t0
  call cpu_time(t0)
  if (item.le.20) timetab_start(item)= t0
end subroutine 

subroutine stop_time(item)
  integer item 
!  real, save ::  timetab_start(20)
!  real, save ::  timetab_stop(20)
!  real ::  timetab_start(20)
!  real ::  timetab_stop(20)
  real       :: t0
  call cpu_time(t0)
  if (item.le.20)  then 
      timetab_stop(item)= t0
      timetab(item) =  timetab(item) + timetab_stop(item)-timetab_start(item)
  endif
end subroutine 
end module timing 
