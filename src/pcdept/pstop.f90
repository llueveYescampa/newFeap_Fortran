subroutine pstop(errStatus)
implicit none
  integer :: errStatus

  if (errStatus == 0) then
    print *, 'Done!!'
  else
    print *, 'Upss... Something went wrong. Error code: ', errStatus
  endif
  call exit(errStatus)
end

