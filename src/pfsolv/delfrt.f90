subroutine delfrt()
implicit none


  character fname*12
  logical   lfil
  integer   n

! Delete all files from frontal solver

  do n = 1,9999

!   Set name

    if(n.lt.10) then
      write(fname,'(a10,i1)') 'Frontal.00',n
    else if(n.lt.100) then
      write(fname,'( a9,i2)') 'Frontal.0',n
    else
      write(fname,'( a8,i3)') 'Frontal.',n
    end if

!   Check and delete

    inquire(file=fname,exist=lfil)
    if(lfil) then
      open (4,file=fname,status='old')
      close(4,status='delete')
    else
      return
    end if
  end do

end
