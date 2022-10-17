subroutine pintio(y,n0)
implicit  none
      character y*80
      integer   n0

!   Purpose:  Input control from current active unit - into a character
!             array for free format processing by 'acheck' into field
!             widths of 'n0'

!   Inputs:
!      n0       - Field width for formatting

!   Outputs:
!      y        - Character array formatted to field widths of 'n0'
!                 N.B. Internal I/O can be made to extract data

    character x*80
    include 'iofile.h'
    if(ioRead.gt.0) then 
       read(ioRead,'(a)',err=100,end=100) x
    else if(ioRead.lt.0) then 
       read(*,'(a)',err=100,end=100) x
    end if 
    
    call acheck(x,y,n0,80)
    return

100 call pperror('pintio',x)

end

