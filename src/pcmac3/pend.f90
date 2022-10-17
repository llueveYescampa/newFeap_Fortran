subroutine pend(subnam)
implicit  none
character subnam*6

!  Purpose: End of file found

!  Inputs:
!     subnam    - Name of subroutine where EOF encountered


   include 'iofile.h'

   if(ioRead.gt.0) then
     write(ioWrite,'(a,a6,a)') &
       ' ** ERROR in ', subnam,' ** end of file encountered'
   end if  
   if(ioRead.lt.0) then
     write(  *,'(a,a6,a)') &
       ' ** ERROR in ', subnam,' ** end of file encountered'
   end if  
   stop

end
