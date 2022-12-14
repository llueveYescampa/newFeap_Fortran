subroutine pend(subnam)
implicit none
  character ::subnam*6

!  Purpose: End of file found

!  Inputs:
!     subnam    - Name of subroutine where EOF encountered


   include 'iofile.h'

   write(iow,'(a,a6,a)') &
     ' ** ERROR in ', subnam,' ** end of file encountered'
   if(ior.lt.0) then
     write(  *,'(a,a6,a)') &
       ' ** ERROR in ', subnam,' ** end of file encountered'
   end if  
   call pstop(-21) ! stop
end
