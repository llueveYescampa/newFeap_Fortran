subroutine pctime(etime)
implicit none

! Purpose: Calls time routine for Micro Soft compilers
!          N.B. Replace !gettim! by appropriate timing routine 
!               for each compiler

! Outputs:
!    etime     - Character array with timing information

   character etime*10
   character mydate*8
   character mytime*10
   character myzone*5
   integer   myvalues(8)
   
   call date_and_time( mydate, mytime, myzone, myvalues)    
   write(etime,'(i3,a1,i2,a1,i2)') myvalues(5),':',myvalues(6),':',myvalues(7)

end


