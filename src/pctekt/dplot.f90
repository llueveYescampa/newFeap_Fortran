subroutine dplot(x,y,ipen)
   double precision x,y
   integer ipen

   character*4  str4,stro
   integer jx1, jy1
   
   common /tekt1/ stro

!  tektronix 4012 device or emulator
!  pen command motions: ipen = 1, move to position x,y
!                       ipen = 2, drawline from to x,y

   jx1 = x*700
   jy1 = y*770
   call brk4(jx1,jy1,str4)
   if(ipen .eq. 2) then
     write(*,'(2a1,2a4)') char(13),char(29),stro,str4
   end if
   stro = str4
end

