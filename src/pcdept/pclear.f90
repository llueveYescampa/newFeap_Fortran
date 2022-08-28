! Esta rutina solo la usa la version DOS. No la version PCTEKT

subroutine pclear()  !     Purpose: Clear PC screen and home cursor on monitor
implicit none

   character*1 cha
   
   cha = '[2J'
   write(*,'(1x,2a1)') char(27), cha
      
end
