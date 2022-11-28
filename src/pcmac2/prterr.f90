subroutine prterr()
implicit none

!  Purpose:  Output error estimation values.

   double precision    elrind,ecrind

   include 'iofile.h'
   include 'errind.h'

!  Output error indicator values

   elrind = 0.0d0
   ecrind = 0.0d0
   if(elproj.ne.0.0d0) then
     elrind = sqrt((efem-elproj)/elproj)
   end if  
   if(ecproj.ne.0.0d0) then
     ecrind = sqrt((efem-ecproj)/ecproj)
     eerror = sqrt(eerror/ecproj)
   end if  
     
   write(iow,2000) efem,elproj,ecproj,elrind,ecrind,eerror
   if(ior.lt.0) then
     write(*,2000) efem,elproj,ecproj,elrind,ecrind,eerror
   end if

2000  format(/'   Finite Element Stress Measure              =',e15.8/ &
              '   L u m p e d  Projected Stress Measure      =',e15.8/ &
              '   Consistent   Projected Stress Measure      =',e15.8/ &
              '   L u m p e d  Error Indicator               =',e15.8/ &
              '   Consistent   Error Indicator               =',e15.8/ &
              '   Direct Error Indicator                     =',e15.8/)

end
