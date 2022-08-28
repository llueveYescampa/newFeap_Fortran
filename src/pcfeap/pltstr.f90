subroutine pltstr(dt,st,numnp)
implicit none
integer numnp
double precision   dt(numnp),st(numnp,*)

!  Stress projections computed by dividing by 'lumped' weightings

   integer  ii,j
   double precision   sig(6)

   double precision eerror,elproj,ecproj,efem,enerr,ebar
   common /errind/  eerror,elproj,ecproj,efem,enerr,ebar

   elproj = 0.0
   do ii = 1,numnp
     if(dt(ii).ne.0.0) then
       do j = 1,4
         sig(j)   = st(ii,j)/dt(ii)
         elproj   = elproj + sig(j)*st(ii,j)
         st(ii,j) = sig(j)
       end do  

!      Compute the principal stress values

       call pstres(sig,sig(4),sig(5),sig(6))
       if(st(ii,5).ne.0.0) then
         dt(ii) = st(ii,5)/dt(ii)
       else
         dt(ii) = sig(6)
       endif
       st(ii,5) = sig(4)
       st(ii,6) = sig(5)
       st(ii,7) = (sig(4)-sig(5)) * 0.5
     end if
   end do  
end
