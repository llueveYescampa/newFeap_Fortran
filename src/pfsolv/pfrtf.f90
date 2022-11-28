subroutine pfrtf(b,dr,nfrt,k,jj,eq)
implicit none
integer               nfrt,k,jj
double precision b(*),dr(*),eq(*)

!  Purpose: Forward elimination macro for front program

!  Inputs:
!     b(*)     - Front residual
!     nfrt     - Size of current front
!     k        - Location of diagonal
!     jj       - Equation number in global equations
!     eq       - Current equation

!  Outputs
!     b(*)     - Reduced front equation


   integer  ii,km,kp
   double precision r

   dr(jj) = dr(jj) + b(k)
   r = dr(jj)/eq(k)
   km = k - 1
   kp = k + 1

   if(km.gt.0) then
     do ii = 1,km
       b(ii) = b(ii) - eq(ii)*r
     end do  
   end if

   if(kp.le.nfrt) then
     do ii = kp,nfrt
       b(ii-1) = b(ii) - eq(ii)*r
     end do  
   end if

   b(nfrt) = 0.0

end
