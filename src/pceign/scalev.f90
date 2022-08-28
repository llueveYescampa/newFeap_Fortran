subroutine scalev(v,nn)
implicit none
integer nn

!  Purpose: Scale vector to have maximum element of 1.0

!  Inputs:
!     v(*)     - Vector to scale
!     nn       - Length of vector

!  Outputs:
!     v(*)     - Scaled vector

   integer n
   double precision  v(*),vmax

   vmax = abs(v(1))
   do n = 1,nn
     vmax = max(vmax,abs(v(n)))
   end do  

   do n = 1,nn
     v(n) = v(n)/vmax
   end do
      
end
