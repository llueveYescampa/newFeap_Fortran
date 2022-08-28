double precision function dot (a,b,n)
implicit none
integer  n
double precision   a(*),b(*)

!  Purpose: Compute dot product of two vectors

!  Inputs:
!     a(*)     - Vector 1
!     b(*)     - Vector 2
!     n        - Length of vectors

!  Outputs:
!     dot      - Dot product value

   integer k

!  dot product function 

   dot = 0.0d0
   do k=1,n
     dot = dot + a(k)*b(k)
   end do  

end
