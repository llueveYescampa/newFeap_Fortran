subroutine addvec(a,b,nn)
implicit none
integer  nn
double precision   a(*),b(*)

!  Purpose: Add two vectors

!  Inputs:
!     a(*)     - Vector 1
!     b(*)     - Vector 2
!     nn       - Length of vectors

!  Outputs:
!     a(*)     - Sum of vector 1 and 2

   integer  n

   do n = 1,nn
     a(n) = a(n) + b(n)
   end do  

end
