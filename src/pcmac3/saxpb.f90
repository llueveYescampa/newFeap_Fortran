subroutine saxpb (a,b,x,n,c)
implicit none
  integer          :: n
  double precision ::  a(*),b(*),x,c(*)

!  Purpose:  Scalar times vector plus vector

!  Inputs:
!     a(*)      - Vector to be multiplied by scalar
!     b(*)      - Vector to add
!     x         - Scalar value multiplying 'a'
!     n         - Length of vectors

!  Outputs:
!     c(*)      - Result of operation

!  integer  k

! vector times scalar added to second vector 

!   do k=1,n
!     c(k) = a(k)*x + b(k)
!   end do  

   c(1:n) = x * a(1:n) + b(1:n)

end
