subroutine colbac(u,s,d,jj)
implicit none
double precision u(*),s(*),d      
integer  jj

!  Purpose: Projection to standard eigenproblem Backsubstitution macro

!  Inputs:
!     u(*)     - Upper triangular factor of matrix
!     d        - Diagonal to divide into solution
!     jj       - Length of column

!  Outputs:
!     s(*)     - Back substituion solution

   double precision dd
   integer  j

   dd = s(jj+1)
   do j = 1,jj
     s(j) = s(j) - dd*u(j)
   end do  
   s(jj) = s(jj)/d

end
