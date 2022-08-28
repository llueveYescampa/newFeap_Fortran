subroutine datest(au,jh,daval)
implicit none
integer  jh
double precision   au(jh),daval

!  Purpose: Test for rank of matrix with zero diagonals
   
!  Inputs:
!     au(*)      - Column for j-th equation
!     jh         - Length of column
   
!  Outputs:
!     daval      - Sum of absolute values of terms, if zero
!                  entire equation assumed to be zero and skipped.
   
   
   integer  j
   
   daval = 0.0d0
   do j = 1,jh
     daval = daval + abs(au(j))
   end do   

end
