subroutine inaccl(id,dr,xm,a,nneq)
implicit  none
integer    nneq,id(nneq)
double precision    dr(*),xm(*),a(nneq)

!  Purpose:  Set initial acceleration
!     id(*)     - Equation numbers for each dof
!     dr(*)     - Residual for initial state
!     xm(*)     - Diagonal mass matrix
!     nneq      - Number of terms in id

!  Outputs:
!     a(*)      - Initial accelerations

   integer   j, n

   do n = 1,nneq
     j = id(n)
     if(j .gt. 0) then
       a(n) = dr(j)/xm(j)
     end if  
   end do  

end
