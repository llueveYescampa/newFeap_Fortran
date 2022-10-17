subroutine rsolve(b,dr,a,ipd,ipr,maxf,nv,neq,nev,engy,ifl)
implicit none
integer  ipd,ipr,maxf,nv,neq,nev,ifl
double precision   b(*),dr(neq,*),a(*),engy

!  Purpose: Resolution for profile solution

!  Inputs: 
!     b(*)      - Not used for profile solution
!     dr(neq,*) - Right hand side of equations
!     a(*)      - Profile storage for matrix
!     ipd       - Precision of !double precision! variables
!     ipr       - Precision of !real*4! variables
!     maxf      - Not used for profile solution
!     nv        - Not used for profile solution
!     neq       - Number of equations in matrix
!     nev       - Number of right hand side vectors to solve
!     ifl       - Not used for profile solution

!  Outputs:
!     dr(neq,*) - Solutions of equations
!     engy      - Energy of solution

   integer  ne,n12

   include 'maxa.h'      
   
   include 'ddata.h'

   n12 = neq*ipd - ipd - ipr + 1
   do ne = 1,nev
     call dasol(a(neq+1),a(neq+1),a,dr(1,ne),im(n12),neq, engy)
   end do  

end
