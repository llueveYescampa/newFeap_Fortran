subroutine rsolve(b,dr,m,ipd,ipr,maxf,nv,neq,nev,aengy,ifl)
implicit none
  integer          :: m(*),ipd,ipr,maxf,nv,neq,nev,ifl
  double precision :: b(maxf,*),dr(neq,*),aengy

!  Purpose: Resolution for frontal solution

!  Inputs:
!     b(*)      - Frontal RHS/solution
!     dr(neq,*) - Residual
!     m(*)      - Buffer for fronal equations
!     ipd       - Precision of 'double precision' variables
!     ipr       - Precision of 'real*4' variables
!     maxf      - Maximum front estimate
!     nv        - Number of frontal buffer blocks
!     neq       - Number of equations to solve
!     nev       - Number of RHS
!     ifl       - File name number for frontal equation storage on disk
     
!  Outputs:
!     dr(neq,*) - Residual
!     aengy     - Energy of solution


   include 'temfl2.h'
   include 'temfl3.h'

   call pfrtfw(b,dr,m,ipd,ibuf,maxf,nv,neq,nev,ifl)
   call pfrtbk(b,dr,m,ipd,ibuf,maxf,nv,neq,nev,aengy,ifl)

end
