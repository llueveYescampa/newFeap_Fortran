subroutine pfrtb(b,dr,nfrt,k,jj,eq,aengy)
implicit none
integer               nfrt,k,jj
double precision b(*),dr(*),eq(*),aengy

!  Purpose: Backsubstitution macro for frontal program

!  Inputs:
!     b(*)     - Front reduced equations
!     nfrt     - Size of current front
!     k        - Location of diagonal
!     jj       - Equation number in global equations
!     eq       - Current equation

!  Outputs
!     b(*)     - Remaining front reduced equations
!     dr(*)    - Solution increment
!     aengy    - Energy increment

   integer kk
   double precision dot

!  Expand b array

   kk = nfrt
   do while (.true.)
     if(kk.le.k) then
       exit
     end if  
     b(kk) = b(kk-1)
     kk = kk - 1
   end do
   
   b(k) = 0.0

!  Extract pivot and solve, also compute energy

   aengy = aengy + dr(jj)**2/eq(k)
   dr(jj)=(dr(jj)-dot(eq,b,nfrt))/eq(k)
   b(k) = dr(jj)

end
