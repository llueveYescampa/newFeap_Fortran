subroutine dredu(al,au,ad,jh,flg,dj)
implicit none
integer  jh
double precision al(jh),au(jh),ad(jh),dj
logical  flg

!  Purpose: Reduce diagonal element in triangular decomposition

!  Inputs:
!     al(*)    - Loser row of equation
!     au(*)    - Upper column of equation
!     ad(*)    - Reciprocal diagonals for equations
!     jh       - length of row/column 
!     flg      - Flag, equations unsymmetric if true

!  Outputs:
!     dj       - Reduced diagonal


   integer  j
   double precision  ud

!  Computation of column for unsymmetric matrices

   if(flg) then
     do j = 1,jh
       au(j) = au(j)*ad(j)
       dj    = dj - al(j)*au(j)
       al(j) = al(j)*ad(j)
     end do  

!  Computation of column for symmetric matrices

   else
     do j = 1,jh
       ud    = au(j)*ad(j)
       dj    = dj - au(j)*ud
       au(j) = ud
     end do  
   end if

end
