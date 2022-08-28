subroutine sprojb(b,v,t,h,neq,nev)
implicit none
integer neq,nev
double precision b(1),v(neq,1),t(1),h(1)

!  Purpose: Compute subspace projection of 'b' to form 'h'

!  Inputs:
!     b(*)     - Diagonal !mass! matrix
!     v(neq,*) - Iteration subspace (eigen) vectors
!     t(*)     - Working vector
!     neq      - Number of equations in b(*)
!     nev      - Size of subspace 

!  Outputs:
!     h(*)     - Subspace projection of b matrix

   integer  i, j, k
   double precision dot

!  Compute 'z' and 'b' projection to form 'h'

   do j = 1,nev
!    Compute 'z' for a lumped mass

     do i = 1,neq
       t(i) = v(i,j)*b(i)
     end do  
!    Project 'z' and 'v' vectors to form 'h'

     k = j*(j+1)/2
     do i = j,nev
       h(k) = dot(t,v(1,i),neq)
       k    = k + i
     end do  
     do i = 1,neq
       v(i,j) = t(i)
     end do  
   end do  

end
