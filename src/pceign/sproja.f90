subroutine sproja(v,z,g,neq,nev)
implicit none
  integer          :: neq,nev
  double precision :: v(neq,*),z(neq,*),g(*)

!  Compute subspace projection of 'aa' to form 'g'

!  Inputs:
!     aa(*)    - Tangent matrix to project
!     z(neq,*) - Iteration subspace (eigen) vectors
!     neq      - Number of equations in b(*)
!     nev      - Size of subspace 

!  Outputs:
!     g(*)     - Subspace projection of aa matrix


   integer  i, j, k, ma

   double precision   dot,engy
   
   include 'adata.h'   
   include 'frdata.h'
   include 'iofild.h'
   include 'nfrta.h'

!  Forward reduce eigenvector estimates

   ma = maxf*nev + 1

!  Copy vectors 'v' into 'z'

   do i = 1,nev
     do j = 1,neq
       z(j,i) = v(j,i)
     end do  
   end do  

!  Solve equations

  !call rsolve(aa,z,aa(ma),ipd,ipr,maxf,nvb,neq,nev,engy,4)
   call rsolve(   z,aa(ma),ipd,             neq,nev,engy  )

!  Compute projection of stiffness

   k = 0
   do j = 1,nev
     do i = 1,j
       k    = k + 1
       g(k) = dot(v(1,i),z(1,j),neq)
     end do
   end do  

end
