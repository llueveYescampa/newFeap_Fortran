subroutine sproja(v,z,g,neq,nev)
implicit none
integer neq,nev
double precision   v(neq,*),z(neq,*),g(*)

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
   
   integer         maxa
include 'maxa.h'      
   double precision aa
   common /adata/   aa(maxa)

   integer         maxf
   common /frdata/ maxf

   integer         iodr,iodw,ipd,ipr,ipi
   common /iofild/ iodr,iodw,ipd,ipr,ipi

   double precision dimx,dimn
   integer                    nvb,npl
   common /nfrta/   dimx,dimn,nvb,npl

!  Forward reduce eigenvector estimates

   ma = maxf*nev + 1

!  Copy vectors 'v' into 'z'

   do i = 1,nev
     do j = 1,neq
       z(j,i) = v(j,i)
     end do  
   end do  

!  Solve equations

   call rsolve(aa,z,aa(ma),ipd,ipr,maxf,nvb,neq,nev,engy,4)

!  Compute projection of stiffness

   k = 0
   do j = 1,nev
     do i = 1,j
       k    = k + 1
       g(k) = dot(v(1,i),z(1,j),neq)
     end do
   end do  

end
