subroutine chlbac(u,s,nn)
implicit none
integer nn
double precision u(*),s(nn,nn)

!  Purpose: Compute eigenvalues of general linear problem

!  Inputs:
!     u(*)     - Cholesky factor of !mass! subspace matrix
!     s(nn,*)  - Eigenvectors of standard problems
!     nn       - Size of subspace probleml

!  Outputs:
!     s(nn,*)  - Eigenvectors of generalized problem

   integer  i,j,jd
   
!  Compute eigenvalues of general linear problem by backsubstitution

   j  = nn
   jd = nn*(nn+1)/2
   do i = 1,nn
     s(nn,i) = s(nn,i)/u(jd)
   end do  

   do while (.true.)
     jd = jd - j
     j  = j - 1
     if(j.le.0) then 
        exit  ! return
     else
       do i = 1,nn
         call colbac(u(jd+1),s(1,i),u(jd),j)
       end do  
       cycle
     end if
   end do
end
