subroutine chlfwd(u,g,s,nn)
implicit none
integer nn
double precision   u(*),g(*),s(nn,nn)

!  Purpose: Forward solution using Choleski factors to project
!           onto a standard eigenproblem

!  Inputs:
!     u(*)     - Cholesky factor of !mass! subspace matrix
!     g(*)     - Subspace stiffness matrix
!     s(nn,*)  - Working storage

!  Outputs:
!     g(*)     - Standard eigenproblem matrix

   integer  i,id,im, j,jd
   double precision   dot

   s(1,1) = g(1)/u(1)
   if(nn .ne. 1) then
     id = 1
     do  i = 2,nn
       s(1,i) = g(id+1)/u(1)
       im = i - 1
       jd = 0
       do j = 1,im
         s(i,j) = (g(id+j) -dot(u(id+1),s(1,j),im))/u(id+i)
         if(j .gt. 1) then
           s(j,i) = (g(id+j) -dot(u(jd+1),s(1,i),j-1))/u(jd+j)
         end if  
         jd = jd + j
       end do  
       id     = id + i
       s(i,i) = (g(id) - dot(u(id-im),s(1,i),im))/u(id)
     end do  
   end if  

!  Complete projection

   g(1) = s(1,1)/u(1)
   if(nn.eq.1) then
     return
   end if  
   jd = 2
   do j = 2,nn
     g(jd) = s(j,1)/u(1)
     id    = 2
     do i = 2,j
       im       = i - 1
       g(jd+im) = (s(j,i) - dot(u(id),g(jd),im))/u(id+im)
       id       = id + i
     end do  
     jd = jd + j
   end do  

end
