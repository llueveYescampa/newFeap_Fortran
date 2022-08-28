subroutine chlfac(a,nn)
implicit none
integer nn
double precision  a(*)
      
!  Purpose: Choleski factorization of a symmetric, positive definite matrix

!  Inputs
!     a(*)     - Matrix to factor
!     nn       - Size of matrix

!  Outputs
!     a(*)     - Factored matrix

   integer  i,id, j,jd,jm
   double precision   dot

   a(1) = sqrt(a(1))
   if(nn.eq.1) then 
     return
   end if  
   jd = 1
   do j = 2,nn
     jm = j - 1
     id = 0
     do i = 1,jm
       if(i-1.gt.0) then 
         a(jd+i) = a(jd+i) - dot(a(id+1),a(jd+1),i-1)
       end if  
       id      = id + i
       a(jd+i) = a(jd+i)/a(id)
     end do  
     a(jd+j) = sqrt(a(jd+j) - dot(a(jd+1),a(jd+1),jm))
     jd      = jd + j
   end do  

end
