subroutine matmult(a,b,c,m,l,n,test)
implicit none
integer m,n,l,test
double precision a(m,l),b(l,n),c(m,n)

!   
!  IF TEST == 0
!     c = a * b
!  IF TEST == 1
!     c = Transpose[a]*b
!
   integer i,j,k
   double precision d
   
   do j=1,n
     do i=1,m
       c(i,j)= 0.0
     end do
     do k=1,l
       do i=1,m
         if (test .eq. 0) then
             d = a(i,k)
         else
             d = a(k,i)
         end if    
         c(i,j) = c(i,j) + d*b(k,j)
       end do
     end do  
   end do
end
