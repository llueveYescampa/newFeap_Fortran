subroutine dasol(al,au,ad,b,jp,neq,energy)
implicit none
  double precision  :: al(*),au(*),ad(*),b(*),energy
  integer           :: jp(*),neq

!  Purpose:
!     Solution of symmetric equations stored in profile form
!     coefficient matrix must be decomposed into its triangular
!     factors using datri before using dasol.

!  Inputs:
!     al(*)   - Lower part of factored matrix
!     au(*)   - Upper part of factored matrix
!     ad(*)   - Reciprocal diagonal of factored matrix
!     b(*)    - Right hand side of equations
!     jp(*)   - Pointer to end of row/columns in profile
!     neq     - Number of equations

!  Outputs:
!     b(*)    - Solution increment
!     energy  - Engergy of incremental solution


   integer  is, j,jh,jr
   double precision   bd,dot

   include 'iofile.h'

!  Find first nonzero entry in right hand side
   energy = 0.0d0

   do is = 1,neq
     if(b(is).ne.0.0d0) then
       exit
     end if  
   end do  

   if (is .gt. neq) then
     write(iow,'(a)') ' **WARNING** Zero right-hand-side vector'
     if(ior.lt.0) then
       write(*,'(a)') ' **WARNING** Zero right-hand-side vector'
     end if  
     return
   end if  
   
   if(is.lt.neq) then

!  Reduce right hand side

     do j = is+1,neq
       jr = jp(j-1)
       jh = jp(j) - jr
       if(jh.gt.0) then
         b(j) = b(j) - dot(al(jr+1),b(j-jh),jh)
       end if
     end do  
   end if

!  Multiply by inverse of diagonal elements

   do j = is,neq
     bd = b(j)
     b(j) = b(j)*ad(j)
     energy = energy + bd*b(j)
   end do  

!  Backsubstitution

   if(neq.gt.1) then
     do j = neq,2,-1
       jr = jp(j-1)
       jh = jp(j) - jr
       if(jh.gt.0) then
         call saxpb(au(jr+1),b(j-jh),-b(j),jh, b(j-jh))
       end if
     end do  
   end if

end
