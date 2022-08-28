subroutine pload(id,fn,b,nn,xm,ac)
implicit  none
integer   id(*),nn
double precision    fn(nn,2),b(*),xm(*),ac(*)

!  Purpose: Form load vector in compact form

!  Inputs:

!     id(*)     - Equation numbers for dof
!     fn(nn,2)  - Forces for active dofs at t_n and t_n+1
!     nn        - Size of id and fn arrays (total dof's)
!     xm(*)     - Diagonal mass array
!     ac(*)     - Acceleration vector

!  Outputs:
!     b(*)      - Residual vector

   integer   j, n

   logical        fl    ,pfr
   common /fdata/ fl(11),pfr

   double precision beta,gamm,theta
   integer                         nop,nt
   common /tbeta/  beta,gamm,theta,nop,nt

   fl(11) = .false.
   call pconsd(b,nn,0.0d0)
   do n = 1,nn
     j = id(n)
     if(j.gt.0) then
       b(j) =  b(j) + theta*fn(n,2) + (1. - theta)*fn(n,1) 
       if(fl(9)) b(j) =  b(j) - xm(j)*ac(n)
     end if
   end do  

end
