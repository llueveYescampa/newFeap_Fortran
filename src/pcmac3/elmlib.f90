subroutine elmlib(d,u,x,ix,t,s,p,i,j,k,iel,isw)
implicit  none
integer ix(*),i, j, k, iel, isw
double precision d(*),u(*),x(*),t(*),s(*),p(*)

!  Purpose: Element library driver

!  Inputs:
!     d(*)      - Material set parameters
!     u(*)      - Element nodal solution parameters
!     x(*)      - Element nodal coordinates
!     ix(*)     - Element node numbers
!     t(*)      - Element nodal temperatures
!     i         - Number dof/node
!     j         - Spatial dimension of mesh
!     k         - Size of element arrays
!     iel       - Element type number
!     isw       - Switch for action to be taken by element routine

!  Outputs:
!     s(k,*)    - Element array
!     p(*)      - Element vector

   include 'iofile.h'

   if(isw.ge.3.and.isw.le.6) then
     call pconsd(p,k,0.0d0)
     call pconsd(s,k*k,0.0d0)
   end if
   select case (iel)
   case (1)
     call elmt01(d,u,x,ix,t,s,p,i,j,k,isw)
   case (2)
     call elmt02(d,u,x,ix,t,s,p,i,j,k,isw)
   case (3)
     call elmt03(d,u,x,ix,t,s,p,i,j,k,isw)
   case (4)
     call elmt04(d,u,x,ix,t,s,p,i,j,k,isw)
   case (5)
     call elmt05(d,u,x,ix,t,s,p,i,j,k,isw)
   case (6)
     call elmt06(d,u,x,ix,t,s,p,i,j,k,isw)
   case (7)
     call elmt07(d,u,x,s,p,i,j,k,isw)
   case default
     write(ioWrite,'(a,i3,a)') &
      '  **ERROR** Element type number',iel,' found.' 
     if(ioRead.lt.0) then 
       write(*,'(a,i3,a)') &
      '  **ERROR** Element type number',iel,' found.' 
     end if  
     stop
   end select
      
end
