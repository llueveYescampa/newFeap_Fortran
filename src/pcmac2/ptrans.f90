subroutine ptrans(ia,angl,ul,p,s,nel,nen,ndf,nst,isw)
implicit  none
integer   ia(2),nel,nen,ndf,nst,isw
double precision    angl(*),ul(ndf,nen,*),p(ndf,*),s(nst,nst)

!  Purpose: Subroutine to make two-dimesional rotations

!  Inputs:
!     ia(2)       - Dof's to rotate
!     angl(*)     - Angle of sloping boundary condition at element nodes
!     nel         - Number of nodes on element
!     nen         - Maximum number of nodes on element
!     ndf         - Number dof/node
!     nst         - Dimension of element arrays
!     ist         - Switch to control actions taken

!  Outputs:
!     ul(ndf,*,*) - Rotated element solution vector/rates (isw = 1)
!     p(ndf,*)    - Rotated element vector (isw = 2)
!     s(nst,*)    - Rotated element matrix (isw = 2)


   integer i, j, i1, ij1,ij2
   double precision ang,cs,sn,tm


!  Recover dof to be rotated

   ij1 = ia(1)
   ij2 = ia(2)
   i1  = 0

!  Transform displacement quantities to element coordinates

   if(isw.eq.1) then
     do i = 1,nel
       if(angl(i).ne.0.0) then
         ang = angl(i)*0.01745329252d0
         cs  = cos(ang)
         sn  = sin(ang)
         do j = 1,4
           tm          = cs*ul(ij1,i,j) - sn*ul(ij2,i,j)
           ul(ij2,i,j) = sn*ul(ij1,i,j) + cs*ul(ij2,i,j)
           ul(ij1,i,j) = tm
         end do  
       end if
       i1 = i1 + ndf
     end do

!  Transform element arrays to global coordinates

   else
     do i = 1,nel
       if(angl(i).ne.0.0) then
         ang = angl(i)*0.01745329252d0
         cs  = cos(ang)
         sn  = sin(ang)

!        Transform load vector

         tm       = cs*p(ij1,i) + sn*p(ij2,i)
         p(ij2,i) =-sn*p(ij1,i) + cs*p(ij2,i)
         p(ij1,i) = tm

!        Postmultiply s by transformation

         do j = 1,nst
           tm         = s(j,i1+ij1)*cs + s(j,i1+ij2)*sn
           s(j,i1+ij2)=-s(j,i1+ij1)*sn + s(j,i1+ij2)*cs
           s(j,i1+ij1)= tm
         end do  

!        Premultiply s by transformation

         do j = 1,nst
           tm         = cs*s(i1+ij1,j) + sn*s(i1+ij2,j)
           s(i1+ij2,j)=-sn*s(i1+ij1,j) + cs*s(i1+ij2,j)
           s(i1+ij1,j)= tm
         end do
       end if
       i1 = i1 + ndf
     end do
   end if
end
