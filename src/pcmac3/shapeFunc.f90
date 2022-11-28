subroutine shapeFunc(ss,tt,xl,shp,xsj,ndm,nel,ix,flg)
implicit none
integer  ix(*),ndm,nel
double precision ss,tt,xl(ndm,1),shp(3,*),xsj
logical  flg
      



!  Purpose: Shape function routine for two dimensional elements

!  Inputs:
!    ss,tt     - Natural coordinates for evaluation
!    xl(ndm,*) - Nodal coordinates for element
!    ndm       - Spatial dimension of mesh
!    nel       - Number of nodes on element
!    ix(*)     - Element nodal connections
!    flg       - Flag, compute globagl derivatives if true

!  Outputs:
!    shp(3,*)  - Shape functions and derivatives
!    xsj       - Jacobian determinant at point


   integer  i, j, k 
   double precision tp,s(4),t(4),xs(2,2),sx(2,2)  

   data s/-.5d0,.5d0,.5d0,-.5d0/
   data t/-.5d0,-.5d0,.5d0,.5d0/

!  form 4-node quadrilateral shape functions
   do i = 1,4
     shp(3,i) = (0.5+s(i)*ss)*(0.5+t(i)*tt)
     shp(1,i) = s(i)*(0.5+t(i)*tt)
     shp(2,i) = t(i)*(0.5+s(i)*ss)
   end do  
   
   if(nel.eq.3) then
!  form triangle by adding third and fourth together

     do i = 1,3
       shp(i,3) = shp(i,3)+shp(i,4)
     end do  
   end if

!  add quadratic terms if necessary

   if(nel.gt.4) then 
     call shap2(ss,tt,shp,ix,nel)
   end if  

!  construct jacobian and its inverse

   do i = 1,2
     do j = 1,2
       xs(i,j) = 0.0
       do k = 1,nel
         xs(i,j) = xs(i,j) + xl(i,k)*shp(j,k)
       end do
     end do
   end do  
   xsj = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
   if(flg) then 
     return
   end if  
   if(xsj.le.0.0d0) then 
     xsj = 1.0
   end if  
   sx(1,1) = xs(2,2)/xsj
   sx(2,2) = xs(1,1)/xsj
   sx(1,2) =-xs(1,2)/xsj
   sx(2,1) =-xs(2,1)/xsj

!  form global derivatives

   do i = 1,nel
     tp        = shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
     shp(2,i)  = shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
     shp(1,i) = tp
   end do  

end
