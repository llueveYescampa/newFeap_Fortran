subroutine shap2(s,t,shp,ix,nel)
implicit  none
integer   ix(*), nel
double precision    s,t,shp(3,*)
      
!  Purpose: Add quadratic shape functions as necessary

!  Inputs:
!    s,t     - Natural coordinate point for evaluation
!    ix(*)   - Element connection list
!    nel     - Number of nodes on element

!  Outputs:
!    shp(3,*) - Quadratic shape functions generated


   integer   i, j, k, l
   double precision s2,t2
   
   s2 = (1.-s*s)*0.5
   t2 = (1.-t*t)*0.5
   do i=5,9
     do j = 1,3
       shp(j,i) = 0.0d0
     end do
   end do  

!  midside nodes (serendipity)

   if(ix(5).ne.0) then 
     shp(1,5) = -s*(1.-t)
     shp(2,5) = -s2
     shp(3,5) = s2*(1.-t)
   end if  
   

   do while (.true.)
   
     if(nel.lt.6) then
       exit
     else if(ix(6).ne.0) then
       shp(1,6) = t2
       shp(2,6) = -t*(1.+s)
       shp(3,6) = t2*(1.+s)
     end if  
     
     
     if(nel.lt.7) then 
       exit
     else if(ix(7).ne.0) then 
       shp(1,7) = -s*(1.+t)
       shp(2,7) = s2
       shp(3,7) = s2*(1.+t)
     end if  
     
     if(nel.lt.8) then 
       exit
     else if(ix(8).ne.0) then   
       shp(1,8) = -t2
       shp(2,8) = -t*(1.-s)
       shp(3,8) = t2*(1.-s)
     end if  
     
!    interior node (lagrangian)
     if(nel.lt.9) then 
       exit
     else if(ix(9).ne.0) then 
       shp(1,9) = -4.*s*t2
       shp(2,9) = -4.*t*s2
       shp(3,9) = 4.*s2*t2
     end if  

!    correct edge nodes for interior node (lagrangian)
     do j= 1,3
       do i = 1,4
         shp(j,i) = shp(j,i) - 0.25*shp(j,9)
       end do  
       do i = 5,8
         if(ix(i).ne.0) shp(j,i) = shp(j,i) - .5*shp(j,9)
       end do  
     end do  
     exit
   end do !! while

!  correct corner nodes for presense of midside nodes
   do i = 1,4
     k = mod(i+2,4) + 5
     l = i + 4
     do j = 1,3
       shp(j,i) = shp(j,i) - 0.5*(shp(j,k)+shp(j,l))
     end do
   end do  

end
