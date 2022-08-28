double precision function coord(xl,shp,ndm,nel)
implicit    none
integer     ndm,nel
double precision  xl(ndm,*),shp(3,*)
      
   integer i
   
   coord = 0.d0
   do i = 1,nel
     coord = coord + shp(3,i)*xl(1,i)
   end do ! i

end
