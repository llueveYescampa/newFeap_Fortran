subroutine sblk(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne, & 
                 ndm,nel1,nodinc,ntyp,nm,ma,prt)
implicit  none
integer nr,ns,ixl(*),ni,ne,ndm,nel1,nodinc,ntyp,nm,ma,ix(nel1,*)
double precision  xl(3,*),shp(3,*),x(ndm,*),dr,ds
logical prt

!  Inputs:
!     nr         - Number increments in 'r' direction
!     ns         - Number increments in 's' direction
!     xl(ndm,*)  - Master element nodal coordinates
!     ixl(*)     - Master element nodal connections list
!     shp(3,*)   - Shape function storage
!     dr         - 'r' coordinate increment
!     ds         - 's' coordinate increment
!     ni         - Initial node number
!     ne         - Initial element number
!     ndm        - Spatial dimension of mesh
!     nel1       - Dimension of ix array
!     nodinc     - Increment to node numbers at end of 'r' line
!     ntyp       - Block type
!     nm         - Number of nodes on master element
!     ma         - Material set number for block
!     prt        - Flag, output results if true

!  Outputs:
!     x(ndm,*)   - Nodal coordinates of mesh
!     ix(nel1,*) - Element nodal connection lists


   integer i,inc, j, k, l, m,me,mct, n
   double precision  r,s,xsj

   include 'cdata.h'

   include 'iofile.h'

   n = ni
   mct = 0
   s = -1.0
   do j = 1,ns
     r = -1.0
       do i = 1,nr
       call shapeFunc(r,s,xl,shp,xsj,3,nm,ixl,.true.)
       do l = 1,ndm
         x(l,n) = 0.0
         do k = 1,9
           m = ixl(k)
           if(m.gt.0) x(l,n) = x(l,n) + shp(3,m)*xl(l,m)
         end do  
       end do
       
       if(prt) then
          mct = mct + 1
          if(mod(mct,50).eq.1) then
            call prthed(ioWrite)
            write(ioWrite,2000) (k,k=1,ndm)
            if(ioRead.lt.0) then
              write(*,2000) (k,k=1,ndm)
            end if  
          end if
          write(ioWrite,'(i10,3f13.4)') n,(x(k,n),k=1,ndm)
          if(ioRead.lt.0) then
            write(*,'(i10,3f13.4)') n,(x(k,n),k=1,ndm)
          end if  
       end if
       n = n + 1
       r = r + dr
     end do
     n = n + nodinc
     s = s + ds
   end do
   if(ne.le.0) then
     return
   end if  
   me = ne - 1
   n = ni
   inc = 1
   if(ntyp.ge.8) inc = 2
   do j = 1,ns-1,inc
     do i = 1,nr-1,inc
       n = n + 1
       me = me + 1
       ix(nel1,me) = ma
       if(ntyp.eq.0) then
          ix(1,me)    = n - 1
          ix(2,me)    = n
          if(ndm.ne.1) then
             ix(3,me)    = n + nr + nodinc
             ix(4,me)    = n + nr - 1 + nodinc
          endif
       else if(ntyp.ge.8) then
          ix(1,me) = n-1
          ix(5,me) = n
          ix(2,me) = n+1
          ix(8,me) = nr+nodinc + n-1
          if(ntyp.gt.8) ix(9,me) = nr+nodinc + n
          ix(6,me) = nr+nodinc + n+1
          ix(4,me) = 2*(nr+nodinc) + n-1
          ix(7,me) = 2*(nr+nodinc) + n
          ix(3,me) = 2*(nr+nodinc) + n+1
          n = n+1
       end if
     end do
     n = n + (inc-1) * nr + nodinc + 1
   end do  

2000 format('    N o d a l   C o o r d i n a t e s'//6x,'node',3(i7,' coord'))
end
