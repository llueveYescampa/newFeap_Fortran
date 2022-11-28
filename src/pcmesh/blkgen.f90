subroutine blkgen(ndm,nel,nel1,x,ix,prt)
implicit none
  integer          :: ndm,nel,nel1,ix(nel1,*)
  double precision :: x(ndm,*)
  logical          :: prt
      
!  Purpose: Generate a block of nodes/elements

!  Inputs:
!     ndm        - Spatial dimension of mesh
!     nel        - Number of nodes/element
!     nel1       - Dimension of ix array
!     prt        - Flag, output results if true

!  Outputs:
!     x(ndm,*)   - Nodal coordinates of mesh
!     ix(nen1,*) - Element nodal connection lists


   character yyy*80
   integer iocheck
   integer i,j,k,l,ma,n,ne,nf,ng,ni,nm,nn,nr,ns,ntyp,nodinc,ixl(9)
   double precision r, s, t,shp(3,9),dr,ds,xl(3,9)

   include 'cdata.h'
   include 'iofile.h'

   do while (.true.)
     if(ior.lt.0) then
       write(*,'(a,/3x,a,$)')' Input: nn,nr,ns,ni,ne,ma,nodinc,ntyp','>'
     end if  
     call pintio(yyy,5)
     read(yyy,'(8i5)',IOSTAT=iocheck) nn,nr,ns,ni,ne,ma,nodinc,ntyp
     if (iocheck .ne. 0) then
       call myPerror('blkgen',yyy)
       cycle 
     end if
     exit
   end do
   nodinc = max(nodinc,0)
   nr = max(nr,1)
   ns = max(ns,1)
   ni = max(ni,1)
   ma = max(ma,1)
   if(prt) then
     call prthed(iow)
     write(iow,2000) nr,ns,ni,ne,ma,nodinc,ntyp
     if(ne.eq.0) then
       write(iow,'(a)') ' **WARNING** No elements are generated '
     end if  
     write(iow,2002) (i,i=1,ndm)
     if(ior.lt.0) then
       write(*,2000) nr,ns,ni,ne,ma,nodinc,ntyp
       if(ne.eq.0) then
         write(*,'(a)') ' **WARNING** No elements are generated '
       end if  
       write(*,2002) (i,i=1,ndm)
     end if
   end if
   
   do n = 1,9
     do j = 1,ndm
       xl(j,n) = 0.0
       ixl(n) = 0
     end do
   end do  
   nm = 0
   
   do n = 1,nn
     do while (.true.)
       if(ior.lt.0) then
         write(*,'(a,/3x,a,$)') ' Input: node, x-1, x-2, x-3','>'
       end if  
       call pintio(yyy,10)
       read(yyy,'(i10,3f10.0)',IOSTAT=iocheck) l,r,s,t
       if (iocheck .ne. 0) then
         call myPerror('blkgen',yyy)
         cycle
       end if
       exit
     end do

     if(l.eq.0) then
       l = n
     end if  
     nm = max(nm,l)
     ixl(l) = l
     xl(1,l) = r
     xl(2,l) = s
     xl(3,l) = t
     if(prt) then
       write(iow,'(i9,1p3e12.3)') l,(xl(i,l),i=1,ndm)
     end if  
   end do
   
   if(prt.and.ior.lt.0) then
     write(*,'(i9,1p3e12.3)') l,(xl(i,l),i=1,ndm)
   end if  
   dr = 2.d0/nr
   ds = 2.d0/ns
   if (ntyp.eq.0) then
      nf = ne + nr*ns - 1
   else if (ntyp.gt.7) then
      nf = ne + (nr*ns)/4 - 1
   else
      nf = ne + 2*nr*ns - 1
   end if
   if(nf.gt.numel.and.ne.gt.0) then
     write(iow,2007) nf,numel
     if(ior.lt.0) then
       write(*,2007) nf,numel
     end if  
     return
   end if  
   nr = nr + 1
   ns = ns + 1
   if(ndm.eq.1) then
     ns = 1
   end if  
   ng = nr*ns + ni -1
   if(ng.gt.numnp) then
     write(iow,2006) ng,numnp
     if(ior.lt.0) then
       write(*,2006) ng,numnp
     end if  
     return
   end if  

!  Form block

   call sblk(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm, &
        nel1,nodinc,ntyp,nm,ma,prt)

!  Print lists if wanted

   if(.not.prt) then
     return
   end if  

!  Print element lists

   if(ne.le.0) return
   do n = ne,nf,50
     call prthed(iow)
     write(iow,2003) (i,i=1,nel)
     if(ior.lt.0) then
       write(*,2003) (i,i=1,nel)
     end if  
     j = min(nf,n+49)
     do i = n,j
       write(iow,'(2i6,8i8/(13x,8i8))') i,ma,(ix(k,i),k=1,nel)
       if(ior.lt.0) then
         write(*,'(2i6,8i8/(13x,8i8))') i,ma,(ix(k,i),k=1,nel)
       end if  
     end do  
   end do  
   return
   
! Formats      

2000  format('    N o d e   G e n e r a t i o n s'//                   &
      9x,'Number of r-increments:',i5/9x,'Number of s-increments:',i5/ &
      9x,'First node number     :',i5/9x,'First element number  :',i5/ &
      9x,'Element material type :',i5/9x,'Node line increment   :',i5/ &
      9x,'Block type (0-9)       :',i5/1x)
     
2002  format(5x,'node',3(i6,' coord'))

2003  format('    E l e m e n t   C o n n e c t i o n s'//  &
       3x,'elmt  matl',8(i3,'-node')/(13x,8(i3,'-node')))

2006  format(' **ERROR** insufficient storage for nodes'/ &
         10x,'final node =',i5,5x,'numnp =',i5)
     
2007  format(' **ERROR** insufficient storage for elements'/ &
         10x,'final element =',i5,5x,'numel =',i5)

end
