subroutine genvec(ndm,x,cd,prt,err,prtz)
implicit none
  integer          ::  ndm
  double precision :: x(ndm,*)
  character        :: cd*12
  logical          :: prt,err,prtz

!  Purpose: Generate real data arrays by linear interpolation

!  Inputs:
!     ndm        - Dimension of x array
!     cd         - Character label for outputs
!     prt        - Flag, output results if true
!     prtz       - Flag, output zero results if true

!  Outputs:
!     x(ndm,*)   - Nodal coordinates of mesh
!     err        - Flag, true if errors occurred


   character yyy*80
   integer   i, j, l,li,lg, mct, n,ng
   integer iocheck
   logical test
   

   double precision    xl(6)

   include 'cdata.h'

   include 'iofile.h'
   
   mct = 0
   n = 0
   ng = 0
   do while (.true.)
     l = n
     lg = ng
     do while (.true.)
       if(ior.lt.0) then
         write(*,'(a,a12,a,/3x,a,$)') &
           ' Input ',cd,' values: node, inc, value(i),i=1,nval','>'
       end if  
       call pintio(yyy,10)
       read(yyy,'(2i10,6f10.0)',IOSTAT=iocheck) n,ng,(xl(i),i=1,ndm)
       if (iocheck .ne. 0) then
         call myPerror('genvec',yyy)
         cycle
       else
         exit
       end if
     end do
     if(n.gt.numnp) then
       write(iow,3001) n,cd
     end if  
     if(ior.lt.0.and.n.gt.numnp) then
       write(*,3001) n,cd
     end if
     if(n.le.0 .or. n.gt.numnp) then
       if(.not.prt) then
         return
       end if  
       do j = 1,numnp
         if(.not. prtz) then
           test= .true.
           do l = 1,ndm
             if(x(l,j).ne.0.0) then
               test = .false.
               exit
             end if  
           end do
           if (test) then
             cycle
           end if
         end if
         mct = mct - 1
         if(mct.le.0) then
           mct = 50
           write(iow,'(a,a12//6x,a,9(i7,a6))') &
             '    N o d a l: ',cd,'node',(l,cd,l=1,ndm)
           if(ior.lt.0) then
             write(*,'(a,a12//6x,a,9(i7,a6))') &
             '    N o d a l: ',cd,'node',(l,cd,l=1,ndm)
           end if  
         end if
         if(x(1,j).eq.-999.0) then
           write(iow,'(i10,a)') j,' has not been input or generated'
         else  
           write(iow,'(i10,9f13.4)') j,(x(l,j),l=1,ndm)
         end if  
         if(ior.lt.0) then
           if(x(1,j).eq.-999.0) then
             write(*,'(i10,a)') j,' has not been input or generated'
           else 
             write(*,'(i10,9f13.4)') j,(x(l,j),l=1,ndm)
           end if  
         end if
       end do
       exit
     else
       do i = 1,ndm
         x(i,n) = xl(i)
       end do
       if(lg .ne. 0) then
         lg = sign(lg,n-l)
         li =(abs(n-l+lg)-1)/abs(lg)
         do i = 1,ndm
           xl(i) = (x(i,n)-x(i,l))/li
         end do
         do while (.true.)
           l = l + lg
           if((n-l)*lg.le.0) then
             test = .false.
             exit
           end if  
           if(l.le.0.or.l.gt.numnp) then
             test = .true.
             exit
           end if  
           do i = 1,ndm
              x(i,l) = x(i,l-lg) + xl(i)
           end do
         end do
         if (test) then
           write(iow,'(a,i5,a,a12)') &
             ' **ERROR** attempt to generate node',l,' in ',cd
           if(ior.lt.0) then
             write(*,'(a,i5,a,a12)') &
             ' **ERROR** attempt to generate node',l,' in ',cd
           end if  
           err = .true.
         end if  
       end if  
     end if
   end do


3001  format(' **ERROR** attempt to input node',i5,', terminate'    &
       ,' input of nodes in ',a12)
       
end
