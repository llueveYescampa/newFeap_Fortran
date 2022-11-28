subroutine pbcin(iii,idl,id,numnp,ndf,prt)
implicit none
  integer :: iii,idl(*),ndf,id(ndf,*),numnp
  logical :: prt

!  Purpose: Input restraint conditions for each node

!  Inputs:
!     iii        - Flag to test change in b.c. inputs
!     idl(*)     - Working vector
!     numnp      - Number nodes in mesh
!     ndf        - Number dof/node
!     prt        - Flag, output results if true

!  Outputs:
!     id(ndf,*)  - Boundary condition codes


   character yyy*80
   integer   i, l,lg, n,ng
   integer iocheck
   logical test

   include 'iofile.h'

!  Read in restraint conditions for each node

   iii = 1
   n = 0
   ng = 0
   test = .true.
   do while(test)
     test = .false.
     l = n
     lg = ng
     do while (.true.)
       if(ior.lt.0) then 
         write(*,'(a/3x,a,$)') ' Input: node, inc, b.c. codes(i),i=1,ndf','>'
       end if  
       call pintio(yyy,10)
       read(yyy,'(8i10)',IOSTAT=iocheck) n,ng,(idl(i),i=1,ndf)
       if (iocheck .ne.0) then
         call myPerror('pbcin ',yyy)
         cycle
       else  
         if(n.gt.0.and.n.le.numnp) then
           do i = 1,ndf
             id(i,n) = idl(i)
             if(l.ne.0.and.idl(i).eq.0.and.id(i,l).lt.0) then
               id(i,n) = -1
             end if  
           end do  
           lg = sign(lg,n-l)
           do while (.true.)
             l = l + lg
             if((n-l)*lg.le.0) then
               test = .true.
               exit
             else  
               do i = 1,ndf
                 if(id(i,l-lg).lt.0) then 
                   id(i,l) = -1
                 end if  
               end do  
             end if  
           end do
         end if
       end if
       exit
     end do
   end do

!  Output nodes with nonzero codes

   if(prt) then
     call prthed(iow)
     write(iow,2000) (i,i=1,ndf)
     do n = 1,numnp
       do l = 1,ndf
         if(id(l,n).ne.0) then
           write(iow,'(i10,8i8)') n,(id(i,n),i=1,ndf)
           exit
         end if
       end do  
     end do  
   end if
   return
      
!  Formats

2000  format('    N o d a l   B. C.'//6x,'node',8(i3,' b.c.')/1x)

end
