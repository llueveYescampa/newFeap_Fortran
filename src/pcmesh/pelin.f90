subroutine pelin(idl,ix,nen1,nen,numnp,numel,error,prt)
implicit  none
integer   idl(*),nen1,ix(nen1,*),nen,numnp,numel
logical   error,prt

!  Purpose: Element data input

!  Inputs:
!     idl(*)     - Working vector
!     nen1       - Dimension of ix array
!     nen        - Number nodes/element
!     numnp      - Number of nodes in mesh
!     numel      - Number of elements in mesh
!     prt        - Flag, output results if true

!  Outputs:
!     ix(nen1,*) - Element nodal connection lists
!     error      - Flag, true if error occurs


   !character yyy*80
   integer i, j, k, l,lk,lx, n,nx,iocheck
   logical test208

   include 'iofile.h'

   include 'ydata.h'

!  Element data input
   nx =0
   l = 0
   do i = 1,numel,50
     if(prt) then
       call prthed(ioWrite)
       write(ioWrite,2001) (k,k=1,nen)
       if(ioRead.lt.0) then
         write(  *,2001) (k,k=1,nen)
       end if  
     end if  
     j = min(numel,i+49)
     
     do n = i,j
       test208=.false.
       if (l .lt. n) then
         do while (.true.)
           if(ioRead.lt.0) then
             write(*,'(a/3x,a,$)')' Input: elm, mat, ix(i),i=1,nen, inc','>'
           end if  
           call pintio(yyy,5)
           read(yyy,'(16i5)',IOSTAT=iocheck) l,lk,(idl(k),k=1,nen),lx
           if (iocheck .ne. 0) then
              call pperror('PELIN ',yyy)
              cycle
           else
             exit
           end if
         end do
         
         if(l.eq.0) then
           l = numel+1
         end if  
         if(lx.eq.0) then
           lx=1
         end if  
         if (l .lt. n) then
           write(ioWrite,'(a,i5,a,i5)') &
           ' **ERROR** element',l,' appears after element',n
           if(ioRead.lt.0) then
             write(*,'(a,i5,a,i5)') &
           ' **ERROR** element',l,' appears after element',n
           end if  
           error = .true.
           cycle
         end if
       end if
       if (l .eq. n) then
         nx = lx
         do k = 1,nen
           if(idl(k).gt.numnp.or.idl(k).lt.0) then
             test208=.true.
             exit
           end if
           ix(k,l) = idl(k)
         end do  
         if (.not. test208) then
           ix(nen1,l) = lk
         end if  
       else if (l .gt. n) then
         ix(nen1,n) = ix(nen1,n-1)
         do k = 1,nen
           ix(k,n) = ix(k,n-1) + nx
           if(ix(k,n-1).eq.0) then
             ix(k,n) = 0
           end if  
           if(ix(k,n).gt.numnp.or.ix(k,n).lt.0) then
             test208=.true.
             exit
           end if  
         end do  
       end if
       if (test208) then
         write(ioWrite,'(a,i5,a)') ' **ERROR** element',n,' has illegal nodes'
         error = .true.
       else
         if(prt) then 
           write(ioWrite,'(2i6,8i8/(13x,8i8))')n,ix(nen1,n),(ix(k,n),k=1,nen)
           if(ioRead.lt.0) then
             write(*,'(2i6,8i8/(13x,8i8))')n,ix(nen1,n),(ix(k,n),k=1,nen)
           end if
         end if  
       end if
     end do
   end do  
   
!  Formats

2001  format('    E l e m e n t   C o n n e c t i o n s'// &
         3x,'elmt  matl',8(i3,'-node')/(13x,8(i3,'-node')))

end
