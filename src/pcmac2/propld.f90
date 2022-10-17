function propld(t,j)
implicit  none
integer   j
double precision    t

!  Purpose: Proportional load table (j load records, maximum 10)

!  Inputs:
!     t      - Time
!     j      - Number of proportional loads to input

!  Outputs:
!     propld - Proportional loading for time 't' (when j = 0)


   !character yyy*80
   integer   i,l,m, nprop
   integer iocheck
   double precision    tmin,tmax,propld

   include 'iofile.h'

   include 'ydata.h'

   include 'prlod.h'

   if(j.gt.0) then
!    Input table of proportional loads
     write(ioWrite,'(30x,a)')'P r o p o r t i o n a l   L o a d   T a b l e'
     if(ioRead.lt.0) then
       write( *,'(30x,a)')'P r o p o r t i o n a l   L o a d   T a b l e'
       write( *,'(a,/,a,$)') &
             ' Input: type, exponent, tmin, tmax, a(i),i=1,4',' >'
     end if
     do i=1,j
       do while (.true.)
         call pintio(yyy,10)
         read(yyy,'(2i10,6f10.0)',IOSTAT=iocheck) ik(i),iexp(i),(a(m,i),m=1,6)
         if (iocheck .eq. 0) then
!          Set a default ramp table if a type !0! input
           if(ik(i).eq.0) then
             a(2,i) = 1.e+6
             a(4,i) = 1.
           end if
           write(ioWrite,2001) i,ik(i),(a(m,i),m=1,6),iexp(i)
           if(ioRead.lt.0) then 
             write(*,2001) i,ik(i),(a(m,i),m=1,6),iexp(i)
           end if
           exit
         else
!          Error message
           call pperror('PROPLD',yyy)
           cycle
         end if
       end do   
     end do  
     nprop = j
   end if

!  Compute value at time t

   propld = 0.0
   do i = 1,nprop
     tmin = a(1,i)
     tmax = a(2,i)
     if(t.lt.tmin.or.t.gt.tmax) then
       cycle
     end if  
     l = max(iexp(i),1)
     propld = a(3,i)+a(4,i)*t+ a(5,i)*(sin(a(6,i)*t+tmin))**l+ propld
   end do

!  Formats

2001  format(/,' number    type      tmin',10x,'tmax',/i3,i10,7x,g10.4, &
       4x,g10.4,/6x,'a(1)',10x,'a(2)',10x,'a(3)',10x,'a(4)',10x,        &
       'exp',/4x,4(g10.4,4x),i5/)

end
