subroutine psetm(na,nl,typ,afl)
implicit none
  integer   :: na,nl
  character :: typ*1
  logical   :: afl

!  Purpose:  Allocate array storage pointer from blank common

!  Inputs:
!     nl        - length of array requested
!     typ       - type of array requested ('i', 'r', 'd')

!  Outputs:
!     na        - Pointer to first word in array
!     afl       - Flag, error if true

   integer          :: np=0,ns=0
   double precision :: amx

   include 'iofild.h'
   include 'iofile.h'
   include 'psize.h'

!  Set data management pointers for arrays

!  alignment for the blank common d(1),r(1),i(xx): xx = size
!  precisions: 'd' = ipd ; 'r' = ipr ; 'i' = ipi

!  N.B. ipd >= ipr >= ipi

   if( typ .eq. 'd' ) then
     np = ipd
     ns = 0
!   else if( typ .eq. 'r' ) then
!     np = ipr
!     ns = ipd/ipr
   else if( typ .eq. 'i' ) then
     np = ipi
!     ns = (ipd+ipr)/ipi
     ns = (ipd)/ipi
   end if
   na = ne
   ne = na + nl*np + mod(ipd - mod(nl*np,ipd),ipd)
   na   = (na + np - 1)/np - ns
   afl = .false.
   amx = maxm
   amx = ne/amx
   if(amx.gt.0.90) then
     write(*,'(a,i6,a,i6,a,f6.3)') &
     '  **Memory warning** Used =', ne,' Avail =',  maxm, ' % =', amx
   end if
     
!   if(ne.le.maxm) then
!     return
!   end if  
   
   if(ne.gt.maxm) then
     write(iow,'(2x,a,/10x,a,i6/10x,a,i6/)')             &
     '**ERROR** Insufficient storage in blank common',   &
     'Required  =', ne,                                  &
     'Available =', maxm
     
     if(ior.lt.0) then
       write(*,'(2x,a,/10x,a,i6/10x,a,i6/)')               &
       '**ERROR** Insufficient storage in blank common',   &
       'Required  =', ne,                                  &
       'Available =', maxm
     end if  
     call pstop(-68) ! stop
   end if  
   !print *, "na: ", na, " allocated space: ", nl, " type: ", typ
   
end


!   if( typ == 'd' ) then
!     na = neDouble
!     
!     neDouble = neDouble + nl
!     neInt = 2*neDouble - 3
!   else  ! assuming typ .eq. 'i'
!     na = neInt
!     neInt = neInt + nl
!     if (mod(neInt,2) == 1) then
!       neDouble = (neInt+3)/2
!     else
!       neDouble = (neInt+4)/2
!     endif
!   endif

