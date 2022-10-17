subroutine psetm(na,nl,typ)
implicit  none
integer   na,nl
character typ*1

!  Purpose:  Allocate array storage pointer from blank common

!  Inputs:
!     nl        - length of array requested
!     typ       - type of array requested ('i', 'r', 'd')

!  Outputs:
!     na        - Pointer to first word in array

   integer   np,ns
   double precision    amx
  
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
   else if( typ .eq. 'r' ) then
     np = ipr
     ns = ipd/ipr
!   else if( typ .eq. 'i' ) then
   else  ! assuming typ .eq. 'i'
     np = ipi
     ns = (ipd+ipr)/ipi
   end if
   na = ne
   ne = na + nl*np + mod(ipd - mod(nl*np,ipd),ipd)
   na   = (na + np - 1)/np - ns
   amx = maxm
   amx = ne/amx
   if(amx.gt.0.90) then
     write(*,'(a,i6,a,i6,a,f6.3)') &
     '  **Memory warning** Used =', ne,' Avail =',  maxm, ' % =', amx
   end if  
   
   
!   if(ne.le.maxm) then
!     return
!   end if  
   
   if(ne >= maxm) then
     write(ioWrite,'(2x,a,/10x,a,i6/10x,a,i6/)')         &
     '**ERROR** Insufficient storage in blank common',   &
     'Required  =', ne,                                  &
     'Available =', maxm
     
     if(ioRead.lt.0) then
       write(*,'(2x,a,/10x,a,i6/10x,a,i6/)')               &
       '**ERROR** Insufficient storage in blank common',   &
       'Required  =', ne,                                  &
       'Available =', maxm
     end if  
     stop
   end if  
   
   return
end