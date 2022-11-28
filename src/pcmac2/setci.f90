subroutine setci(ioRead)
implicit none
  integer :: ioRead

!  Purpose: Compute integration constants 'c1' to 'c5' for current 'dt'

!  Inputs:
!     ioRead      - Data input logical unit number

!  Outputs:
!     none     - Output through common block

   include 'tbeta.h'
   include 'tdata.h'

   if(dt.le.0.0) then
     write(*,'(a)') ' **ERROR** Input DT as nonzero number.'
     if(ioRead.gt.0) call pstop(-18) ! stop
     return
   end if

!  Compute integration constants 'c1' to 'c5' for current 'dt'

   if(nop.eq.1) then
     c1 = 1.d0/(beta*dt)
     c2 = c1
     c3 = 1.d0/beta
     c6 = beta
   else if(nop.eq.2) then
     c5 = 2.d0/(gamm*dt)
     c1 = c5/dt
     c2 = c5*beta
     c3 = dt*beta
     c4 = 1.d0/gamm
     c6 = gamm
   else
     c1 = 1.d0/(beta*dt*dt)
     c2 = gamm/(dt*beta)
     c3 = 1.d0 - 1.d0/(beta+beta)
     c4 = 1.d0 - gamm/beta
     c5 = (1.d0 - gamm/(beta+beta))*dt
     c6 = dt*c1
   end if
end
