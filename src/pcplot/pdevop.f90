subroutine pdevop()
implicit none

! Purpose: Open graphics workstation


   integer*2 status,vopnwk,vclrwk,vslcol
   integer*2 icl
  
!  Open kernel system
  
   status = vopnwk()
   if(status.lt.0) then
     write(*,'(a)') ' Graphics Device Driver not installed rectly.'
     return
   end if
   status = vclrwk()
   icl    = 2
   status = vslcol(icl)
   call dplot(0.0000d0,0.0000d0,1)
   call dplot(1.4545d0,0.0000d0,2)
   call dplot(1.4545d0,1.0090d0,2)
   call dplot(0.0000d0,1.0090d0,2)
   call dplot(0.0000d0,0.0000d0,2)
end
