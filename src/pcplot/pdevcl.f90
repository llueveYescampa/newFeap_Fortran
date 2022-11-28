subroutine pdevcl()
implicit none

!  Purpose: Close the plotting device

   integer*2 status,vclrwk,vencur,vclswk,vrqstr
!  integer*2 ixy(2)   ! yo la elimine

   character*1 xxx

!  status = vrqstr(1,0,ixy,xxx)     ! yo la modifique

   status = vrqstr(xxx)
   status = vclrwk()
   status = vencur()
   status = vclswk()
end
