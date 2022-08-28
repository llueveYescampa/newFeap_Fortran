subroutine prthed(iow)
implicit none
integer  iow

!  Purpose: Output a header to printed outputs

!  Inputs:
!     iow      - Logical unit number for outputs

   character     head*4
   common/bdata/ head(20)

   write(iow,'(/1x,20a4//1x)') head
end
