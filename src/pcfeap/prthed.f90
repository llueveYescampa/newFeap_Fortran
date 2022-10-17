subroutine prthed(iow)
implicit none
integer  iow

!  Purpose: Output a header to printed outputs

!  Inputs:
!     iow      - Logical unit number for outputs

   include 'bdata.h'

   write(iow,'(/1x,20a4//1x)') head
end
