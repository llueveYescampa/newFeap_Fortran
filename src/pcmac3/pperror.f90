subroutine pperror(subnam,yy)
implicit  none
character subnam*6, yy*80

! Purpose: Read error encountered

! Inputs:
!    subnam    - Name of subroutine where ERROR encountered
!    yy        - String of data read


   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite
   
   if(ioRead.gt.0) then
     write(ioWrite,'(a,a6,a,/1x,a80)') &
       ' **ERROR in ', subnam, '** Reinput last record:', yy  
     stop
   else
     write(*,'(a,a6,a,/1x)') &
       ' **ERROR in ', subnam, '** Reinput last record:'
   end if
end
