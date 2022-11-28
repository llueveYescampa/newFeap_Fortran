subroutine pdisk(disk,files)
implicit none
  character*1 :: disk,files(*)

!  Purpose: Set disk name character string

!  Inputs:
!     disk        - Disk label for files (e.g., c for DOS hard)

!  Output:
!     files(*)    - File name with disk label added

   integer i

   i = 1
   if((files(1).ne.' ') .and. (files(2).eq.':')) then
     i = 3
   end if  
   files(i) = disk
end
