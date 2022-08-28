subroutine pdefil(tfile,n1,n2)
implicit   none
integer n1,n2
character  tfile(n2)*12

!  Purpose: Destroy temporary files used to solve problems

!  Inputs:
!    tfile(*)  - Filenames to delete
!    n1        - First number of array to delete
!    n2        - Last  number of array to delete


   logical    lfil
   integer    n
   do n = n1,n2
     inquire(file=tfile(n),exist=lfil)
     if(lfil) then
       open (4,file=tfile(n),status='old')
       close(4,status='delete')
     end if
   end do  

!  Frontal files only

   call delfrt()

end
