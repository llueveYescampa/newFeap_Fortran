subroutine phelp(wd,nwd,wrd,isw)
implicit none
  integer   :: nwd,isw
  character :: wd(nwd)*4,wrd*5

!  Purpose: Help file for macro command list

!  Inputs:
!     wd(*)     - List of FEAP command names
!     nwd       - Number of command names
!     wrd       - Request name for help
!     isw       - Switch: =1 for solution commands; otherwise mesh.


   include 'iofile.h'

   if(ior.gt.0) then
     return
   end if  
   if(isw.eq.1) then 
     write(*,2000)
   end if  
   if(isw.ne.1) then 
     write(*,2001) wrd
   end if
   write(*,'(/8(3x,a4))') wd
   write(*,'(/a,a5,a)') ' Terminate ', wrd,' execution with an !end! command.'

2000  format(//' The following macro commands are available'//        &
              ' use of loop must terminate with a matching next'//    &
              ' multiple loop-next pairs may occur up to depth of 8')
2001  format(//' The following ',a5,'commands are available:')

end
