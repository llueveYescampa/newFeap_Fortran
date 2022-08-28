subroutine pconsi(iv,nn,ic)
implicit  none
integer   nn,ic
integer   iv(nn)

!  Purpose: Set integer array to constant value

!  Inputs:
!     nn        - Length of array
!     ic        - Value to assign

!  Outputs:
!     iv(*)     - Array set to 'ic'

   integer n

   do n = 1,nn
     iv(n) = ic
   end do  
end
