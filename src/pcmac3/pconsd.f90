subroutine pconsd(v,nn,cc)
implicit none
integer  nn
double precision  v(nn),cc

! Purpose: Set double precision array to constant value

! Inputs:
!    nn        - Length of array
!    cc        - Value to assign

! Outputs:
!    v(*)      - Array set to 'cc'

  integer   n
  do n = 1,nn
    v(n) = cc
  end do  

  !v(1:nn) = cc

end
