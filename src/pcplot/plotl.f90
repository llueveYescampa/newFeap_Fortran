subroutine plotl(x1,x2,x3,ipen)
implicit  none
double precision    x1,x2,x3
integer   ipen

! Line drawing command

! Inputs:
!    x1        - X-mesh coordinate to plot
!    x2        - Y-mesh coordinate to plot
!    x3        - Z-mesh coordinate to plot
!    ipen      - Pen action: = 2 draw; = 3 move

! Outputs:
!    none

  double precision s1,s2
  include 'pdata1.h'

! Compute the normal coordinates

  s1 = x3
  s1 = max(0.0d0,min(1.45d0,myScale*(x1 + x1 - sx(1)) + 0.725d0))
  s2 = max(0.0d0,min(1.00d0,myScale*(x2 + x2 - sx(2)) + 0.500d0))
  call dplot(s1,s2,ipen)
  
end
