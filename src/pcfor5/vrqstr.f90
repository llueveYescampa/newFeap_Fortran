!integer*2 function vrqstr(i1,i2,ixy,xx)  !.... input character to quit plot
integer*2 function vrqstr(xx)  !.... input character to quit plot
!      integer   i1,i2
!      integer*2 ixy
      character*1 xx
      
      read(*,'(a1)') xx
      vrqstr = 0
end
