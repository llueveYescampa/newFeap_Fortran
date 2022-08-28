subroutine dplot(x,y,ipen)
implicit none

!     Purpose: Plot line routine

!     Inputs:
!        x         - X-normalized coordinate to plot
!        y         - Y-normalized coordinate to plot
!        ipen      - Pen action: = 2 draw; = 3 move

!     Outputs:
!        none


      integer*2  status,vpline
      integer    ipen
      double precision   x,y

      integer*2       ixy
      common /pdata2/ ixy(4)
      
!     Pen command motions  (ipen = 2, pendown)
      ixy(3) = 22000*x
      ixy(4) = 22000*y
      if(ipen .eq.2 ) then
        status = vpline(2,ixy)
      end if
      ixy(1) = ixy(3)
      ixy(2) = ixy(4)

end
