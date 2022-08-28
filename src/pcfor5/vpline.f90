      include 'fgraph.fi'
      integer*2 function vpline(npt,ixy)
      implicit integer*2 (a-z)

      include 'fgraph.fd'

      integer*2 ixy(2,*), ir1,ir2,irc
      integer   npt
      
      record /xycoord/ xy
      common /vgraph/ idxl,idyl

!     Draw line

      ir1 = idxl
      ir2 = idyl
      irc = 22200
      call moveto( int2(ixy(1,1)/ir1), int2((irc - ixy (2,1))/ir2), xy )
      do n = 2,npt
        vpline = lineto( int2(ixy(1,n)/ir1),int2((irc - ixy(2,n))/ir2) )
      end do  

      end

