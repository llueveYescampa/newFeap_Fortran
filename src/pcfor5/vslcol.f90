      include 'fgraph.fi'
      integer*2 function vslcol(it)

      implicit integer*2 (a-z)

      include 'fgraph.fd'

      integer*2 ipal(7)
      integer         jfill,jplot
      logical                     lfil
      common /instl2/ jfill,jplot,lfil

!     Set line color

      data ipal/ 15, 4, 2, 1, 14, 3, 5/
      if(it.gt.0 .and. it.le.7 ) then
        icll = ipal(it)
        if(jfill.lt.2) then 
          icll = 1
        end if  
      else
        icll = 0
      end if
      
      vslcol = setcolor( icll )

      end

