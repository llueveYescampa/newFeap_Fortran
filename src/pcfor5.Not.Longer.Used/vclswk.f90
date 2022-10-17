      include 'fgraph.fi'
      integer*2 function vclswk()

      implicit integer*2 (a-z)

      include 'fgraph.fd'

      record /rccoord / s

!     Home cursor - text mode

      call     settextposition( 1 , 1 , s )
      vclswk = settextcursor ( #0607 )
      vclswk = setvideomode  ( $TEXTC80 )
      vclswk = displaycursor ( $GCURSORON )

      end

