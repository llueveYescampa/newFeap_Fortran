      include 'fgraph.fi'
      
      integer*2 function vclrwk()
      implicit integer*2 (a-z)

      include 'fgraph.fd'

!     Clear workstation

      call     clearscreen( $GCLEARSCREEN )
      vclrwk = 0

      end

