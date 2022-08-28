      include 'fgraph.fi'
      integer*2 function vopnwk()
      implicit integer*2 (a-z)

      include 'fgraph.fd'

      record    /videoconfig/ myscreen

      integer            jfill,jplot
      logical                        lfil
      common    /instl2/ jfill,jplot,lfil

      common    /vgraph/ idxl,idyl

      save


!     Open workstation, home cursor, set up scaling

      call     getvideoconfig(myscreen)
      status = setvideomode ( $MAXRESMODE )
      call     getvideoconfig( myscreen )
      ixln   = myscreen.numxpixels - 1
      iyln   = myscreen.numypixels - 1
      idxl   = 32640.0/(ixln+1) + 0.5
      idyl   = 22480.0/(iyln+1) + 0.5
      if(myscreen.numcolors .le. 4) then
        jfill = 1
      else
        jfill = 2
      end if
      call     clearscreen( $GCLEARSCREEN )
      vopnwk = 0

      end

