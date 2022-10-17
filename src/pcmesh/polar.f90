subroutine polar(x,ndm,prt)
implicit none
integer  ndm
double precision   x(ndm,*)
logical  prt

!     Purpose: Convert polar to cartesian coordinates

!     Inputs:
!        x(ndm,*)   - Polar form of nodal coordinates
!        ndm        - Spatial dimension of mesh
!        prt        - Flag, output results if true

!     Outputs:
!        x(ndm,*)   - Cartesian form of nodal coordinates


      !character yyy*80
      integer   i,inc, mct, n,ni,ne
      integer iocheck
      double precision  x0,y0, r,th

      include 'cdata.h'
      include 'ydata.h'
      include 'iofile.h'

      if(ndm.eq.1) then 
        return
      end if

      mct = 0
      th = atan(1.0d0)/45.0d0
      do while (.true.)
        if(ioRead.lt.0) then
          write(*,'(a/3x,a,$)') &
            ' Input: 1-node, 2-node, inc, x1(cent), x2(cent)','>'
        end if  
        call pintio(yyy,10)
        read(yyy,'(3i10,2f10.0)',IOSTAT=iocheck) ni,ne,inc,x0,y0
        if(iocheck .ne. 0) then
          call pperror('POLAR ',yyy)
          cycle
        else
          if(ni.le.0) then
            return
          end if  
          if(ni.gt.numnp.or.ne.gt.numnp) then
!           Error
            write(ioWrite,'(a,i6,a,i6)') &
                ' **ERROR** attempt to convert nodes ni= ',ni,' - ne= ',ne
            if(ioRead.lt.0) then
              write(*,'(a,i6,a,i6)') &
                ' **ERROR** attempt to convert nodes ni= ',ni,' - ne= ',ne
            end if  
            stop
          end if  
          inc = sign(max(abs(inc),1),ne-ni)
          if(ne.eq.0) then 
            ne = ni
          end if  
          n = ni
          do while (.true.)
            r = x(1,n)
            x(1,n) = x0 + r*cos(x(2,n)*th)
            x(2,n) = y0 + r*sin(x(2,n)*th)
            if(mct.le.0) then
              if(prt) then
                call prthed(ioWrite)
                write(ioWrite,2000) x0,y0,(i,i=1,ndm)
              end if  
              if(ioRead.lt.0.and.prt) then
                write(*,2000) x0,y0,(i,i=1,ndm)
              end if  
              mct = 50
            end if
            if(prt) then 
              write(ioWrite,'(i10,6f13.4)') n,(x(i,n),i=1,ndm)
            end if  
            if(ioRead.lt.0.and.prt) then 
              write(*,'(i10,6f13.4)') n,(x(i,n),i=1,ndm)
            end if  
            mct = mct - 1
            n = n + inc
            if((ne-n)*inc.ge.0) then
              cycle
            end if  
            if(mod(ne-ni,inc).eq.0) then
              exit
            end if  
            ni = ne
            n = ne
          end do
        end if
      end do

!     Formats

2000  format('    P o l a r   t o   C a r t e s i a n   C o o r d s.'/   &
       8x,'Center: x0 = ',e12.4,' y0 = ',e12.4/6x,'node',6(i7,'-coord'))

end
