subroutine elmt04(d,u,x,ix,t,s,p,ndf,ndm,nst,isw)
implicit none
integer ndf,ndm,nst,isw, ix(*)
double precision  d(*),u(ndf,*),x(ndm,*),t(*),s(nst,*),p(*)
      
!     Purpose: Any dimensional truss element routine

!     Inputs:
!        d(*)      - Material set parameters
!        ul(ndf,*) - Solution variables for element
!        xl(ndf,*) - Nodal coordinates for element
!        ix(*)     - Nodal connections for element
!        tl(*)     - Nodal temperatures for element
!        ndf       - Number dof/node
!        ndm       - Spatial dimension of FEM mesh
!        nst       - Size of element matrix
!        isw       - Switch parameter to control action taken by routine

!     Outputs:
!        s(nst,*)  - Element matrix
!        p(*)      - Element vector

      

      integer         i,i1,ii,j,j1,jj
      double precision  db(3),dx(3),xx(3),xl,eps,sig,ad,sm
      !character       yyy*80
      integer iocheck
      logical test

      include 'cdata.h'

      include 'ydata.h'

      include 'eldata.h'
      
      include 'hdata.h'

      include 'iofile.h'

      select case (isw)
      case (1)
!       Input material properties
        test = .true.
        do while (test)
          if(ioRead.lt.0) then 
             write(*,3000) 
          end if   
          call pintio(yyy,10)
          read(yyy,1000,IOSTAT=iocheck) (d(i),i=1,6)
          if (iocheck .eq. 0) then
            d( 7) = d(1)*d(2)
            d( 8) = d(4)*d(2)
            d( 9) = d(5)*d(2)
            d(10) = d(6)*d(2)
            d(11) = d(3)*d(2)
            call pconsd(xx,3,0.0d0)
            if(d(4).gt.0.0d0) then
              write(ioWrite,2000) (d(i),i=1,6)
              if(ioRead.lt.0) write(*,2000) (d(i),i=1,6)
            else
              write(ioWrite,2001) (d(i),i=1,3)
              if(ioRead.lt.0) write(*,2001) (d(i),i=1,3)
            end if
            nh1 = 3
            test = .false.
            return
          else   
            call pperror('PCELM4',yyy)
          end if 
        end do  
      case (3,4,5,6)
!     Compute element arrays
        xl  = 0.0
        eps = 0.0
        do i = 1,ndm
          dx(i) = x(i,2) - x(i,1)
          xl = xl + dx(i)**2
          eps = eps + dx(i)*(u(i,2)-u(i,1))
          xx(i) = (x(i,2) + x(i,1))/2.
        end do ! i
        eps = eps/xl
        if(mod(isw,3).eq.0) then
          call modl04(d,eps, sig,ad)
!         Form a residual
          sig = sig/sqrt(xl)
          do i = 1,ndf
            p(i) = dx(i)*sig
            p(i+ndf) = -p(i)
          end do ! i
!         Compute tangent stiffness
          if(isw.eq.3) then
            xl = xl*sqrt(xl)
            do i = 1,ndm
              db(i) = ad*dx(i)
              dx(i) = dx(i)/xl
            end do ! i
            i1 = 0
            do ii = 1,2
              j1 = i1
              do jj = ii,2
                do i = 1,ndm
                  do j = 1,ndm
                    s(i+i1,j+j1) = db(i)*dx(j)
                  end do ! j
                end do ! i
                j1 = j1 + ndf
                do j = 1,ndm
                  dx(j) = -dx(j)
                end do ! j
              end do ! jj
              i1 = i1 + ndf
            end do ! ii
            do i = 1,ndm
              do j = 1,ndm
                s(i+ndf,j) = s(j,i+ndf)
              end do ! j
            end do ! i
          end if
!         Output stress and strain in element
        else if(isw.eq.4) then
          call modl04(d,eps, sig,ad)
          mct = mct - 1
          if(mct.le.0) then
            call prthed(ioWrite)
            write(ioWrite,2002) 
            if(ioRead.lt.0) write(*,2002)
            mct = 50
          end if
          write(ioWrite,2003) n,ma,xx,sig,eps
          if(ioRead.lt.0) then 
            write(*,2003) n,ma,xx,sig,eps
          end if  
!         Compute element lumped mass matrix
        else if(isw.eq.5) then
          sm = d(11)*sqrt(xl)/2.0d0
          do i = 1,ndm
            p(i    ) = sm
            p(i+ndf) = sm
          end do ! i
        end if
      case (2,7,8)
!        Check for errors
         return
      end select

!     Formats

1000  format(8f10.0)

2000  format(5x,'T r u s s    E l e m e n t 4'//         &
         10x,'Modulus  =',e12.5/10x,'Area     =',e12.5/  &
         10x,'Density  =',e12.5/10x,'Yield    =',e12.5/  &
         10x,'Iso. Hard=',e12.5/10x,'Kin. Hard=',e12.5/)

2001  format(5x,'T r u s s    E l e m e n t 4'//         &
         10x,'Modulus  =',e12.5/10x,'Area     =',e12.5/  &
         10x,'Density  =',e12.5/)

2002  format(5x,'T r u s s    E l e m e n t 4'//' elem mate',          &
        4x,'1-coord',4x,'2-coord',4x,'3-coord',5x,'force',7x,'strain')

2003  format(2i5,3f11.4,2e13.5)

3000  format(' Input  El/Plas: E, A, rho, Y, H-iso, H-Kin'/   &
             '        Elastic: E, A, rho'/3x,'>',$)
     

end
