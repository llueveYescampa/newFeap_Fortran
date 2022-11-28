subroutine elmt02(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
implicit none
  integer          ::  ix(*),ndf,ndm,nst,isw
  double precision :: d(*),ul(ndf,*),xl(ndm,*),tl(*),s(nst,*),p(*)
      

!   Purpose: Plane/Axisymmetric Linear Element -- bbar formulation

!   Inputs:
!      d(*)      - Material set parameters
!      ul(ndf,*) - Solution variables for element
!      xl(ndf,*) - Nodal coordinates for element
!      ix(*)     - Nodal connections for element
!      tl(*)     - Nodal temperatures for element
!      ndf       - Number dof/node
!      ndm       - Spatial dimension of FEM mesh
!      nst       - Size of element matrix
!      isw       - Switch parameter to control action taken by routine

!   Outputs:
!      s(nst,*)  - Element matrix
!      p(*)      - Element vector


    logical flg

    integer i,i1,ii,j,j1,jj, ib,ityp,l,lint
    double precision rr,zz,xr0,xz0,type,vol,yld,dot
    double precision eps(4),sig(6),bbar(4,2,4),bbd(4,2),sg(4),tg(4), &
                    wg(4),shp3(4,4),shp(3,4,4),xsj(4),ang,siga(6)
    
    include 'adata.h'
    include 'cdata.h'


   include 'eldata.h'
   include 'elcom2.h'

   include 'iofile.h'


!   Go to correct array processor

    select case (isw)
    case(1)
!     Input material properties

      call matl02(d,ityp,ib)
      d(5) = ityp
      d(6) = ib
      l    = 2
      call pgauss(l,lint,sg,tg,wg)
      return
    case(2)
!     Check for element errors

      call ckisop(ix,xl,shp,ndm)
      return
    case(3,4,6,8)
!     Compute tangent stiffness and residual force vector
      
      type = d(5)
      ib   = d(6)
      
!     Compute volumetric integrals
      
      call pconsd(g,8,0.0d0)
      do l = 1,lint
        call shapeFunc(sg(l),tg(l),xl,shp(1,1,l),xsj(l),ndm,4,ix,.false.)
        call gvc02(shp(1,1,l),shp3(1,l),xsj(l),wg(l),xl,type,ndm)
      end do ! l
      vol  = xsj(1) + xsj(2) + xsj(3) + xsj(4)
      do i = 1,4
        g(1,i) = g(1,i)/vol
        g(2,i) = g(2,i)/vol
      end do ! i
      if(isw.eq.4) then
!       Compute stresses at center of element

        call pconsd(sig,6,0.0d0)
        rr = 0.0d0
        zz = 0.0d0
        do l = 1,lint
        
!         Compute stress, strain, and material moduli
        
          call strn02(shp(1,1,l),xl,ul,type,xr0,xz0,ndm,ndf,eps)
          call modl02(d,ul,eps,siga,1.d0,ndf,ib)
          do i = 1,5
            sig(i) = sig(i) + 0.25d0*siga(i)
          end do ! i
          rr = rr + 0.25d0*xr0
          zz = zz + 0.25d0*xz0
        end do ! l
        yld = sig(5)
        call pstres(sig,sig(5),sig(6),ang)
        
!       Output stresses
        
        mct = mct - 2
        if(mct.le.0) then
          call prthed(iow)
          write(iow,2001)
          if(ior.lt.0) then
            write(*,2001)
          end if  
          mct = 50
        end if
        write(iow,2002) n,ma,sig,rr,zz,yld,ang
        if(ior.lt.0) then
          write(*,2002) n,ma,sig,rr,zz,yld,ang
        end if  
        return
      end if  
      if(isw.eq.8) then 
!     Stress computations for nodes
      
        call stcn02(ix,d,xl,ul,shp,aa,aa(numnp+1),ndf,ndm,numnp,sg,tg,sig,&
                   eps,lint,type,ib)
        return
      end if  
      flg  = isw .eq. 3
      do l = 1,lint
      
!       Compute stress, strain, and material moduli
      
        call strn02(shp(1,1,l),xl,ul,type,xr0,xz0,ndm,ndf,eps)
        call modl02(d,ul,eps,sig,xsj(l),ndf,ib)
        i1 = 0
        do i = 1,4
      
!         Compute the internal stress divergence term
      
          p(i1+1) =                                                            &
          p(i1+1) - shp(1,i,l)*sig(1) - shp(2,i,l)*sig(2) - shp3(i,l)*sig(4)
          p(i1+2) = p(i1+2) - shp(2,i,l)*sig(3) - shp(1,i,l)*sig(2)
      
!         Compute stiffness
      
      
          if(flg) then     ! isw = 3 Nota de edgar
      
!           Compute b-bar matrix
      
            call bmat02(shp3(i,l),shp(1,i,l),g(1,i),bbar(1,1,i),ib)
            do ii = 1,2
              do jj = 1,4
                bbd(jj,ii) = dot(bbar(1,ii,i),ad(1,jj),4)
              end do ! jj
            end do ! ii
            j1 = 0
            do j  = 1,i
              do ii = 1,2
                do jj = 1,2
                  s(ii+i1,jj+j1) = s(ii+i1,jj+j1)                              &
                                   + dot(bbd(1,ii),bbar(1,jj,j),4)
                end do ! jj
              end do ! ii
      
              j1 = j1 + ndf
            end do ! j
          end if
          i1 = i1 + ndf
        end do ! i
      end do ! l
      
!     Form lower part by symmetry
      
      if(flg) then      ! isw = 3 Nota de edgar
        do i = 1,nst
          do j = 1,i
            s(j,i) = s(i,j)
          end do ! j
        end do ! i
      end if
      return
    case(5)
!   Compute lumped mass matrix

      do l = 1,lint
        call shapeFunc(sg(l),tg(l),xl,shp,xsj(1),ndm,4,ix,.false.)
      
!       Compute radius and multiply into jacobian for axisymmetry
      
        if(d(5).ne.0.0d0) then
          rr = 0.0d0
          do j = 1,4
            rr = rr + shp(3,j,1)*xl(1,j)
          end do ! j
          xsj(1) = xsj(1)*rr
        end if
        xsj(1) = wg(l)*xsj(1)*d(4)
      
!       For each node j compute db = rho*shape*dv
      
        j1 = 1
        do j = 1,4
          p(j1  ) = p(j1) + shp(3,j,1)*xsj(1)
          p(j1+1) = p(j1)
          j1 = j1 + ndf
        end do ! j
      end do ! l
      return
    case(7)
!     Dummy call

      return
    end select

!   Formats for input-output


2001  format('  Element Stresses'//'  elmt  matl  11-stress  12-stress', &
       '  22-stress  33-stress   1-stress   2-stress'/'  1-coord',       &
       '  2-coord  ',33x,'yield    angle')

2002  format(2i6,1p6e11.3/0p2f9.3,30x,1p1e11.3,0pf8.2/1x)

end
