subroutine elmt05(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
implicit none
integer                   ix(*),    ndf,ndm,nst,isw
double precision  d(*),ul(ndf,*),xl(ndm,*),tl(*),s(nst,*),p(nst)

!  Purpose: Geometrical nonlinear axisymmetric shell: (c) w.wagner

!  Inputs:
!     d(*)      - Material set parameters
!     ul(ndf,*) - Solution variables for element
!     xl(ndf,*) - Nodal coordinates for element
!     ix(*)     - Nodal connections for element
!     tl(*)     - Nodal temperatures for element
!     ndf       - Number dof/node
!     ndm       - Spatial dimension of FEM mesh
!     nst       - Size of element matrix
!     isw       - Switch parameter to control action taken by routine

!  Outputs:
!     s(nst,*)  - Element matrix
!     p(*)      - Element vector



   logical     bs
   integer     i,i1,ii,j,j1,jj,k,lin
   double precision  sn,cs,sl,f,r,z,recr,dv,wks
   double precision  shp(2),dd(5,5),btd(5,3)
   double precision  bm(5,3),sig(5),vl(3,2),dot

   include 'bdata.h'

   include 'cdata.h'

   include 'eldata.h'

   include 'iofile.h'

   ! include 'ydata.h'

!  Input material properties

   if(isw.eq.1) then
     call matl05(d)
     call pconsd(dd,25,0.0d0)
     call pconsd(bm,15,0.0d0)
   else

!  Length, angle, radius, jacobian

     sn = xl(2,2)-xl(2,1)
     cs = xl(1,2)-xl(1,1)
     sl = sqrt(cs*cs + sn*sn)
     bs = d(13).eq.0.0d0

!    Check mesh for error

     if(isw.eq.2) then
       if(sl.le.0.0d0) then
         write(ioWrite,3001) n
         if(ioRead.lt.0) write(*,3001) n
       end if

!    Form lumped mass matrix; length, volume*density

     else if(isw.eq.5) then
       dv    = d(11)*sl/8.0d0
       if(bs) then
         p(1)     = 4.0d0*dv
         p(ndf+1) = p(1)
       else
         p(1)     = dv*(3.d0*xl(1,1) + xl(1,2))
         p(ndf+1) = dv*(3.d0*xl(1,2) + xl(1,1))
       end if
       p(2)     = p(1)
       p(3)     = p(1)*d(12)
       p(ndf+2) = p(ndf+1)
       p(ndf+3) = p(ndf+1)*d(12)

!    Form stiffness/residual

     else
       sn   = sn/sl
       cs   = cs/sl
       if(bs) then
         r    = 1.0d0
         recr = 0.0d0
       else
         r    = 0.5d0*(xl(1,1) + xl(1,2))
         recr = 1.0d0/r
       end if
       dv   = sl*r

!      Shape function derivatives

       shp(2) =  1.0d0/sl
       shp(1) = -shp(2)

!      Local  displacements

       do k = 1,2
         vl(1,k) = cs*ul(1,k) + sn*ul(2,k)
         vl(2,k) =-sn*ul(1,k) + cs*ul(2,k)
         vl(3,k) = ul(3,k)
       end do ! k

!      Derivative w,s

       lin = int(d(5))
       if(lin.eq.0) then
         wks = 0.0d0
       else
         wks = (vl(2,2) - vl(2,1))/sl
       end if

!      Stresses, strains, D-matrix

       call modl05(sig,vl,dd,d,sn,cs,sl,recr,wks)
       if(mod(isw,3).eq.0) then

!        Load vector in local coordinates (reference system)

         i = ndf + 2
         f = d(4)*sl/8.0d0*mydm
         if(bs) then
           p(2) = f*4.0d0
           p(i) = p(2)
         else
           p(2) = f*(3.d0*xl(1,1) + xl(1,2))
           p(i) = f*(3.d0*xl(1,2) + xl(1,1))
         end if

!        K-sigma tangent matrix

         if(lin.ne.0) then
           s(2,2) =   r*sig(1)/sl
           s(i,i) =   s(2,2)
           s(2,i) = - s(2,2)
         end if

!        Multiply stress and moduli by jacobian

         do k = 1,5
           sig(k) = sig(k)*dv
           do j = 1,5
             dd(j,k) = dd(j,k)*dv
           end do ! j
         end do ! k
         i1=0
         do ii = 1,2

!          Residual G = P - Bt*S 

           call bmat05(bm,shp(ii),sn,cs,recr,wks)
           do i = 1,3
             p(i1+i) = p(i1+i) - dot(bm(1,i),sig,5)
           end do ! i

!          Tangent stiffness matrix

           if(isw.eq.3) then
             do i = 1,3
               do k = 1,5
                 btd(k,i) = dot(bm(1,i),dd(1,k),5)
                 end do ! k
               end do ! i
             j1 = i1
             do jj = ii,2
               call bmat05(bm,shp(jj),sn,cs,recr,wks)
               do i = 1,3
                 do j = 1,3
                   s(i1+i,j1+j)=s(i1+i,j1+j)+dot(btd(1,i),bm(1,j),5)
                 end do ! j
               end do ! i
               j1 = j1 + ndf
             end do ! jj
           end if
           i1 = i1 + ndf
         end do ! ii

!        Lower part stiffness matrix and transform to global frame

         if(isw.eq.3) then
           do i = 1,3
             do j = 1,3
               s(i+ndf,j) = s(j,i+ndf)
             end do ! j
           end do ! i
           call tran05(s,cs,sn,nst,ndf,1)
         end if
         call tran05(p,cs,sn,nst,ndf,2)

!      Output stresses (N, M, Q)

       else if(isw.eq.4) then
         mct = mct - 1
         if(mct.le.0) then
            write(ioWrite,2001) head
            if(ioRead.lt.0) write(*,2001) head
            mct = 50
         end if
         r = 0.5*(xl(1,1) + xl(1,2))
         z = 0.5*(xl(2,1) + xl(2,2))
         write(ioWrite,'(i5,i4,0p2f8.3,5(1x,1p1e10.3))') n,ma,r,z,sig
         if(ioRead.lt.0) then 
           write(*,'(i5,i4,0p2f8.3,5(1x,1p1e10.3))') n,ma,r,z,sig
         end if  
       end if
     end if
   end if

!  Format statements

2001  format(1x,20a4//2x,'E L E M E N T   S T R E S S E S'//  &
        '  El  Mat  1-Coor  2-Coor   **NS**    **NPHI**'      &
        '    **MS**    **MPHI**    **QS**'/1x)
     
3001  format(' *ERROR* Element',i5,' has zero length.')

end
