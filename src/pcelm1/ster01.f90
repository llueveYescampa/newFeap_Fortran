subroutine ster01(ix,d,xl,ul,tl,shp,st,ndf,ndm,numnp,numel, &
                 sg,tg,wg,sig,eps,lint,ityp)
implicit  none
integer     ix(*),numnp,numel,         ndf,ndm,lint,ityp
double precision    d(*),xl(ndm,*),ul(ndf,*),tl(*),shp(3,4), &
                   st(numnp,*),sg(16),tg(16),wg(16),sig(6),eps(4)

!  Purpose: Compute error estimator values

!  Inputs
!     ix(*)     - Nodes connected to element
!     d(*)      - Material set parameters
!     xl(ndf,*) - Nodal coordinates for element
!     ul(ndf,*) - Solution variables for element
!     ix(*)     - Nodal connections for element
!     tl(*)     - Nodal temperatures for element
!     shp(3,*)  - Shape function and derivatives
!     st(*,*)   - Projected nodal stresses
!     ndf       - Number dof/node
!     ndm       - Spatial dimension of FEM mesh
!     nel       - Number of nodes on element
!     numnp     - Number of nodes in mesh
!     numel     - Number of elements in mesh
!     sg(*)     - Quadrature points
!     tg(*)     - Quadrature points
!     wg(*)     - quadrature weights
!     sig(4)    - Stresses at quadrature points
!     eps(4)    - Stresses at quadrature points
!     lint      - number of quadrature points
!     ityp      - Problem type

!  Outputs:
!     none      - done through common blocks


   integer   i,ii, j,ll
   double precision  gr,xsj, psi,psis, xx,yy
   double precision  sigp(4),dsig(4)


   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   double precision          dm
   integer            n,ma,mct,iel,nel
   common /eldata/ dm,n,ma,mct,iel,nel

   double precision          eerror,elproj,ecproj,efem,enerr,ebar
   common /errind/ eerror,elproj,ecproj,efem,enerr,ebar

!  Stress error computations

   psis = 0.0
   psi  = 0.0
   gr   = (1.+d(17))/d(16)
   do ii = 1,lint
     call shapeFunc(sg(ii),tg(ii),xl,shp,xsj,ndm,nel,ix,.false.)
     call stre01(d,xl,ul,tl,shp,eps,sig,xx,yy,ndm,ndf,nel,ityp)
     xsj = xsj*wg(ii)
     do i = 1,4
       sigp(i)= 0.0d0
     end do  
     do i = 1,nel
       ll = abs(ix(i))
       if(ll.ne.0) then
         do j = 1,4
           sigp(j) = sigp(j) + shp(3,i)*st(ll,j)
         end do  
       end if
     end do  

!    Compute the integral of the stress squares for error indicator use

     do i = 1,4
       dsig(i) = sigp(i)-sig(i)
       efem    = efem   + sig(i)*sig(i)*xsj
       ecproj  = ecproj + sigp(i)*sigp(i)*xsj
       psis    = psis   + (dsig(i)**2)*xsj
       psi     = psi    +  gr*(dsig(i)**2)*xsj
     end do  
     psi = psi +  gr*dsig(2)**2*xsj +                 &
          ((dsig(1)+dsig(3)+dsig(4))**2)*d(18)*xsj
   end do

   eerror = eerror + psis
   if(elproj.ne.0.0d0) then
     psi    = sqrt(abs(psi))/ebar
     psis   = 20.0*sqrt(abs(psis)/elproj*numel)
     if(mct.eq.0) then
       write(ioWrite,2000)
       if(ioRead.lt.0) write(*,2000)
       mct = 50
     end if
     mct = mct - 1
     write(ioWrite,'(i8,1p2e12.4)') n,psis,psi
     if(ioRead.lt.0) then 
       write(*,'(i8,1p2e12.4)') n,psis,psi
     end if  
   end if

2000  format('   M e s h   R e f i n e m e n t s   f o r   5%', &
             '   E r r o r'//'    elmt   h-sigma    h-energy'/)

end
