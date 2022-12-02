subroutine  elmt07(d,u,x,     s,p,ndf,ndm,nst,isw)
implicit none
  integer          :: ndf,ndm,nst,isw
  double precision :: d(*),u(ndf,*),x(ndm,*),s(nst,*),p(*)
      
! Purpose: Any dimensional truss element routine
! Inputs:
!    d(*)      - Material set parameters
!    u(ndf,*) - Solution variables for element
!    x(ndm,*) - Nodal coordinates for element
!    ndf       - Number dof/node
!    ndm       - Spatial dimension of FEM mesh
!    nst       - Size of element matrix
!    isw       - Switch parameter to control action taken by routine
!    dx()       - direction cosines
! Outputs:
!    s(nst,*)  - Element matrix
!    p(*)      - Element vector

  integer          :: i,j,k
  double precision :: t(9),up(6), temp(9),temp1(9)
  double precision :: kp(9),axial,sm
  double precision :: dx(4),xx(3),xl,eps
  double precision :: t1,t2,t3=0.0d0,sig  
  character        :: yyy*80
  
  integer iocheck
  logical test

  include 'cdata.h'
  include 'eldata.h'
  include 'hdata.h'  
  include 'iofile.h'
  
  xl  = 0.0
  eps = 0.0
  do i = 1,ndm
    dx(i) = x(i,2) - x(i,1)
    xl = xl + dx(i)*dx(i)
    xx(i) = (x(i,2) + x(i,1))*0.5
    eps = eps + dx(i)*(u(i,2) - u(i,1))
  end do 
  eps = eps/xl
  xl = sqrt(xl)
  
  select case (isw)
  case (1)
!   input material properties
    test = .true.
    do while (test)
      if(ior.lt.0) then 
        write(*,3000) 
      end if  
      call pintio(yyy,10)
      read(yyy,'(8f10.0)',IOSTAT=iocheck) (d(i),i=1,6)
      if (iocheck .eq. 0) then
        d( 7) = d(1)*d(2)
        d( 8) = d(4)*d(2)
        d( 9) = d(5)*d(2)
        d(10) = d(6)*d(2)
        d(11) = d(3)*d(2)
        call pconsd(xx,3,0.0d0)
        if(d(4).gt.0.0d0) then
          write(iow,2000) (d(i),i=1,6)
          if(ior.lt.0) then 
            write(*,2000) (d(i),i=1,6)
          end if  
        else
          write(iow,2001) (d(i),i=1,3)
          if(ior.lt.0) then 
            write(*,2001) (d(i),i=1,3)
          end if  
        end if
        nh1 = 3
        test = .false.
        return
      else 
          call myPerror('elmt07',yyy)
      end if
    end do
  case (2,7,8)
    return
  case (3,6)
!   compute element arrays
    if (ndm .eq. 3) then
      dx(4) = sqrt(dx(1)*dx(1) + dx(3)*dx(3))/xl
    end if
    do  i = 1,ndm
      dx(i) = dx(i)/xl
    end do
    do i=1,2*ndm
      up(i)=0.0
    end do  
    if (ndm .eq. 2) then
!     Rotation matrix for 2D case
      t(1) =  dx(1)
      t(2) = -dx(2)
      t(3) =  dx(2)
      t(4) =  dx(1)
    else
      if(dx(4) .ne. 0.0) then
!     Rotation matrix for 3D case
        t(1) = dx(1)
        t(2)= -dx(1)*dx(2)/dx(4)
        t(3)= -dx(3)/dx(4)
        t(4) = dx(2)
        t(5) = dx(4)
        t(6) = 0.0
        t(7) = dx(3)
        t(8)= -dx(2)*dx(3)/dx(4)
        t(9) = dx(1)/dx(4)
      else
!     Rotation matrix for 3D case
!     case of vertical member
        t(1) = 0.0
        t(2)= -dx(2)
        t(3)= 0.0
        t(4) = dx(2)
        t(5) = 0.0
        t(6) = 0.0
        t(7) = 0.0
        t(8)= 0.0
        t(9) = 1
      end if
    end if
    call matmult(t, u(1,1), up(1), ndm, ndm, 1,0)
    call matmult(t, u(1,2), up(ndm+1), ndm, ndm, 1,0)
    if (ndm .eq. 3) then
      t1 = (up(4)-up(1))/xl
      t2 = (up(5)-up(2))/xl
      t3 = (up(6)-up(3))/xl
      axial = (-up(1)*(2+t1) - up(2)*(t2) - up(3)*(t3) +          &
           up(4)*(2+t1) + up(5)*(t2) + up(6)* (t3))*0.5*d(7)/xl
      temp1(1)=(1+t1)*axial
      temp1(2)=t2*axial
      temp1(3)=t3*axial
    else
      t1 = (up(3)-up(1))/xl
      t2 = (up(4)-up(2))/xl
      axial = (-up(1)*(2+t1) - up(2)*(t2) +                      &
               up(3)*(2+t1) + up(4)*(t2))*0.5*d(7)/xl
      temp1(1)=(1+t1)*axial
      temp1(2)=t2*axial
    end if
    call matmult(t, temp1, p, ndf, ndf, 1,1) 
    do i=1,ndf
      p(i+ndf) = -p(i)
    end do
    
    if (isw .eq. 3) then
      if (ndm .eq. 3) then
        kp(1) = (1+t1)*(1+t1)*d(7)/xl
        kp(2) = (1+t1)*t2*d(7)/xl
        kp(3) = (1+t1)*t3*d(7)/xl
        kp(4) = (1+t1)*t2*d(7)/xl
        kp(5) = t2*t2*d(7)/xl
        kp(6) = t2*t3*d(7)/xl
        kp(7) = (1+t1)*t3*d(7)/xl
        kp(8) = t2*t3*d(7)/xl
        kp(9) = t3*t3*d(7)/xl
!       including the geometric stiffness matrix 
        kp(1) = kp(1)+ axial/xl
        kp(5) = kp(5)+ axial/xl
        kp(9) = kp(9)+ axial/xl
      else
        kp(1) = (1+t1)*(1+t1)*d(7)/xl
        kp(2) = (1+t1)*t2*d(7)/xl
        kp(3) = (1+t1)*t2*d(7)/xl
        kp(4) = t2*t2*d(7)/xl
!       including the geometric stiffness matrix 
        kp(1) = kp(1)+ axial/xl
        kp(4) = kp(4)+ axial/xl
      end if
      call matmult(t, kp, temp, ndm, ndm, ndm,1)
      call matmult (temp,t, temp1, ndm, ndm, ndm,0)
      k=0
      do j=1,ndm
        do i=1,ndm
          k = k + 1
          s(i,j)=temp1(k)
          s(i+ndm,j+ndm) = temp1(k)
          s(i+ndm,j) =-temp1(k)
          s(i,j+ndm) =-temp1(k)
        end do
      end do  
    end if
!   output stress and strain in element
  case (4)
      mct = mct - 1
      if(mct.le.0) then
        call prthed(iow)
        write(iow,2002) 
        if(ior.lt.0) then 
          write(*,2002)
        end if  
        mct = 50
      end if
      sig = eps*d(7)
      write(iow,'(2i5,3f11.4,2e13.5)') n,ma,xx,sig,eps
      if(ior.lt.0) then 
        write(*,'(2i5,3f11.4,2e13.5)') n,ma,xx,sig,eps
      end if  
!   compute element lumped mass matrix
  case (5)
    sm = d(11)*xl*0.5
    do i = 1,ndm
      p(i)     = sm
      p(i+ndf) = sm
    end do  
  end select
   
!.... formats
2000 format(5x,'T r u s s    E l e m e n t 7'//             &
        10x,'Modulus  =',e12.5/10x,'Area     =',e12.5/      &
        10x,'Density  =',e12.5/10x,'Yield    =',e12.5/      &
        10x,'Iso. Hard=',e12.5/10x,'Kin. Hard=',e12.5/)
2001 format(5x,'T r u s s    E l e m e n t 7'//             &
        10x,'Modulus  =',e12.5/10x,'Area     =',e12.5/      &
        10x,'Density  =',e12.5/)
2002 format(5x,'T r u s s    E l e m e n t 7'//' elem mate',          &
       4x,'1-coord',4x,'2-coord',4x,'3-coord',5x,'force',7x,'strain')
3000 format(' Input  El/Plas: E, A, rho, Y, H-iso, H-Kin'/           &
            '        Elastic: E, A, rho'/3x,'>',$)
end
