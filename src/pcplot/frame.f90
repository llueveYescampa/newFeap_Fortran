subroutine frame(x,ndm,numnp)
implicit none
  integer          :: ndm,numnp
  double precision :: x(ndm,*)
      
! Purpose: Compute scaling for plot area

! Inputs:
!    x(ndm,*)  - Nodal coordinates for mesh
!    ndm       - Spatial dimension of mesh
!    numnp     - Number of nodes in mesh

! Outputs:
!    none      - Passed through common blocks


  logical   iflg
  integer   i,ii,ij, n
  double precision  xmn(2),xmx(2),xmin(3),xmax(3),fact,xcen

  include 'pdata1.h'

! Determine window coordinates
  if(ndm.eq.1) then
    dx(2) = 0.0
    sx(2) = 0.0
  end if
  ii = min(ndm,3)
  ij = min(ndm,2)

! Find the minimum and maximum coordinate of input nodes
  iflg = .true.
  do n = 1,numnp
    if(x(1,n).ne. -999.) then
      if(iflg) then
        do i = 1,ii
          xmin(i) = x(i,n)
          xmax(i) = x(i,n)
        end do  
        iflg = .false.
      else
        do i = 1,ii
          xmin(i) = min(xmin(i),x(i,n))
          xmax(i) = max(xmax(i),x(i,n))
        end do  
      end if
    end if
  end do  
  scale  = max(xmax(1)-xmin(1),xmax(2)-xmin(2))

! Plot region determination
  do i = 1,ij
    xmn(i) = min(xmin(i),xmax(i))
    xmx(i) = max(xmin(i),xmax(i))
    dx(i) = xmx(i) - xmn(i)
    sx(i) = xmx(i) + xmn(i)
  end do  

! Rescale window
  if(dx(1).gt.1.45*dx(2)) then
    xmn(2) = (sx(2) - dx(1))*0.5
    xmx(2) = (sx(2) + dx(1))*0.5
    fact   = 0.58
  else
    xmn(1) = (sx(1) - dx(2))*0.5
    xmx(1) = (sx(1) + dx(2))*0.5
    fact   = 0.40
  end if
  do i = 1,ij
    xmin(i) = max(xmin(i),xmn(i)) - scale * 0.01
    xmax(i) = min(xmax(i),xmx(i)) + scale * 0.01
  end do  

! Default values
  scale  = max(xmax(1)-xmin(1),xmax(2)-xmin(2))

! Reset values for deformed plotting
  do i = 1,ij
    xcen = xmax(i)+xmin(i)
    xmax(i) = (xcen + 1.1*scale)*0.5
    xmin(i) = (xcen - 1.1*scale)*0.5
  end do  
  scale = fact/scale

end
