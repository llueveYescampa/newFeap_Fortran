subroutine pdefm(x,b,c,ndm,ndf,numnp, dr)
implicit  none
integer ndf,ndm,numnp
double precision x(ndm,*),b(ndf,*),c,dr(ndm,*)

! Purpose: Compute deformed position of two-dimensional meshes

! Inputs:
!    x(ndm,*)  - Nodal coordinates of mesh
!    b(ndf,*)  - Solution vector
!    c         - Scale factor for added displacements
!    ndm       - Spatial dimension of mesh
!    ndf       - Number dof/node
!    numnp     - Number of nodes in mesh

! Outputs:
!    dr(ndf,*) - Deformed nodal coordinates

  integer   i, n
  
  call pconsd(dr,ndm*numnp,0.0d0)
  do n = 1,numnp
    if(x(1,n) .ne. -999.) then
      do i = 1,ndm
        dr(i,n) = x(i,n) + c*b(i,n)
      end do  
    end if
  end do  
end
