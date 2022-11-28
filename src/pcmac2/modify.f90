subroutine modify(b,s,dul,nst)
implicit none
  integer          :: nst
  double precision :: b(*),s(nst,*),dul(*)

!  Purpose: Modify for non-zero displacement boundary conditions

!  Inputs:
!     s(nst,*) - Element tangent matrix
!     dul(*)   - Element vector of specified displacements
!     nst      - Dimension of element vectors
!  Outputs:
!     b(*)     - Residual modified for specified displacements

   integer  i,j

   do i = 1,nst
     do j = 1,nst
       b(i) = b(i) - s(i,j)*dul(j)
     end do
   end do  

end
