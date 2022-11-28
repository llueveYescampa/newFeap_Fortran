subroutine pangl(ix,nen,angl,angg,nrot)
implicit none
  integer :: nen,ix(nen),nrot
  double precision :: angl(nen),angg(*)

!  Purpose: Set up table of rotation angles

!  Inputs:
!     ix(*)     - Element nodal connections
!     nen       - Number of nodes on element
!     angg(*)   - Angle for sloping b.c. of each node point

!  Outputs:
!     angl(*)   - Angle for each element node
!     nrot      - Number of nodes requiring rotation


   integer   ii,n

   nrot = 0
   do n = 1,nen
     angl(n) = 0.0
     ii = abs(ix(n))
     if (ii.gt.0) then
       angl(n) = angg(ii)
       if (angg(ii).ne.0.0) then 
         nrot = nrot + 1
       end if  
     end if
   end do  

end
