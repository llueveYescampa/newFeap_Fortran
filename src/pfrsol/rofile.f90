subroutine profil (jd,idl,id,ix,ndf,nen1)
implicit  none
integer   jd(*),idl(*),id(*),ndf,nen1,ix(nen1,*)

!  Purpose: Compute front profile of global arrays

!  Inputs:
!     jd(*)      - Working array
!     idl(*)     - Element elimination order
!     id(ndf,*)  - Boundary condition array
!     ix(nen1,*) - Element nodal connection list
!     ndf        - Number dof/node
!     nen1       - Dimension of ix array

!  Outputs:
!     id(ndf,*)  - Equation numbers for each dof
!     ix(nen1,*) - Element nodal connection list tagged for frontal
!                  eliminations.

   integer   j, n , nneq
   
   include 'cdata.h'

   include 'iofile.h'

   include 'frdata.h'

!  Set up the equation numbers

   neq = 0
   nneq = ndf*numnp
   do n = 1,nneq
     j = id(n)
     if(j.eq.0) then
       neq = neq + 1
       id(n) = neq
     else
       id(n) = 0
     end if
   end do  

!  Set default element elimination order
!  N.B. Can replace with optimized ordering.

   do n = 1,numel
     idl(n) = n
   end do  

!  Compute front width

   call prefrt(jd,idl,ix,maxf,ndf,nen,nen1,numel,numnp)

!  Max front width must produce a triangle which fits in adata
!  plus some buffer space for reduced equations.

   if(maxf.gt.150) then
     write(*,'(a,i8)') ' *ERROR* Front requires too much storage =',maxf
     if(ioRead.lt.0) then
       write(ioWrite,'(a,i8)') ' *ERROR* Front requires too much storage =',maxf
     end if  
   end if

end
