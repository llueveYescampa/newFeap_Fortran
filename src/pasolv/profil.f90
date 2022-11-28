subroutine profil (jd,idl,id,ix,ndf,nen1)
implicit none
  integer :: ndf,nen1
  integer :: jd(*),idl(*),id(*),ix(nen1,*)

!  Purpose: Compute profile of global arrays

!  Inputs:
!     idl(ndf,*) - Temporary storage
!     id(*)      - Boundary condition codes
!     ix(nen1,*) - Element nodal connections
!     ndf        - Number dof/node
!     nen1       - Dimension of ix array

!  Outputs:
!     id(ndf,*)  - Equation numbers for each dof
!     jd(*)      - Profile for equations

   integer  i,ii, j,jj, mm, n,nad,nneq


   include 'cdata.h'

   include 'frdata.h'

   include 'iofile.h'

!  Set up equation numbers

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

!  Compute column heights

   call pconsi(jd,neq,0)
   do n = 1,numel
     mm  = 0
     nad = 0
     do i = 1,nen
       ii = iabs(ix(i,n))
       if(ii.gt.0) then
         do j = 1,ndf
           jj = id((ii-1)*ndf+j)
           if(jj.gt.0) then
             if(mm.eq.0) mm = jj
             mm = min(mm,jj)
             nad = nad + 1
             idl(nad) = jj
           end if
         end do  
       end if
     end do  

     if(nad.gt.0) then
       do i = 1,nad
         ii = idl(i)
         jj = jd(ii)
         jd(ii) = max(jj,ii-mm)
       end do  
     end if
   end do  

!  Compute diagonal pointers for profile

   nad = 0
   jd(1) = 0
   if(neq.gt.1) then
     do n = 2,neq
       jd(n) = jd(n) + jd(n-1)
     end do  
     nad = jd(neq)
   end if

!  Set element search order to sequential

   do n = 1,numel
     idl(n) = n
   end do  

!  Equation summary

   maxf = 0
   mm   = 0
   if(neq.gt.0) then
     mm = (nad+neq)/neq
   end if  
   write(iow,2001) neq,numnp,mm,numel,nad,nummat
   if(ior.lt.0) then
     write(*,2001) neq,numnp,mm,numel,nad,nummat
   end if  

2001  format(/'   E q u a t i o n    P r o b l e m   S u m m a r y:'//  &
         5x,'Number of equations  =',i7,5x,'Number nodes      =',i5/    &
         5x,'Average col. height  =',i7,5x,'Number elements   =',i5/    &
         5x,'No. terms in profile =',i7,5x,'Number materials  =',i5/)

end
