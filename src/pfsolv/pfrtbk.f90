subroutine pfrtbk(b,dr,m,ipd,ibuf,maxf,nv,neq,nev,aengy,ifl)
implicit none
  integer          :: m(*),ipd,ibuf,maxf,nv,neq,nev,ifl
  double precision :: b(maxf,*),dr(neq,*),aengy

!  Purpose: Backsubstitution for frontal solution

!  Inputs:
!     b(*)     - Frontal right hand side/solution vector
!     dr(*)    - Forward reduced solution vector
!     m(*)     - Buffer storage array
!     ipd      - Real*8 precision
!     ibuf     - Buffer length
!     maxf     - Maximum front allowed
!     nv       - Number of buffer blocks
!     neq      - Number of equations to solve
!     nev      - Number of right hand sides
!     ifl      - File name number for frontal equation storage on disk

!  Outputs:
!     dr(*)    - Solution incrrement
!     aengy    - Energy of solution increment


   integer  i, jj, k, n,np,nfrt

   include 'nfrta.h'

   call pconsd(b,maxf*nev,0.0d0)
   aengy = 0.0

!  Recover block

   np = npl
   do n = nv,1,-1
     if(nv.gt.1) then
       call pbuff(m,ibuf,np,n,1,ifl)
     end if  
     do while (.true.)
       nfrt = m(np)
       np = np - 8 - nfrt*ipd
       k  = m(np+2)
       jj = m(np+3)
       do i = 1,nev
         call pfrtb(b(1,i),dr(1,i),nfrt,k,jj,m(np+5),aengy)
       end do  
       if(np.gt.0) then
         cycle
       else
         exit
       end if
     end do
   end do  

end
