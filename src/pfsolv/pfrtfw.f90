subroutine pfrtfw(b,dr,m,ipd,ibuf,maxf,nv,neq,nev,ifl)
implicit none
integer m(*),ipd,ibuf,maxf,nv,neq,nev,ifl
double precision  b(maxf,1),dr(neq,1)

!  Purpose: Forward solution for frontal resolutions

!  Inputs:
!     b(*)     - Frontal right hand side/solution vector
!     dr(*)    - Residual vector
!     m(*)     - Buffer storage array
!     ipd      - Double precision
!     ibuf     - Buffer length
!     maxf     - Maximum front allowed
!     nv       - Number of buffer blocks
!     neq      - Number of equations to solve
!     nev      - Number of right hand sides
!     ifl      - File name number for frontal equation storage on disk

!  Outputs:
!     dr(*)    - Forward reduced solution


   integer  i,ilast, jj, k, n,np,nfrt

   call pconsd(b,maxf*nev,0.0d0)

   do n = 1,nv
     call pbuff(m,ibuf,ilast,n,1,ifl)
     np = 1
     do while (.true.)
       nfrt = m(np)
       k    = m(np+1)
       jj   = m(np+2)
       do i = 1,nev
         call pfrtf(b(1,i),dr(1,i),nfrt,k,jj,m(np+4))
       end do  
       
!      Align last word of buffer for backsubstitution reads
       
       np = np + 8 + m(np)*ipd
       
       if(np.lt.ilast) then
         cycle
       else
         exit
       end if
     end do  
   end do  

end
