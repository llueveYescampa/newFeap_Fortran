subroutine pfrtd(a,b,r,ig,lg,solv,eq,nfrt,k)
implicit none
  integer          :: ig(*),lg(*),nfrt,k
  double precision :: a(*),b(*),r,eq(*)
  logical          :: solv

!  Purpose: Triangular decomposition for frontal program

!  Inputs:
!     a(*)      - Frontal matrix before reduction
!     b(*)      - RHS of frontal equ;ations
!     r         - Residual for current equation
!     ig(*)     - Pointer to diagonals in a(*)
!     lg(*)     - Global equation numbers for each frontal equation
!     solv      - Flag, solve equations if true; otherwise just factor

!  Outputs:
!     eq(*)     - Eliminated equation
!     nfrt      - Size of front
!     k         - Location of diagonal in eq(*)


   integer  j,jj, kk,kl,km,kp, l
   double precision pivot,term

   include 'nfrta.h'

!  Extract equation

   do j=1,nfrt
     l = max(j,k)
     l=ig(l)-l+min(j,k)
     eq(j)=a(l)
   end do  

   pivot=eq(k)

   if(pivot.eq.0.0d0) then
     call pconsd(eq,nfrt,0.0d0)
     write(*,'(a/a)')' WARNING -- Zero pivot, check boundary codes.', &
                     '            Pivot set to 1.0 and solution continued.'
     pivot = 1.0
     eq(k) = 1.0
   else
     dimx = max(dimx,abs(pivot))
     if(dimn.eq.0.0d0) then
       dimn = abs(pivot)
     end if  
     dimn = min(dimn,abs(pivot))
   end if
   km = k - 1
   kp = k + 1
   kk = 1
   if(km.gt.0) then
     do jj = 1,km
       if(eq(jj).ne.0.0d0) then
         term = -eq(jj)/pivot
         if(solv) b(jj) = b(jj) + term*r
         call saxpb(eq,a(kk),term,jj,a(kk))
       end if
       kk = ig(jj) + 1
     end do  
   end if
   kk = ig(k) + 1
   if(kp.le.nfrt) then
     do jj = kp,nfrt
       lg(jj-1) = lg(jj)
       kl = kk - jj
       term = -eq(jj)/pivot
       if(solv) then
         b(jj-1) = b(jj) + term*r
       end if  
       call saxpb(eq(kp),a(kk+k),term,jj-k,a(kl+k))
       if(km.gt.0) then
         kl = kl + 1
         call saxpb(eq,a(kk),term,km, a(kl))
       end if
       kk = ig(jj) + 1
     end do  
   end if
   call pconsd(a(kk-nfrt),nfrt,0.0d0)
   b(nfrt)  = 0.0
   lg(nfrt) = 0

end
