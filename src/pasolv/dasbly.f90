subroutine dasbly(s,p,ld,jp,ns,alfl,aufl,bfl,b,al,au,ad)
implicit none
integer  ns,ld(ns),jp(*)
double precision   s(ns,ns),p(ns),b(*),al(*),au(*),ad(*)
logical  alfl,aufl,bfl
      

!  Purpose:  Assemble symmetric or unsymmetric arrays for 'dasol'

!  Inputs:
!     s(ns,ns)   - element matrix
!     p(ns)      - Element vector
!     ld(ns)     - Equation numbers for elements
!     jp(*)      - Pointer to end of row/columns in profile
!     alfl       - Flag, unsymmetric if true
!     aufl       - Flag, assemble tangent if true
!     bfl        - Flag, assembel vector if true

!  Outputs:
!     b(*)       - Assembled vector
!     al(*)      - Assembled lower part of matrix
!     au(*)      - Assembled upper part of matrix
!     ad(*)      - Assembled matrix diagonal

   integer  i,ii, j,jc,jj


!  Loop through rows to perform assembly
   do i = 1,ns
     ii = ld(i)
     if(ii.gt.0) then
       if(aufl) then

!       Loop through columns to perform assembly

         do j = 1,ns
           if(ld(j).eq.ii) then
             ad(ii) = ad(ii) + s(i,j)
           else if(ld(j).gt.ii) then
             jc = ld(j)
             jj = ii + jp(jc) - jc + 1
             au(jj) = au(jj) + s(i,j)
             if(alfl) then
               al(jj) = al(jj) + s(j,i)
             end if  
           end if
         end do  
       end if
       if(bfl) then
         b(ii)  = b(ii)  + p(i)
       end if  
     end if
   end do  

end
