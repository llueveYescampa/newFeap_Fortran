subroutine pfrtas(a,s,ld,ig,lg,ie,nfrt,nst)
implicit none
integer  ld(*),ig(*),lg(*),ie, nfrt, nst
double precision   a(*),s(nst,*)

!  Purpose: Assembly and elimination determination for frontal program

!  Inputs:
!     s(nst,*) - Element matrix to assemble
!     ld(*)    - Equation numbers for element (negative tag to elim)
!     ig(*)    - Pointers to end of columns in frontal matrix
!     lg(*)    - Front/global equation numbers
!     nfrt     - Size of front before eliminations
!     nst      - Size of element arrays

!  Outputs:
!     a(*)     - Assembled frontal matrix
!     ld(*)    - Locations to eliminate
!     ie       - Number of equations to eliminate after assembly step

   integer  i,ii, j, k, l, m
   logical  test
!  Convert ld to front order

   do j=1,nst
     ii=abs(ld(j))
     i = 0
!    Check if ii is already in list

     test = .true.
     if(ii.ne.0) then
       if(nfrt.ne.0) then
         do i=1,nfrt
           if(ii.eq.abs(lg(i))) then
             test = .false.
             exit
           end if
         end do  
       end if
       
       if (test) then
!        Assign ii to next available entry and increase front width
         nfrt = nfrt + 1
         i = nfrt
       end if

!      Replace destination value by new value

       lg(i)=ld(j)
!      Set ld for assembly
       ld(j)=i
     end if
   end do  

!  Assemble element into front

   do j = 1,nst
     k = abs(ld(j))
     if(k.gt.0) then
       l = ig(k) - k
       do i = 1,nst
         m = abs(ld(i))
         if(m.gt.0 .and. m.le.k) then
           a(m+l) = a(m+l) + s(i,j)
         end if
       end do  
     end if
   end do  

!  Set up equations to be eliminated

   ie=0
   do i = nfrt,1,-1
     if(lg(i).lt.0) then
       ie = ie + 1
       ld(ie) = i
     end if
   end do  

end
