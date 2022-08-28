subroutine just(y,k,n0)
implicit    none
integer     k,n0
character*1 y(k)

!  Purpose: Complete parser alignment of input data

!  Inputs:
!     y(*)        - Input string
!     k           - Length of input string
!     n0          - Field width for parsed data

!  Outputs:
!     y(*)        - Adjusted data to field widths of n0

   integer     i,j,kl,l,n1,n2
   character*1 yi,ze,ni,mi,pl,dt,sp
   logical test

   data ze,ni,mi,pl,dt,sp/'0','9','-','+','.',' '/

   test = .true.
   n1 = n0 - 1
   n2 = n1 - 1
   do i = 1,k,n0
     do j = i,i+n1
       if(y(j).ne.sp) then 
         test = .false.
         exit
       end if  
     end do  
     if (test) then
       y(i+n1) = ze
     end if  
     if(y(i+n1).ne.sp) then 
       exit
     end if  
     yi = y(i)
     if((yi.ge.ze.and.yi.le.ni).or.(yi.eq.mi).or.(yi.eq.pl).or.(yi.eq.dt)) then
       do j = i+n2,i,-1
         if(y(j).ne.sp) then 
           exit
         end if
       end do  
       kl = n1 + i - j
       do l = j,i,-1
         y(l+kl) = y(l)
         y(l) = sp
       end do  
     end if
   end do  

end
