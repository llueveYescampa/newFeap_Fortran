subroutine acheck(x,y,n0,nl)
implicit none
integer  n0,nl
character*1 x(nl),y(nl)
      
!  Purpose: Data parser

!  Inputs:
!     x        - Character input data
!     n0       - Field widths for parse
!     nl       - Number of characters to check

!  Outputs:
!     y        - Parsed array


   integer  i,ii,il, k

   ii = nl
   do while (x(ii).eq.' ')
      ii = ii -1
   end do

   do i = 1,nl
      y(i) = ' '
   end do   

   k = 0
   il= 0
   do i = 1,ii
     if(x(i).eq.',') then
        k  = k + n0
        if(k .gt. nl-n0) then
          exit 
        end if  
        il = k - i
     else
        y(i+il) = x(i)
     end if
   end do
   k  = k + n0

   call just(y,k,n0)

   do i = n0,nl,n0
     if(y(i).eq.' ') then
       y(i) = '0'
     end if  
   end do
end
