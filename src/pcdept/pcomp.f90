logical function pcomp(a,b)
implicit  none
character a(4)*1,b(4)*1

! Edgar note: For C++ version use a string comparation functions.

!  Purpose: Determine match between alphanumeric data:
!           N.B. inc = parameter for difference between 'a' and 'A'

!  Inputs:
!     a(*)      - 4-character array to test
!     b(*)      - 4-character comparison array

!  Outputs:
!     pcomp     - Flag, true if match made.  
!                 N.B.  converts all characters to lower case before
!                       test performed

   integer   i,ia,ib,inc

   data inc /32/
   pcomp = .false.

   do i = 1,4
      ia = ichar(a(i))
      ib = ichar(b(i))
      if(ia .ne. ib .and. ia+inc .ne. ib .and. ia .ne. ib+inc ) then
         return
      end if  
   end do   
   
   pcomp = .true.

end

