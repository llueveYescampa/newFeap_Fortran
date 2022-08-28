subroutine pcdate(edate)
implicit  none

!     Purpose:  Put date on output file

!     Outputs:
!        edate     - character array with date information

   character edate*12, mon(12)*3
   character mydate*8
   character mytime*10
   character myzone*5
   integer   myvalues(8)

   data mon/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug', &
            'Sep','Oct','Nov','Dec'/

    call date_and_time( mydate, mytime, myzone, myvalues)    
    write(edate,'(a3,1x,i2,a2,i4)') mon(myvalues(2)),myvalues(3),', ',myvalues(1)
   
      

end
