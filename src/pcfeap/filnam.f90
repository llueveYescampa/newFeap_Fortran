subroutine filnam()
implicit none

!  Purpose:  Set filenames for execution

   logical   linp,lout,lres,lsav
   integer   iop,irs,isv,ifill ! i,
   character y*1,disknm*1,pdisknm*1,wd(2)*6
   
   character*12 finp,fout,fres,fsav,pinp,pout,pres,psav,edate
   
   integer iocheck, jfill,jplot
   logical test, lfil

   include 'iofild.h'
   include 'temfl1.h'
   include 'vdata.h'

   data wd/'new   ','exists'/


!  output version data to screen
   call pcdate(edate)
   write(*,2000) versn,edate

!  look to see if any problem has been run
   inquire(file='feap.nam',exist=lfil)
   if(lfil) then
      test = .false.
      open(3,file='feap.nam',status='old')
      read(3,'(4a12,a1,i2)') pinp,pout,pres,psav,disknm,jfill
      finp = pinp
      fout = pout
      fres = pres
      fsav = psav
   else
      test = .true.
!     default installation parameters
      pinp   = 'NONE'
      disknm = 'c'
      jfill  = 2
      
!     set scratch disk name
      write(*,2007) disknm
      read (*,'(a1)') pdisknm
      if(pdisknm.ne.' ') then
        disknm = pdisknm
      end if  
!     set graphics monitor type
      write(*,2008) jfill
      read (*,'(i1)') ifill
      if(ifill.ne.0) then 
        jfill = ifill
      end if  
   end if
   
!  name file for input data
   do while(.true.)
     if (test) then
       write(*,2000) versn,edate
       do while (.true.)
         write(*,2001) pinp
         read (*,'(a12)',IOSTAT=iocheck) finp
         if (iocheck .ne. 0) then
           write(*, '(a)') ' *** ERROR on read *** reinput'
           cycle
         end if
         if(finp.eq.' ') then 
           finp = pinp
         end if  
!        check if the input files exists
         inquire(file=finp,exist=linp)
         if(.not.linp) then
            write(*, '(a)') ' *** ERROR - - Specified input file does not exist'
            cycle
         else
            pout = finp
            pres = finp
            psav = finp
            call pdisk('O',pout)
            call pdisk('R',pres)
            call pdisk('S',psav)
            exit
         end if
       end do
       pinp = finp
       
!      name file for output data
       do while (.true.)
         write(*,2002) pout
         read (*,'(a12)',IOSTAT=iocheck) fout
         if (iocheck .ne. 0) then
           write(*, '(a)') ' *** ERROR on read *** reinput'
           cycle
         end if
         if(fout.eq.' ') then 
           fout = pout
         end if  
         pout = fout
         exit
       end do
       
!      name file for restart read data
       
       do while (.true.)
         write(*,2003) pres
         read (*,'(a12)',IOSTAT=iocheck) fres
         if (iocheck .ne. 0) then
           write(*, '(a)') ' *** ERROR on read *** reinput'
           cycle
         end if
         if(fres.eq.' ') then 
           fres = pres
         end if  
         pres = fres
         exit
       end do
       
!      name file for restart save data
       do while (.true.)
         write(*,2004) psav
         read (*,'(a12)',IOSTAT=iocheck) fsav
         if (iocheck .ne. 0) then
           write(*, '(a)') ' *** ERROR on read *** reinput'
           cycle
         end if
         if(fsav.eq.' ') then 
           fsav = psav
         end if  
         psav = fsav
         exit
       end do  
     end if
!    check file status and input if necessary
     inquire(file=finp,exist=linp)
     if(.not.linp) then
       test=.true.
       cycle
     end if  
     inquire(file=fout,exist=lout)
     iop = 1
     if(lout) then 
       iop = 2
     end if  
     inquire(file=fres,exist=lres)
     irs = 1
     if(lres) then 
       irs = 2
     end if  
     inquire(file=fsav,exist=lsav)
     isv = 1
     if(lsav) then 
       isv = 2
     end if  
     write(*,2005) finp,wd(2),fout,wd(iop),fres,wd(irs),fsav,wd(isv)
     read(*,'(a1)') y
     if(y.eq.'s' .or.  y.eq.'S') then
       call pstop(-162) ! stop
     end if  
     if(y.ne.'y' .and. y.ne.'Y') then
       test=.true.
       cycle
     end if  
     
!    save a copy of the current filenames
     
     if(.not.lfil) then 
       open(3,file='feap.nam',status='new')
     end if  
     rewind 3
     write(3,'(4a12,a1,i2)') finp,fout,fres,fsav,disknm,jfill
     close(3)
     
!    Erase the output file if it exists
     
     if(lout) then
       open(3,file=fout,status='old')
       close(3,status='delete')
     end if
     write(*,2006)
     
!    Open the files for input and output
     
     open(unit=iodr,file=finp,status='old')
     open(unit=iodw,file=fout,status='new')
     
!    Set the scratch disk names and locations
     
     tfile(1) = 'frnt.tem'
     tfile(2) = 'hist.tem'
     tfile(3) = 'mesh.tem'
     tfile(4) = 'stre.tem'
     tfile(5) = fres
     tfile(6) = fsav
     !do i = 1,4
     !  call pdisk(disknm,tfile(i))
     !end do  
     exit
   end do
   return


!  Format statements

2000  format(////6x,                                                     &
      ' F I N I T E   E L E M E N T   A N A L Y S I S   P R O G R A M'// &
       13x,'VERSION: ',3a12/34x,a12)
2001  format(/13x,'I n p u t    F i l e n a m e s'//     &   
            15x,'Input data   (default: ',a12,') :',$)
            
2002  format(15x,'Output data  (default: ',a12,') :',$)

2003  format(15x,'Restart read (default: ',a12,') :',$)

2004  format(15x,'Restart save (default: ',a12,') :',$)

2005  format( /13x,'Files are set as follows :'//                  &
        32x,'Filename',9x,'Status'/                                &
        15x,'Input    (read ) : ',a12,5x,a6/                       &
        15x,'Output   (write) : ',a12,5x,a6/                       &
        15x,'Restart  (read ) : ',a12,5x,a6/                       &
        15x,'Restart  (write) : ',a12,5x,a6//                      &
        13x,'Caution, existing write files will be overwritten.'// &
        13x,'Are filenames correct? ( y or n, s=stop) : ',$)
     
2006  format(/12x,'R U N N I N G   P C F E A P   P R O B L E M   N O W')

2007  format(/13x,'I N S T A L L A T I O N   P A R A M E T E R S'//  &
             15x,'Disk Name For Scratch Files:'/                     &
             17x,'             (default = ',a1,9x,') :',$)
     
2008  format(15x,'Set: Graphics Terminal Type:'/              &
             17x,'1=CGA; 2=EGA (default = ',i2,8x,') :',$)

end
