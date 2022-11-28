subroutine pmacio (jct,lct,ct,wd,nwd,nlp,ll)
implicit none
  integer          :: jct(*),nwd,nlp,ll
  character*4      :: lct(*),wd(nwd)
  double precision ::  ct(3,*)

!  Purpose: Command language instruction input subprogram

!  Inputs:
!     wd(*)      - Command names
!     nwd        - Number of commands
!     nlp        - Number of !loop! command in list

!  Outputs
!     jct(*)     - Number of command to execute
!     lct(*)     - Command option name
!     ct(3,*)    - Command parameters
!     ll         - Number of commands to execute


   logical    pcomp
   integer    i, j, l
   integer iocheck, test
   character  clab1*4,clab2*4,tary*10,yyy*80

   include 'iofile.h'

!  initiate the read of macro statements

   if(ior.gt.0) then
     call prthed(iow)
     write(iow,2001)
   end if

!  read macro cards

   ll = 1
   jct(1) = nlp
   ct(1,1) = 1.0
   do while (.true.)
     do while (.true.)
       if(ior.lt.0) then
        call pctime(tary)
        write(*,2002) tary,ll
       end if
       ll = ll + 1
       call pintio(yyy,15)
       read(yyy,'(2(a4,11x),3f15.0)',IOSTAT=iocheck)clab1,clab2,(ct(i,ll),i=1,3)
       if (iocheck .ne. 0) then
         call myPerror('pmacio',yyy)
         cycle
       end if
       if(ior.lt.0.and.pcomp(clab1,'help')) then
         call phelp(wd,nwd,'MACRO',1)
         ll = ll - 1
         cycle
       end if
       if(ior.gt.0.and.pcomp(clab1,'end ')) then
         exit
       end if
       if(ior.lt.0) then
         if(pcomp(clab1,'exit')) then
           ll = -1
         end if
         if(pcomp(clab1,'q   ').or.pcomp(clab1,'quit')) then
           ll = -2
         end if
         if(ll.lt.0) then
           return
         end if
       end if
     
!      set execution flag
     
       lct(ll) = clab2
       test=0
       do j = 1,nwd
         if(pcomp(clab1,wd(j))) then
           jct(ll) = j
           if(ior.gt.0) then
             write(iow,'(7x,a4,1x,a4,1x,3g12.5)') clab1,clab2,(ct(l,ll),l=1,3)
             test=100
             exit !! go to 100  or cycle
           end if
           ll = ll + 1
           test=150
           exit !! go to 150  or out oh the while
         end if
       end do
       if (test .eq. 100) then
         cycle
       else if(test .eq. 150) then
         exit
       else
         call myPerror('pmacio',yyy)
         ll = ll - 1
         cycle
       end if
     
     exit
     end do
     jct(ll)= nlp+1
     
!    set loop markers
     
     j = 0
     test = 0
     do l = 2,ll-1
       if(jct(l).eq.nlp) then
         j = j + 1
       end if
       if(j.gt.8) then
         test  = 400
         exit
       end if
       if(jct(l).eq.nlp+1) then
         j = j - 1
       end if
       if(j.lt.0) then
         test = 400
         exit
       end if
     end do
     
     if(test .eq. 400) then
!      error messages
       write(iow,3000)
       if(ior.gt.0) then
         call pstop(-129) ! stop
       end if
       if(ior.lt.0) then
         write(*,3000)
       end if
       cycle
     end if
     
     test=0
     
     if (j.ne.0) then
       if (ior.gt.0) then
         test = 400
       else if (ior.lt.0) then
         ll = ll - 1
       end if
       if (test .eq. 0) then
         cycle
       end if  
     end if
     
     if(test .eq. 400) then
!      error messages
       write(iow,3000)
       if(ior.gt.0) then
         call pstop(-154) ! stop
       end if
       if(ior.lt.0) then
         write(*,3000)
       end if
       cycle
     end if
     
     
     
     do l = 1,ll-1
       if(jct(l).ne.nlp) then
         cycle
       end if
       j = 1
       test=0
       do i = l+1,ll
         if(jct(i).eq.nlp) then
           j = j + 1
         else if(jct(i).eq.nlp+1) then
           j = j - 1
         end if
         if(j.eq.0) then
           ct(2,i) = l
           ct(2,l) = i
           test=220
           exit
         end if
       end do
       
       if (test .eq. 220) then
         cycle
       else
         test = 400
         exit
       end if
     end do
     
     if(test .eq. 400) then
!      error messages
       write(iow,3000)
       if(ior.gt.0) then
         call pstop(-196) ! stop
       end if
       if(ior.lt.0) then
         write(*,3000)
       end if
       cycle
     end if

   exit
   end do

   return



2001  format('  M a c r o   I n s t r u c t i o n s'//         &
      '  Macro Statement  Variable 1  Variable 2  Variable 3')
     
2002  format(' Input MACRO: !exit! = stop with restart',        &
         '; !quit! = quit.'/3x,'Time = ',a10,' Macro',i3,'> ',$)
     
3000  format(' **ERROR**  Wrong loop/next order, or > 8 loops')
end
