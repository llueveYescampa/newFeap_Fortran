subroutine phstio(iu,irec,hh,nh,isw,tfile,itrec)
implicit  none
integer          iu,irec,    nh,isw,      itrec
double precision          hh(nh)
character                          tfile*12

!  Purpose: I/O control for history terms on disk
!  esta subroutina la modifique yo para que pudiera correr en el compilador fortran
!  el 8096 parece que es una limitacion el el tamano del record en la version DOS ???.
!  OS/2 tambien tiene una limitacion. Un record  de 103784 fallo, pero el limite
! exacto no lo conozco todavia

!  Inputs:
!     iu         - Logical unit number for I/O
!     irec       - Record number
!     hh(*)      - History terms to write (isw = 2 or 22 only)
!     nh         - Number of history terms
!     isw        - Switch: =1 or 11 for read; = 2 or 22 for write
!     tfile      - File name for history terms
!     itrec      - Record length of data

!  Outputs:
!     hh(*)      - History terms read (isw = 1 or 11 only)


   integer jtrec,i
   jtrec=itrec/nh
   
!  direct access read/write

   write(*,*) 'RECL =',tfile,' ',itrec
   
   if(isw.lt.10) then
     if(itrec .le. 8096) then
       open(iu,file=tfile,access='direct',recl=itrec,status='unknown')
       select case (isw)
       case (1)
         read (iu,rec=irec) hh
       case (2)
         write(iu,rec=irec) hh
       end select
     else
       open(iu,file=tfile,access='direct',recl=jtrec,status='unknown')
       select case (isw)
       case(1)
         do i=1,nh
           read (iu,rec=i) hh(i)
         end do  
       case(2)
         do i=1,nh
           write(iu,rec=i) hh(i)
         end do
       end select
     end if  
     close(iu)  
   else
!    sequential access read/write

     if(isw.eq.11) then
       read (iu) hh
     else if(isw.eq.22) then
       write(iu) hh
     end if  
   end if
      
end
