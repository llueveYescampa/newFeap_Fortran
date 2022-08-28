subroutine pbuff(m,ibuf,ilast,nv,is,ifl)
implicit none
integer  ibuf,ilast,nv,is,ifl,m(ibuf)

!  Purpose: Input/output routine for frontal program

!  Inputs:
!     m(*)     - Buffer storage
!     ibuf     - Buffer length (in integer words)
!     ilast    - Uset to go backward through data files
!     nv       - Record number for I/O
!     is       - Switch: =1 for read; =2 for write
!     ifl      - File name number for frontal equation storage on disk


   integer         i,j
   integer iocheck
   character*12 fname

   double precision         dimx,dimn
   integer                  nvp,npl
   common /nfrta/ dimx,dimn,nvp,npl

   character*12    tfile
   common /temfl1/ tfile(6)

   integer         itrec   ,nw1,nw2
   common /temfl2/ itrec(4),nw1,nw2


   if(is.eq.2) then
     nv = nv + 1
   end if  
   if(nv.lt.10) then
     write(fname,'(a10,i1)') 'Frontal.00',nv
   else if(nv.lt.100) then
     write(fname,'( a9,i2)') 'Frontal.0',nv
   else
     write(fname,'( a8,i3)') 'Frontal.',nv
   end if
   open(ifl,file=fname,form='unformatted',status='unknown')

!  Read record 'nv' from the file

   if(is.eq.1) then

     do i = 1,ibuf,8000
       read(ifl,IOSTAT=iocheck)  ilast,(m(j),j=i,min(ibuf,i+7999))
       if (iocheck .ne. 0) then
         call pend('PBUFF ')
       end if
     end do

!  Write record 'nv' from the file

   else if(is.eq.2) then

     if(nv.eq.1) then
       npl = ilast
     end if  
     do i = 1,ibuf,8000
       write(ifl,IOSTAT=iocheck) ilast,(m(j),j=i,min(ibuf,i+7999))
       if (iocheck .ne. 0) then
         call pend('PBUFF ')
       end if
     end do
     ilast = 0

   end if
   
   close(ifl)
   return
end
