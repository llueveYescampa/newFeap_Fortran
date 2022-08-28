subroutine pmatin(d,x,ix,idl,ie,nummat,ndm,ndf,prt)
implicit  none
integer   ix(*),idl(*),ie(9,*),nummat,ndm,ndf
double precision d(18,*),x(ndm,*)
logical   prt

!  Purpose: Input material set data

!  Inputs:
!     idl(*)     - Working vector
!     nummat     - Number of material sets
!     ndm        - Spatial dimension of mesh
!     ndf        - Number dof/node
!     prt        - Flag, output results if true

!  Outputs:
!     d(*)       - Material set parameters
!     x(ndm,*)   - Nodal coordinates of mesh
!     ix(nen1,*) - Element nodal connection lists
!     ie(9,*)    - Material set assembly data


   integer         maxa
include 'maxa.h'      

   integer iocheck
   logical test
   integer   i
   character yyy*80


   double precision ad
   common /adata/   ad(maxa)

   double precision dq
   integer             n,ma,mct,iel,nel
   common /eldata/  dq,n,ma,mct,iel,nel

   integer         nh1,nh2
   common /hdata/  nh1,nh2

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

!  Material data input

   if(prt) then
     call prthed(ioWrite)
     write(ioWrite,'(a)') '    M a t e r i a l    P r o p e r t i e s'
     if(ioRead.lt.0) then
       write(*,'(a)') '    M a t e r i a l    P r o p e r t i e s'
     end if  
   end if
   do n = 1,nummat
     do while(.true.)
       if(ioRead.lt.0) then
         write(*,'(a/3x,a,$)') ' Input: matl. no., elmt type','>'
       end if  
       call pintio(yyy,10)
       read(yyy,'(8i10)',IOSTAT=iocheck) ma,iel,(idl(i),i=1,ndf)
       if (iocheck .ne. 0) then
         call pperror('PMESH ',yyy)
         cycle
       else
         exit
       end if
     end do  
     if(ma.le.0) then
       return
     end if  

!    Set all zero inputs

     do i = 1,ndf
       ie(i,ma) = idl(i)
     end do
     test = .true.
     do i = 1,ndf
       if(idl(i).ne.0) then
         test = .false.
         exit
       end if  
     end do
     if (test) then
       do i = 1,ndf
         ie(i,ma) = i
       end do
     end if  
     ie(7,ma) = iel
     mct = 0
     nh1 = 0
     if(prt) then
       write(ioWrite,2000) ma,iel,(i,ie(i,ma),i=1,ndf)
       if(ioRead.lt.0) then
         write(*,2000) ma,iel,(i,ie(i,ma),i=1,ndf)
       end if  
     end if
     call elmlib(d(1,ma),ad,x,ix,dq,ad,ad,ndf,ndm,ndf,iel,1)
     if(nh1.eq.0 .and. mct.ne.0) then
       nh1 = mct
     end if  
     ie(8,ma) = nh1
   end do  

!  Formats

2000  format(/5x,'Material Set',i3,' for Element Type',i2,5x,//   &
         10x,'degree of freedom assignments    local    global' / &
         42x, 'number',4x,'number'/(36x,2i10))
end
