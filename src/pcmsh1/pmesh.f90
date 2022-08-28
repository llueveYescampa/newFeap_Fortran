subroutine pmesh(idl,ie,d,id,x,ix,f,t,ndf,ndm,nen1,iii,prt)
implicit none
integer  ndf,ndm,nen1,iii,idl(6),ie(9,*),id(ndf,*),ix(nen1,*)
double precision   d(18,*),x(ndm,*),f(ndf,*),t(*)
logical  prt      

!  Purpose: Data input routine for mesh description

!  Inputs:
!     idl(*)     - Working vector
!     ndf        - Number dof/node
!     ndm        - Spatial dimension of mesh
!     nen1       - Dimension of ix array
!     iii        - Flag to initialize arrays
!     prt        - Flag, output results if true

!  Outputs:
!     ie(9,*)    - Material set assembly data
!     d(*)       - Material set parameters
!     id(ndf,*)  - Boundary condition codes
!     x(ndm,*)   - Nodal coordinates of mesh
!     ix(nen1,*) - Element nodal connection lists
!     f(ndf,*)   - Nodal force/displacements
!     t(*)       - Nodal temperatures


   logical     error,pcomp
   integer     i, list, iocheck

   character*4 wd(13),cc,yyy*80,cds*12,tmp*12,fds*12,an*12

   integer        numnp,numel,nummat,nen,neq
   common /cdata/ numnp,numel,nummat,nen,neq

   double precision dq
   integer             n,ma,mct,iel,nel
   common /eldata/  dq,n,ma,mct,iel,nel

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   integer        n11a,n11b,n11c,ia
   common /mdat2/ n11a,n11b,n11c,ia(2,11)

   integer         maxa
include 'maxa.h'      
   double precision dm
   common           dm(maxa)

   data wd/'coor','elem','mate','boun','forc','temp','prin','nopr', &
           'bloc','pola','angl','ebou','end '/

   data  an/'  angles    '/
   data cds/' coordinates'/
   data tmp/' temperature'/

   data fds/' force/displ'/
   data list/13/
   
!  Initialize arrays

   error = .false.
   if(iii.ge.0) then
     prt = .true.
     call pconsd(dm(n11b),numnp,0.0d0)
     call pconsd(   f,numnp*ndf,0.0d0)
     call pconsi(id,numnp*ndf,0)
     if(iii.eq.0) then
       call pconsd( t,numnp,0.0d0)
       call pconsd( x,numnp*ndm,-999.0d0)
     end if  
   end if
   
   do while (.true.)
     if(ioRead.lt.0) then
       write(*,'(a,$)') '     Mesh 1 > '
     end if  
     call pintio(yyy,10)
     read(yyy,'(a4)',IOSTAT=iocheck) cc
     if (iocheck .eq. -1) then
!      End of file encountered
       call pend('pmesh ')
     else if(iocheck .ne. 0) then
!       Error was found
       call pperror('PMESH ',yyy)
       cycle
     end if
     if((ioRead.lt.0) .and. pcomp(cc,'help')) then
        call phelp(wd,list,'MESH ',0)
        cycle
     end if
     do i = 1,list
       if(pcomp(cc,wd(i))) then
         exit
       end if  
     end do
     if (i .gt. list) then
       cycle
     end if
     
   
     select case (i)
     case ( 1)
!      Nodal coordinate data input
       call genvec(ndm,x,cds,prt,error,.true.)
     case ( 2)
!      Element data input
       call pelin(idl,ix,nen1,nen,numnp,numel,error,prt)
     case ( 3)
!      Material data input
       call pmatin(d,x,ix,idl,ie,nummat,ndm,ndf,prt)
     case ( 4)
!      Read in restraint conditions for each node
       call pbcin(iii,idl,id,numnp,ndf,prt)
     case ( 5)
!      Force/displ data input
       call genvec(ndf,f,fds,prt,error,.false.)
     case ( 6)
!      Temperature data input
       call genvec(1,t,tmp,prt,error,.false.)
     case ( 7,8)
!      Set print flag
       prt = i.eq.7
     case ( 9)
!      Generate block of nodes and 4-node elements
       if(iii.lt.0) then
         write(ioWrite,3000)
       end if  
       call blkgen(ndm,nen,nen1,x,ix,prt)
     case (10)
!      Convert polar/cylindrical to cartesian coordinates
       call polar(x,ndm,prt)
     case (11)
!      Set boundary angles
       call genvec(1,dm(n11b),an,prt,error,.false.)
     case (12)
!      Set edge boundary values
       call pedges(iii,x,id,ndm,ndf,numnp,prt)
     case (13)
!    End of mesh
       if(error) then
         stop
       end if  
       return
     end select
   end do
      
3000  format(' **WARNING** Element connections necessary to use '  &
            ,'block in macro program')

end
