subroutine pmacr (ld,s,ie,d,id,x,ix,f,t,jd,f0,b,dr,ndf,ndm,nen1,nst,prt)
implicit none
  integer          :: ld(*),ie(*),id(*),ix(*),jd(*),ndf,ndm,nen1,nst
  double precision :: s(*),d(*),x(*),f(*),t(*),f0(*),b(*),dr(*)
  logical          :: prt


!  Purpose: Command language instruction driver

!  Inputs:
!     ld(*)     - Element equation number assembly vector
!     s(*)      - Element vector/matrix storage
!     ie(*)     - Material set assembly data
!     d(*)      - Material set parameters
!     id(*)     - List of equation numbers for each dof
!     x(*)      - Nodal coordinates for mesh
!     ix(*)     - Element nodal connection lists
!     f(*)      - Nodal specified force/displacmeents
!     t(*)      - Nodal temperatures
!     jd(*)     - Pointer to end of row/columns for profiles
!     f0(*)     - Base solution
!     b(*)      - Nodal solution array
!     dr(*)     - Residual/reactions
!     ndf       - Number dof/node
!     ndm       - Spatial dimension of FEM mesh
!     nen1      - DImension of ix array
!     nst       - Dimension of element arrays
!     prt       - Flag, output results if true

!  Outputs:
!     none


   logical   lxst
   integer   i, j, k, ll, nm1,nm2,nlp,nneq
   double precision  ct(3,100)
   character wd(24)*4,tary*10,lwd*3
   character*4 :: lct(100)
   ! common /ldatb/ lct(100)


   include 'cdata.h'
   include 'ddata.h'
   include 'fdata.h'
   include 'hdatb.h'
   include 'iofile.h'
   include 'ldata.h'
   include 'ndata.h'
   include 'rdata.h'
   include 'tbeta.h'
   include 'tdata.h'
   include 'temfl1.h'
   include 'temfl2.h'
   include 'prlod.h'

   data wd/'stre','utan','tang','form','mass','reac','chec','disp',&
           'solv','mesh','rest',                                   &
           'tol ','dt  ','loop','next','prop','data','time','beta',&
           'newf','subs','eigv','iden',                            &
           'plot'/


!  nmi = no. macro commands in 'pmacri'; nlp = loop number

   data nm1,nm2/11,12/
   
   nlp = nm1 + 3

!  set initial values of parameters

   call pinitc(engy,rnmax,shift,tol,dt,prop,ttim,npld)
   nw1  = nm1
   nw2  = nm2 + nw1
   nneq = ndf*numnp
   call psetm (nt,nneq*2,'d',fl(1))
   do j = 0,nneq*2-1
     dm(nt+j) = 0.0d0
   end do

!  input the macro commands

   do while (.true.)
     call pmacio (jct,lct,ct,wd,nw2+1,nlp,ll)
     if(ll.le.0) then
       exit
     end if  
!    execute macro instruction program
     
     lv = 0
     l = 1
     do while (l.le.ll)
       j = jct(l)
       i = l - 1
       call pctime (tary)
       if(l.ne.1.and.l.ne.ll) then
        write(iow,2001) i,wd(j),lct(l),(ct(k,l),k=1,3),tary
        if(ior.lt.0) then
          write(*,2001) i,wd(j),lct(l),(ct(k,l),k=1,3),tary
        end if  
       end if
       if(j.le.nw1) then
         call pmacr1(id,ie,ix,ld,d,s,x,f,f0,t,jd,b,dr,lct,ct, &
                     ndf,ndm,nen1,nst,nneq,prt,j)
       end if
       if(j.ge.nw1+1.and.j.le.nw2) then
         call pmacr2(id,x,f,f0,jd,b,dr,lct,ct,ndf,nneq,j-nw1)
       end if  
     
!      plot macro call
     
       !if(j.eq.nw2+1) then
       !  call pplotf(x,ix,b,lct(l),ct(1,l),ndf,ndm,nen1)
       !end if
       l = l + 1
     end do
     if (ior.lt.0) then
       cycle
     end if
     exit
   end do

   call pctime(tary)
   write(iow,'(a,15x,a,a10)') ' *End of macro execution*','time=',tary
   if(ior.lt.0) then
     write(*,'(a,15x,a,a10)') ' *End of macro execution*','time=',tary
   end if  
   if(.not.fl(4)) then
     close(4,status='delete')
   end if  

!  save restart information

   if(ll.lt.-1.or.fl(7)) then
     close(3,status='delete')
     return
   end if  
   
   inquire(file=tfile(6),exist=lxst)
   lwd = 'new'
   if(lxst) lwd = 'old'
   open (7,file=tfile(6),form='unformatted',status=lwd)
   rewind 7
   write(7) numnp,numel,nummat,ndm,ndf,nhi,nhf,nrec
   write(7) ttim,(b(i),i=1,3*nneq)
   write(7) (dm(i),i=nt,nt+nneq)
   if(fl(9)) then
     write(7) (dm(i),i=nv,nv+4*neq)
   end if  
   if(nrec.gt.0) then
     do j = 1,nrec
       call phstio(3,j,dm(nhi),nhf-nhi+1,1,tfile(2),itrec(2))
       call phstio(7,j,dm(nhi),nhf-nhi+1,22,tfile(6),0)
     end do  
    call pdefil(tfile,2,2)
   end if
   close(7)

!  formats

2001  format(' *Macro ',i3,' *',2(a4,1x), &
           'V1=',g10.3,' V2=',g10.3,' V3=',g10.3/40x,'time=',a10)
end
