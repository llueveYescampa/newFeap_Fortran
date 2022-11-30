subroutine pmacr2(id,x,f,f0,jd,b,dr,lct,ct,ndf,nneq,j)
implicit none
  integer          :: id(*),    jd(*), ndf,nneq,j
  double precision :: x(*),f(*),f0(*),b(*),dr(*),ct(3,*)
  character*4      :: lct(*)

!  Purpose: Command language instruction executions

!  Inputs:
!     id(*)     - List of equation numbers for each dof
!     x(*)      - Nodal coordinates for mesh
!     f(*)      - Nodal specified force/displacmeents
!     f0(*)     - Base solution
!     jd(*)     - Pointer to end of row/columns for profiles
!     b(*)      - Nodal solution array
!     dr(*)     - Residual/reactions
!     lct       - Command option to execute
!     ct(3,*)   - Command parameters
!     ndf       - Number dof/node
!     nneq      - Size of id, f, f0, and b arrays
!     j         - Command number to execute

!  Outputs:
!     various   - Depends on command number

   logical   pcomp,sfl
   integer   iocheck
   integer   i,mg,mh,mp,mdp,mdt, n,n1,n2,ndm,nlm
   integer  :: md = 0, mv = 0, mf = 0, mq = 0
   
   double precision propld,xtl
   character*4 ctl(2)
   character yyy*80
   !double precision uu(4000)
    double precision, dimension(:), allocatable :: uu

   include 'cdata.h'
   include 'fdata.h'
   include 'hdatb.h'
   include 'iofile.h'
   include 'iofild.h'
   include 'ldata.h'
   include 'ndata.h'
   include 'prlod.h'
   include 'psize.h'
   include 'rdata.h'
   include 'tbeta.h'
   include 'tdata.h'
   include 'temfl1.h'
   include 'ddata.h'

   allocate ( uu(4000) )  ! real dimension og uu should be used here

!  transfer to correct process

   select case (j)
   case ( 1)
!    set solution tolerance
     tol = ct(1,l)
   case ( 2)
!    set time increment
     dt = ct(1,l)
     if(fl(9)) then
       call setci(ior)
     end if
   case ( 3)
!    set loop start indicators
     lv = lv + 1
     lvs(lv) = l
     lve(lv) = int(ct(2,l))
     ct(1,lve(lv)) = 1.
   case ( 4)
!    loop terminator control
     n = int(ct(2,l))
     ct(1,l) = ct(1,l) + 1.0
     if(ct(1,l).gt.ct(1,n)) then
       lv = lv - 1
     end if
     if(ct(1,l).le.ct(1,n)) then
       l = n
     end if
   case ( 5)
!    input proportional load table
     npld = int(ct(1,l))
     prop = propld (ttim,npld)
   case ( 6)
!    data command
     do ! while (.true.)
       if(ior.lt.0) then
         write(*,'(a,a4,a,$)') ' Input ',lct(l),' macro >'
       end if
       call pintio(yyy,10)
       read(yyy,'(a4,6x,a4,6x,f15.0)',IOSTAT=iocheck) (ctl(i),i=1,2),xtl
       if (iocheck .ne. 0) then
          call myPerror('pmacr2',yyy)
          cycle
       end if
       if(pcomp(lct(l),ctl(1))) then
         if(pcomp(ctl(1),'tol ')) then
           tol = xtl
         end if
         if(pcomp(ctl(1),'dt  ')) then
           dt  = xtl
         end if
         exit ! retur_n
       else
!        error diagnostics
         if(ior.gt.0) then
           write(iow,'(a)')' **ERROR** Macro label mismatch on data command.'
           call pstop (-110) ! stop
         end if
         write(*,'(a)')' **ERROR** Macro label mismatch on data command.'
       end if
     end do
   case ( 7)
!    increment time - initialize force / solution vectors
     ttim = ttim + dt
     do i = 0,nneq-1
       dm(nt+i) = dm(nt+nneq+i)
     end do
     if(npld.gt.0) then
       prop = propld(ttim,0)
     end if
     write(iow,'(a,1pe12.5/,a,1pe12.5)') &
     '   Computing solution for time',ttim,'   Proportional load value is ',prop
     if(ior.lt.0) then
       write(  *,'(a,1pe12.5/,a,1pe12.5)') &
     '   Computing solution for time',ttim,'   Proportional load value is ',prop
     end if
     engy = 0.0
     rnmax = 0.0
!    update history on the disk
     if(.not.hfl) then
       hout = .true.
       call formfe(b,dr,.false.,.false.,.false.,.false.,6,1,numel,1)
       hout = .false.
     end if
!    update dynamic vectors for time step
     if(fl(9)) then
       call setci(ior)
       call update(id,dm(nt),b,dm(nv),dr,nneq,fl(9),1)
     end if
!    zero displacement increment for next time step
     call pconsd(b(nneq+1),nneq+nneq,0.0d0)
     fl(10) = .true.
   case ( 8)
!  input integration parameters and initialize vectors
     call param(lct(l),ct(1,l))
     if( .not. fl(9)) then ! retur_n
       call psetm(nv,nneq*4,'d',fl(9))
       call pconsd(dm(nv),nneq*4,0.0d0)
       nw = nv + nneq
!      set initial condition for transient solution

       n1 = nw + nneq - 1
       do i = 1,nneq
         dm(n1+i) = b(i)
       end do
       fl(9) = .true.
     endif     
   case ( 9)
!    update the current force vector f0
!    call raxpb(f, f0, prop, nneq, f0)
     call saxpb(f, f0, prop, nneq, f0)
   case (10)
!    subspace eigencomputations
     mf = int(ct(1,l))
     n2 = int(ct(2,l))
     mf = min(neq,max(1,mf))
     mq = min(mf+mf,mf+8,neq)
     if(n2.gt.0) then
       mq = min(mf+n2,neq)
     end if
     if(fl(2)) then
       call numass(dm(nl),neq,mq)
     end if
     if(mq.lt.mf) then
       write(iow,'(a,i4,a)') &
         ' **WARNING** Subspace set to',mq,' number of mass terms.'
       if (ior.lt.0) then
         write(*,'(a,i4,a)') &
           ' **WARNING** Subspace set to',mq,' number of mass terms.'
       end if
     end if
     mf = min(mf,mq)
     
     
     if( mf > 0) then  ! mf.le.0
       sfl = pcomp(lct(l),'prin')
!      establish addresses for eigen solutions
       md =  1
       mv = md + mq
       mg = mv + mq*neq
       mh = mg + mq*(mq+1)/2
       mdp= mh + mq*(mq+1)/2
       mdt= mdp+ mq
       mp = mdt+ mq
       nlm= max(mp+mq*mq,neq*mq+neq)
       if(nlm.gt.maxm/ipd) then
         if(ior.lt.0) then
           write(*,'(a,/5x,a,i7,a,i7)') &
           ' **ERROR** Subspace memory too large','Need =',nlm, &
           ' : Available =',maxm/ipd
           return
         else
           write(iow,'(a,/5x,a,i7,a,i7)') &
           ' **ERROR** Subspace memory too large','Need =',nlm, &
           ' : Available =',maxm/ipd
           call pstop (-209) ! stop
         end if
       end if
       nlm = ne/ipd + 1
       open(7,file=tfile(3),form='unformatted',status='new')
       write(7) (dm(i),i=1,nlm)
       close (7)
       nlm = nl - 1
       n1  = neq*ipd - ipd !  - ipr
       do i = 1,neq
         dm(i)    = dm(i+nlm)
         im(n1+i) = jd(i)
       end do
       nlm = neq + (neq + ipd -1)/ipd + 1
       call subsp(dm(1),uu(mv),uu(mg),uu(mh),uu(md),uu(mdp), &
          uu(mdt),uu(mp),dm(nlm),mf,mq,neq,shift,tol,sfl,25)
!      restore the solution to continue macro executions
       nlm = ne/ipd + 1
       open(7,file=tfile(3),form='unformatted',status='old')
       rewind 7
       read(7) (dm(i),i=1,nlm)
       close(7,status='delete')
     end if     
   case (11)
!    print eigenvectors
     if(mf.gt.0) then
       n1   = int(ct(1,l))
       n1   = max(1,min(mq,n1))
       n2   = mv + neq*(n1-1) -1

!    expand and move eigenvector for prints

       call pconsd(dm,nneq,0.0d0)
       do i = 1,nneq
         n = id(i)
         if(n.ne.0) then
           dm(i) = uu(n2+n)
         end if
       end do
       call prtdis(x,dm,ttim,uu(md+n1-1),ndm,ndf,1,numnp,1)
     else
       write(iow,'(a)')' **ERROR** Need eigensolution.'
       if(ior.lt.0) then
         write(*,'(a)')' **ERROR** Need eigensolution.'
       end if
     end if
   case (12)
!  set identity matrix for mass
     if(fl(5)) then
       call psetm(nl,neq,'d',fl(5))
     end if
     fl(2) = .true.
     call pconsd(dm(nl),neq,1.d0)
   end select
!  formats
   deallocate (uu)  
end

