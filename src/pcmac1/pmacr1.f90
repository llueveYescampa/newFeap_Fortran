subroutine pmacr1(id,ie,ix,ld,d,s,x,f,f0,t,jd,b,dr,lct,ct, &
           ndf,ndm,nen1,nst,nneq,prt,j)
implicit none
  integer          :: id(*),ie(*),ix(*),ld(*),jd(*),ndf,ndm,nen1,nst,nneq,j
  double precision :: d(*),s(*),x(*),f(*),f0(*),t(*),b(*),dr(*),ct(3,*)
  character*4      :: lct(*)
  logical          :: prt


!  Purpose: Command language instruction executions

!  Inputs:
!     id(*)     - List of equation numbers for each dof
!     ie(*)     - Material set assembly data
!     ix(*)     - Element nodal connection lists
!     ld(*)     - Element equation number assembly vector
!     d(*)      - Material set parameters
!     s(*)      - Element vector/matrix storage
!     x(*)      - Nodal coordinates for mesh
!     f(*)      - Nodal specified force/displacmeents
!     f0(*)     - Base solution
!     t(*)      - Nodal temperatures
!     jd(*)     - Pointer to end of row/columns for profiles
!     b(*)      - Nodal solution array
!     dr(*)     - Residual/reactions
!     lct       - Command option to execute
!     ct(3,*)   - Command parameters
!     ndf       - Number dof/node
!     ndm       - Spatial dimension of FEM mesh
!     nen1      - DImension of ix array
!     nst       - Dimension of element arrays
!     nneq      - Size of id, f, f0, and b arrays
!     prt       - Flag, output results if true
!     j         - Command number to execute

!  Outputs:
!     various   - Depends on command number


   logical   pcomp,sflg,tflg
   integer   ietyp,ml1,na,nal,n1,n2,n3,n4,i
   integer   nnpo,nnlo,nnmo,ndmo,ndfo,nhio,nhfo,nrco
   double precision    shft,rnorm,dot
   double precision :: enold != 0.0d0

   character tt*10

   include 'adata.h'   
   include 'cdata.h'
   include 'errind.h'
   include 'fdata.h'
   include 'frdata.h'
   include 'hdatb.h'
   include 'iofild.h'
   include 'iofile.h'
   include 'ldata.h'
   include 'mdat2.h'
   include 'ndata.h'
   include 'prlod.h'
   include 'rdata.h'
   include 'tbeta.h'
   include 'tdata.h'
   include 'temfl1.h'
   include 'temfl2.h'
   include 'ddata.h'

!  transfer to correct process

   n1 = 1
   n3 = 1
   select case (j)
   case ( 1)
!    print stress values
     n1 = int(ct(1,l))
     n2 = int(ct(2,l))
     n3 = int(ct(3,l))
     n3 = max(n3,1)
     n4 = numnp + 1
     if (pcomp(lct(l),'node')) then
       n1 = max(1,min(numnp,n1))
       n2 = max(n1,min(numnp,n2))
       enerr = 0.0
       if(.not.fl(11)) then
         call pconsd(aa,8*numnp,0.0d0)
         call formfe(b,dr,.false.,.false.,.false.,.false.,8,1,numel,1)
         call pltstr(aa,aa(n4),numnp)
       end if
       call prtstr(aa,aa(n4),numnp,n1,n2,n3)
       fl(11) = .true.
     else if (pcomp(lct(l),'erro')) then
       n1 = max(n1,1)
       n2 = 8*numnp
       !call pconsd(aa(n2+1),n2-numnp,0.0d0)
       do i = 1,n2
         aa(n2+i) = 0.0
       end do  
       enerr = 0.0
       do i = 1,n1
         call pconsd(aa,n2,0.0d0)
         call formfe(b,dr,.false.,.false.,.false.,.false.,8,1,numel,1)
         call pltstr(aa,aa(n4),numnp)
         call addvec(aa(n2+1),aa(n4),n2-numnp)
       end do  
       fl(11) = .true.
       eerror = 0.0
       efem   = 0.0
       ebar   = 0.05*sqrt(enerr/numel)
       ietyp  = 1
       call formfe(b,dr,.false.,.false.,.false.,.false.,7,1,numel,1)
       call prterr
     else
       if(pcomp(lct(l),'all ')) then
         n2 = numel
       else
         n1 = max(1,min(numel,n1))
         n2 = max(n1,min(numel,n2))
       endif
       call formfe(b,dr,.false.,.false.,.false.,.false.,4,n1,n2,n3)
     end if
     return
   case (2,3,9)
!    solution step
     
     shft = c1
     sflg = fl(9)
     if(j.eq.9) then
       if(.not.fl(8)) then
         return
       end if  
       fl(7) = .false.
       tflg  = .false.
     
!    form tangent stiffness
     
     else
     
!     unsymmetric tangent (must use profile solver)
     
       if(j.eq.2) then
         fl(6) = .true.
       end if  
     
!     symmetric tangent
     
       if(j.eq.3) then
         fl(6) = .false.
       end if  
       if(ct(1,l).ne. 0.0d0) then
         fl(8) = .true.
         fl(7) = .false.
!         call raxpb(f, f0, prop, nneq, dm(nt+nneq))
         call saxpb(f, f0, prop, nneq, dm(nt+nneq))
         call pload(id,dm(nt),dr,nneq,dm(nl),dm(nw))
       end if
       shift= 0.
       tflg = .true.
       if(.not.fl(9).and.ct(2,l).ne. 0.0d0) then
         if(fl(2)) then
           sflg = .true.
           shift= ct(2,l)
           shft = -shift
           if(ior.lt.0) then
             write(*,'(a,1pe11.4,a)')'   Shift of',shift,' applied with mass'
           end if  
           write(iow,'(a,1pe11.4,a)')'   Shift of',shift,' applied with mass'
         else
           if(ior.lt.0) then
             write(*,'(a)') '   Shift requested but no mass matrix exists.'
           end if  
           write(iow,'(a)') '   Shift requested but no mass matrix exists.'
           if(ior.gt.0) then
             call pstop(-172) !stop
           end if  
           return
         end if
       end if
     end if
     
!    call the solve routine to assemble and solve the tangent matrix
     
     na = maxf + 1
     nal= (maxf*(maxf+1))/2 + na
     call psolve(b,aa(na),dr,dm(nl),s,ld,jd,nst,1  ,tflg,fl(8),sflg,shft,rnorm,engy)
     call pctime(tt)
     if(fl(8)) then
       fl(8) = .false.
       if(tflg) then
         write(iow,'(a,1pe15.7,12x,a,a10)') '   | R(i) | =',rnorm,'time=',tt
         if(ior.lt.0) then
           write(*,'(a,1pe15.7,12x,a,a10)') '   | R(i) | =',rnorm,'time=',tt
         end if  
       end if
       
       if (rnmax.eq.0.0d0) then
         rnmax = abs(engy)
         if(ct(3,l).le.0.0) then
           ct(3,l) = 0.6
         end if  
         enold = rnmax*0.9999
       end if
       write(iow,2004) rnmax,engy,tol
       if(ior.lt.0) then
         write(*,2004) rnmax,engy,tol
       end if  
       if(abs(engy).le.tol*rnmax) then
         ct(1,lve(lv)) = ct(1,lvs(lv))
         l = lve(lv) - 1
       else if(pcomp(lct(l),'line')) then
!    line search
         if(abs(engy).gt.ct(3,l)*enold) then
          ml1 = 1 + nneq
          call serchl(dm(nt),id,engy,aa,b,dr,ct(3,l),aa(ml1),nneq)
         end if
       end if
       call update(id,dm(nt),b,dm(nv),dr,nneq,fl(9),2)
     else
       write(iow,'(40x,a,a10)') 'time=',tt
       if(ior.lt.0) then
         write(*,'(40x,a,a10)') 'time=',tt
       end if  
     end if
     return
   case (4)
!    form out of balance force for time step/iteration
     
     if(fl(8)) then
       return
     end if
     
!     call raxpb(f, f0, prop, nneq, dm(nt+nneq))
     call saxpb(f, f0, prop, nneq, dm(nt+nneq))
     call pload(id,dm(nt),dr,nneq,dm(nl),dm(nw))
     call formfe(b,dr,.false.,.true.,.false.,.false.,6,1,numel,1)
     if(pcomp(lct(l),'acce').and.ttim.eq.0.0 .and. nop.eq.3) then
       call inaccl(id,dr,dm(nl),dm(nw),nneq)
     end if
     rnorm = sqrt(dot(dr,dr,neq))
     write(iow,'(a,1pe15.7)') '   | R(i) | = ',rnorm
     if(ior.lt.0) then
       write(*,'(a,1pe15.7)') '   | R(i) | = ',rnorm
     end if  
     fl(8) = .true.
     return
   case (5)
!    form a lumped mass approximation
     
     if(fl(5)) then
       call psetm(nl,neq,'d',fl(5))
     end if  
     call pconsd(dm(nl),neq,0.0d0)
     fl(2) = .true.
     call formfe(b,dm(nl),.false.,.true.,.false.,.false.,5,1,numel,1)
     return
   case (6,8)
!    compute reactions and print
     
     if(pcomp(lct(l),'all ')) then
       n2 = numnp
     else
       n1 = int(ct(1,l))
       n2 = int(ct(2,l))
       n3 = int(ct(3,l))
       n1 = max(1,min(numnp,n1))
       n2 = max(n1,min(numnp,n2))
       n3 = max(1,n3)
     end if
     if(j.eq.6) then
       call pconsd(dr,nneq,0.0d0)
       call formfe(b,dr,.false.,.true.,.false.,.true.,6,1,numel,1)
       call prtrea(dr,ndf,numnp,n1,n2,n3)
     else
       if(nop.le.2) then
         call prtdis(x,dm(nw+nneq),ttim,prop,ndm,ndf,n1,n2,n3)
       else
         call prtdis(x,b,ttim,prop,ndm,ndf,n1,n2,n3)
       end if
     end if
     return
   case (7)
!    check mesh for input errors
     
     call formfe(b,dr,.false.,.false.,.false.,.false.,2,1,numel,1)
     return
   case (10)
!    modify mesh data (cannot change profile of stiffness/mass)
     
     i = -1
     call pmesh(ld,ie,d,id,x,ix,f,t,ndf,ndm,nen1,i,prt)
     if (i.gt.0) then
!      error diagnostics
       write(iow,'(a)')' **ERROR** Attempt to change profile during mesh.'
       if(ior.gt.0) then
         call pstop(-293) ! stop
       else  
         write(*,'(a)')' **ERROR** Attempt to change profile during mesh.'
       end if  
     end if  
     return
   case (11)
!    restart previously run problem
     
     open (7,file=tfile(5),form='unformatted',status='old')
     read(7) nnpo,nnlo,nnmo,ndmo,ndfo,nhio,nhfo,nrco
     if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat) &
       .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf).and.(nrco.eq.nrec)   &
       .and.(nhfo-nhio.eq.nhf-nhi) ) then
        read(7) ttim,(b(i),i=1,3*nneq)
        read(7) (dm(i),i=nt,nt+nneq)
        if(fl(9)) read(7) (dm(i),i=nv,nv+4*neq)
        if(nrec.gt.0) then
          do j = 1,nrec
            call phstio(7,j,dm(nhi),nhf-nhi+1,11,tfile(5),0)
            call phstio(3,j,dm(nhi),nhf-nhi+1,2,tfile(2),itrec(2))
          end do  
        end if
        close(7)
     else
        if(ior.gt.0) then
          write(iow,'(a)')' **ERROR** File incompatible with current problem.'
        else if(ior.lt.0) then
          write(  *,'(a)')' **ERROR** File incompatible with current problem.'
        end if  
     end if
     return
   end select

!  formats

2004  format('   Energy Convergence Test'/'   E(1)=',1pe21.14, &
       ', E(i)=',1pe21.14,', Tol.=',1pe11.4)
     
end
