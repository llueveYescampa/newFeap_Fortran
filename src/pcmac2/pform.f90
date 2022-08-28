subroutine pform(ul,xl,tl,ld,p,s,ie,d,id,x,ix,t,idl,u,b,ndf,ndm,nen1,nst)
implicit  none
integer ndf,ndm,nen1,nst,ld(ndf,*),ie(9,*), id(ndf,*), ix(nen1,*),idl(*)
double precision ul(ndf,*),xl(ndm,*),tl(*),p(*),s(nst,*),d(18,*),x(ndm,*),t(*),&
                 u(ndf,*),b(*)

!  Purpose: Compute element arrays and assemble global arrays

!  Inputs:
!     ul(ndf,*) - Element nodal solution parameters
!     xl(ndm,*) - Element nodal coordinates
!     tl(*)     - Element nodal temperatures
!     ld(*)     - Element equation numbers for each dof
!     p(*)      - Element vector
!     s(nst,*)  - Element matrix
!     ie(9,*)   - Material set assembly data
!     d(*)      - Material set parameters
!     id(*)     - Equation numbers for each dof
!     x(ndm,*)  - Nodal coordinates for mesh
!     ix(nen1,*)- Element nodal connection lists
!     t(*)      - Nodal temperatures for mesh
!     idl(*)    - Working vector
!     u(*)      - Nodal displacements and increments
!     ndf       - Number dof/node
!     ndm       - Spatial dimension of mesh
!     nen1      - Dimension of ix array
!     nst       - Dimension of element arrays

!  Outputs:
!     b(*)      - Residual/reactions or diagonal mass


   logical   efl
   integer   i,ii,iid,ijj,isx,ixi, j,jj, k
   integer   ne2,ne3,ne4, nu,nt0,nt1,nrot,numnp2
   double precision tm
   double precision dun,un


   integer        numnp,numel,nummat,nen,neq
   common /cdata/ numnp,numel,nummat,nen,neq

   double precision dq
   integer             n,ma,mct,iel,nel
   common /eldata/  dq,n,ma,mct,iel,nel

   logical        fl    ,pfr
   common /fdata/ fl(11),pfr

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   integer        n11a,n11b,n11c,ia
   common /mdat2/ n11a,n11b,n11c,ia(2,11)

   integer        nv,nw,nl
   common /ndata/ nv,nw,nl

   integer        nh1,nh2
   common /hdata/ nh1,nh2

   integer        nhi,nhf,ihbuff,irec,jrec,nrec
   logical                                      hfl,hout
   common /hdatb/ nhi,nhf,ihbuff,irec,jrec,nrec,hfl,hout

   double precision prop,ap
   integer                       iexp    ,ik    ,npld
   common /prlod/  prop,ap(6,10),iexp(10),ik(10),npld

   double precision beta,gamm,theta
   integer                          nop,nt
   common /tbeta/   beta,gamm,theta,nop,nt

   character*12    tfile
   common /temfl1/ tfile(6)

   integer         itrec   ,nw1,nw2
   common /temfl2/ itrec(4),nw1,nw2

   integer        isw,nn1,nn2,nn3
   logical                        afl,bfl,cfl,dfl
   common /xdata/ isw,nn1,nn2,nn3,afl,bfl,cfl,dfl

   integer         maxa
include 'maxa.h'      
   double precision dm
   integer                   im
   common           dm(maxa),im(maxa)

!  Set up local arrays before calling element library

   iel = 0
   tm  = 1. - theta
   efl = .false.
   if(.not.dfl.and.isw.eq.6) then 
     efl = .true.
   end if  
   if(bfl.and.isw.eq.3) then
     efl = .true.
   end if  
   if(isw.ne.3.or.nn1.eq.1) then
     irec = 0
   end if  
   ne2 = nen + nen
   ne3 = nen + ne2
   ne4 = nen + ne3
   nt0 = nt  - 1
   nt1 = nt0 + ndf*numnp
   numnp2 = numnp + numnp
   do nu = 1,numel
     n = idl(nu)
     if( (n.ge.nn1 .and. n.le.nn2) .and. (mod(n-nn1,nn3).eq.0)) then
     
!      Set up history terms
     
       ma  = ix(nen1,n)
       nh1 = ix(nen+1,n)
       nh2 = nh1
     
!      write(ioWrite,*) 'HFL,HOUT:',hfl,hout
!      write(ioWrite,*) 'IREC,JREC:',irec,ix(nen+2,n)
!      write(ioWrite,*) 'NHI,NHF:',nhi,nhf,nhf-nhi+1,itrec(2)
!      write(*,*) 'HFL,HOUT:',hfl,hout
     
       if(.not.hfl) then
         jrec= ix(nen+2,n)
         if(jrec.ne.irec) then
           if(hout .and. irec.ne.0) then
             call phstio(3,irec,dm(nhi),nhf-nhi+1,2,tfile(2),itrec(2))
           end if
           call phstio(3,jrec,dm(nhi),nhf-nhi+1,1,tfile(2),itrec(2))
           irec = jrec
         end if
       end if
       call pconsd(ul,5*nen*ndf,0.0d0)
       call pconsd(xl,nen*ndm,0.0d0)
       call pconsd(tl,nen,0.0d0)
       call pconsd(dm(n11a),nen,0.0d0)
       call pconsi(ld,nst,0)
       un = 0.0
       dun= 0.0
       call pangl(ix(1,n),nen,dm(n11a),dm(n11b),nrot)
       do i = 1,nen
         ixi= ix(i,n)
         ii = abs(ixi)
         if(ii.gt.0) then
           iid = ii*ndf - ndf
           nel = i
           tl(i) = t(ii)
           do j = 1,ndm
             xl(j,i) = x(j,ii)
           end do  
           do j = 1,ndf
             jj = ie(j,ma)
             if(jj.le.0) then
               cycle
             end if  
             ijj     = iid + jj
             k       = id(jj,ii)
             ul(j,i) = u(jj,ii)
             ul(j,i+nen) = u(jj,ii+numnp)
             ul(j,i+ne2) = u(jj,ii+numnp2)
             if(fl(9)) ul(j,i+ne3) = dm(nv+ijj)
             if(k.le.0) then
               ul(j,i+ne4) = theta*dm(nt1+ijj)+tm*dm(nt0+ijj)-u(jj,ii)
               dun = max(dun,abs(ul(j,i+ne4)))
             end if
             un = max(un,abs(ul(j,i)))
             if(dfl) then
               ld(j,i) = ijj
             else
               ld(j,i) = sign(k,ixi)
             end if
           end do  
         end if
       end do  
     
!      Form element array
     
       dq  = prop
       if(ie(7,ma).ne.iel) mct = 0
       iel = ie(7,ma)
       isx = isw
       if(efl .and. dun.gt.0.0000001d0*un .and. .not.afl) then
         isx = 3
       end if  
       if(nrot.gt.0) then
         call ptrans(ia(1,iel),dm(n11a),ul,p,s,nel,nen,ndf,nst,1)
       end if  
       call elmlib(d(1,ma),ul,xl,ix(1,n),tl,s,p,ndf,ndm,nst,iel,isx)
       if(nrot.gt.0) then
         call ptrans(ia(1,iel),dm(n11a),ul,p,s,nel,nen,ndf,nst,2)
       end if
     
!      Modify for non-zero displacement boundary conditions
     
       if(efl) then
         call modify(p,s,ul(1,ne4+1),nst)
       end if  
     
!      Assemble a vector if needed
     
       if(bfl) then
         do i = 1,nst
           j = abs(ld(i,1))
           if(j.ne.0) then
             b(j) = b(j) + p(i)
           end if
         end do  
       end if
     end if
   end do 

!  Put last history state on disk

   if(hout) then
     call phstio(3,jrec,dm(nhi),nhf-nhi+1,2,tfile(2),itrec(2))
   end if  

end
