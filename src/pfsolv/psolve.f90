subroutine psolve(u,a,b,dr,m,xm,s,ld,ig,idl,nst,nrs,afac,solv, &
                        dyn,c1,ipd,rnorm,aengy,ifl)
implicit  none
logical   afac,solv,dyn
integer   m(*),ld(*),ig(*),idl(*),nst,nrs,ipd,ifl
double precision u(*),a(*),b(*),dr(*),xm(*),s(*),c1,aengy,rnorm

!  Purpose:  Frontal assembly and solution of equations
!            N.B. interface for psolve same as for profile solutions

!  Inputs:
!     u(*)      - Current nodal solution
!     a(*)      - Storage for frontal tangent matrix
!     b(*)      - Frontal right hand side/solution vector
!     dr(*)     - Residual
!     m(*)      - Buffer for front equations eliminated
!     xm(*)     - Mass matrix (diagonal only)
!     s(nst,*)  - Element tangent matrix
!     ld(nst)   - Element equation numbers
!     ig(*)     - Pointers to end of frontal matrix row/columns
!     idl(*)    - Element ordering for frontal solution
!     nst       - Size of element arrays
!     nrs       - Number of right hand sides to solve
!     afac      - Flag, factor matrix if true
!     solv      - Flag, perform forward/backward solution if true
!     dyn       - Flag, transient solution if true
!     c1        - Factor to add mass matrix to tangent
!     ipd       - Precision of !real! variables

!  Outputs:
!     dr(*)     - Solution increment
!     rnorm     - Residual norm
!     aengy     - Energy of solution increment
!     ifl       - File name number for frontal equation storage on disk


   integer   maxa
!   parameter(maxa = 100000)
include 'maxa.h'      


   logical   fa
   integer   i,ie,ii, jj, k, n,ne,np,nal,nfrt,ihsiz
   double precision      r,rxm


   integer         numnp,numel,nummat,nen,neq
   common /cdata/  numnp,numel,nummat,nen,neq

   logical         fl    ,pfr
   common /fdata/  fl(11),pfr

   integer         maxf
   common /frdata/ maxf

   integer         ior,iow
   common /iofile/ ior,iow

   double precision dimx,dimn
   integer                    nv,npl
   common /nfrta/   dimx,dimn,nv,npl

   character*12    tfile
   common /temfl1/ tfile(6)

   integer         itrec   ,nw1,nw2
   common /temfl2/ itrec(4),nw1,nw2

   integer         ihfac,ibuf
   common /temfl3/ ihfac,ibuf

   data fa/.false./

!  Record length factors (ihsiz = # int words; ihfac = # bytes/int wd)

   ipd   = 2
   ihfac = 4
   ihsiz = maxa
   if(afac) then
!    Set control data and zero
     
     if(fl(4)) then
       nal= 1 + (maxf*(maxf+3))/2
       if((nal+maxf)*ipd.gt.ihsiz) stop 'front too large'
       ibuf = (min(ihsiz - nal,(maxf+2)*neq))*ipd
       itrec(1) = ibuf*ihfac + ihfac
!      open (4,file=tfile(1),status='new',access='direct',form='unformatted',recl=itrec(1))
!      close(4)
       fl(4) = fa
     end if
     call pconsd(a,maxf*(maxf+1)/2,0.0d0)
     call pconsd(b,maxf,0.0d0)
     ig(1)=1
     ig(maxf+1)=0
     
     do n=2,maxf
       ig(n)=ig(n-1)+n
       ig(maxf+n)=0
     end do  
      
!    Begin loop through elements to perform front solution
     
     np  = 0
     nv  = 0
     nfrt= 0
     dimx  = 0.0
     dimn  = 0.0
     rnorm = 0.0
     
!    Frontal elimination program
     
     if(ior.lt.0) then
       write(*,'(a/a)') '   Solution status',' '
     end if  
     
     do n = 1,numel
     
!      Pick next element
     
       ne = idl(n)
     
!      Compute element arrays
     
       call formfe(u,dr,.true.,solv,fa,fa,3,ne,ne,1)
       if(ior.lt.0 .and. mod(n,20).eq.0) then
         write(*,'(a,i4,a,i4)') &
         '+  ->', n, ' elements completed.  Current front width is', nfrt
       end if  
     
!      Assemble element and determine eliminations
     
       call pfrtas(a,s,ld,ig,ig(maxf+1),ie,nfrt,nst)
     
       if(ie.ne.0) then
     
!        Eliminate equations
     
         do i=1,ie
           k=ld(i)
           jj=-ig(maxf+k)
           if(solv) then
             rnorm = rnorm + dr(jj)**2
             r=b(k)+dr(jj)
             dr(jj) = r
           end if
           if(dyn) then
            ii = ig(k)
            a(ii) = a(ii) + c1*xm(jj)
           end if
           if(np+8+nfrt*ipd.gt.ibuf) then
            call pbuff(m,ibuf,np,nv,2,ifl)
           end if
           m(np+1) = nfrt
           m(np+2) = k
           m(np+3) = jj
           call pfrtd(a,b,r,ig,ig(maxf+1),solv,m(np+5),nfrt,k)
     
!          Align last word of buffer for backsubstitution reads
     
           np = np + 8 + nfrt*ipd
           m(np) = nfrt
           nfrt = nfrt - 1
         end do  
       end if
     
!    End of forward elimination and triangularization
     
     end do
     
     rnorm = sqrt(rnorm)
     rxm   = 0.0
     if(dimn.ne.0.0d0) then
       rxm = dimx/dimn
     end if  
     if(ior.lt.0) then
       write(*,2002) dimx,dimn,rxm
     end if  
     write(iow,2002) dimx,dimn,rxm
     
!    Clear buffer one last time
     
     if(np.gt.0) then
       call pbuff(m,ibuf,np,nv,2,ifl)
     end if  
   end if  

!  Forward and back substitutions for solution

   if(solv .and. (.not. afac)) then
     call pfrtfw(b,dr,m,ipd,ibuf,maxf,nv,neq,nrs,ifl)
   end if  
   
   if(solv) then
     call pfrtbk(b,dr,m,ipd,ibuf,maxf,nv,neq,nrs,aengy,ifl)
   end if  

2002  format(5x,'Condition check: D-max',1p1e11.4,'; D-min',1p1e11.4, &
               '; Ratio',1p1e11.4)

      end
