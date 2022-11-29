!subroutine psolve(u,a,b,dr,m,xm,s,ld,ig,idl,nst,nrs,afac,solv,dyn,c1,ipd,rnorm,aengy,ifl)
subroutine psolve(u,a     ,dr,xm    ,s,ld,ig,nst,nrs,afac,solv ,dyn ,c1  ,ipd,rnorm,aengy)
           
implicit none
  integer          :: ld(*),ig(*),nst,nrs,ipd
  double precision :: u(*),a(*),dr(*),xm(*),s(nst,*),c1,rnorm,aengy
  logical          :: afac,solv,dyn
      

!  Purpose:  Active column assembly and solution of equations
!            N.B. interface for psolve same as for frontal solutions

!  Inputs:
!     u(*)      - Current nodal solution
!     a(*)      - Storage for tangent matrix
!     b(*)      - Not used for profile solution
!     dr(*)     - Residual
!     m(*)      - Not used for profile solution
!     xm(*)     - Mass matrix (diagonal only)
!     s(nst,*)  - Element tangent matrix
!     ld(nst)   - Element equation numbers
!     ig(*)     - Pointers to end of profile row/columns
!     idl(*)    - Not used for profile solution
!     nst       - Size of element arrays
!     nrs       -
!     afac      - Flag, factor matrix if true
!     solv      - Flag, perform forward/backward solution if true
!     dyn       - Flag, transient solution if true
!     c1        - Factor to add mass matrix to tangent
!     ipd       - Precision of !real! variables

!  Outputs:
!     dr(*)     - Solution increment
!     rnorm     - Residual norm
!     aengy     - Energy of solution increment
!     ifl       - Not used for profile solution


   include 'maxa.h'

   logical   afl,fa
   integer   n,ne,nep
   integer   ibuf,ihsiz,ihfac,iz
   double precision    dot


   include 'cdata.h'
   include 'fdata.h'

   include 'iofile.h'

   include 'temfl1.h'
   include 'temfl2.h'

   data   fa /.false./

!  Form and assemble matrix

!  Record length factors (may be too long for some machines)

   ihfac = 8 ! record length factor (8)
   ihsiz = maxArray

   if(afac) then
     if(fl(6)) then
       nep = neq + ig(neq)
       afl = .true.
     else
       nep = neq
       afl = .false.
     end if
     if(fl(3).or.fl(4)) then
       ibuf = ig(neq)+neq
       if(fl(3)) then
         ibuf  = ibuf + ig(neq)
         fl(3) = fa
       end if
       if(ibuf.gt.ihsiz) stop 'profile too large'
       iz       = ibuf
       itrec(1) = iz*ihfac
!       open (4,file=tfile(1),status='new',access='direct',form='unformatted',recl=itrec(1))
!       close(4)
       fl(4) = fa
     end if

     call pconsd(a,ibuf,0.0d0)

!    Modify tangent form lumped mass effects

     if(dyn) then
       do n = 1,neq
         a(n)  =  c1*xm(n)
       end do  
     end if

     if(ior.lt.0) then 
       write(*,'(a/a)') '+  Solution status',' '
     end if  

!    Compute and assemble element arrays

     do n = 1,numel
       ne = n
       call formfe(u,dr,.true.,solv,fa,fa,3,ne,ne,1)
       if(ior.lt.0 .and. mod(n,50).eq.0) then
         write(*,'(a,i4,a)') '+  ->',n,' Elements completed.'
       end if  
       call dasbly(s,s,ld,ig,nst,afl,afac,fa,dr,a(nep+1),a(neq+1),a)
     end do  

     rnorm = sqrt(dot(dr,dr,neq))

     call datri(a(nep+1),a(neq+1),a,ig,neq,afl)
     call phstio(4,1,a,ibuf,2,tfile(1),itrec(1))

   end if

   if(solv) then
     if(.not.afac) then
       call phstio(4,1,a,ibuf,1,tfile(1),itrec(1))
     end if  
     do n = 1,nrs
       ne = (n-1)*neq + 1
       call dasol(a(nep+1),a(neq+1),a,dr(ne),ig,neq, aengy)
     end do  
   end if

end
