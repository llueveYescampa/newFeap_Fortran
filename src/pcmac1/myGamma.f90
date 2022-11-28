double precision function myGamma(fn,id,u,dr,du,t,s,nqe,nneq)
implicit none
  integer          :: id(*), nqe,nneq
  double precision :: fn(nneq,2),u(*),dr(*),du(*),t(*),s

!  Purpose:  Compute energy for line search step

!  Inputs:
!     fn(*)     - Force vector for active dof's
!     id(*)     - Equation numbers for each dof
!     u(*)      - Currenet solution
!     dr(*)     - Residual
!     du(*)     - Current increment to solution
!     t(*)      - Working space
!     nqe       - Number equations
!     nneq      - Total dof in problem

!  Outputs:
!     myGamma     - Value of energy for step


   integer j, n, nneq2
   double precision  ctem,db,dot


   include 'cdata.h'
   include 'fdata.h'
   include 'mdata.h'
   include 'ndata.h'
   include 'prlod.h'
   include 'tbeta.h'
   include 'tdata.h'
   include 'ddata.h'
   
   logical   fa,tr
   data fa /.false./
   data tr /.true./

!  get a search displacement

   nneq2 = nneq + nneq
   ctem  = theta
   theta = s*theta
   do n = 1,nneq
     j = id(n)
     if(j.gt.0) then
       db         = s*du(j)
       t(n)       = u(n)      + db
       t(n+nneq ) = u(n+nneq) + db
       t(n+nneq2) = db
     else
       db = theta*fn(n,2) + (1.0-theta)*fn(n,1)
       t(n+nneq2) = db - u(n)
       t(n+nneq ) = u(n+nneq) - u(n) + db
       t(n)       = db
     end if
   end do  

!  compute a residual

   call pload(id,fn,dr,nneq,dm(nl),dm(nw))
   call formfe(t,dr,fa,tr,fa,fa,6,1,numel,1)
   theta = ctem

!  update the residual for lumped mass inertial effects

   if(fl(9)) then
     ctem = s*c1
     do n = 1,nneq
       j = id(n)
       if(j.gt.0) then
         dr(j) = dr(j) - dm(nl+j-1)*(dm(nw+n-1) + ctem*du(j))
       end if
     end do  
   end if

!  compute the value of gamma

   myGamma = dot (du,dr,nqe)

end
