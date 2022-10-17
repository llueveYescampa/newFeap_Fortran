subroutine pplotf(x,ix,b,lci,ct,ndf,ndm,nen1)
implicit none
integer  ix(*),ndf,ndm,nen1
character lci*4
double precision x(ndm,*),b(*),ct(2)

!  Purpose: Plot control subroutine for feap

!  Inputs:
!     x(ndm,*)   - Nodal coordinates for mesh
!     ix(nen1,*) - Element nodal connection lists
!     b(*)       - Current solution vecotr
!     lci        - Character plot option
!     ct(2)      - Parameters for plot
!     ndf        - Number dof/node
!     ndm        - Spatial dimension of mesh
!     nen1       - Dimension of ix array

!  Outputs:
!     none       - Graphics to plot area

   logical   pcomp,oflg
   ! integer*2 status,vslcol,coli
   integer   i,ic,n1
   double precision    c

   include 'maxa.h'      

   include 'adata.h'

   include 'cdata.h'
   double precision dr(1)
   integer                il(1)
   equivalence(aa(1), dr(1))
   equivalence(aa(1), il(1))

   double precision, dimension (:), allocatable :: dm
   allocate ( dm(maxa) )

!  Open kernel system and plot mesh or outline of parts
   call pdevop()
   call frame(x,ndm,numnp)
!  Plot mesh or outline of parts
   oflg = .not.pcomp(lci,'outl')
   if(ct(1).ne.0.0.or.pcomp(lci,'eigv')) then
     ic = 2
     if(ct(1).eq.0.0) ct(1) = 1.0
   else
     ic = 1
   end if
   c  = 0.0
   n1 = 2*ndm*numnp
   do i = ic,1,-1
     ! coli   = 8 - 2*i
     ! status = vslcol(coli)
     if(pcomp(lci,'eigv')) then
       call pdefm(x,dm,c,ndm,ndf,numnp, dr)
     else
       call pdefm(x, b,c,ndm,ndf,numnp, dr)
     end if
     call pline(dr,ix,il(n1),numnp,numel,ndm,nen1,nen,ct(2),oflg)
     c  = ct(1)
   end do  
!  Close plot
   call pdevcl()
   deallocate (dm) 
end
