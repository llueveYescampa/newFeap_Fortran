subroutine prefrt(il,idl,ix,maxf,ndf,nen,nen1,numel,numnp)
implicit  none
integer   maxf,ndf,nen,nen1,numel,numnp,il(*),idl(*),ix(nen1,*)

!  Purpose: Prefrontal routine to flag last occurance of nodes

!  Inputs:
!     il(*)      - Working vector
!     idl(*)     - Elemenet elimination order
!     ix(nen1,*) - Element nodal connections
!     ndf        - Number dof/node
!     nen        - Maximum number of nodes/element
!     nen1       - Dimension of ix array
!     numel      - Number of elements in mesh
!     numnp      - Number of nodes in mesh

!  Outputs:
!     ix(nen1,*) - Element nodal connections (tag last occurances)
!     maxf       - Maximum front estimate


   integer   i,ii, jj, n,nu,nowf

   integer         ior,iow
   common /iofile/ ior,iow

!  Preset check array

   do n=1,numnp
     il(n)=0
   end do  

!  Set last occurance of nodes

   do nu=numel,1,-1
     n = idl(nu)
     do i=nen,1,-1
       ii=abs(ix(i,n))
       if((ii.ne.0).and.(il(ii).eq.0)) then
         il(ii)=n
         ii=-ii
       end if
       ix(i,n)=ii
     end do
   end do  

!  Get estimate to maximum frontwith

   maxf=0
   nowf=0
   do nu=1,numel
     n = idl(nu)
     do i=1,nen
       ii=ix(i,n)
       if(ii.ne.0) then
         jj=abs(ii)
         if(il(jj).ne.0) nowf=nowf+ndf
         maxf=max(maxf,nowf)
         il(jj)=0
       end if
     end do  
     do i = 1,nen
       if(ix(i,n).lt.0) then
         nowf = max(0,nowf-ndf)
       end if
     end do
   end do  
   if(ior.lt.0) then
     write(*,'(a,i4,a)')'   -> Maximum Front Estimate =',maxf,' d.o.f.'
   end if  

end
