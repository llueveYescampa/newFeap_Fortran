subroutine pline(x,ix,ic,numnp,numel,ndm,nen1,nen,ct,isw)
implicit none
integer numnp,numel,ndm,nen1,nen,ix(nen1,*),ic(numnp,*)
logical isw
double precision    x(ndm,*),ct

!  Purpose: Plot mesh or outline

!  Inputs:
!    x(ndm,*)   - Nodal coordinates of mesh
!    ix(nen1,*) - Element connection list
!    ic(numnp,*)- Working array
!    numnp      - Number of nodes in mesh
!    numel      - Number of elements in mesh
!    ndm        - Spatial dimension of mesh
!    nen1       - Dimension of ix array
!    nen        - Number of nodes/element, maximum
!    ct         - Color selection for plots
!    isw        - Flag, plot mesh if true, otherwise outline

!  Outputs:
!    none

   logical   ifl,iend, test
   double precision    x3
   integer   i,ii,ij, j,jj, k, n,n1,n2,ni
   
   integer   iplt(8)
   data iplt/5,2,6,3,7,4,8,1/
   
!  Initialize connection array
   call pconsi(ic,numnp*4,0)
!  Loop through elements to set up list
   do n = 1,numel
     jj = abs(ct)
     ii = abs(ix(nen1,n))
     if(jj.eq.0 .or. jj.eq.ii) then
       i = 1
       ii = abs(ix(i,n))
       do ij = 1,8
         j = iplt(ij)
         if(j.le.nen.and.ix(j,n).ne.0) then
           jj = abs(ix(j,n))
           if(jj.ne.ii) then
             n1 = min(ii,jj)
             n2 = max(ii,jj)
             do k = 1,4
               if(ic(n1,k).eq.0) then
                 ic(n1,k) = n2
                 ii=jj
                 exit
               else if(ic(n1,k).eq.n2) then
                 ic(n1,k) = -n2
                 ii = jj
                 exit
               end if
             end do
           end if
         end if
       end do
     end if
   end do
!  change signs to permit mesh plot
   if(isw) then
     do n = 1,numnp
       do i = 1,4
         ic(n,i) = abs(ic(n,i))
       end do
     end do
   end if
!  plot outline of part with continuous lines
   x3 = 0.0
   do ni = 1,numnp
     iend = .true.
     do n = 1,numnp
       ifl = .true.
       n1 = n
       test = .true.
       do while (test)
         test = .false.
         do i = 1,4
           if(ic(n1,i) .gt. 0) then
             iend = .false.
             if(ndm.ge.3) then
               x3 = x(3,n1)
             end if  
             if(ifl) then
               call plotl(x(1,n1),x(2,n1),x3,1)
             end if  
             ifl = .false.
             n2 = ic(n1,i)
             ic(n1,i) = -n2
             if(ndm.ge.3) then 
               x3 = x(3,n2)
             end if  
             call plotl(x(1,n2),x(2,n2),x3,2)
             n1 = n2
             test = .true.
             exit
           else if (ic(n1,i) .eq. 0) then 
             exit ! go to 303
           else  
             cycle 
           end if
         end do  
       end do
     end do  
     if(iend) then
        return
      end if
   end do   
   return
   end
