subroutine subsp(b,v,g,h,d,dp,dtol,p,z,nf,nev,neq,shift,tol,prt,its)
implicit none
  logical          :: prt
  integer          :: nf,nev,neq,its
  double precision :: b(*),v(neq,*),g(*),h(*),d(*),dp(*),dtol(*),p(nev,*),&
                      z(neq,*),shift,tol
                 
!  Purpose: Subspace iteration to extract lowest nf eigenpairs

!  Inputs:
!     v(*)     - Diagonal !mass! matrix
!     g(*)     - Subspace projection of stiffness matrix storage
!     h(*)     - Subspace projection of !mass! matrix storage
!     dp(*)    - Storage for previous eigenvalues
!     dtol(*)  - Storage for tolerance of computed eigenvalues
!     p(nev,*) - Storage for subspace eigenvectors
!     z(neq,*) - Storage for iteration vectors
!     nf       - Number of eigenvalues to converge
!     nev      - Size of subspace problem
!     neq      - Size of finite element equations
!     shift    - Value of shift 
!     tol      - Solution tolerance on eigenvalues
!     prt      - Flag, output each iteration result
!     its      - Maximum number of subspace iterations to perform

!  Outputs:
!     v(neq,*) - Eigenvectors
!     d(*)     - Eigenvalaues


   logical          :: conv
   integer          :: i,it,itt=0,itlim, j, k, n,nmas
   double precision :: dm,tolmx,told

   include 'iofild.h'
   include 'iofile.h'

!  Compute initial iteration vectors

   call pconsd(v,nev*neq,0.0d0)
   nmas = 0

!  Count number of nonzero masses

   do n = 1,neq
     if(b(n) .ne. 0.0d0) then
       nmas = nmas + 1
     end if
   end do

   nmas = nmas/nev
   i = 0
   j = 1
   do n = 1,neq
     dm = b(n)
     if(dm.ne.0.0d0) then
       v(n,j) = dm
       i = i + 1
       if(mod(i,nmas).eq.0) then
         j = j + 1
       end if  
       j = min(j,nev)
     end if
   end do     
   do i = 1,nev
     dp(i)   = 0.0d0
     dtol(i) = 1.0d0
     call scalev(v(1,i),neq)
   end do

!  Compute new vectors and project 'a' onto 'g'

   told  = tol
   conv  = .false.
   itlim = its
   if(nev.eq.nf) then
     itlim = 1
   end if  
   
!   do it = 1,itlim
   it = 1
   do while(it .le. itlim .and. (.not. conv))
     itt = it

!    Project 'b' matrix to form 'h' and compute 'z' vectors

     call sprojb(b,v,z,h,neq,nev)

!    Project 'a' matrix to form 'g'

     call sproja(v,z,g,neq,nev)

!    Solve reduced eigenproblem 'g*p = h*p*d'

     call geig(g,h,d,p,v,nev,prt)

!    Check for convergence

     tolmx = 0.0d0
     do n = 1,nev
       if(d(n).ne.0.0d0) then 
         dtol(n) = abs((d(n)-dp(n))/d(n))
       end if  
       dp(n) = d(n)
       if(n.le.nf) then 
         tolmx = max(tolmx,dtol(n))
       end if  
     end do  

     if(prt) then
       write(iow,'(/5x,a,i4/(4d20.8))') &
          'Current reciprocal shifted eigenvalues, iteration',it,(d(n),n=1,nev)
       if(ior.lt.0) then
         write(*,'(/5x,a,i4/(4d20.8))') &
          'Current reciprocal shifted eigenvalues, iteration',it,(d(n),n=1,nev)
       end if  
       if(itlim.gt.1) then 
         write(iow,'( 5x,a/(4d20.8))') 'Current residuals', (dtol(n),n=1,nev)
       end if  
       if(ior.lt.0 .and. itlim .gt. 1) then
         write(*,'( 5x,a/(4d20.8))') 'Current residuals', (dtol(n),n=1,nev)
       end if
     else
       if(ior.lt.0) then 
         write(*,'(a,i3,a,1p1e11.4)') '+  Iteration',it,' Max tol =',tolmx
       end if  
     end if

!    Tolerance check

     do n = 1,nf
       if(dtol(n).gt.told) then 
         exit
       end if  
     end do  
     if (n .gt. nf) then    ! completo el ciclo
       conv = .true.
     end if  

!    Divide eigenvectors by eigenvalue to prevent overflows

     do i = 1,nev
       dm = d(i)
       if(p(i,i).lt.-0.00001d0) then 
         dm = -dm
       end if  
       do j = 1,nev
         p(j,i) = p(j,i)/dm
       end do  
     end do  

!    Compute new iteration vector 'u' from 'z'

     do i = 1,neq
       do j = 1,nev
         v(i,j) = 0.0d0
         do k = 1,nev
           v(i,j) = v(i,j) + z(i,k)*p(k,j)
         end do  
       end do  
     end do  
     
     it=it+1
   end do  

!  Scale vectors to have maximum element of 1.0

   do n = 1,nev
     d(n)  = 1.0/d(n) + shift
     dp(n) = sqrt(abs(d(n)))
     call scalev(v(1,n),neq)
   end do  

   write(iow,'(/5x,a,i4/(4d20.8))') &
       'Solution  for  eigenvalues, iteration',itt,(d(n),n=1,nev)
   write(iow,'( 5x,a/(4d20.8))') 'Square root of eigenvalues', (dp(n),n=1,nev)
   if(itt.gt.1) then 
     write(iow,'( 5x,a/(4d20.8))') 'Current residuals', (dtol(n),n=1,nev)
   end if  
   if(ior.lt.0) then
     write(*,'(/5x,a,i4/(4d20.8))') &
       'Solution  for  eigenvalues, iteration',itt,(d(n),n=1,nev)
     write(*,'( 5x,a/(4d20.8))')  'Square root of eigenvalues', (dp(n),n=1,nev)
     if(itt.gt.1) then 
       write(*,'( 5x,a/(4d20.8))') 'Current residuals', (dtol(n),n=1,nev)
     end if  
   end if
end
