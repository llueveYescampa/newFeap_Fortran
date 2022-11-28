subroutine datri(al,au,ad,jp,neq,flg)
implicit none
  integer          :: jp(*),neq
  double precision :: al(*),au(*),ad(*)
  logical          :: flg

!  Purpose: Triangular decomposition of a matrix stored in profile form

!  Inputs:
!     al(*)    - Lower part of matrix
!     au(*)    - Upper part of matrix
!     ad(*)    - Diagonal of matrix
!     jp(*)    - Pointer to end of row/columns
!     neq      - Number of equations to factor
!     flg      - Flag, equations unsymmetric if true

!  Outputs:
!     al(*)    - Lower part of factored matrix
!     au(*)    - Upper part of factored matrix
!     ad(*)    - Reciprocal diagonal of factored matrix


   integer  i,id,ie,ih,is,idh,ifig, j,jd,jh,jr,jrh
   double precision   dd,dimx,dimn,dfig,tol,daval,dot


   include 'iofile.h'

!  N.B.  tol should be set to approximate half-word precision.

   data tol/0.5d-07/

!  Set initial values for conditioning check

   dimx = 0.0d0
   dimn = 0.0d0
   do j = 1,neq
     dimn = max(dimn,abs(ad(j)))
   end do  
   dfig = 0.0d0

!  Loop through columns to perform triangular decomposition

   jd = 1
   do j = 1,neq
     jr = jd + 1
     jd = jp(j)
     jh = jd - jr
     if(jh.gt.0) then
       is = j - jh
       ie = j - 1

!      If diagonal is zero compute a norm for singularity test

       if(ad(j).eq.0.0d0) then
         call datest(au(jr),jh,daval)
       end if  
       do i = is,ie
         jr = jr + 1
         id = jp(i)
         ih = min(id-jp(i-1),i-is+1)
         if(ih.gt.0) then
           jrh = jr - ih
           idh = id - ih + 1
           au(jr) = au(jr) - dot(au(jrh),al(idh),ih)
           if(flg) al(jr) = al(jr) - dot(al(jrh),au(idh),ih)
         end if
       end do  
     end if

!    Reduce diagonal

     dd = ad(j)
     if(jh.ge.0) then
       jr = jd - jh
       jrh = j - jh - 1
       call dredu(al(jr),au(jr),ad(jrh),jh+1,flg  ,ad(j))

!      Check for possible errors and print warnings

       if(abs(ad(j)).lt.tol*abs(dd)) then
         write(iow,2000) j
       end if  
       if(dd.lt.0.d0.and.ad(j).gt.0.d0) then
         write(iow,2001) j
       end if  
       if(dd.gt.0.d0.and.ad(j).lt.0.d0) then
         write(iow,2001) j
       end if  
       if(ad(j) .eq.  0.d0) then
         write(iow,'(a,i5)') &
              ' **WARNING** Reduced diagonal is zero for equation',j
       end if  
       if(dd.eq.0.d0.and.jh.gt.0) then
         if(abs(ad(j)).lt.tol*daval) then
           write(iow,2003) j
         end if  
       end if
       if(ior.lt.0) then
         if(abs(ad(j)).lt.tol*abs(dd))  then
           write(*,2000) j
         end if  
         if(dd.lt.0.d0.and.ad(j).gt.0.d0) then
           write(*,2001) j
         end if  
         if(dd.gt.0.d0.and.ad(j).lt.0.d0) then
           write(*,2001) j
         end if  
         if(ad(j) .eq.  0.d0) then
           write(*,'(a,i5)') &
              ' **WARNING** Reduced diagonal is zero for equation',j
         end if  
         if(dd.eq.0.d0.and.jh.gt.0) then
           if(abs(ad(j)).lt.tol*daval) then
             write(*,2003) j
           end if  
         end if
       end if
     end if

!    Store reciprocal of diagonal, compute condition checks

     if(ad(j).ne.0.d0) then
       dimx  = max(dimx,abs(ad(j)))
       dimn  = min(dimn,abs(ad(j)))
       dfig  = max(dfig,abs(dd/ad(j)))
       ad(j) = 1.d0/ad(j)
     end if
     if(ior.lt.0 .and. mod(j,50).eq.0) then
       write(*,'(a,i5,a,i5,a)') '+  ->',j,' Equations of',neq,' reduced'
        
     end if
   end do

!  Print conditioning information

   dd = 0.0d0
   if(dimn.ne.0.0d0) then
     dd = dimx/dimn
   end if  
   ifig = log10(dfig) + 0.6
   write(iow,2004) dimx,dimn,dd,ifig
   if(ior.lt.0) then
     write(*,2004) dimx,dimn,dd,ifig
   end if  

!  Formats

2000  format(' **WARNING** Loss of at least 7 digits in reducing', &
       ' diagonal of equation',i5)
     
2001  format(' **WARNING** Sign of diagonal changed when reducing', &
       ' equation',i5)
     

2003  format(' **WARNING** Rank failure for zero unreduced',         &
       ' diagonal in equation',i5)
     
2004  format(' Condition check: D-max',1pe11.4,'; D-min',1pe11.4,     &
       '; Ratio',1pe11.4/' Maximum no. diagonal digits lost:',i3)

      end
