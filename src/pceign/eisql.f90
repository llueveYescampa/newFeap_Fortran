subroutine eisql(a,d,e,z,n,ierr)
implicit none
integer  n, ierr
double precision   a(*),d(*),e(*),z(n,n)

!  Purpose: Solve standard eigen problem

!  Inputs:
!     a(*)     - Matrix of standard problem
!     e(*)     - Working vector
!     n        - Size of eigen problem

!  Inputs:
!     d(*)     - Eigenvalues
!     z(n,*)   - Eigenvalues of standard problem
!     ierr     - Error counter: zero if no errors

   integer  i,ii, j,jp1, k, l,l1, m,mml, n2
   double precision   b, c, f, g, h, hh, p, r, s
   double precision   scale,machep
   logical test

   data machep/0.222044605d-15/
 
!  Eispac ql algorithm: Adapted from 'tred2' and 'tql2'

   n2 = 0
   do i = 1,n
     do j = 1,i
       n2     = n2 + 1
       z(i,j) = a(n2)
     end do
   end do  

   if(n.ne.1) then
     n2 = n + 2
     do ii = 2,n
       i     = n2 - ii
       l     = i - 1
       h     = 0.0d0
       scale = 0.0d0
       if(l.lt.2) then
         e(i) = z(i,l)
         d(i) = h
         cycle
       end if  
       do k = 1,l
         scale = scale + abs(z(i,k))
       end do  
       if(scale.eq.0.0d0) then 
         e(i) = z(i,l)
         d(i) = h
         cycle
       end if  
       do k = 1,l
         z(i,k) = z(i,k)/scale
         h      = h + z(i,k)*z(i,k)
       end do  
       f      = z(i,l)
       g      = -sign(sqrt(h),f)
       e(i)   = scale*g
       h      = h - f*g
       z(i,l) = f - g
       f      = 0.0d0
       do j = 1,l
         z(j,i) = z(i,j)/h
         g      = 0.0d0
         do k = 1,j
           g = g + z(j,k)*z(i,k)
         end do  
         jp1 = j + 1
         if(l.ge.jp1) then 
           do k = jp1,l
             g = g + z(k,j)*z(i,k)
           end do  
         end if  
         e(j) = g/h
         f    = f + e(j)*z(i,j)
       end do  
       hh = f/(h+h)
       do j = 1,l
         f    = z(i,j)
         g    = e(j) - hh*f
         e(j) = g
         do k = 1,j
           z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
         end do  
       end do  
       d(i) = h
     end do
   end if

!  Set transformation array for ql

   d(1)   = z(1,1)
   z(1,1) = 1.0d0
   e(1)   = 0.0d0
   ierr   = 0
   
   if(n.eq.1) then 
     return
   end if  
   
   do i = 2,n
     l = i - 1
     if(d(i).ne.0.0d0) then
       do j = 1,l
         g = 0.0d0
         do k = 1,l
           g = g + z(i,k)*z(k,j)
         end do  
         do k = 1,l
           z(k,j) = z(k,j) - g*z(k,i)
         end do  
       end do  
     end if
     d(i)   = z(i,i)
     z(i,i) = 1.0d0
     do j = 1,l
       z(i,j) = 0.0d0
       z(j,i) = 0.0d0
     end do
   end do  

!  Begin 'ql' algorithm on tridagonal matrix now stored in 'd' and 'e

   do i = 2,n
     e(i-1) = e(i)
   end do  
   f    = 0.0d0
   b    = 0.0d0
   e(n) = 0.0d0
   test = .true.
   
!   do l = 1,n
   l=1
   do while (l.le.n)
     if (test) then
       j = 0
       h = machep*(abs(d(l)) + abs(e(l)))
       if(b.lt.h) then
         b = h
       end if  
       do m = l,n
         if(abs(e(m)).le.b) then
           exit
         end if  
       end do  
       if(m.eq.l) then
         d(l) = d(l) + f
         l=l+1
         cycle
       end if  
     end if

     if(j.eq.30) then 
       ierr = l
       return
     end if  
     j    = j + 1
     l1   = l + 1
     g    = d(l)
     p    = (d(l1)-g)/(e(l)+e(l))
     r    = sqrt(p*p+1.0d0)
     d(l) = e(l)/(p+sign(r,p))
     h    = g - d(l)
     do i = l1,n
       d(i) = d(i) - h
     end do  
     f   = f + h
     p   = d(m)
     c   = 1.0d0
     s   = 0.0d0
     mml = m - l
     do ii = 1,mml
       i = m - ii
       g = c*e(i)
       h = c*p
       if(abs(p).ge.abs(e(i))) then
         c      = e(i)/p
         r      = sqrt(c*c+1.0d0)
         e(i+1) = s*p*r
         s      = c/r
         c      = 1.0d0/r
       else
         c      = p/e(i)
         r      = sqrt(c*c+1.0d0)
         e(i+1) = s*e(i)*r
         s      = 1.0d0/r
         c      = c*s
       end if
       p      = c*d(i) - s*g
       d(i+1) = h + s*(c*g + s*d(i))
       do k = 1,n
         h        = z(k,i+1)
         z(k,i+1) = s*z(k,i) + c*h
         z(k,i  ) = c*z(k,i) - s*h
       end do  
     end do  
     e(l) = s*p
     d(l) = c*p
     if(abs(e(l)).gt.b) then
       test = .false.
       cycle
     end if
     test = .true.
     d(l) = d(l) + f
     l=l+1
   end do  

   do ii = 2,n
     i = ii - 1
     k = i
     p = d(i)
     do j = ii,n
       if(abs(d(j)).gt.abs(p)) then
         k = j
         p = d(j)
       end if
     end do  
     if(k.ne.i) then
       d(k) = d(i)
       d(i) = p
       do j = 1,n
         p      = z(j,i)
         z(j,i) = z(j,k)
         z(j,k) = p
       end do  
     end if
   end do  

end
