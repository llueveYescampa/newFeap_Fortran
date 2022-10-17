subroutine pian03(d,ul,s,p,nst,ndf,isw)
implicit  none
integer                    nst,ndf,isw
double precision  d(*), ul(ndf,*),s(nst,*),p(*)

!  Pian-Sumihara stiffness matrix developed explicitly

   integer   i,i1,j,j1
   double precision    r1,r2,vol,d11,d12,d33,hx,hy,h44,h45,h55
   double precision    bd11,bd12,bd13,bd21,bd22,bd23

   double precision    rr(5),x1(4),x2(4),y1(4),y2(4)
   double precision    ax(4),bx(4),cx(4),ay(4),by(4),cy(4)

   include 'elcom3.h'

!  1.) set up stress interpolants for the 4-5 term

   r1 = xj1/xj0
   r2 = xj2/xj0
   a1(1) = xs*xs
   a1(2) = ys*ys
   a1(3) = xs*ys
   a2(1) = xt*xt
   a2(2) = yt*yt
   a2(3) = xt*yt

!  2.) set up shape function coefficients - jacobian weighted

   ax(1) = -yt + ys
   ax(2) =  yt + ys
   ax(3) = -ax(1)
   ax(4) = -ax(2)
   bx(1) = -yh - ys
   bx(2) = -bx(1)
   bx(3) =  yh - ys
   bx(4) = -bx(3)
   cx(1) =  yt + yh
   cx(2) = -yt + yh
   cx(3) = -cx(2)
   cx(4) = -cx(1)
   ay(1) =  xt - xs
   ay(2) = -xt - xs
   ay(3) = -ay(1)
   ay(4) = -ay(2)
   by(1) =  xh + xs
   by(2) = -by(1)
   by(3) = -xh + xs
   by(4) = -by(3)
   cy(1) = -xt - xh
   cy(2) =  xt - xh
   cy(3) = -cy(2)
   cy(4) = -cy(1)

!  3.) compute volume and stabilization h-array

   vol = 4.*xj0
   d11 = d(1)/vol
   d12 = d(2)/vol
   d33 = d(3)/vol
   hy  = vol*3.
   hx  = hy*d(5)
   h44 = hx*(1.-r2*r2/3.)*(a1(1)+a1(2))**2
   h55 = hx*(1.-r1*r1/3.)*(a2(1)+a2(2))**2
   h45 =-(r1*r2/3.)*(hx*(xs*xt+ys*yt)**2+d(6)*hy*(ys*xt-xs*yt)**2)

!  4.) Invert stabilization h-array

   hx  = h44*h55 - h45*h45
   hy  = h55/hx
   h55 = h44/hx
   h45 =-h45/hx
   h44 = hy

!  5.) Compute the current stress parameters

   call pconsd(rr,5,0.0d0)
   do j = 1,4
     hx = cx(j) - r2*ax(j)
     hy = cy(j) - r2*ay(j)
     x1(j) = a1(1)*hx + a1(3)*hy
     x2(j) = a1(2)*hy + a1(3)*hx
     hx = bx(j) - r1*ax(j)
     hy = by(j) - r1*ay(j)
     y1(j) = a2(1)*hx + a2(3)*hy
     y2(j) = a2(2)*hy + a2(3)*hx
     rr(1) = rr(1) + ax(j)*ul(1,j)
     rr(2) = rr(2) + ay(j)*ul(2,j)
     rr(3) = rr(3) + ay(j)*ul(1,j) + ax(j)*ul(2,j)

!    (stabilization terms)

     rr(4) = rr(4) + x1(j)*ul(1,j) + x2(j)*ul(2,j)
     rr(5) = rr(5) + y1(j)*ul(1,j) + y2(j)*ul(2,j)
   end do ! j
   beta(1) = d11*rr(1) + d12*rr(2)
   beta(2) = d12*rr(1) + d11*rr(2)
   beta(3) = d33*rr(3)
   beta(4) = (h44*rr(4) + h45*rr(5))*3.
   beta(5) = (h45*rr(4) + h55*rr(5))*3.

!  6.) Form stiffness matrix for 1-pt (constant) terms

   if(isw.eq.3) then
     d11 = d11*d(11)
     d12 = d12*d(11)
     d33 = d33*d(11)
     i1 = 1
     do i = 1,2
       bd11 = ax(i)*d11
       bd12 = ax(i)*d12
       bd13 = ay(i)*d33
       bd21 = ay(i)*d12
       bd22 = ay(i)*d11
       bd23 = ax(i)*d33
       j1   = i1
       do j = i,2
         s(i1  ,j1  ) = bd11*ax(j) + bd13*ay(j)
         s(i1  ,j1+1) = bd12*ay(j) + bd13*ax(j)
         s(i1+1,j1  ) = bd21*ax(j) + bd23*ay(j)
         s(i1+1,j1+1) = bd22*ay(j) + bd23*ax(j)
         j1 = j1 + ndf
       end do ! j
       i1 = i1 + ndf
     end do ! i

!    7.) Copy other parts from computed terms

     i1 = ndf + ndf
     do i = 1,i1
       do j = i,i1
         s(i   ,j+i1) =-s(i,j)
         s(i+i1,j+i1) = s(i,j)
       end do ! j
     end do ! i
     j1 = i1 + ndf
     s(2    ,i1+1) =-s(2,1)
     s(ndf+2,j1+1) =-s(ndf+2,ndf+1)
     s(ndf+1,i1+1) =-s(1,ndf+1)
     s(ndf+1,i1+2) =-s(2,ndf+1)
     s(ndf+2,i1+1) =-s(1,ndf+2)
     s(ndf+2,i1+2) =-s(2,ndf+2)

!    8.) Add stabilization matrix

     h44 = h44*d(11)
     h45 = h45*d(11)
     h55 = h55*d(11)
     j1 = 1
     do j = 1,4
       bd11 = h44*x1(j) + h45*y1(j)
       bd12 = h44*x2(j) + h45*y2(j)
       bd21 = h45*x1(j) + h55*y1(j)
       bd22 = h45*x2(j) + h55*y2(j)
       i1 = 1
       do i = 1,j
         s(i1  ,j1  ) = s(i1  ,j1  ) + x1(i)*bd11 + y1(i)*bd21
         s(i1  ,j1+1) = s(i1  ,j1+1) + x1(i)*bd12 + y1(i)*bd22
         s(i1+1,j1  ) = s(i1+1,j1  ) + x2(i)*bd11 + y2(i)*bd21
         s(i1+1,j1+1) = s(i1+1,j1+1) + x2(i)*bd12 + y2(i)*bd22
         i1 = i1 + ndf
       end do ! i
       j1 = j1 + ndf
     end do ! j
   end if

!  9.) compute the residual force vector

   if(mod(isw,3).eq.0) then

!    a.) compute the constant part

     do i = 1,5
       beta(i) = beta(i)*d(11)
     end do ! i
     do i = 1,2
       p(2*i-1) = -(ax(i)*beta(1) + ay(i)*beta(3))
       p(2*i+3) = -p(2*i-1)
       p(2*i  ) = -(ay(i)*beta(2) + ax(i)*beta(3))
       p(2*i+4) = -p(2*i  )
     end do ! i

!    b.) compute the stabilization part

     do i = 1,4
       p(2*i-1) = p(2*i-1) - (x1(i)*beta(4) + y1(i)*beta(5))/3.0
       p(2*i  ) = p(2*i  ) - (x2(i)*beta(4) + y2(i)*beta(5))/3.0
     end do ! i
   end if

end
