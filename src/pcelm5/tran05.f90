subroutine tran05(s,cs,sn,nst,ndf,itype)
implicit none
  integer          :: nst,ndf,itype
  double precision :: s(nst,*),cs,sn

!  Itype: 1  Transform matrix s(nst,nst)
!         2  Transform vector s(nst,1)

   integer i,j,nn
   double precision t      

   if(itype.eq.1) then
     do i = 1,nst,ndf
       do j = 1,nst
         t        = s(j,i)*cs - s(j,i+1)*sn
         s(j,i+1) = s(j,i)*sn + s(j,i+1)*cs
         s(j,i  ) = t
       end do ! j
     end do ! i
     nn = nst
   else
     nn = 1
   end if
   do i = 1,nst,ndf
     do j = 1,nn
       t        = s(i,j)*cs - s(i+1,j)*sn
       s(i+1,j) = s(i,j)*sn + s(i+1,j)*cs
       s(i  ,j) = t
     end do ! j
   end do ! i
end
