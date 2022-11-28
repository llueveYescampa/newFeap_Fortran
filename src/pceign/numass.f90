subroutine numass(b,neq,mq)
implicit none
  double precision :: b(*)
  integer          ::  neq,mq

!  Purpose: Count number of non-zero entries in mass

!  Inputs
!     b(*)     - Diaonal mass
!     neq      - Length of b array

!  Outputs
!     nq       - Number of nonzero entries in b array

   integer  n,nn   
   include 'iofile.h'

   nn = 0
   do n = 1,neq
     if(b(n).ne.0.0d0) then 
       nn = nn + 1
     end if  
   end do  

   if(nn.lt.mq) then 
     write(iow,'(1x,a,i4,a)') &
     'Subspace reduced to',nn,' by number of nonzero lumped mass terms'
   end if  
   if(ior.lt.0.and.nn.lt.mq) then 
     write(*,'(1x,a,i4,a)')   &
     'Subspace reduced to',nn,' by number of nonzero lumped mass terms'
   end if  
   mq = min(mq,nn)

end
