subroutine param(ctl,ct)
implicit none
character        :: ctl*4
double precision :: ct(3)

!  Purpose: Set appropriate time integration parameters

!  Inputs:
!     ctl       - Character string of integration type
!     ct(3)     - Integration type parameters

!  Outputs:
!     none      - Output through common block


   logical   pcomp

   include 'iofile.h'
   include 'tbeta.h'

!  Set integration parameters

   beta = ct(1)
   gamm = ct(2)

!  SS11 algorithm

   if(pcomp(ctl,'ss11')) then
     nop = 1
     if(beta.le.0.0) beta = 0.50
     write(iow,2001) beta
     if(ior.lt.0) write(*,2001) beta
     theta = beta

!  SS22 algorithm

   else if(pcomp(ctl,'ss22')) then
     nop = 2
     if(beta.le.0.0) beta = 0.50
     if(gamm.le.0.0) gamm = 0.50
     write(iow,2002) beta,gamm
     theta = gamm
     if(ior.lt.0) then
       write(*,2002) beta,gamm
     end if  

!  Newmark algorithm

   else
     nop = 3
     if(beta.le.0.0) beta = 0.25
     if(gamm.le.0.0) gamm = 0.50
     write(iow,2003) beta,gamm
     if(ior.lt.0) write(*,2003) beta,gamm
     theta = beta
   end if

2001  format(' SS-11 Method Parameter: theta = ',f9.4)
2002  format(' SS-22 Method Parameters'/' theta-1 = ',f9.4,' ;  theta-2 = ',f9.4)             
2003  format(' Newmark Method Parameters'/' beta = ',f9.4,' ;  gamma = ',f9.4)
             

end
