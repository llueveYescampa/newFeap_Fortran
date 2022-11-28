subroutine pinitc(engy,rnmx,shift,tol,dt,prop,ttim,npld)
implicit none
  double precision :: engy,rnmx,shift,tol,dt,prop,ttim
  integer          :: npld

!  Purpose: Set initial values of solution parameters

!  Outputs:
!     engy       - Solution step energy         ( 0.0d0)
!     rnmx       - First iterate residual norm  ( 0.0d0)
!     whift      - Shift to eigenproblem        ( 0.0d0)
!     tol        - Solution tolerance           (1.d-12)
!     dt         - Solution time step size      ( 0.0d0)
!     prop       - Proportional loading         ( 1.0d0)
!     ttim       - Solution time                ( 0.0d0)
!     npld       - Number of proportional loads (  0   )

   integer i

   include 'fdata.h'
   include 'tbeta.h'

!  set initial values of parameters

   nop   = 3
   theta = 1.0
   npld  = 0
   engy  = 0.0
   rnmx  = 0.0
   shift = 0.0
   tol   = 1.e-12
   dt    = 0.0
   prop  = 1.0
   ttim  = 0.0
   do i = 1,7
     fl(i+4) = .false.
     fl(i)   = .true.
   end do
      
end
