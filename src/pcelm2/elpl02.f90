subroutine elpl02(d,eps,epln,alph,ep,sig,ul,ndf,ib)
implicit none
integer ndf,ib
double precision d(*),eps(4),epln(4),alph(4),ep,sig(5),ul(ndf,*)

!  Purpose: Elasto-plastic model with isotropic / kinematic hardening

!  Inputs: - for t-n
!     d      - array of material constants
!     eps    - strains (at t-n+1)
!     sig    - stresses
!     alph   - back stress
!     ep     - effective plastic strain

!  Outputs: - at t-n+1
!     ep     - effective plastic strain
!     sig    - stress tensor
!     alph   - back stress tensor
!     ad     - !tangent! matrix


   integer i,j
   double precision oneg,twog,bulk,elam,treps,tt,beta,gam1,gam2,gamn,press, &
                    radius,psitr,psi(4),en(4)


   include 'eldata.h'

   include 'elcom2.h'


!  Set parameters

   data tt/.816496580927726D+00/

   oneg = d(2)
   twog = d(2) + d(2)
   bulk = d(1)
   elam = bulk - twog/3.d0

!  Compute the trial deviator stress

   treps  = (eps(1) + eps(2) + eps(3))/3.0d0
   do i = 1,3
     eps(i) = eps(i) - treps
   end do ! i
   do i = 1,4
     sig(i) = twog*(eps(i) - epln(i))
     psi(i) =       sig(i) - alph(i)
   end do ! i

!  Set up elastic tangent

   do i = 1,3
     do j = 1,3
       ad(i,j) = elam
     end do ! j
     ad(i,i) = elam + twog
   end do ! i
   ad(4,4) = oneg

!  Compute the yield state - J2d

   if(d(11).gt.0.0d0) then
     radius = tt*(d(11) + d(12)*ep)
     psitr  = sqrt(psi(1)**2+psi(2)**2+psi(3)**2+2.d0*psi(4)**2)
     sig(5) = psitr/tt/d(11)

!    Compute plasticity solution state

     if (psitr.gt.radius) then
       beta   = 1.d0/(1.d0 + (d(12) + d(13))/(3.d0*oneg))
       gamn   = beta*(psitr - radius)/twog
       gamn   = gamn*(1.d0 -1.d-10)
       ep     = ep + tt*gamn

!      Gam1 ensures stress is slightly outside yield surface.

       gam1   = gamn*twog
       gam2   = (d(13)+d(13))*gamn/3.d0
       sig(5) = (psitr - gam1 - gam2)/tt/d(11)
       do i = 1,4
         en(i)   = psi(i)/psitr
         sig(i)  =  sig(i)  - gam1*en(i)
         alph(i) =  alph(i) + gam2*en(i)
         epln(i) =  epln(i) + gamn*en(i)
       end do ! i

!      Plastic modification for tangent

       gam1 = gam1/psitr*twog
       gam2 = gam1/3.0d0

!      I-dev part

       do i = 1,3
         do j = 1,3
           ad(i,j) = ad(i,j) + gam2
         end do ! j
         ad(i,i) = ad(i,i) - gam1
       end do ! i
       ad(4,4) = ad(4,4) - 0.5*gam1

!      n x n part

       gam1 = gam1 - twog*beta
       do i = 1,4
         gam2    = gam1*en(i)
         do j = 1,4
           ad(i,j) = ad(i,j) + gam2*en(j)
         end do ! j
       end do ! i
     end if
   end if

!  Compute trace of strain and add pressure term

   if(ib.eq.0) then
     treps = 0.0d0
     do i = 1,4
       treps = treps + g(1,i)*ul(1,i) + g(2,i)*ul(2,i)
     end do ! i
   else
     treps = treps*3.d0
   end if
   press = bulk*treps
   do i = 1,3
     sig(i) = sig(i) + press
   end do ! i

end
