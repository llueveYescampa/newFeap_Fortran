subroutine matl02(d,it,ib)
implicit none
double precision  d(15)
integer it,ib

!  Purpose:  Input material set parameters for elements

!  Outputs:
!        d(*)  -  Constants for constitutive equation.
!        it    -  Geometry type (0=plane, 1=axisymmetric)
!        ib    -  Formulation type (ib = 0 b-bar,
!                   ib = 1 normal b matrix).
!        nh1   -  number of history variables needed for each
!                   stress point (i.e., at gauss points in element).


   integer       nh
   integer iocheck
   double precision        ee,xnu
   character*24  wa(2)

!  Parameter specification for FEAP materials

   integer         nh1,nh2
   common /hdata/  nh1,nh2

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   character       yyy*80
   common /ydata/  yyy


   data wa/' P l a n e   S t r a i n',' A x i s y m m e t r i c'/
   
!  Input material parameters
   do while (.true.)
     if(ioRead.lt.0) then 
       write(*,3000)
     end if
     call pintio(yyy,10)
     read(yyy,'(i10,f10.0)',IOSTAT=iocheck) it,d(4)
     if (iocheck .eq. 0) then
       it  = max(0,min(1,it))
       ib  = 0
       nh  = 9
       nh1 = 36
       if(ioRead.lt.0) then
         write(*,3001)
       end if  
       call pintio(yyy,10)
       read(yyy,'(8f10.0)',IOSTAT=iocheck) ee,xnu,d(11),d(12),d(13)
       if (iocheck .eq. 0) then
         d(1)    = ee/(1. - 2.*xnu)/3.0d0
         d(2)    = ee/(1.+xnu)/2.
         if(d(11).ne.0.0d0) then
           write(ioWrite,2000) wa(it+1),ee,xnu,d(4),d(11),d(12),d(13)
           if(ioRead.lt.0) then 
             write(*,2000) wa(it+1),ee,xnu,d(4),d(11),d(12),d(13)
           end if
         else
           write(ioWrite,2001) wa(it+1),ee,xnu,d(4)
           if(ioRead.lt.0) then
             write(*,2001) wa(it+1),ee,xnu,d(4)
           end if  
         end if
         exit
       else
         call pperror('PCELM2',yyy)
         cycle
       end if  
       exit 
     else
       call pperror('PCELM2',yyy)
       cycle
     end if
   end do
   return

!     Formats


2000  format(2x,'E l a s t i c / P l a s t i c   M a t e r i a l'/ &
             2x,'v o n  M i s e s  Y i e l d - ',a24/              &
        10x,'Youngs Modulus (E)',e17.5/                            &
        10x,'Poisson Ratio (nu)',e17.5/                            &
        10x,'Mass Density (rho)',e17.5/                            &
        10x,'Yield Stress      ',e17.5/                            &
         8x,'Linear hardening moduli'/                             &
        10x,'Isotropic Hardening',e16.5/                           &
        10x,'Kinematic Hardening',e16.5/)

2001  format(2x,'E l a s t i c   M a t e r i a l'/ &
             2x,a24,'   A n a l y s i s'/          &
        10x,'Youngs Modulus (E)',e17.5/            &
        10x,'Poisson Ratio (nu)',e17.5/            &
        10x,'Mass Density (rho)',e17.5/)

3000  format(' Input: it, rho'/                                 &
        4x,'it = 0: Plane'/4x,'it = 1: Axisymmetric'/' >',$)

3001  format(' Input (Elastic/Plastic) : E, nu, Y, H-iso, H-kine'/ &
             '       (Linear  Elastic) : E, nu'/3x,'>',$)

end
