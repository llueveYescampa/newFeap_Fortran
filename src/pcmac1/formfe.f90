subroutine formfe(u,b,af,bf,cf,df,is,ne1,ne2,ne3)
implicit  none
integer                           is,ne1,ne2,ne3
double precision         u(*),b(*)
logical              af,bf,cf,df

!  Purpose: Form finite element arrays as required

!  Inputs:
!     u(*)      - Current solution state
!     af        - Flag, assemble matrix if true
!     bf        - Flag, assemble vector if true
!     cf        - Flag, matrix unsymmetric if true
!     df        - Flag, assemble vector in full form if true
!     is        - Switch to control actions to be taken
!     ne1       - First element to process
!     ne2       - Last  element to process
!     ne3       - Increment to ne1

!  Outputs:
!     b(*)      - Element vector

   include 'maxa.h'      

   include 'mdata.h'

   include 'mdat2.h'

   include 'sdata.h'
   
   include 'xdata.h'

   include 'ddata.h'

!  form appropriate f.e. array

   afl = af
   bfl = bf
   cfl = cf
   dfl = df
   isw = is
   nn1 = ne1
   nn2 = ne2
   nn3 = ne3
   call pform(dm(nn),dm(n0),dm(n1),im(n2),dm(n3),dm(n4),im(n5),dm(n6), &
                im(n7),dm(n8),im(n9),dm(n11),im(n11c),u,b,ndf,ndm,nen1,nst)
end
