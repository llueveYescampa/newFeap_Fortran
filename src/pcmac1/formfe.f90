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



   integer        nn,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13
   common /mdata/ nn,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13

   integer        n11a,n11b,n11c,ia
   common /mdat2/ n11a,n11b,n11c,ia(2,11)

   integer        ndf,ndm,nen1,nst
   common /sdata/ ndf,ndm,nen1,nst

   integer        isw,nn1,nn2,nn3
   logical                        afl,bfl,cfl,dfl
   common /xdata/ isw,nn1,nn2,nn3,afl,bfl,cfl,dfl
   integer         maxa
include 'maxa.h'      

   double precision dm
   integer                  m
   common          dm(maxa),m(maxa)

!  form appropriate f.e. array

   afl = af
   bfl = bf
   cfl = cf
   dfl = df
   isw = is
   nn1 = ne1
   nn2 = ne2
   nn3 = ne3
   call pform(dm(nn),dm(n0),dm(n1),m(n2),dm(n3),dm(n4),m(n5),dm(n6), &
                m(n7),dm(n8),m(n9),dm(n11),m(n11c),u,b,ndf,ndm,nen1,nst)
end
