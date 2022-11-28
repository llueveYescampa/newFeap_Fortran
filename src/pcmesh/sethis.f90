subroutine sethis(ie,ix,idl,ipd,nen,nen1,numel)
implicit none
   integer :: ipd,nen,nen1,numel,ie(9,*),ix(nen1,*),idl(*)

!  Purpose: Set up history addresses in ix array

!  Inputs:
!     ie(9,*)    - Material set assembly data
!     ix(nen1,*) - Element nodal connection lists
!     idl(*)     - Working vector
!     ipd        - Precision of 'double precision' variables
!     nen        - Number nodes/element
!     nen1       - Dimension of ix array
!     numel      - Number of elements in mesh

!  Outputs:
!     ix(nen1,*) - Element nodal connection lists with history
!                  requirements added.


   integer   ihsiz,ihfac,ihmin, n,nu,nh0,nhinc


   include 'hdatb.h'
   include 'psize.h'
   include 'temfl1.h'
   include 'temfl2.h'
   include 'maxa.h'
   include 'ddata.h'


   nh0 = (ne+ipd-1)/ipd
   nhi = nh0

!  Record length factors

   ihfac  = 8
   ihsiz  = maxa
   ihbuff = ihsiz + 1 - nh0
   hfl    = .true.
   hout   = .false.
   irec   = 0
   nrec   = 1

!  Determine buffer size needed for history terms

   ihmin  = 0
   do nu = 1,numel
     n     = idl(nu)
     ihmin = ihmin + ie(8,ix(nen1,n))
   end do  

!  Set buffer length and record length to minimum possible

   ihbuff = min(ihbuff,ihmin)
   nhf    = nhi + ihbuff - 1

!  Integer values

   ihbuff = ihbuff*ihfac

!  Set history area

   do nu = 1,numel
     n     = idl(nu)
     nhinc = ie(8,ix(nen1,n))
     if(hfl .and. nhinc .ne. 0) then
       itrec(2) = ihbuff
       open(3,file=tfile(2),access='direct', status='new', &
              form='unformatted',recl=itrec(2))
       close(3)
       hfl = .false.
       call pconsd(dm(nhi),ihbuff/ihfac,0.0d0)
     endif
     
!    Integer values
     
     if(nh0+nhinc .gt. ihsiz) then
       call phstio(3,nrec,dm(nhi),ihsiz-nhi+1,2,tfile(2),itrec(2))
       nrec = nrec + 1
       nh0 = nhi
     end if
     ix(nen+1,n) = nh0
     ix(nen+2,n) = nrec
     nh0 = nh0 + nhinc
   end do

!  Check for errors and finish initialization

   if(nrec.gt.numel) then
     write(*,'(a)')' **ERROR** Insufficient storage for history terms'
     call pdefil(tfile,2,2)
     call pstop(-93) ! stop
   else if(nh0.gt.nhi) then
     call phstio(3,nrec,dm(nhi),nhf-nhi+1,2,tfile(2),itrec(2))
   else
     nrec = nrec - 1
   end if
   ne   = (nhf+1)*ipd + 1
   close(3)

end
