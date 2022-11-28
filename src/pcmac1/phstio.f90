subroutine phstio(iu,irec,hh,nh,isw,tfile,itrec)
implicit none
  integer          :: iu,irec,nh,isw,itrec
  double precision :: hh(*)
  character*12     :: tfile

!  Purpose: I/O control for history terms on disk

!  Inputs:
!     iu         - Logical unit number for I/O
!     irec       - Record number
!     hh(*)      - History terms to write (isw = 2 or 22 only)
!     nh         - Number of history terms
!     isw        - Switch: =1 or 11 for read; = 2 or 22 for write
!     tfile      - File name for history terms
!     itrec      - Record length of data
!  Outputs:
!     hh(*)      - History terms read (isw = 1 or 11 only)

!  direct access read/write

  !print *, 'RECL =',tfile,' ',itrec
   
  if(isw < 10) then
    open(iu,file=tfile,access='direct',recl=itrec,status='unknown')
    select case (isw)
    case (1)
      read (iu,rec=irec) hh(nh)
    case (2)
      write(iu,rec=irec) hh(nh)
    end select
    close(iu)  
  else ! isw >= 10; e.g. isw == 11 or isw == 22
!    sequential access read/write
    select case (isw)
    case (11)
      read (iu) hh(nh)
    case (22)
      write(iu) hh(nh)
    end select
  end if      
end
