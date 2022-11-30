!|---------------------------------------------------------------------+
!|                                                                     |
!|......      F E A P    - -  A Finite Element Analysis Program for    |
!|                                     Mini and Mainframe Computers    |
!|                                                                     |
!|......    p c F E A P  - -  A Finite Element Analysis Program for    |
!|                                     Personal Computers              |
!|                                                                     |
!|.... A (PC)  Finite Element Analysis Program for solution of general |
!|.... problem classes using the finite element method.   Problem size |
!|.... is controlled by the dimension of blank common and value of max |
!|.... set below.                                                      |
!|                                                                     |
!|.... Programmed by:                                                  |
!|                R. L. Taylor                                         |
!|                Department of Civil Engineering                      |
!|                University of California                             |
!|                Berkeley, California 94720                           |
!|                                                                     |
!|.... Mini Version 2.04 of FEAP - June    1992                        |
!|                Last Updated   - May     1997                        |
!|                Use Skyline Solution system in file: pasolv.for      |
!|                User must supply graphics interface                  |
!|                                                                     |
!|         Uses: Fortran 77 Compiler System                            |
!|                                                                     |
!|.... PC Version 2.04 of pcFEAP - June    1992                        |
!|                Last Updated   - May     1997                        |
!|                Use Frontal Solution system in file: pfsolv.for      |
!|                  or                                                 |
!|                Use Skyline Solution system in file: pasolv.for      |
!|                                                                     |
!|         Uses: Microsoft (R) Fortran Version 3.3x                    |
!|               Graphics interfaced using files: pcdigl.for, pcio.asm |
!|                                                                     |
!|         Uses: Microsoft (R) Fortran Version 5.0 or PowerStation 1.0 |
!|               Graphics interfaced using file: pcfor5.for            |
!|                                                                     |
!|.... (C) Copyright -  Robert  L.  Taylor  -  1985 to 1997  -         |
!|                                                                     |
!|---------------------------------------------------------------------+
program         pcfeap
implicit        none

   include 'cdata.h'
   include 'iofild.h'
   include 'psize.h'
   include 'temfl1.h'
   include 'vdata.h'

!  Parameters control program capacity:
!                 mmax      = mesh size
!                 maxArray  = equation size (in several routines)

   !integer        mmax
   include 'maxa.h'
   !parameter     (mmax = 3*maxArray)

   integer        m
   common         m(3*maxArray)

!  Set version data

   versn(1) = '   pcFEAP   '
   versn(2) = ' -- 2.04 -- '
   versn(3) = '  05/20/97  '

!  Reserve memory size; set default input/output units
!      N.B. maxm is defined in number of INTEGER words.

   maxm = 3*maxArray ! mmax
   ne   = 1
   iodr = 15
   iodw = 16

!  Set precision values: ipd = double; ipr = real; ipi = integer
!      N.B. Precisions are set as multiples of the integer value.
!
!        a. For short integers (integer*2) the values are:
!               ipd = 4, ipr = 2, ipi = 1
!
!        b. For normal integer (integer*4) the values are:
!               ipd = 2, ipr = 1, ipi = 1
!
!  integer*2

!  ipd  = 4
!  ipr  = 2

!  integer

   ipd  = 2
   ! ipr  = 0
   ipi  = 1

!  Start system if necessary

   call pstart

!  Open files; erase scratch files if they exist; start execution

   call filnam
   call pdefil(tfile,1,4)
   call pcontr

!  Close input and output files; destroy temporary disk files

   close(iodr)
   close(iodw)
   call pdefil(tfile,1,4)
   call pstop(0)
end
