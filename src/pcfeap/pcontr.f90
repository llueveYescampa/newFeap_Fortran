subroutine pcontr()
implicit  none

!  Control program for feap

   logical   prt,pcomp
   integer   i,iii, n14,nad
   integer iocheck
   character titl(20)*4! ,yyy*80 ! why not use ydata.h for the yyy variable?

   include 'maxa.h'     
   include 'bdata.h'
   include 'cdata.h'
   include 'ddata.h'
   include 'iofile.h'
   include 'iofild.h'
   include 'mdata.h'
   include 'mdat2.h'
   include 'sdata.h'
   include 'vdata.h'
   include 'ydata.h'


!   Set parameters for rotation dof

    do i = 1,4
      ia(1,i) = 1
      ia(2,i) = 2
    end do  

!   Set files back to default values

    do while(.true.)
      ioRead = iodr
      ioWrite = iodw
      
!     Read a card and compare first 4 columns with macro list
      
      read(ioRead,'(20a4)',IOSTAT=iocheck) (titl(i),i=1,20)
      if (iocheck .eq. -1) then
!        EOF was found
         call pend('pcontr')
      else if (iocheck .ne. 0) then
!        Error was found
         call pperror('PCONTR',yyy)
         return
      end if   
      if(pcomp(titl(1),'feap')) then
!       Read and print control information
        do i = 1,20
           head(i) = titl(i)
        end do
        call pintio(yyy,10)
!        read(yyy,'(8i10)',err=600) numnp,numel,nummat,ndm,ndf,nen,nad
         read(yyy,'(8i10)',IOSTAT=iocheck) numnp,numel,nummat,ndm,ndf,nen,nad
         if (iocheck.ne.0) then
!          Error was found
           call pperror('PCONTR',yyy)
           return
         end if
        write(ioWrite,2000) head,versn,numnp,numel,nummat,ndm,ndf,nen,nad
        
!       Set pointers for allocation of data arrays
        
        nen1 = nen + 4
        nst  = nen*ndf + nad
        call psetm(nn, numnp*max(ndm,ndf,2),'d')
        call psetm(nn, 5*nen*ndf,           'd')
        call psetm(n0, nen*ndm,             'd')
        call psetm(n1, nen,                 'd')
        call psetm(n2, nst,                 'i')
        call psetm(n3, nst,                 'd')
        call psetm(n4, nst*nst,             'd')
        call psetm(n5, nummat*9,            'i')
        call psetm(n6, nummat*18,           'd')
        call psetm(n7, ndf*numnp,           'i')
        call psetm(n8, ndm*numnp,           'd')
        call psetm(n9, nen1*numel,          'i')
        call psetm(n10,numnp*ndf,           'd')
        call psetm(n11,numnp,               'd')
        call psetm(n11a,nen,                'd')
        call psetm(n11b,numnp,              'd')
        call psetm(n11c,max(numel,nst),     'i')
        call psetm(n12,ndf*numnp,           'i')
        
!       Call mesh input subroutine to read and print all mesh data
        
        iii = 0
        prt = .true.
        call pmesh(im(n2),im(n5),dm(n6),im(n7),dm(n8),im(n9),dm(n10),dm(n11), &
                   ndf,ndm,nen1,iii,prt)
        cycle
      else if(pcomp(titl(1),'inte') .or. pcomp(titl(1),'macr')) then
        if(pcomp(titl(1),'inte')) then
!         Set files for interactive macro execution
          ioRead = -iodr
        end if 
!       Compute profile

        call profil(im(n12),im(n11c),im(n7),im(n9),ndf,nen1)
        call psetm(n13,numnp*ndf,   'd')
        call psetm(n14,3*numnp*ndf, 'd')

!       Set up stress history addresses

        call sethis(im(n5),im(n9),im(n11c),ipd,nen,nen1,numel)

!       Zero the initial force and solution vectors

        call pconsd(dm(n13),  ndf*numnp,0.0d0)
        call pconsd(dm(n14),3*ndf*numnp,0.0d0)

!       Macro module for establishing solution algorithm

        call pmacr(im(n2),dm(n4),im(n5),dm(n6),im(n7),dm(n8),im(n9),dm(n10), &
              dm(n11),im(n12),dm(n13),dm(n14),dm(1),ndf,ndm,nen1,nst,prt)
        cycle
      else if(pcomp(titl(1),'stop')) then 
        return
      end if
    end do  
    
!   Input/output formats
2000  format(1x,20a4//5x,'VERSION :',3a12//        &
      5x,'Number of nodal points       =',i6/      &
      5x,'Number of elements           =',i6/      &
      5x,'Number of material sets      =',i6/      &
      5x,'Dimension of coordinate space=',i6/      &
      5x,'Degree of freedoms/node      =',i6/      &
      5x,'Nodes per element (maximum)  =',i6/      &
      5x,'Extra d.o.f. to element      =',i6)
end
