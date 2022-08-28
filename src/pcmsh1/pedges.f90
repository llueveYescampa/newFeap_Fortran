subroutine pedges(iii,x,id,ndm,ndf,numnp,prt)
implicit none
integer  iii,ndm,ndf,numnp,id(ndf,numnp)
double precision   x(ndm,numnp)
logical  prt
      
!     Purpose: Set boundary edge for restraints

!     Inputs:
!        iii        - Flag to test change in b.c. inputs
!        x(ndm,*)   - Nodal coordinates of mesh
!        ndm        - Spatial dimension of mesh
!        ndf        - Number dof/node
!        numnp      - Number of nodes in mesh
!        prt        - Flag, output results if true

!     Outputs:
!        id(ndf,*)  - Boundary condition codes


      character yyy*80
      integer   i, j, n,idl(6)
      integer iocheck
      double precision x0,dx

      integer         ioRead,ioWrite
      common /iofile/ ioRead,ioWrite

!     Read input of boundary edge for restraints

      iii = 1
      do while(.true.) 
        if(ioRead.lt.0) then
          write(*,'(a,/3x,a,$)') ' Input: i-dir, x(i), b.code(j),j=1,ndf','>'
        end if  
        call pintio(yyy,10)
        read(yyy,'(i10,f10.0,6i10)',IOSTAT=iocheck) i,x0,idl
        if (iocheck .ne. 0) then
          call pperror('PEDGES',yyy)
          cycle ! go to 100
        else
          if(i.le.0.or.i.gt.ndm) then 
            exit ! go to 4
          end if  
          dx = 1.e-04
          do n = 1,numnp
            dx = max( dx , x(i,n)*1.e-04 )
          end do  
          do n = 1,numnp
            if(x(i,n).ne.-999. .and. abs(x(i,n)-x0).le.dx) then
              do j = 1,ndf
                id(j,n) = max(abs(id(j,n)),abs(idl(j)))
              end do  
            end if
          end do  
          cycle ! go to 100
        end if
      end do
      
      if(prt) then
        call prthed(ioWrite)
        write(ioWrite,2000) (i,i=1,ndf)
        do n = 1,numnp
          do i = 1,ndf
            if(id(i,n).ne.0) then
              write(ioWrite,'(i10,8i8)') n,(id(j,n),j=1,ndf)
              exit 
            end if
          end do  
        end do  
      end if

2000  format(/'    E d g e    N o d a l    B. C.'/        &
             /6x,'node',8(i3,' b.c.')/1x)
     
end
