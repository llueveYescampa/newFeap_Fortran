subroutine update(id,fn,u,v,du,nneq,fdyn,isw)
implicit  none
integer   id(*),nneq,isw
double precision    fn(nneq,2),u(nneq,*),v(nneq,*),du(*)
logical   fdyn

!   Purpose: Update displacements (and velocities and accelerations)

!   Inputs:
!     id(*)      - Equation numbers for each dof
!     fn(*)      - Forced solution values at restrained dof
!     u(*)       - Solution vector and incrmeents: isw = 1
!     du(*)      - Solution increment: isw = 2
!     nneq       - Number of total parameters in u and v
!     fdyn       - Flag, transient solution if true
!     isw        - Switch

!   Outputs:
!     u(*)       - Initialization for   t_n+1: isw = 1
!                  Solution iterates at t_n+1: isw = 2
!
!     v(*)       - Initialized rates for t_n+1: isw = 1
!                  Rate iterates at      t_n+1: isw = 2


    integer   j,n
    double precision ur1,ur2,dot


    include 'iofile.h'

   include 'tbeta.h'

   include 'tdata.h'

!   update solution vectors to begin a step

    if(isw.eq.1) then

      ur1 = sqrt(dot(v(1,1),v(1,1),nneq))
      ur2 = sqrt(dot(v(1,2),v(1,2),nneq))
      write(ioWrite,2000) ur1,ur2
      if(ioRead.lt.0) then
        write(*,2000) ur1,ur2
      end if  

      do n = 1,nneq

!       SS11 algorithm u=u-bar; v1=v-bar; v2=a=abar; v3=u;

        if(nop.eq.1) then
          u(n,1) = v(n,3)
          v(n,1) = 0.0d0
          v(n,2) = 0.0d0

!       SS22 algorithm u=u-bar; v1=v-bar; v2=a=abar; v3=u; v4=v

        elseif(nop.eq.2) then
          u(n,1) = v(n,3) + c3*v(n,4)
          v(n,1) = v(n,4)
          v(n,2) = 0.0d0
          v(n,3) = v(n,3) + dt*v(n,4)

!       Newmark algorithm

        else
          ur2    = - c6*v(n,1) + c3*v(n,2)
          v(n,1) =   c4*v(n,1) + c5*v(n,2)
          v(n,2) =   ur2
        end if
      end do  
    else

!     Update displacement and its increments within time step

      ur2 = 1. - theta
      do n = 1,nneq
        j = id(n)
        if (j.gt.0) then

!         For active degrees-of-freedom compute values from solution

          u(n,1) = du(j) + u(n,1)
          u(n,2) = du(j) + u(n,2)
          u(n,3) = du(j)
        else

!         For fixed degrees-of-freedom compute values from forced inputs

          ur1    = theta*fn(n,2) + ur2*fn(n,1)
          u(n,3) = ur1    - u(n,1)
          u(n,2) = u(n,2) + u(n,3)
          u(n,1) = ur1
        end if
      end do  

!     For time dependent solutions update rate terms

      if(fdyn) then

        do n = 1,nneq
          v(n,1) = v(n,1) + c2*u(n,3)
          v(n,2) = v(n,2) + c1*u(n,3)

!         SS11 algorithm

          if(nop.eq.1) then
            v(n,3) = v(n,3) + c3*u(n,3)

!         SS22 algorithm

          else if(nop.eq.2) then
            v(n,3) = v(n,3) + c4*u(n,3)
            v(n,4) = v(n,4) + c5*u(n,3)
          end if
        end do  
      end if
    end if

2000  format('   N o r m s   f o r   D y n a m i c s'/  &
         10x,'Velocity:',e13.5,' Acceleration:',e13.5)

end
