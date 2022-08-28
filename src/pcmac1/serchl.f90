subroutine serchl(fn,id,g0,rsd,u,d,stol,t,nneq)
implicit  none
integer id(*),nneq
double precision fn(nneq,2),g0,rsd(*),u(*),d(*),stol,t(*)
      

!  linear line search for nonlinear problems

   integer j, linmax
   double precision s,step,gamma,g,ga,gb,sa,sb
   double precision gama
   
   integer         numnp,numel,nummat,nen,neq
   common /cdata/  numnp,numel,nummat,nen,neq

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

!  compute step size for line search in direction d

   linmax = 10
   sb     = 0.0
   sa     = 1.0
   s      = 1.0
   g      = gama(fn,id,u,rsd,d,t,s,neq,nneq)

!  find bracket on zero

   if(g*g0.gt.0.0d0) then
     write(ioWrite,'(4x,a)') 'No line search, end points both positive.'
     if(ioRead.lt.0) then
       write(*,'(4x,a)') 'No line search, end points both positive.'
     end if  
   else
     j    = 0
     gb   = g0
     ga   = g
     sb   = 0.0d0
     sa   = 1.0d0
     do while (.true.)
       j    = j + 1
       step = sa - ga*(sa-sb)/(ga-gb)
       g    = gama(fn,id,u,rsd,d,t,step,neq,nneq)
       gb   = 0.5d0*gb
       if (g*ga.lt.0.0d0) then
         sb = sa
         gb = ga
       end if
       sa = step
       ga = g
       write(ioWrite,'(4x,a,i2,a,e12.5,a,e12.5)') &
             'Iter =',j,' Step Size =',step,' Energy =',g
       if(ioRead.lt.0) then
         write(*,'(4x,a,i2,a,e12.5,a,e12.5)') &
             'Iter =',j,' Step Size =',step,' Energy =',g
       end if  
       if (j.lt.linmax) then
         if(abs(g)     .gt. stol*abs(g0)      ) then
           cycle
         end if  
         if(abs(sb-sa) .gt. stol*0.5d0*(sa+sb)) then
           cycle
         end if
       end if
       exit
     end do
     do j = 1, neq
       d(j) = step*d(j)
     end do  
   end if

end
