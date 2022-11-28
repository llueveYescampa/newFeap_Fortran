subroutine wprojm(a,nn,ia)
implicit none
  double precision  ::  a(*)
  integer  :: ia,nn

  integer ::  i, j, k, n
  character :: ah(2)*1

  include 'iofile.h'

  data ah(1),ah(2) /'g','h'/

  write(iow,'(1x,a,a1)') 'Matrix ',ah(ia) 
  if(ior.lt.0) then 
    write(*,'(1x,a,a1)') 'Matrix ',ah(ia) 
  end if  

  i = 1
  do n = 1,nn
    j = i + n - 1
    write(iow,'(1p8d10.2)') (a(k),k=i,j)
    if(ior.lt.0) then
      write(*,'(1p8d10.2)') (a(k),k=i,j)
    end if  
    i = i + n
  end do  

end
