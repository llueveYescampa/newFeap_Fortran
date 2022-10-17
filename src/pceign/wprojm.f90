subroutine wprojm(a,nn,ia)
implicit  none

integer   i,ia, j, k, n,nn

character ah(2)*1
double precision    a(*)

include 'iofile.h'

data ah(1),ah(2) /'g','h'/

write(ioWrite,'(1x,a,a1)') 'Matrix ',ah(ia) 
if(ioRead.lt.0) then 
  write(*,'(1x,a,a1)') 'Matrix ',ah(ia) 
end if  

i = 1
do n = 1,nn
  j = i + n - 1
  write(ioWrite,'(1p8d10.2)') (a(k),k=i,j)
  if(ioRead.lt.0) then
    write(*,'(1p8d10.2)') (a(k),k=i,j)
  end if  
  i = i + n
end do  

end
