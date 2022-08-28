integer*2 function vclrwk()  !.... clear tektronix 4012 device

      write(*,'(2a1)') char(27),char(12)
      vclrwk = 0
      
end
      
!      write(*,1000) char(27),char(12)
!1000  format(2a1)

