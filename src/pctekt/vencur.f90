integer*2 function vencur()   !.... enter text mode with cursor

 write(*,'(7a1)') char(31),char(27),char(50),char(27),char(91),char(50),char(74)
 vencur = 0
  
end

