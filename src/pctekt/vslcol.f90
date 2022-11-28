integer*2 function vslcol(icol)   !.... set line type

  integer iln
  
  iln = 0
  if(icol.eq.4) then
     iln = 3
  end if   
  write(*,'(2a1)') char(27),char(iln+96)
  vslcol = 0
end
