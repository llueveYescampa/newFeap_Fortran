subroutine brk4(jx1,jy1, str)   !.... convert integers for tektronix 4012
      integer jx1, jy1
      character*1 str(4)
      
      str(1) = char(32 + jy1/32)
      str(2) = char(96 + mod(jy1,32))
      str(3) = char(32 + jx1/32)
      str(4) = char(64 + mod(jx1,32))

end
