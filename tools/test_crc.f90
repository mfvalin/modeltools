program testcrc
  implicit none
  integer *4, dimension(32) :: iarray
  integer *8, dimension(64) :: larray
  integer *4 :: i,mask
  integer *8 :: lmask
  integer *4, external :: f_update_crc_ne
  integer *4 :: crc16, crc24, crc32

  crc16=0
  crc24=0
  crc32=0
  mask=1
  do i=1,32
  iarray(i)=mask
  mask=ishft(mask,1)
  enddo
  lmask=1
  do i=1,64
  larray(i)=lmask
  lmask=lmask+lmask
  enddo
  crc16=f_update_crc_ne(crc16,16,iarray,4,32,0)
  crc24=f_update_crc_ne(crc24,24,iarray,4,32,0)
  crc32=f_update_crc_ne(crc32,32,iarray,4,32,0)
  print 102,iarray(1),iarray(32)
  print *,"crc16=0000f39a, crc24=006b7cac, crc32=e3b736e4 is the correct answer"
  print 101,' crc16=',crc16,'crc24=',crc24,'crc32=',crc32
  crc16=f_update_crc_ne(crc16,16,larray,8,64,2)
  crc24=f_update_crc_ne(crc24,24,larray,8,64,2)
  crc32=f_update_crc_ne(crc32,32,larray,8,64,2)
  print 103,larray(1),larray(64)
  print *,"crc16=00000fc6, crc24=0061f67c, crc32=5560ac91 is the correct answer"
  print 101,' crc16=',crc16,'crc24=',crc24,'crc32=',crc32
  101 format(3(A,Z8.8,2X))
  102 format(8Z10.8)
  103 format(4Z20.16)
  stop
end
