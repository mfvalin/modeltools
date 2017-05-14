#define MARKEOS 3
#define MARKEOI 2
#define MARKEOF 1
#define MARKFRM 21845
#define MBUFSZ 90

      program analyze
      integer :: STATUS
      status = 0
      do while(status == 0)
        call TTTPRCS( STATUS )
      enddo
      end program

      SUBROUTINE TTTPRCS( STATUS )
      implicit none
      INTEGER  STATUS
      integer, save :: MCUNIT=0
      integer, dimension(1024) :: BPAIRS
      integer :: NPAIRS, OPCODE, option, opval
      integer, external :: TTTBPR

      STATUS = TTTBPR(MCUNIT, BPAIRS, 2)
      if(status .ne. 0) return

      IF(ISHFT(BPAIRS(1), -15) .EQ. 0) THEN   ! 4 byte instrunction
        print 100, 'I4B',BPAIRS(1),BPAIRS(2)
100     format(A4,2Z10.8)
        return
      endif
      IF(ISHFT(BPAIRS(1), -14) .EQ. 2) THEN   ! illegal instruction
        status = -1
        return
      endif
      NPAIRS = ISHFT(IAND(255 , BPAIRS(1)+1),-1) -1
      STATUS = TTTBPR(MCUNIT, BPAIRS(3), NPAIRS)
      OPCODE = IAND(ISHFT(BPAIRS(1),-8), 63)       ! multibyte opcode
      if(OPCODE == 35) then                        ! optn
        option = IAND(ISHFT(BPAIRS(2),-8), 255)
        opval  = BPAIRS(3) - 32768
        print 102, 'OPT',OPCODE,(NPAIRS+1)*2,option,opval,opval
102     format(A4,2I4,2I10,Z10.8)
      else
        print 101, 'IMB',OPCODE,(NPAIRS+1)*2,BPAIRS(2:min(NPAIRS+2,16))
101     format(A4,2I4,16Z10.8)
      endif
      return
      
      end SUBROUTINE TTTPRCS

      integer FUNCTION TTTBPR( MCUNIT, BPAIRS, NPAIRS)
      implicit none
      INTEGER  MCUNIT, NPAIRS, BPAIRS(NPAIRS)
      EXTERNAL WAREAD, WACLOS, DFLERMS, POPPAR
!
      INTEGER, save :: BUFFER(MBUFSZ), TWOBYT(2*MBUFSZ)
      INTEGER, save :: bufindex = 0
      INTEGER, save :: ninbuf = 0
      INTEGER, save :: ADDRESS = 1
      integer :: demande, servi, NPB, status, i
      integer, save :: unit = 0
      integer, external :: fnom

      if(MCUNIT == 0) then
        status = fnom(MCUNIT,'metacod','RND', 0)
      endif

      TTTBPR  = 0
      DEMANDE = NPAIRS
      SERVI   = 0

 10   continue

      NPB = MIN(ninbuf,DEMANDE)
      DO I = 1,NPB
         BPAIRS(SERVI+I) = TWOBYT(bufindex+I-1)
         BPAIRS(SERVI+I) = IAND(bpairs(servi+i), 65535)    ! on fait le mask de 16 bits
      ENDDO

      SERVI       = SERVI    + NPB
      bufindex    = bufindex + NPB
      DEMANDE     = DEMANDE  - NPB
      ninbuf      = ninbuf   - NPB
      if(DEMANDE <= 0) return   ! done

      CALL WAREAD(MCUNIT, BUFFER, ADDRESS, MBUFSZ)                 ! read metacode buffer
      ADDRESS = ADDRESS + MBUFSZ
      CALL IIPAK(TWOBYT, BUFFER, MBUFSZ*2, 1, -16, 0, 2)        ! unpack into 16 bit tokens

      bufindex = 3
      ninbuf   = TWOBYT(1)
      if(TWOBYT(2) .ne. MARKFRM) then
        TTTBPR = -1
        return
      endif
      if(TWOBYT(2) == MARKEOI) then
        TTTBPR = MARKEOI
        return
      endif
      if(TWOBYT(2) == MARKEOF) then
        TTTBPR = MARKEOF
        return
      endif
      if(TWOBYT(2) == MARKEOS) then
        TTTBPR = MARKEOS
        return
      endif
      IF(DEMANDE > 0)  GO TO 10

      RETURN
      END
