      program test
      call varspec2
      end
 	SUBROUTINE varspec2
#if ! defined(VARSPEC2_CODE)
*
*     CALL SEQUENCE...
*      varspec  IN SP2OUT SP1OUT XOUT INPUT OUTPUT [ -a VAL1 -b VAL2 \
*                                                    -c VAL3 -d VAL4 -j VAL5 -k VAL6 ]
*
*
*     DESCRIPTION...
*      VARSPEC - COMPUTE SPECTRUM FROM GRIDED VALUES
*
*     AUTHOR - B.Denis/J.Cote,  May 2002
*        REV - B.Denis/B.Dugas, Feb 2004
*        REV - B.Denis/B.Dugas, Jun 2005
*        REV - B.Denis          Dec. 2005 : 1-detrend routine modified to allow 1D detrending
*                                           2-bugfix for the input array maping for the DFT case. 
*                                           3-REAL*8 consistency in SPECT_* routine
*     
*     CONTACT: BERTRAND DENIS, PH.D.
*              METEOROLOGICAL SERVICE OF CANADA
*              EMAIL:bertrand.denis@ec.gc.ca
*              TEL:514-421-7264
*-------------------------------------------------------------------------------------
*
*     LAST REVISION: $Header: /data/cvs/RMNLIB/utils/rdiag/lspgm/varspec.ptn,v 4.3 2002/12/06 16:10:35 dugas Exp $
*
*
*     PURPOSE - COMPUTES SPECTRAL VARIANCE ARRAY AND POWER SPECTRUM FROM
*               GRIDED DATA. DEPENDING ON THE NATURE OF THE PERIODICITY OF
*               THE 2D GRID, THE USER MAY HAVE THE CHOICE OF 2 TRIGOMETRIC
*               TRANSFORMS:   THE DCT (DISCRETE COSINE  TRANSFORM) 
*                          OR THE DFT (DISCRETE FOURIER TRANSFORM).
*   
*               REFERENCE:'SPECTRAL DECOMPOSTION OF TWO-DIMENSIONAL ATMOSPHERIC 
*                          FIELDS ON LIMITED-AREA DOMAINS USING 
*                          THE DISCRETE COSINE TRANSFORM (DCT)'
*                          BY B. DENIS, J. COTE AND R. LAPRISE 
*                          IN MONTHLY WEATHER REVIEW, 2002, VOL.130,NO.7
*
*     INPUT FILE...
*      IN = INPUT ARRAYS (TYPE 'GRID' OR 'SUBA') FILE
*                                      
*     OUTPUT FILES...                        
*      SP2OUT = 2D ARRAY CONTAINING THE VARIANCE ASSOCIATED
*               WITH EACH SPECTRAL COEFFICIENT. 
*      SP1OUT = 1D ARRAY (TYPE 'ZONL') CONTAINING THE VARIANCE COMPUTED
*               BY SUMMING VARIANCES FROM SP1OUT OVER ANNULUS WAVEBANDS
*               THIS GIVES THE STANDARD POWER SPECTRUM.
*      XOUT   = WAVELENGTH ASSOCIATED TO EACH COMPONENT.
*                                      
*      NOTE THAT THE FIRST ELEMENT OF SP1OUT AND SP2OUT IS
*      THE SQUARE OF THE FIELD AVERAGE.             
*                                      
*                                      
*     PARAMETRES...                   
*      CAS   (a) = TYPE OF TRANSFORM.            
*                                      
*       1 => DCT: RECOMMENDED FOR NON-PERIODIC DATA. THIS IS THE
*                 SO-CALLED DISCRETE COSINE TRANSFORM.
*                                      
*       2 => DFT: STANDARD PERIODIC DISCRETE FOURIER TRANSFORM
*                 NOT RECOMMENDED FOR NON-PERIODIC DATA: 
*                 USE CAS=1 INSTEAD.
*                 NB: IT IS POSSBILE TO APPLY DETRENDING TO MAKE
*                 THE DATA PERIODIC BEFORE THE APPLICATION OF THE DFT.
*                 SEE THE 'DET' OPTION BELOW.
*                                      
*       NB: ( CAS=1 ) THE DIMENSIONS OF THE INPUT DATA MUST BE *EVEN*
*
*           ( CAS=2 ) THE DIMENSIONS OF THE INPUT DATA MUST BE
*                     *ODD* AND THE FIELD WIL BE CONSIDERED PERIODIC,
*                     I.E. IN(1,J) = IN(NI+1,J) AND IN(I,1) = (I,NJ+1).            
*                                      
*                                      
*      LOGV  (b) = TO GET OR NOT THE LOG BASE 10 OF VARIANCE COEFF. IN SP2OUT.
*                                      
*             0=>  LOG BASE 10 IS *NOT* COMPUTED
*             1=>  LOG BASE 10 IS COMPUTED (DEFAULT)
*                                      
*      CLUST (c) = ONLY USED WHEN CAS=1 TO CLUSTER HALF AND FULL
*                  WAVENUMBERS TOGETHER TO GIVE SIMILAR POWER SPECTRA AS
*                  PRODUCED BY CAS=2
*                                      
*             0=>  NO CLUSTERING. A NI X NI GRID WOULD GIVE SPECTRUM
*                  GOING FROM 0 TO NI-1 ELEMENT.        
*                                      
*             1=>  CLUSTER WAVENUMBERS 1 & 1 1/2, 2 & 2 1/2 & 3 & 3 1/2 ....
*                  - IN THIS CASE VARIANCE RELATED TO WAVENUMBER 1/2 WOULD
*                    BE MISSING.                
*                                      
*             2=>  CLUSTER WAVENUMBERS 1/2 & 1,  1 1/2 @ 2,  2 1/2 & 3  ....
*                  - IN THIS CASE VARIANCE RELATED TO THE LAST WAVENUMBER WOULD
*                    BE MISSING. *** THIS IS THE DEFAULT *** 
*
*             3=>  CLUSTER :  0.5*N1/2 + N + 0.5*N3/2               
*                  THIS IS A MOVING AVERAGE.                                        
*                                      
*      DX    (d) = GRID SPACING IN KILOMETER. THIS IS OPTIONAL: IT USED TO
*                  PRODUCE FILE XOUT THAT CONTAINS A VECTOR OF THE 
*                  WAVELENGTH ASSOCIATED TO EACH COMPONANT. THIS IS USEFUL
*                  FOR PLOTING PURPOSES.
*
*      DET   (j) = FOR DETRENDING THE FIELDS BEFORE APPLYING THE DFT 
*                  WHEN THE DFT IS USED.
*
*            (0) = NO DETRENDING APPLYIED
*
*            (1) = DETRENDING APPLYIED (DEFAULT FOR THE DFT; NOT AVAILABLE FOR DCT)
*
*      AVG1D (k) = TO COMPUTE SPECTRA BY DOING 2D SPECTRAL TRANSFORMS
*                  OR 1D SPECTRAL TRANSFORMS (IN X OR IN Y).
*                  THIS LAST OPTION IS TO MIMIC IN-SITU MEASURES TAKEN BY 
*                  AIRFLIGHTS. FIRSTLY, 1D SPECTRA ARE COMPUTED FOR EVERY ROW 
*                  OR COLUMN.SECONDLY, THE AVERAGE OF THOSE 1-D SPECTRA ARE
*                   COMPUTED. THIS METHODOLOGY IS INSPIRED FROM THIS PAPER:
*                  'EVALUATING MESOSCALE NWP MODELS USING KINETIC ENERGY SPECTRA'
*                  BY W.C. SKAMAROCK. MWR, DECEMBRE 2004.
*
*            (0)= DO 2D SPECTRAL TRANSFORM (DEFAULT).
*                 (DO NOT AVERAGE 1D SPECTRA)
*
*            (1)= DO AVERAGE OF 1D WEST-EAST SPECTRA 
*
*            (2)= DO AVERAGE OF 1D SOUTH-NORTH SPECTRA
*
*     READ(5,5000) CAS,LOGV,CLUST,DX,DET,AVG1D
*5000 FORMAT(10X,3I5,E5.0,2I5)
*
*                                      
*     EXAMPLE OF INPUT CARD...              
*                                      
* VARSPEC     1    1    2  45.   1    0
*
*        0123456789012345678901234567890123456789012345678901234567890
*        1         2         3         4         5      
*
*
*     EXIT CONDITIONS...
*       0 NO PROBLEMS
*       1 PROBLEM WITH CAS VALUE
*       2 PROBLEM WITH LOGV VALUE
*       3 PROBLEM WITH CLUST VALUE
*       4 PROBLEM WITH DET VALUE
*       5 PROBLEM WITH DET VALUE
*       6 UNABLE TO READ FIRST FIELD
*       7 MORE THAN ONE ARRAY SIZE IN INPUT FILE
*       8 EVEN DIMENSIONS WITH STANDARD PERIODIC FFTs
*       9 ODD DIMENSIONS WITH DISCRETE COSINE TRANSFORMS
*       10 UNABLE TO READ INPUT DIRECTIVES

*    +     (IN,      SP2OUT,      SP1OUT,      XOUT,      INPUT,      OUTPUT,
*    +TAPE1=IN,TAPE2=SP2OUT,TAPE3=SP1OUT,TAPE4=XOUT,TAPE5=INPUT,TAPE6=OUTPUT)
*     ---------------------------------------------------------------

*     $Log: varspec.ptn,v $
*     Revision 4.3  2002/12/06 16:10:35  dugas
*     Forcer les sortie a 32 bits (sans compaction)
*
*     Revision 4.2  2002/02/02 00:32:18  dugas
*     Ajouter l'option CLUST=3 (B.Denis)
*
*     Revision 4.1  2000/06/08 17:46:38  armnrbd
*     Modifier le nom des variables d'entree.
*     Sauver les indicateurs de temps et de niveau a la sortie.
*
*     Revision 4.0  2000/05/31 16:04:06  armnrbd
*     Version initiale.
*

*----------------------------------------------------------------------------
      IMPLICIT      none

      INTEGER       HEAD

*BD<< taille_entete doit normalement etre definit au pre-processing.
*BD   Ici on le fixe temporairement
*BD==      PARAMETER   ( HEAD = taille_entete )
      PARAMETER   ( HEAD = 28 )
      CHARACTER*4   ZONL,LIGNE*80
      LOGICAL       OK,REGROUPE,SKIP,MOY
      LOGICAL       LDCT,LDFT
      LOGICAL       DO_DETREND,DO_AVG1D,WEST_EAST
      LOGICAL       LAST_WAVENUMBER
      REAL*8        DX

      INTEGER       NI,NJ,IJ,IJT,NID,NJD,NFACT,MMAX,NMAX
      INTEGER       KBUF(HEAD),MAXW,LEN,KPAK,NFF
      INTEGER       NIP,NJP,NIT,NJT,IWAY,IAXE,NX,NY,I,J,L
      INTEGER       NCOEF,NWAVE,NRECS,IO,MAXI,MAXL,IER
      INTEGER       NWDS,CAS,LOGV,CLUST,DET,AVG1D


      REAL*8        AX(500*500),AY(500*500)
      REAL*8        U1(1), U2(1), U3(1)   ! dummies to pass as unused arguments to TRANS1D

!       POINTER     ( PA,         A(1) ),( PB,   B(1) )
!       POINTER     ( PR1,       R1(1) ),( PR2, R2(1) )
!       POINTER     ( PR3,       R3(1) )
!       POINTER     ( PB,         B(1) )
!       POINTER     ( PVA, VARIANCE(1) )
!       POINTER     ( PSP,  SPECTRE(1) )
!       POINTER     ( PIB,     IBUF(8) )
!       POINTER     ( PJB,     JBUF(8) )
      REAL*8,  dimension(:), pointer :: A,B,VARIANCE,SPECTRE
      integer, dimension(:), pointer :: IBUF,JBUF

      CHARACTER     NOMPRG*256
      COMMON       /PROGNAM/ NOMPRG

      INTEGER       GETSIZ
      LOGICAL       RPBLOC,SETIO64
      EXTERNAL      RPBLOC,JCLPNT,GETLIGN,BURNF,XIT
      EXTERNAL      GETFLD2,NGFFT,TRANS1D,SPECT1,SPECT2
      EXTERNAL      PUTFLD2,GETSIZ,SETIO64
      EXTERNAL      CMPLBL

      DATA          ZONL /'ZONL'/
*===================================================================

      NOMPRG =
     +'$SOURCE: /DATA/CVS/RMNLIB/UTILS/RDIAG/LSPGM/VARSPEC.PTN,V $'

*-----------------------------------------------------------------------

      NFF = 6
      CALL JCLPNT( NFF, 1,2,3,4, 5,6 )

***    SETUP FOR 64-BIT I/O.

      OK = SETIO64(.TRUE.)


***    SET THE DEFAULTS FIRST

      CAS   = 1           ! DCT
      LOGV  = 1           ! VARIANCE 2D FIELD IS EXPRESSED IN LOG 
      CLUST = 2           ! CLUSTER WAVENUMBERS 1/2 & 1,  1 1/2 @ 2, ...
      DX    = 1.          ! GRID SPACING = 1 KILOMETER.
      DET   = 1           ! APPLY DETRENDING FOR THE DFT CASE
      AVG1D = 0           ! SPECTRA ARE NOT AVERAGES OF 1D SPECTRA
                          ! BUT ARE 2D SPECTRA.



***    READ INPUT AND OUTPUT GRID PARAMETERS 

      IF (RPBLOC( ' ',LIGNE))                                  THEN

          OK = RPBLOC('A',LIGNE)
          IF (OK) READ(LIGNE,0005,ERR=900,IOSTAT=IO) CAS
          OK = RPBLOC('B',LIGNE)
          IF (OK) READ(LIGNE,0005,ERR=900,IOSTAT=IO) LOGV
          OK = RPBLOC('C',LIGNE)
          IF (OK) READ(LIGNE,0005,ERR=900,IOSTAT=IO) CLUST
          OK = RPBLOC('D',LIGNE)
          IF (OK) READ(LIGNE,0010,ERR=900,IOSTAT=IO) DX
          OK = RPBLOC('J',LIGNE)
          IF (OK) READ(LIGNE,0005,ERR=900,IOSTAT=IO) DET
          OK = RPBLOC('K',LIGNE)
          IF (OK) READ(LIGNE,0005,ERR=900,IOSTAT=IO) AVG1D

      ELSE

          CALL GETLIGN( 5, LIGNE,80, OK )

          IF (.NOT.OK)                                         THEN
              GOTO 901
          ELSE
              READ( LIGNE, 5000, ERR=900,END=901,IOSTAT=IO ) 
     +                           CAS,LOGV,CLUST,DX,DET,AVG1D
              CALL BURNF
          END IF  

      END IF


      WRITE(6,6001) CAS 

      LDCT   =  ( CAS .EQ. 1 )  ! LDCT  => DCT
      LDFT    = ( CAS .EQ. 2 )  ! LDFT  => DFT

      IF ( LDCT )                                              THEN
          WRITE(6,*) 'USING DISCRETE COSINE TRANSFORM  (DCT)'

      ELSE IF ( LDFT  )                                        THEN 
          WRITE(6,*) 'USING DISCRETE FOURIER TRANSFORM (DFT)'

      ELSE 
          WRITE(6,*) 'PROBLEM TH THE CHOICE OF TRANSFORM. ',
     +               'CAS MUST BE 1 (DCT) OR 2 (DFT).'
          CALL                                    XIT(' VARSPEC',-1 )

      END IF

      WRITE(6,6002) LOGV

      IF (LOGV.EQ.0)                                           THEN
          WRITE(6,*) 'NOT COMPUTING THE LOG BASE 10 FOR THE 2D OUTPUT'
      ELSE IF (LOGV.EQ.1)                                      THEN
          WRITE(6,*) 'COMPUTING THE LOG BASE 10 FOR THE 2D OUTPUT'
      ELSE
          WRITE(6,*) 'PROBLEM WITH THE CHOICE OF LOGV. MUST BE 0 OR 1'
          CALL                                     XIT(' VARSPEC',-2 )
      END IF
      
      IF (CAS.EQ.1)                                            THEN

          WRITE(6,6003) CLUST

***       DEFAULT IS CLUST = 2


          IF (CLUST.EQ.0)                                      THEN
              WRITE(6,*) 'NO CLUSTERING OF HALF AND FULL WAVENUMBERS'
              REGROUPE = .FALSE.
              SKIP=      .FALSE.
              MOY=       .FALSE.
          ELSE IF (CLUST.EQ.1)                                 THEN
              WRITE(6,*) 'CLUSTERING WAVE # 1 & 1 1/2, 2 & 2 1/2, ..'
              REGROUPE = .TRUE.
              SKIP =     .TRUE.
              MOY=      .FALSE.
          ELSE IF (CLUST.EQ.2)                                 THEN
              WRITE(6,*) 'CLUSTERING WAVE # 1/2 & 1,  1 1/2 @ 2, ..'
              REGROUPE = .TRUE.
              SKIP =    .FALSE.
              MOY=      .FALSE.
          ELSE IF (CLUST.EQ.3)                                 THEN
              WRITE(6,*) 'CLUSTERING WAVE 0.5*N1/2 + N + N3/2, ..'
              REGROUPE = .TRUE.
              SKIP =    .FALSE.
              MOY=      .TRUE.
          ELSE 
              WRITE(6,*) 'PROBLEM WITH THE CHOICE OF CLUST. ',
     +                   'MUST BE 0,1OR 3'
              CALL                                 XIT(' VARSPEC',-3 )

          END IF

      END IF

***    WRITE OUT THE GRID SPACING GIVEN FOR PRODUCING 
***    THE WAVELENTGH FILE.

      WRITE(6,6004) DX

***   CHECK AND SET THE DETRENDING OPTION 
      IF (CAS.EQ.2) THEN

              WRITE(6,6005) DET

              IF (DET.EQ.0) THEN

                 DO_DETREND=.FALSE.
                 WRITE(6,*)'DETRENDING IS NOT APPLIED BEFORE THE DFT'

              ELSE IF (DET.EQ.1) THEN
                    
                 DO_DETREND=.TRUE.
                 WRITE(6,*)'DETRENDING IS APPLIED BEFORE THE DFT'
                       
              ELSE
                     
              WRITE(6,*)'PROBLEM WITH THE CHOICE OF THE DETRENDING OPT'
              CALL                                 XIT(' VARSPEC',-4 )
                      
              ENDIF
                    
         ENDIF

      WRITE(6,6006) AVG1D

      IF (AVG1D.EQ.0) THEN

         DO_AVG1D=.FALSE.
         WRITE(6,*) 'FULL 2D SPECTRA ARE COMPUTED'
         WRITE(6,*) '==> NO AVERAGE OF 1D SPECTRA WILL BE COMPUTED'

         WRITE(6,*) 'TOTO'

         ELSEIF (AVG1D.EQ.1) THEN
            
            DO_AVG1D=.TRUE.
            WEST_EAST=.TRUE.
        WRITE(6,*)'    AVERAGE OF 1D WEST-EAST SPECTRA WILL BE COMPUTED'

         ELSEIF (AVG1D.EQ.2) THEN
            
            DO_AVG1D=.TRUE.
            WEST_EAST=.FALSE.
      WRITE(6,*)'    AVERAGE OF 1D SOUTH-NORTH SPECTRA WILL BE COMPUTED'


         ELSE
               WRITE(6,*)'AVG1D OPTION VALUE IS NOT CORRECT.'
               CALL                                XIT(' VARSPEC',-5 )
         ENDIF

***    ALLOCATE WORK ARRAYS. help 

      WRITE(6,*) 'TOTO 1'

      MAXW = GETSIZ( 1, KBUF,LEN,NWDS,KPAK )

      WRITE(6,*) 'TOTO 2'

      IF (MAXW.LE.0) CALL                          XIT(' VARSPEC',-6 )

      WRITE(6,*) 'TOTO 3'

      NI = KBUF(5)
      NJ = KBUF(6)

      MAXI = MAX( NI+2 , NJ+2 )
      MAXL =     (NI+2)*(NJ+2)

!       print*,'JFC: HPALLOC 1 start'
!       CALL HPALLOC( PA, MAXL*7+LEN*2+MAXI, IER, 8 )
!       print*,'JFC: HPALLOC 1 end'
!       PB  = LOC(        A(MAXL+1) )
!       PR1 = LOC(        B(MAXL+1) )
!       PR2 = LOC(       R1(MAXL+1) )
!       PR3 = LOC(       R2(MAXL+1) )
!       PVA = LOC(       R3(MAXL+1) )
!       PSP = LOC( VARIANCE(MAXL+1) )
!       PIB = LOC(  SPECTRE(MAXI+1) )
!       PJB = LOC(     IBUF(LEN +1) )
      allocate(        A(MAXL+1),       B(MAXL+1))
      allocate( VARIANCE(MAXL+1), SPECTRE(MAXI+1))
      allocate(     IBUF(LEN +1),    JBUF(LEN +1))
*BD<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*BD<< A MODIFIER MEMOIRE   ---> FIN <---
*BD<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

***    INITIALIZE JBUF

      DO  J=1,HEAD
          JBUF(J) = KBUF(J)
      END DO

*BD<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*BD<< A MODIFIER POUR I/O   ---> DEBUT <---
*BD<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

***    READ GRID FROM FILE 1 
*      =======================

  100 CALL GETFLD2( 1, B, 'GRID'//'SUBA',-1,' ',-1, IBUF,MAXW,OK )

      IF (.NOT.OK)                                             THEN 
          IF (NRECS.EQ.0)                                      THEN
              CALL                                 XIT(' VARSPEC',-6 )
          ELSE
              WRITE(6,6010) NRECS 
              CALL                                 XIT(' VARSPEC', 0 ) 
          END IF 
      END IF 

      IF (NRECS.EQ.0) WRITE(6,6000) IBUF

      CALL CMPLBL( 0,IBUF, 0,KBUF, OK )

      IF (.NOT.OK)                                             THEN
          WRITE(6,6000) IBUF
          CALL                                     XIT(' VARSPEC',-7 )
      END IF


      WRITE (6,*) ' '
      WRITE (6,*) '    INPUT FIELD READ:'
      WRITE (6,*) '    ================='


***    PREPARATION OF THE INPUT ARRAY FOR THE SELECTED TRANSFORM
*      ===========================================================

*           ---------------------------------------------
***         * STANDARD PERIODIC FOURIER TRANSFORM CASE. *
*           ---------------------------------------------

      IF (LDFT)                                                THEN
         
         IF (.NOT.DO_AVG1D) THEN
           IF (MOD( NI,2 ).EQ.0 .OR.
     +         MOD( NJ,2 ).EQ.0 )                               THEN
               
          WRITE(6,*)'**************************************************'
          WRITE(6,*)'PROBLEM WITH NI= ', NI, ' AND/OR NJ=' ,NJ
          WRITE(6,*)'-------'
          WRITE(6,*)'STANDARD PERIODIC FOURIER TRANSFORMS (DFT) REQUIRE'
          WRITE(6,*)'THAT INPUT FIELD DIMENSIONS BE ODD. THE DFT WORKS'
          WRITE(6,*)'ON EVEN DIMENSIONS OF ROWS (COLUMNS)  MADE OF '
          WRITE(6,*)'DICTINCT  ELEMENTS. THE LAST ELEMENT OF THE INPUT'
          WRITE(6,*)'ROWS (COLUMNS) ARE ASSUMED TO BE A REPETITION OF'
          WRITE(6,*)'THE FIRST ELEMENT BY PERIODICITY.'
          WRITE(6,*)'**************************************************'
            CALL                                 XIT(' VARSPEC',-8 )
           ENDIF
  
        ELSE

           IF (WEST_EAST.AND.MOD( NI,2 ).EQ.0) THEN
          WRITE(6,*)'**************************************************'
          WRITE(6,*)'PROBLEM WITH NI= ', NI
          WRITE(6,*)'-------'
          WRITE(6,*)'STANDARD PERIODIC FOURIER TRANSFORMS (DFT) REQUIRE'
          WRITE(6,*)'THAT INPUT FIELD DIMENSIONS BE ODD IN X. THE DFT '
          WRITE(6,*)'WORKS ON EVEN DIMENSIONS OF ROWS MADE OF DICTINCT '
          WRITE(6,*)'ELEMENTS. THE LAST ELEMENT OF THE INPUT ROWS'
          WRITE(6,*)'ARE ASSUMED TO BE A REPETITION OF THE FIRST '
          WRITE(6,*)'ELEMENT BY PERIODICITY.'
          WRITE(6,*)'**************************************************'
            CALL                                 XIT(' VARSPEC',-8 )
           ENDIF

           IF ((.NOT.WEST_EAST).AND.MOD( NJ,2 ).EQ.0) THEN
          WRITE(6,*)'**************************************************'
          WRITE(6,*)'PROBLEM WITH NJ= ', NJ
          WRITE(6,*)'-------'
          WRITE(6,*)'STANDARD PERIODIC FOURIER TRANSFORMS (DFT) REQUIRE'
          WRITE(6,*)'THAT INPUT FIELD DIMENSIONS BE ODD IN Y. THE DFT'
          WRITE(6,*)'WORKS ON EVEN DIMENSIONS OF COLUMNS MADE OF '
          WRITE(6,*)'DISTINCT ELEMENTS. THE LAST ELEMENT OF THE INPUT'
          WRITE(6,*)'COLUMNS ARE ASSUMED TO BE A REPETITION OF THE '
          WRITE(6,*)'FIRST ELEMENT BY PERIODICITY.'
          WRITE(6,*)'**************************************************'
            CALL                                 XIT(' VARSPEC',-8 )
           ENDIF

            
        ENDIF
         
***      APPLY DETRENDING IF ASKED FOR
 
         IF (DO_DETREND) THEN

            IF (.NOT.DO_AVG1D) THEN             ! 2D case

               CALL DETREND(B,NI,NJ,1,1)

            ELSE                               ! 1D case
               IF (WEST_EAST) THEN           

                  CALL DETREND(B,NI,NJ,1,0)

                ELSE                           ! NORTH_SOUTH case

                  CALL DETREND(B,NI,NJ,0,1)

               ENDIF

            ENDIF

         ELSE
          WRITE(6,*)' ***WARNING****: NO DETRENDING IS APPLIED !!'
          WRITE(6,*)'=>DETRENDING IS RECOMMENDED WHEN DFT IS CHOSEN<='
         ENDIF

***    AT THIS POINT THE INPUT GRID IS ASSUMED TO BE PERIODIC SO 
***    THAT THE LAST COLUMN AND LAST ROWS HAVE TO BE REMOVED 
***    BEFORE PASSING IT TO THE  TRANSFORMS WHICH WORKS ON VECTOR 
***    CONTANING DISTINCT DATA. SINCE THE INPUT FIELDS WHERE REQUIRED 
***    TO BE ODD, THEY WILL BECOME EVEN NOW.
         
          NID = NI - 1            ! NUMBER OF DISTINCT POSITION IN X.
          NJD = NJ - 1            ! NUMBER OF DISTINCT POSITION IN X.

***    ISSUE A WARNING IF A FAST TRANSFORM CANNOT BE USED.

          NFACT = NID
          CALL NGFFT( NFACT )
         
          IF ((.NOT.DO_AVG1D).OR.(DO_AVG1D.AND.WEST_EAST)) THEN
           IF (NFACT.NE.NID)                                    THEN 
             WRITE(6,*) 'WARNING: A FAST TRANSFORM CANNOT BE USED IN X:'
             WRITE(6,6110) NID,NFACT
           END IF
          ENDIF

          NFACT = NJD
          CALL NGFFT( NFACT )
         
          IF ((.NOT.DO_AVG1D).OR.(DO_AVG1D.AND.(.NOT.WEST_EAST))) THEN
           IF (NFACT.NE.NJD)                                    THEN 
             WRITE(6,*) 'WARNING: A FAST TRANSFORM CANNOT BE USED IN Y:'
             WRITE(6,6120) NJD,NFACT
           END IF
          ENDIF

***    SET THE MAXIMUM OF INTEGER WAVENUMBERS SUPPORTED BY THE GRID.
      
          MMAX = NID/2
          NMAX = NJD/2

***    FOR THIS STANDARD FFT WE NEED TO PADD THE ARRAYS SUCH THIS:

*                       --------------------
*                       |O O O O O O O O O |
*                       |O O O O O O O O O |
*                       |X X X X X X X O O |
*                       |X X X X X X X O O |
*                       |X X X X X X X O O |
*                       |X X X X X X X O O |
*                       --------------------

          NIP = 2 ! EXTRA SPACES NEEDED BY THE FFT
          NJP = 2

      END IF

*               -------------------------------------
***             * DISCRETE COSINE TRANSFORM CASE.*
*               -------------------------------------

      IF (LDCT)                                               THEN
         
          NIP  = 0
          NJP  = 0
          NID  = NI
          NJD  = NJ
          MMAX = NID/2
          NMAX = NJD/2

***    STOP IF THE DIMENSIONS ARE NOT EVEN 
         IF (.NOT.DO_AVG1D) THEN

          IF (MOD( NID,2 ).NE.0 .OR.
     +        MOD( NJD,2 ).NE.0 )                              THEN

          WRITE(6,*)'**************************************************'
          WRITE(6,*)'PROBLEM WITH NI= ', NI, ' AND/OR NJ=' ,NJ
          WRITE(6,*)'-------'
          WRITE(6,*) ' DISCRETE COSINE TRANSFORMS REQUIRE THAT '
          WRITE(6,*) ' DIMENSIONS BE EVEN '
              CALL                                 XIT(' VARSPEC',-9 )
           ENDIF

          ELSE

           IF (WEST_EAST.AND.MOD( NI,2 ).NE.0) THEN
          WRITE(6,*)'**************************************************'
          WRITE(6,*)'PROBLEM WITH NI= ', NI
          WRITE(6,*)'-------'
          WRITE(6,*) ' DISCRETE COSINE TRANSFORMS REQUIRE THAT '
          WRITE(6,*) ' DIMENSION NI BE EVEN '
              CALL                                 XIT(' VARSPEC',-9 )
           ENDIF

           IF ((.NOT.WEST_EAST).AND.MOD( NJ,2 ).NE.0) THEN
          WRITE(6,*)'**************************************************'
          WRITE(6,*)'PROBLEM WITH NJ= ', NJ
          WRITE(6,*)'-------'
          WRITE(6,*) ' DISCRETE COSINE TRANSFORMS REQUIRE THAT '
          WRITE(6,*) ' DIMENSION NJ BE EVEN '
              CALL                                 XIT(' VARSPEC',-9 )
           ENDIF
         ENDIF

         
***    ISSUE A WARNING IF A FAST TRANSFORM CANNOT BE USED.

          NFACT = NID
          CALL NGFFT( NFACT )
         
        IF ((.NOT.DO_AVG1D).OR.(DO_AVG1D.AND.WEST_EAST)) THEN
          IF (NFACT.NE.NID)                                    THEN 
             WRITE(6,*) 'WARNING: A FAST TRANSFORM CANNOT BE USED IN X:'
             WRITE(6,6130) NID,NFACT
          END IF
        ENDIF

          NFACT = NJD
          CALL NGFFT( NFACT )
         
         IF ((.NOT.DO_AVG1D).OR.(DO_AVG1D.AND.(.NOT.WEST_EAST))) THEN 
          IF (NFACT.NE.NJD) THEN 
             WRITE(6,*) 'WARNING: A FAST TRANSFORM CANNOT BE USED IN Y:'
             WRITE(6,6140) NJD,NFACT
          END IF
         ENDIF

      END IF
      
***    BUILT THE ARRAY "A" FOR THE TRANSFORMS.
***    THIS DOUBLE LOOP IS COMMON TO ALL CASES.

      NIT = NID+NIP
      NJT = NJD+NJP

      IJ  = 0
      IJT = 0

      DO  J=1,NJT
          DO  I=1,NIT
            
              IJT = IJT+1

              IF (I.GT.NID .OR. 
     +            J.GT.NJD )                                   THEN
                  A(IJT) = 0.0
              ELSE

!BD<< * BUGFIX FOR THE MAPPING OF THE INPUT ARRAY IN THE CASE OF DFT
!BD     WHICH HAS ONE ADDITIONAL COLUMN THAN FOR THE DCT CASE.            

                  IF ( LDFT.and. I.EQ.1) IJ=IJ+1
!BD>>
                  IJ     = IJ+1

                  A(IJT) = B(IJ)
              END IF

          END DO
      END DO


      IF(DO_AVG1D) THEN  ! WE NEED TO MAKE A SEPARATE COPIES OF THE INPUT FIELDS


        DO I=1,NIT*NJT
            AX(I)=A(I)
            AY(I)=A(I)
         ENDDO

      ENDIF


***    DO THE TRANSFORMS IN TWO 1-D PASSES.
*      ======================================

***    FIRST PASS ==> ALONG INDEX "I"
*      ----------

      IWAY = -1                   ! REVERSE  TRANSFORM (GRIDPOINT TO FOURIER)
      

      IAXE = 0                    ! TRANSFORM FOLLOWING X
      
      IF (IAXE.EQ.0)                                           THEN
          NX = NID                ! TRANSFORM FOLLOWING X
          NY = NJD
      ELSE 
          NX = NJD                ! TRANSFORM FOLLOWING Y
          NY = NID
      END IF

      IF(DO_AVG1D) THEN
         IF (WEST_EAST) THEN
            CALL TRANS1D( AX,U1,U2,U3,MAXL,NX,NY,CAS,IAXE,IWAY)
         ENDIF
      ELSE                      ! THIS ONE IS FOR THE FULL 2D CASE 
         CALL TRANS1D( A,U1,U2,U3,MAXL,NX,NY,CAS,IAXE,IWAY )
      ENDIF
      
***    SECOND PASS => ALONG INDEX "J"
*      -----------

      IAXE = 1                    ! TRANSFORM FOLLOWING Y

      IF (IAXE.EQ.0)                                           THEN
          NX = NID                ! TRANSFORM FOLLOWING X
          NY = NJD                ! TRANSFORM FOLLOWING X
      ELSE 
          NX = NJD                ! TRANSFORM FOLLOWING Y
          NY = NID                ! TRANSFORM FOLLOWING X
      END IF


      IF(DO_AVG1D) THEN
         IF (.NOT.WEST_EAST) THEN ! DO THE SOUTH-NORTH CASE
            CALL TRANS1D( AY,U1,U2,U3,MAXL,NX,NY,CAS,IAXE,IWAY )
         ENDIF
         
      ELSE ! THIS ONE IS FOR THE FULL 2D CASE
         CALL TRANS1D( A,U1,U2,U3,MAXL,NX,NY,CAS,IAXE,IWAY )
      ENDIF
      

***    COMPUTE THE SPECTRA (2D AND 1D, RESPECTIVELY)
*      ====================

***    SPECT1/SPECT2 COMPUTE THE VARIANCE FROM ARRAY OF 
***    COEFFICIENT "A". THE OUTPUT IS THE 2D VARIANCE ARRAY "VARIANCE"
***    AND THE 1-D POWER SPECTRUM "SPECTRE. 

      IF (LDFT)                                                THEN 
         
         IF (DO_AVG1D) THEN
            CALL SPECT_DFT1D( SPECTRE,VARIANCE, AX,AY,
     +                        MMAX,NMAX,NIT,NJT,MAXI,WEST_EAST )
        ELSE
            CALL SPECT_DFT2D( SPECTRE,VARIANCE,A,
     +                        MMAX,NMAX,NIT,NJT,MAXI )
         ENDIF
         
         
      ELSE IF (LDCT)                                           THEN
         
         IF (DO_AVG1D) THEN
            CALL SPECT_DCT1D( SPECTRE,VARIANCE,AX,AY,
     +                        NIT,NJT,MAXI,REGROUPE,SKIP,MOY,WEST_EAST )
         ELSE
            CALL SPECT_DCT2D( SPECTRE,VARIANCE,A,
     +                        NIT,NJT,MAXI,REGROUPE,SKIP,MOY)
         ENDIF
      END IF

      IF (DO_AVG1D) THEN
         IF (WEST_EAST) THEN 
         WRITE(6,*)'AVERAGE OF ', NJD,' 1D WEST-EST SPECTRA COMPUTED'
         ELSE
         WRITE(6,*)'AVERAGE OF ', NID,' 1D SOUTH-NORTH SPECTRA COMPUTED'
         ENDIF
      ENDIF


***    OUTPUT SECTION
*      ================

      DO  I=1,HEAD
          JBUF(I) = IBUF(I)
      END DO

      JBUF(8) = -32        ! NO PACKING !!!

***    FIRST: * 2D VARIANCE OUPUT FIELD.
*      ---------------------------------

      IF (LDFT)                                                THEN
         
          JBUF(5) = MMAX+1        ! THE "+1" IS FOR THE WAVE COMPONENT "0"
          JBUF(6) = NMAX+1

      ELSE IF (LDCT)                                           THEN
         
          JBUF(5) = 2*MMAX
          JBUF(6) = 2*NMAX

      END IF

***    COMPUTE THE LOG OF THE 2D VARIANCE OUPUT FIELD
***     IF REQUIRED AND WRITE IT.

      IF (LOGV.EQ.0)                                           THEN
          DO  L=1,JBUF(5)*JBUF(6)
              B(L) = VARIANCE(L)
          END DO
      ELSE
          DO  L=1,JBUF(5)*JBUF(6)
              IF (VARIANCE(L).EQ.0.)
     +            VARIANCE(L) =  1.0  ! TO AVOID NAN !
              B(L) = LOG10( VARIANCE(L) )
          END DO
      END IF

      CALL PUTFLD2( 2, B, JBUF,MAXW ) 

      WRITE (6,*) ' '
      WRITE (6,*) '    2D ARRAY OF VARIANCE WRITTEN:'
      WRITE (6,*) '    -----------------------------' 
      WRITE (6,6000) JBUF



***    SECOND: * 1D POWER SPECTRUM AND THE WAVELENGTH FILE
***    ---------------------------------------------------

      READ(ZONL,0004) JBUF(1)
      JBUF(6) = 1

      LAST_WAVENUMBER=.FALSE.   ! SWITCH TO OUTPUT OR NOT THE VALUE
                                ! ASSOCIATED WITH THE LAST WAVENUMBER

      IF (LDFT)                                                THEN

         IF (DO_AVG1D) THEN  
            IF (WEST_EAST) THEN
               NCOEF = MMAX + 1
            ELSE
               NCOEF = NMAX + 1
            ENDIF
         ELSE                       ! FULL 2D CASE

            NCOEF = MIN( MMAX,NMAX ) + 1

         ENDIF

          DO  L=1,NCOEF 
              B(L) = SPECTRE(L)
          END DO
         
          IF (LAST_WAVENUMBER)                                      THEN
              JBUF(5) = NCOEF
           ELSE                     
              JBUF(5) = NCOEF - 1   ! WE OMIT THE LAST ELEMENT
          END IF

          CALL PUTFLD2( 3, B, JBUF,MAXW ) 

          WRITE (6,*) ' '
          WRITE (6,*) '    POWER SPECTRUM WRITTEN:'
          WRITE (6,*) '    -----------------------' 
          WRITE( 6,6000) JBUF


C         ****************************************
C         ****************************************


          B(1) = 1.E38           ! WAVELENTGH OF WAVENUMBER ZERO SHOULD
                                 ! BE INFINITY BUT WE PUT IT TO 1.E38.

          IF (DO_AVG1D) THEN  
             IF (WEST_EAST) THEN
                NWAVE = MMAX
             ELSE
                NWAVE = NMAX
             ENDIF
          ELSE  ! FULL 2D CASE
             NWAVE=MIN( MMAX,NMAX )
          ENDIF
          
          DO  L=1,NWAVE
              B(L+1) = 2*NWAVE*DX/FLOAT( L )
          END DO
         

          CALL PUTFLD2( 4, B, JBUF,MAXW ) 

          WRITE (6,*) ' '
          WRITE (6,*) '    ASSOCIATED WAVELENGTHS WRITTEN:'
          WRITE (6,*) '    -----------------------' 
          WRITE (6,6000) JBUF


C         ****************************************
C         ****************************************


      ELSE IF (LDCT)                                           THEN
         
          IF (REGROUPE)                                        THEN

              IF (SKIP)                                        THEN

                 IF (DO_AVG1D) THEN  
                    IF (WEST_EAST) THEN
                       NCOEF = MMAX
                    ELSE
                       NCOEF = NMAX 
                    ENDIF
                    
                 ELSE           ! FULL 2D CASE
                    NCOEF = MIN( MMAX,NMAX )
                 ENDIF
                 
                 DO  L=1,NCOEF 
                    B(L) = SPECTRE(L)
                 END DO
               
                  JBUF(5) = NCOEF


                  CALL PUTFLD2( 3, B, JBUF,MAXW ) 

                  WRITE (6,*) ' '
                  WRITE (6,*) '    POWER SPECTRUM WRITTEN:'
                  WRITE (6,*) '    ----------------------' 
                  WRITE( 6,6000) JBUF
              
                  B(1) = 1.E38       ! WAVELENTGH OF WAVENUMBER ZERO SHOULD
                                     ! BE INFINITY BUT WE PUT IT TO 1.E38.
               
                  IF (DO_AVG1D) THEN  
                     IF (WEST_EAST) THEN
                        NWAVE = MMAX - 1
                     ELSE
                        NWAVE = NMAX - 1
                     ENDIF
                  ELSE               ! FULL 2D CASE
                     NWAVE = MIN( MMAX,NMAX ) - 1
                     
                  ENDIF

                  DO  L=1,NWAVE
                      B(L+1) = 2*(NWAVE+1)*DX/FLOAT( L )
                  END DO
 
              
                  CALL PUTFLD2( 4, B, JBUF,MAXW ) 

          WRITE (6,*) ' '
          WRITE (6,*) '    ASSOCIATED WAVELENGTHS WRITTEN:'
          WRITE (6,*) '    -----------------------' 
          WRITE (6,6000) JBUF


C         ****************************************
C         ****************************************

              
              ELSE
               
                  IF (DO_AVG1D ) THEN
                    IF (WEST_EAST) THEN
                       NCOEF  = MMAX + 1 ! THE "+ 1" IS FOR COEF. "0"
                    ELSE
                       NCOEF  = NMAX  + 1 ! THE "+ 1" IS FOR COEF. "0"
                    ENDIF
                    
                  ELSE                          ! FULL 2D CASE
                     NCOEF  = MIN( MMAX,NMAX ) + 1 ! THE "+ 1" IS FOR COEF. "0"
                  ENDIF

                  DO  L=1,NCOEF 
                      B(L) = SPECTRE(L)
                  END DO
                              
                  IF (LAST_WAVENUMBER)                         THEN
                     JBUF(5) = NCOEF
                  ELSE 
                      JBUF(5) = NCOEF - 1 ! WE OMIT THE LAST ELEMENT
                  END IF
 
                  CALL PUTFLD2( 3, B, JBUF,MAXW ) 

                  WRITE (6,*) ' '
                  WRITE (6,*) '    POWER SPECTRUM WRITTEN:'
                  WRITE (6,*) '    -----------------------' 
                  WRITE (6,6000) JBUF



                  B(1) = 1.E38          ! WAVELENTGH OF WAVENUMBER ZERO SHOULD
                                        ! BE INFINITY BUT WE PUT IT TO 1.E38.
            
                  IF (DO_AVG1D ) THEN
                     IF (WEST_EAST) THEN
                        NWAVE = MMAX
                     ELSE
                        NWAVE = NMAX
                     ENDIF
                  ELSE                  ! FULL 2D CASE      
                     NWAVE = MIN( MMAX,NMAX )
                  ENDIF

                  DO  L=1,NWAVE
                      B(L+1) = 2*(NWAVE)*DX/FLOAT( L )
                  END DO

                  CALL PUTFLD2( 4, B, JBUF,MAXW ) 

                  WRITE (6,*) ' '
                  WRITE (6,*) '    ASSOCIATED WAVELENGTHS WRITTEN:'
                  WRITE (6,*) '    -----------------------' 
                  WRITE (6,6000) JBUF

C                 ****************************************
C                 ****************************************

              END IF

          ELSE
             IF (DO_AVG1D)                               THEN  

                IF (WEST_EAST)                           THEN
                   NCOEF = 2*MMAX
                ELSE
                   NCOEF = 2*NMAX
                ENDIF

             ELSE               ! FULL 2D CASE
                NCOEF = MIN( 2*MMAX,2*NMAX )

             ENDIF
             
            
              DO  L=1,NCOEF 
                  B(L) = SPECTRE(L)
              END DO
        
              JBUF(5) = NCOEF


              CALL PUTFLD2( 3, B, JBUF,MAXW ) 

              WRITE (6,*) ' '
              WRITE (6,*) '    POWER SPECTRUM WRITTEN:'
              WRITE (6,*) '    -----------------------' 
              WRITE (6,6000) JBUF



C             ****************************************
C             ****************************************

           
              B(1) = 1.E38          ! WAVELENTGH OF WAVENUMBER ZERO SHOULD
                                    ! BE INFINITY BUT WE PUT IT TO 1.E38.
            
              IF (DO_AVG1D)                                    THEN  

                 IF (WEST_EAST)                                THEN
                    NWAVE = 2*MMAX
                 ELSE
                    NWAVE = 2*NMAX
                 ENDIF

              ELSE              ! FULL 2D CASE

                 NWAVE = MIN( 2*MMAX,2*NMAX )

              ENDIF
              
              DO  L=1,NWAVE
                 B(L+1) = 2*NWAVE*DX/FLOAT( L )
              END DO

              CALL PUTFLD2( 4, B, JBUF,MAXW ) 

              WRITE (6,*) ' '
              WRITE (6,*) '    ASSOCIATED WAVELENGTHS WRITTEN:'
              WRITE (6,*) '    -----------------------' 
              WRITE (6,6000) JBUF

C             ****************************************
C             ****************************************
       
          END IF
         
      END IF

      WRITE (6,*) ' '
      WRITE (6,*) ' '     
 
      IF (NRECS.EQ.0) WRITE(6,6000) JBUF

      NRECS = NRECS+1 
      GOTO 100 

*---------------------------------------------------------------------
***    E.O.F. OR ERROR ON UNIT 5.

  900 IF (IO.NE.0)
     +    WRITE(6,6160) IO

  901 CALL                                         XIT(' VARSPEC',-10 )
  
*-----------------------------------------------------------------------
 0004 FORMAT(A4)
 0005 FORMAT(BN,I5)
 0010 FORMAT(BN,E10.0)

 5000 FORMAT(10X,3I5,E5.0,2I5) 

 6000 FORMAT(5X,A4,I12,1X,A4,I10,2I4,I8,I6)
 6001 FORMAT(' CAS  =',I1)
 6002 FORMAT(' LOGV =',I1)
 6003 FORMAT(' CLUST=',I1)
 6004 FORMAT(' DX   =',F5.1,' KILOMETERS')
 6005 FORMAT(' DET  =',I1)
 6006 FORMAT(' AVG1D  =',I1)
 6010 FORMAT('0',I6,' RECORDS PROCESSED.')
 6110 FORMAT(' N = NI-1  =  ', I4,' THE NEAREST FACTORIZABLE N = ',I4)
 6120 FORMAT(' N = NJ-1  =  ', I4,' THE NEAREST FACTORIZABLE N = ',I4)
 6130 FORMAT(' N = NI = ', I4,' THE NEAREST FACTORIZABLE N = ',I4)
 6140 FORMAT(' N = NJ = ', I4,' THE NEAREST FACTORIZABLE N = ',I4)
 6160 FORMAT(' PROBLEM READING INPUT PARAMETRES: IOSTAT=',I5)

#endif
      END
      SUBROUTINE DETREND(B, NI,NJ,IX,IY)

*PURPOSE  - REMOVE LINEAR TENDENCY OF A 2-D FIELD USING THE METHODOLOGY
*           PROPOSED BY ERRICO, R. M., 1985: SPECTRA COMPUTED FROM
*           A LIMITED AREA GRID. MWR, VOL. 113, P.1554-1562.
* -----------------------------------------------------------------------
      IMPLICIT  NONE

      INTEGER I,J,IJ,NI,NJ,INC,IX,IY
      REAL SLOPE
      REAL*8 B(NI*NJ)

      IF (IX.EQ.1) THEN
*           DETRENDING ALONG THE X AXIS
*
         DO 200 J=1,NJ
            INC=NI*(J-1)
            SLOPE=( B(NI+INC) - B(1+INC) ) / (NI-1)
            DO 200 I=1,NI
               B(I+INC)=B(I+INC)-0.5*(2*I-NI-1)*SLOPE
 200        CONTINUE

      ENDIF   

      IF (IY.EQ.1) THEN
*
*           DETRENDING ALONG THE Y AXIS
*

         DO 300 I=1,NI
            SLOPE=( B(I+NI*(NJ-1)) - B(I) ) / (NJ-1)
            DO 300 J=1,NJ
               B(I+NI*(J-1))=B(I+NI*(J-1))-0.5*(2*J-NJ-1)*SLOPE
 300        CONTINUE
            
      ENDIF

      RETURN
      END

      SUBROUTINE SPECT_DFT2D( SPECTRE,VARIANCE,AMP,
     +                        MMAX,NMAX,NIT,NJT,MAXI )

      IMPLICIT  NONE

*     AUTHOR  - B.DENIS  AVRIL 2000
*
*     PURPOSE - COMPUTES SPECTRAL VARIANCE ARRAY AND POWER SPECTRUM FROM             
*               2D SPECTRAL COEFFICIENT ARRAY PRODUCED BY THE DFT
*               STANDARD PERIODIC FOURIER TRANSFORM.
*               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*
*     INPUT VARIABLES...                                                   
*                                                                               
*       AMP  = 2D SPECTRAL COEFFICIENT ARRAY FROM TRANSFORM.
*       NIT  = X DIMENSION OF THE SPECTRAL COEFFICIENT ARRAY.
*       NJT  = Y DIMENSION OF THE SPECTRAL COEFFICIENT ARRAY.                                                                     
*       MMAX = MAX X WAVENUMBER  
*       NMAX = MAX Y WAVENUMBER 
*       MAXI = 1D MEMORY DIMENSION

*      OUTPUT VARIABLES...

*       SPECTRE =  1-D POWER SPECTRUM FIELD.
*       VARIANCE=  2-D POWER SPECTRUM FIELD.

*       NB: AMP, SPECTRE AND VARIANCE ARRAYS START AT INDEX "0" FOR
*           WAVENUMBER "0".IN THE MAIN, THEY START AT INDEX "1" 
*           FOR WAVENUMBER "0".
*-----------------------------------------------------------------------

      INTEGER   NIT,NJT,MMAX,NMAX,MAXI
      REAL*8    AMP(0:NIT-1,0:NJT-1)
      REAL*8    VARIANCE(0:MMAX,0:NMAX), SPECTRE(0:MAXI-1)
      REAL*8    R, A, B

*BD<< bugfix that correct improper results 
*BD== REAL      DKM,DKN
      REAL*8    DKM,DKN

*BD== REAL      RCUT_LOW,RCUT_HIGH
      REAL*8    RCUT_LOW,RCUT_HIGH
*BD>>

      INTEGER   RMTRUNC,RNTRUNC
      INTEGER   I,J,K,N,M,NBWV

***    FIRST: COMPUTE THE VARIANCE OF EACH COEFFICIENT.
***    -----

***    BECAUSE THE TRANSFORM ROUTINE WORKS FROM REAL TO REAL
***    (NO COMPLEX), FOR 2D REAL FIELD EACH VARIANCE (M,N) IS 
***    MADE OF FOR REAL PARTS SINCE:

***    FOR A GIVEN (M,N), 
***    VARIANCE(M,N) = ACC**2 + ACS**2 + ASC**2 + ASS**2

***    ACC  COMES FROM COS(M)*COS(N) TERM
***    ACS     "    "  COS(M)*SIN(M)  "
***    ASC     "    "  SIN(M)*COS(N)  "
***    ASS     "    "  SIN(M)*SIN(N)  "

***    FURTHERMORE, FOR NORMALIZATION REASONS THESE FACTORS ARE USED.

*      ----------------------------------------------------
*     |     "            "            "           "   
*     | 2*VAR{0,3} || 4*VAR{1,3} | 4*VAR{2,3} |   "
*      -----------------------------------------------------
*     | 2*VAR{0,2} || 4*VAR{1,2} | 4*VAR{2,2} |   "
*      -----------------------------------------------------
*     | 2*VAR{0,1} || 4*VAR{1,1} | 4*VAR{2,1} |   "
*      ====================================================
*     | 1*VAR{0,0} || 2*VAR{1,0} | 2*VAR{2,0} |   "
*      -----------------------------------------------------
*    
*-----------------------------------------------------------------------
***    INITIALISATION

c     Debut - addition JFC
      open (unit=1,file="spec_comp_VSP.txt",action="write",status="new")
      DO  J=0,NMAX      
          DO I=0,MMAX
              write (1,'(E15.6)') amp(i,j)
          END DO
      END DO
      close(unit=1)
c     Fin   - addition JFC

      DO  J=0,NMAX      
          DO I=0,MMAX
              VARIANCE(I,J) = 0
          END DO
      END DO


      DO  J=0,NMAX      
          DO I=0,MMAX
              VARIANCE(I,J) = 0
          END DO
      END DO

***    1- TERMS INVOLVING THE FACTOR "4" (UPPER RIGHT CORNER).

      DO  N=1,NMAX
          DO M=1,MMAX
             I = 2*M
             J = 2*N

             VARIANCE(M,N) = 4.0*(  AMP(I,  J  )**2 
     +                            + AMP(I,  J+1)**2
     +                            + AMP(I+1,J  )**2
     +                            + AMP(I+1,J+1)**2)
          END DO
      END DO

***    2- TERMS INVOLVING THE FACTOR "2" IN THE LOWER ROW. 

      DO  M=1,MMAX
          I = 2*M
          VARIANCE(M,0) = 2.0*(  AMP(I,  0)**2
     +                         + AMP(I,  1)**2
     +                         + AMP(I+1,0)**2
     +                         + AMP(I+1,1)**2)
      END DO

***    3- TERMS INVOLVING THE FACTOR "2" IN THE LEFT COLUMN. 

      DO  N=1,NMAX
          J = 2*N
          VARIANCE(0,N) = 2.0*(  AMP(0,  J)**2
     +                         + AMP(0,J+1)**2
     +                         + AMP(1,  J)**2
     +                         + AMP(1,J+1)**2)
      END DO

***    4- ELEMENT 0,0 (THE MEAN).

         VARIANCE(0,0) = (  AMP(0,0)**2
     +                    + AMP(0,1)**2
     +                    + AMP(1,0)**2
     +                    + AMP(1,1)**2)


***    SECOND: COMPUTE THE POWER SPECTRUM
***    ------

***    WE DEFINE THE NUMBER OF WAVEBAND. SHOULD CORRECTLY TAKE
***    CARE OF NON-SQUARE ARRAYS.
 
      NBWV = MIN( MMAX,NMAX )
      DKM  = MMAX/FLOAT( NBWV )   ! DEFINE INCREMENT IN WAVENUMBER K.
      DKN  = NMAX/FLOAT( NBWV )
      
***    WE INITIALIZE SPECTRE.
      
      DO  K=0,MAXI-1
          SPECTRE(K) = 0.0
      END DO
      

***    FOR EACH WAVE COMPONANT WE GO THROUGH THE ARRAY OF VARIANCE.

***    EACH SPECTRAL COMPONANT "KS" OF THE POWER SPECTRUM "SPECTRE" IS 
***    MADE FROM THE CONTRIBUTION OF THE VARIANCES THAT HAS A 
***    WAVENUMBER K SUCH THAT  KS-1 < K <= KS.

***         K       |   KS
***     ------------------      
***     0          -- > 0
***     0 < K <= 1 -- > 1
***     1 < K <= 2 -- > 2
***     2 < K <= 3 -- > 3
***       " " ""   """"  

***    ELEMENT K=0 IS THE MODE 0, I.E. THE SQUARE OF THE MEAN.
***    WE GOT IT FROM VARIANCE(0,0):
      
      SPECTRE(0) = VARIANCE(0,0)

      DO  K=1,NBWV

          RCUT_LOW  = (K - 1.0)/FLOAT( K )
          RCUT_HIGH = (K + 0.0)/FLOAT( K )

          RMTRUNC   = DKM*K
          RNTRUNC   = DKN*K

          DO  N=0,NMAX
              DO  M=0,MMAX
               
                  A = M/FLOAT( RMTRUNC )
                  B = N/FLOAT( RNTRUNC )
               
                  R = SQRT((A**2) + (B**2))

                  IF ((R.GT.RCUT_LOW).AND.(R.LE.RCUT_HIGH))    THEN
                      SPECTRE(K)  = SPECTRE(K) + VARIANCE(M,N)
                  END IF
               
              END DO
          END DO
         
      END DO

      RETURN
*-----------------------------------------------------------------------

      END
      SUBROUTINE SPECT_DFT1D( SPECTRE,VARIANCE,AMPX,AMPY,
     +     MMAX,NMAX,NIT,NJT,MAXI,WEST_EAST )

      IMPLICIT  NONE

*     AUTHOR  - B.DENIS  MAY 2005
*
*     PURPOSE - COMPUTES SPECTRAL VARIANCE ARRAY AND POWER SPECTRUM FROM             
*               1D SPECTRAL COEFFICIENT ARRAY PRODUCED BY THE DFT
*               STANDARD PERIODIC FOURIER TRANSFORM.
*               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*
*     INPUT VARIABLES...                                                   
*                                                                               
*       AMPX  = 1D SPECTRAL COEFFICIENT (ALONG X)  ARRAY FROM TRANSFORM.
*       AMPY  = 1D SPECTRAL COEFFICIENT (ALONG Y) ARRAY FROM TRANSFORM.
*       NIT  = X DIMENSION OF THE SPECTRAL COEFFICIENT ARRAY.
*       NJT  = Y DIMENSION OF THE SPECTRAL COEFFICIENT ARRAY.                                                                     
*       MMAX = MAX X WAVENUMBER  
*       NMAX = MAX Y WAVENUMBER 
*       MAXI = 1D MEMORY DIMENSION

*      OUTPUT VARIABLES...

*       SPECTRE =  1-D POWER SPECTRUM FIELD.
*       VARIANCE=  2-D POWER SPECTRUM FIELD.

*       NB: AMP, SPECTRE AND VARIANCE ARRAYS START AT INDEX "0" FOR
*           WAVENUMBER "0".IN THE MAIN, THEY START AT INDEX "1" 
*           FOR WAVENUMBER "0".
*-----------------------------------------------------------------------

      INTEGER   NIT,NJT,MMAX,NMAX,MAXI
      REAL*8    AMPX(0:NIT-1,0:NJT-1)
      REAL*8    AMPY(0:NIT-1,0:NJT-1)
      REAL*8    VARIANCE(0:MMAX,0:NMAX), SPECTRE(0:MAXI-1)
      REAL*8    R, A, B

*BD<< bugfix to insure proper results  
*BD== REAL      DKM,DKN
      REAL*8    DKM,DKN

*BD== REAL      RCUT_LOW,RCUT_HIGH
      REAL*8    RCUT_LOW,RCUT_HIGH
*BD>>
      
      LOGICAL   WEST_EAST

      INTEGER   RMTRUNC,RNTRUNC
      INTEGER   I,J,K,N,M,NBWV

*-----------------------------------------------------------------------
***    INITIALISATION

      DO  J=0,NMAX      
          DO I=0,MMAX
              VARIANCE(I,J) = 0
          END DO
      END DO

      IF(WEST_EAST) THEN 

***   FOLLOWING X: WE TAKE THE AVERAGE IN Y FOR EACH WAVENUMBER M
         DO M=1,MMAX-1
            DO J=0,NJT-3
               I = 2*M
               VARIANCE(M,0)=VARIANCE(M,0)+ 2.0*((  AMPX(I,  J  )**2 
     +                                            + AMPX(I+1,J  )**2 ))
     +                                         / FLOAT(NJT-2)
            ENDDO
         ENDDO
         
*     THE FIRST ELEMENT (M=0, I.E. THE MEAN) AND THE LAST NEED A SPECIAL 
*     TREATMENT.

         DO J=0,NJT-3
            I = 0
            VARIANCE(0,0)   =VARIANCE(0,0)    + 1.0*( AMPX(I,  J )**2)
     +                                             /  FLOAT(NJT-2)
            I=2*MMAX
            VARIANCE(MMAX,0)=VARIANCE(MMAX,0) + 1.0*( AMPX(I,  J )**2)
     +                                             / FLOAT(NJT-2)
         ENDDO
         
      ELSE                      ! DO THE SOUTH-NORTH CASE.
         
***   FOLLOWING Y: WE TAKE THE AVERAGE IN Y FOR EACH WAVENUMBER N
         
         DO N=1,NMAX-1
            DO I=0,NIT-3
               J = 2*N
               VARIANCE(0,N)=VARIANCE(0,N) + 2.00*( AMPY(I,  J  )**2 
     +                                     +        AMPY(I,  J+1)**2 )
     +                                           / FLOAT(NIT-2)
            ENDDO
         ENDDO

*     THE FIRST ELEMENT (M=0, I.E. THE MEAN) AND THE LAST NEED A SPECIAL 
*     TREATMENT.
         
         DO I=0,NIT-3
            J = 0
            VARIANCE(0,0)=VARIANCE(0,0) + 1.0*(  AMPY(I,  J  )**2 )
     +                                       / FLOAT(NIT-2)
            J=2*NMAX
            VARIANCE(0,NMAX)=VARIANCE(0,NMAX) + 1.0*( AMPY(I, J )**2)
     +                                             / FLOAT(NIT-2)
         ENDDO
         
      ENDIF 
      
***    SECOND: COMPUTE THE POWER SPECTRUM
***    ------

***    WE DEFINE THE NUMBER OF WAVEBAND. SHOULD CORRECTLY TAKE
***    CARE OF NON-SQUARE ARRAYS.
 
      IF(WEST_EAST) THEN
         NBWV = MMAX
         DKM  = MMAX/FLOAT( NBWV ) 
      ELSE
         NBWV = NMAX
         DKN  = NMAX/FLOAT( NBWV )
      ENDIF
***    WE INITIALIZE SPECTRE.
      
      DO  K=0,MAXI-1
          SPECTRE(K) = 0.0
      END DO
      

***    FOR EACH WAVE COMPONANT WE GO THROUGH THE ARRAY OF VARIANCE.

***    EACH SPECTRAL COMPONANT "KS" OF THE POWER SPECTRUM "SPECTRE" IS 
***    MADE FROM THE CONTRIBUTION OF THE VARIANCES THAT HAS A 
***    WAVENUMBER K SUCH THAT  KS-1 < K <= KS.

***         K       |   KS
***     ------------------      
***     0          -- > 0
***     0 < K <= 1 -- > 1
***     1 < K <= 2 -- > 2
***     2 < K <= 3 -- > 3
***       " " ""   """"  

***    ELEMENT K=0 IS THE MODE 0, I.E. THE SQUARE OF THE MEAN.
***    WE GOT IT FROM VARIANCE(0,0):
      
      SPECTRE(0) = VARIANCE(0,0)

      DO  K=1,NBWV

         RCUT_LOW  = (K - 1.0)/FLOAT( K )
         RCUT_HIGH = (K + 0.0)/FLOAT( K )
         
         IF(WEST_EAST) THEN

            RMTRUNC   = DKM*K
            
            DO  M=0,MMAX
               
               R = M/FLOAT( RMTRUNC )
               
               IF ((R.GT.RCUT_LOW).AND.(R.LE.RCUT_HIGH))    THEN
                  SPECTRE(K)  = SPECTRE(K) + VARIANCE(M,0)
               ENDIF
               
            ENDDO
            
         ELSE
            
            RNTRUNC   = DKN*K
            
            DO  N=0,NMAX
               
               R = N/FLOAT( RNTRUNC )
               
               IF ((R.GT.RCUT_LOW).AND.(R.LE.RCUT_HIGH))    THEN
                  SPECTRE(K)  = SPECTRE(K) + VARIANCE(0,N)
               END IF
               
            END DO
         ENDIF
      END DO

      RETURN
*-----------------------------------------------------------------------

      END

      SUBROUTINE SPECT_DCT2D( SPECTRE,VARIANCE,AMP,
     +                        NIT,NJT,MAXI,REGROUPE,SKIP,MOY)

      IMPLICIT  NONE

*     AUTHOR  - B.DENIS  AVRIL 2000
*                                                                               
*     PURPOSE - COMPUTES SPECTRAL VARIANCE ARRAY AND POWER SPECTRUM FROM             
*               2D SPECTRAL COEFFICIENT ARRAY PRODUCED BY THE DCT

*     INPUT VARIABLES...                                                   
*                                                                               
*       AMP  = 2D SPECTRAL COEFFICIENT ARRAY FROM TRANSFORM.
*       NIT  = X DIMENSION OF THE SPECTRAL COEFFICIENT ARRAY.
*       NJT  = Y DIMENSION OF THE SPECTRAL COEFFICIENT ARRAY.
*       MAXI = 1D MEMORY DIMENSION
*   REGROUPE = SWITCH TO GATHER TWO-BY-TWO HALF AND FULL WAVENUMBERS.
*      SKIP  = SWITCH TO SKIP WAVE NUMBER 1/2 WHEN REGROUPE=.TRUE.
*       MOY  = SWITCH TO DO A MOVING AVERAGE.

*     OUTPUT VARIABLES...

*       SPECTRE =  1-D POWER SPECTRUM FIELD.
*       VARIANCE=  2-D POWER SPECTRUM FIELD.

*       NB#1: AMP, SPECTRE AND VARIANCE ARRAYS START AT INDEX "0" FOR
*             WAVENUMBER "0".IN THE MAIN, THEY START AT INDEX "1" 
*             FOR WAVENUMBER "0".

*       NB#2: BECAUSE THE SHIFTED COS/SIN TRANSFORM GIVES FULL *AND*
*             HALF WAVENUMBER COEFFICIENTS, THIS ROUTINE CAN GATHER FULL
*             AND HALF WAVENUMBER VARIANCES OF THE SPECTRUM
*             (SEE SKIP AND REGROUPE SWITCH).
*             THIS SIMULATES SPECTRUM PRODUCED BY THE STANDARD FOURIER
*             TRANSFORM WHERE ONLY FULL WAVE NUMBER ARE PRODUCED.
*-----------------------------------------------------------------------

      INTEGER   NIT,NJT,MAXI
      LOGICAL   SKIP, REGROUPE, MOY
      REAL*8    R, A, B
      REAL*8    AMP(0:NIT-1,0:NJT-1)
      REAL*8    VARIANCE(0:NIT-1,0:NJT-1), SPECTRE(0:MAXI-1)

*BD<< bugfix to insure proper results  
*BD== REAL      DKM,DKN
      REAL*8    DKM,DKN

*BD== REAL      RCUT_LOW,RCUT_HIGH
      REAL*8    RCUT_LOW,RCUT_HIGH
*BD>>
      

      INTEGER   RMTRUNC,RNTRUNC
      INTEGER   I,J,K,NBWV

***    FIRST: COMPUTE THE VARIANCE OF EACH COEFFICIENT FIRST.
      
***    IN TERMS ON OF WAVENUMBERS (FULL AND HALF), THE COEFFICIENTS 
***    ARE ORGANIZED LIKE THIS:

*      ---------------------------------------------------------
*     |     "            "            "           "   
*     | COEF{0, 1 1/2} || COEF{1/2, 1 1/2} | COEF{2, 1 1/2} |   "
*      -----------------------------------------------------
*     | COEF{0, 1    } || COEF{1/2, 1    } | COEF{2, 1    } |   "
*      -----------------------------------------------------
*     | COEF{0, 1/2  } || COEF{1/2, 1/2  } | COEF{1, 1/2  } |   "
*      ====================================================
*     | COEF{0, 0    } || COEF{1/2, 0    } | COEF{1, 0}     |   "
*      ----------------------------------------------------------

***    THESE COEF. ARE IN THE ARRAY "AMP" WITH THIS FORTRAN STRUCTURE.

*      ----------------------------------------------------
*     |     "            "            "           "   
*     | AMP{0,3}       || AMP{1,3}         | AMP{2,3} |   "
*      -----------------------------------------------------
*     | AMP{0,2}       || AMP{1,2}         | AMP{2,2} |   "
*      -----------------------------------------------------
*     | AMP{0,1}       || AMP{1,1}         | AMP{2,1} |   "
*      ====================================================
*     | AMP{0,0}       || AMP{1,0}         | AMP{2,0} |   "
*      -----------------------------------------------------

***    THE VARIANCE IS SIMPLY VARIANCE(I,J)= AMP(I,J)**2.

***    HERE I,J "EVEN" ARE FOR THE FOR THE FULL WAVENUMBERS 
***    AND  I,J  "ODD" ARE FOR THE HALF WAVENUMBERS.

***    FOR NORMALIZATION REASONS THESE FACTORS ARE USED.

*      ----------------------------------------------------
*     |     "            "            "           "   
*     | 1/2*VAR{0,3} || 1/4*VAR{1,3} | 1/4*VAR{2,3} |   "
*      -----------------------------------------------------
*     | 1/2*VAR{0,2} || 1/4*VAR{1,2} | 1/4*VAR{2,2} |   "
*      -----------------------------------------------------
*     | 1/2*VAR{0,1} || 1/4*VAR{1,1} | 1/4*VAR{2,1} |   "
*      ====================================================
*     |   1*VAR{0,0} || 1/2*VAR{1,0} | 1/2*VAR{2,0} |   "
*      -----------------------------------------------------

*-----------------------------------------------------------------------
***    INITIALISATION

      DO  J=0,NJT-1      
          DO  I=0,NIT-1
              VARIANCE(I,J) = 0.0
          END DO
      END DO


***    1- TERMS INVOLVING THE FACTOR "1/4" (UPPER RIGHT CORNER).

      DO  I=1,NIT-1
          DO  J=1,NJT-1 
              VARIANCE(I,J) = 0.25*(AMP(I,J)**2)
          END DO
      END DO

***    2- TERMS INVOLVING THE FACTOR "1/2" IN THE LOWER ROW. 

      DO  I=1,NIT-1
          VARIANCE(I,0) = 0.5*(AMP(I,0)**2)
      END DO

***    3- TERMS INVOLVING THE FACTOR "1/2" IN THE LEFT COLUMN. 

      DO  J=1,NJT-1
          VARIANCE(0,J) = 0.5*(AMP(0,J)**2)
      END DO

***    4- ELEMENT 0,0 (THE MEAN).

      VARIANCE(0,0) = (AMP(0,0)**2)

***    WE DEFINE THE NUMBER OF WAVEBAND (INCLUDING WAVENUMBER ZERO)
***    THIS SHOULD CORRECTLY TAKE CARE OF NON-SQUARE ARRAYS.

      NBWV = MIN(NIT,NJT)
      DKM  = NIT/FLOAT(NBWV)     ! DEFINE INCREMENT IN WAVENUMBER K.
      DKN  = NJT/FLOAT(NBWV)

***    WE INITIALIZE SPECTRE.
      
      DO  K=0,MAXI-1
          SPECTRE(K) = 0.0
      END DO
      
***    FOR EACH WAVE COMPONENT WE GO THROUGH THE ARRAY OF VARIANCE.

***    EACH SPECTRAL COMPONANT "KS" OF THE POWER SPECTRUM "SPECTRE" IS 
***    MADE FROM THE CONTRIBUTION OF THE VARIANCES THAT HAS A 
***    WAVENUMBER K SUCH THAT  KS-1 < K <= KS.

***         K       |   KS
***     ------------------      
***     0          -- > 0
***     0 < K <= 1 -- > 1
***     1 < K <= 2 -- > 2
***     2 < K <= 3 -- > 3
***       " " ""   """"  

***    ELEMENT K=0 IS THE MODE 0, I.E. THE SQUARE OF THE MEAN.
***    WE GOT IT FROM VARIANCE(0,0):

      SPECTRE(0) = VARIANCE(0,0)

      DO  K=1,NBWV-1

          RCUT_LOW  = (K - 1.0)/FLOAT(K)
          RCUT_HIGH = (K + 0.0)/FLOAT(K)

          RMTRUNC   = DKM*K
          RNTRUNC   = DKN*K
         
          DO  J=0,NJT-1
              DO  I=0,NIT-1
               
                  A = I/FLOAT( RMTRUNC )
                  B = J/FLOAT( RNTRUNC )
               
                  R = SQRT(A**2 + B**2)
               
                  IF ((R.GT.RCUT_LOW).AND.(R.LE.RCUT_HIGH))    THEN
                      SPECTRE(K)  = SPECTRE(K) + VARIANCE(I,J)
                  END IF
               
              END DO
          END DO
      END DO
      
      IF(REGROUPE)                                            THEN
         
***   ON REGROUPE LES VARIANCES DEUX PAR DEUX POUR SIMIMULER LE SPECTRE
***   D'UNE TRANSFORMEE DE FOURIER NORMAL, I.E. QUE POUR N POINTS ON AURA N/2 TERMES.
         
         IF(SKIP) THEN
            
***   ON SKIP L'ELEMENT 1/2 ET REGROUPE 1&11/2-->SPECTRE(1)
***   ET L'ELEMENT 2&21/2--> SPECTRE(2), ECT
***   DANS CE CAS L'ELEMENT 1/2 EST ELIMINE DU SPECTRE.
            
            DO K=0,(NBWV/2)-2
               SPECTRE(K+1)=SPECTRE(2*K+2)+ SPECTRE(2*K+3)
            ENDDO
            
         ELSE 
            
            IF(MOY) THEN
               
***   ON REGROUPE 0.5*N1/2 + N + 0.5*N3/2
***   CA FAIT UNE MOYENNE "MOVING".
               
               DO K=0,(NBWV/2)-2
                  SPECTRE(K+1)=0.5*SPECTRE(2*K+1)+ SPECTRE(2*K+2)
     +                 + 0.5*SPECTRE(2*K+3)
               ENDDO
               
               SPECTRE(NBWV/2)=0.5*SPECTRE(NBWV-1)
               
            ELSE
***   ON REGROUPE 1/2&1--> SPECTRE(1), 11/2&2--> SPECTRE(2),
***   LE DERNIER ELEMENT SE RETROUVE SEUL.
               
               DO K=0,(NBWV/2)-2
                  SPECTRE(K+1)=SPECTRE(2*K+1)+ SPECTRE(2*K+2)
               ENDDO
               
               SPECTRE(NBWV/2)=SPECTRE(NBWV-1)
               
            ENDIF
         ENDIF
      ENDIF
      
      RETURN
*-----------------------------------------------------------------------
      
      END
      SUBROUTINE SPECT_DCT1D( SPECTRE,VARIANCE,AMPX,AMPY,
     +                        NIT,NJT,MAXI,REGROUPE,SKIP,MOY,WEST_EAST)

      IMPLICIT  NONE

*     AUTHOR  - B.DENIS  MAY 2005
*                                                                               
*     PURPOSE - COMPUTES SPECTRAL VARIANCE ARRAY AND POWER SPECTRUM FROM             
*               1D SPECTRAL COEFFICIENT ARRAY PRODUCED BY THE DCT
*
*     NB: THIS ROUTINE WAS INSPIRED BY THE ROUTINE SPECT_DCT2D.

*     INPUT VARIABLES...                                                   
*                                                                               
*       AMPX  = 1D SPECTRAL COEFFICIENT ARRAY FROM TRANSFORM ALONG THE 
*               THE WEST-EAST DIRECTION.
*       AMPY  = 1D SPECTRAL COEFFICIENT ARRAY FROM TRANSFORM ALONG THE
*               THE WEST-EAST DIRECTION.
*       NIT  = X DIMENSION OF THE SPECTRAL COEFFICIENT ARRAY.
*       NJT  = Y DIMENSION OF THE SPECTRAL COEFFICIENT ARRAY.
*       MAXI = 1D MEMORY DIMENSION
*   REGROUPE = SWITCH TO GATHER TWO-BY-TWO HALF AND FULL WAVENUMBERS.
*      SKIP  = SWITCH TO SKIP WAVE NUMBER 1/2 WHEN REGROUPE=.TRUE.
*       MOY  = SWITCH TO DO A MOVING AVERAGE.
*   WEST_EAST= SWITCH TO PROCESS WEST-EAST (OR NORTH-SOUTH IF FALSE) SPECTRA.

*     OUTPUT VARIABLES...

*       SPECTRE =  1-D POWER SPECTRUM FIELD.
*       VARIANCE=  2-D POWER SPECTRUM FIELD.

*       NB#1: AMP, SPECTRE AND VARIANCE ARRAYS START AT INDEX "0" FOR
*             WAVENUMBER "0".IN THE MAIN, THEY START AT INDEX "1" 
*             FOR WAVENUMBER "0".

*       NB#2: BECAUSE THE DCT TRANSFORM GIVES FULL *AND*
*             HALF WAVENUMBER COEFFICIENTS, THIS ROUTINE CAN GATHER FULL
*             AND HALF WAVENUMBER VARIANCES OF THE SPECTRUM
*             (SEE SKIP AND REGROUPE SWITCH).
*             THIS SIMULATES SPECTRUM PRODUCED BY THE STANDARD FOURIER
*             TRANSFORM WHERE ONLY FULL WAVE NUMBER ARE PRODUCED.
*-----------------------------------------------------------------------

      INTEGER   NIT,NJT,MAXI
      LOGICAL   SKIP, REGROUPE, MOY,WEST_EAST
      REAL*8    R, A, B
      REAL*8    AMPX(0:NIT-1,0:NJT-1)
      REAL*8    AMPY(0:NIT-1,0:NJT-1)
      REAL*8    VARIANCE(0:NIT-1,0:NJT-1), SPECTRE(0:MAXI-1)

*BD<< bugfix to insure proper results  
*BD== REAL      DKM,DKN
      REAL*8    DKM,DKN

*BD== REAL      RCUT_LOW,RCUT_HIGH
      REAL*8    RCUT_LOW,RCUT_HIGH
*BD>>

      INTEGER   RMTRUNC,RNTRUNC
      INTEGER   I,J,K,NBWV

*-----------------------------------------------------------------------
***    INITIALISATION

      DO  J=0,NJT-1      
          DO  I=0,NIT-1
              VARIANCE(I,J) = 0.0
          END DO
      END DO


      IF (WEST_EAST) THEN
         
***   FOLLOWING X: WE TAKE THE AVERAGE IN Y FOR EACH WAVENUMBER 
***   OF THE SPECTRAL COMPUTED ALONG THE WEST-EAST DIRECTION
         
***   THE FIRST ELEMENTS (MEAN**2) NEED A SPECIAL TREATMENT.
         
         DO  J=0,NJT-1
            VARIANCE(0,0) =  VARIANCE(0,0) + 1.0*(AMPX(0,J)**2)
     +                                          / FLOAT(NJT)
         END DO
         
***   DO THE REST OF THE SPECTRA ELEMENTS
         DO  I=1,NIT-1
            DO  J=0,NJT-1 
               VARIANCE(I,0) = VARIANCE(I,0)+ 0.5*(AMPX(I,J)**2)
     +                                           / FLOAT(NJT)
            END DO
         END DO

      ELSE ! THEN DO THE SOUTH-NORTH DIRECTION 
         
***   FOLLOWING Y: WE TAKE THE AVERAGE IN X FOR EACH WAVENUMBER 
***   OF THE SPECTRAL COMPUTED ALONG THE SOUTH-NORTH DIRECTION
         
         DO  I=0,NIT-1
            VARIANCE(0,0) =  VARIANCE(0,0) + 1.0*(AMPY(I,0)**2)
     +                                          / FLOAT(NIT)
         END DO
         
         DO  J=1,NJT-1 
            DO  I=0,NIT-1
               VARIANCE(0,J) = VARIANCE(0,J)+ 0.5*(AMPY(I,J)**2)
     +                                          / FLOAT(NIT)
            END DO
         END DO
      ENDIF
      
***   NB1:AT THIS POINT THE 2D-ARRAY 'VARIANCE' IS ZERO EVERYWHERE
***       EXCEPT ALONG ITS FIRST ROW FOR THE WEST-EAST CASE OR 
***       ALONG THE ITS FIRST COLUMN FOR THE SOUTH-NORTH CASE.

***   NB2:THIS LAST PART OF THE CODE COMES FROM THE ROUTINE SPECT_DCT2D.
***       SINCE IT WAS GENERAL ENOUGH, IT DID NOT NEED MODIFICATIONS 
***       FOR THESE 1D CASES. 
*** 

***    WE DEFINE THE NUMBER OF WAVEBAND (INCLUDING WAVENUMBER ZERO)
***    THIS SHOULD CORRECTLY TAKE CARE OF NON-SQUARE ARRAYS.

      IF (WEST_EAST) THEN
         NBWV = NIT
         DKM  = NIT/FLOAT(NBWV) 
      ELSE
         NBWV = NJT
         DKN  = NJT/FLOAT(NBWV)
      ENDIF

***    WE INITIALIZE SPECTRE.
      
      DO  K=0,MAXI-1
          SPECTRE(K) = 0.0
      END DO
      
***    FOR EACH WAVE COMPONENT WE GO THROUGH THE ARRAY OF VARIANCE.

***    EACH SPECTRAL COMPONANT "KS" OF THE POWER SPECTRUM "SPECTRE" IS 
***    MADE FROM THE CONTRIBUTION OF THE VARIANCES THAT HAS A 
***    WAVENUMBER K SUCH THAT  KS-1 < K <= KS.

***         K       |   KS
***     ------------------      
***     0          -- > 0
***     0 < K <= 1 -- > 1
***     1 < K <= 2 -- > 2
***     2 < K <= 3 -- > 3
***       " " ""   """"  

***    ELEMENT K=0 IS THE MODE 0, I.E. THE SQUARE OF THE MEAN.
***    WE GOT IT FROM VARIANCE(0,0):

      SPECTRE(0) = VARIANCE(0,0)

      DO  K=1,NBWV-1

          RCUT_LOW  = (K - 1.0)/FLOAT(K)
          RCUT_HIGH = (K + 0.0)/FLOAT(K)

          IF(WEST_EAST) THEN
             
             RMTRUNC   = DKM*K
             
             DO  I=0,NIT-1
                
                R = I/FLOAT( RMTRUNC )
                
                IF ((R.GT.RCUT_LOW).AND.(R.LE.RCUT_HIGH))    THEN
                   SPECTRE(K)  = SPECTRE(K) + VARIANCE(I,0)
                END IF
                
             END DO
          
          ELSE
             
             RNTRUNC   = DKN*K
             
             DO  J=0,NJT-1
                
                R = J/FLOAT( RNTRUNC )
                
                IF ((R.GT.RCUT_LOW).AND.(R.LE.RCUT_HIGH))    THEN
                   SPECTRE(K)  = SPECTRE(K) + VARIANCE(0,J)
                END IF
                
             END DO

          ENDIF
       END DO
     
      IF(REGROUPE)                                            THEN
         
***   ON REGROUPE LES VARIANCES DEUX PAR DEUX POUR SIMIMULER LE SPECTRE
***   D'UNE TRANSFORMEE DE FOURIER NORMAL, I.E. QUE POUR N POINTS ON AURA N/2 TERMES.
         
         IF(SKIP) THEN
            
***   ON SKIP L'ELEMENT 1/2 ET REGROUPE 1&11/2-->SPECTRE(1)
***   ET L'ELEMENT 2&21/2--> SPECTRE(2), ECT
***   DANS CE CAS L'ELEMENT 1/2 EST ELIMINE DU SPECTRE.
            
            DO K=0,(NBWV/2)-2
               SPECTRE(K+1)=SPECTRE(2*K+2)+ SPECTRE(2*K+3)
            ENDDO
            
         ELSE 
            
            IF(MOY) THEN
               
***   ON REGROUPE 0.5*N1/2 + N + 0.5*N3/2
***   CA FAIT UNE MOYENNE "MOVING".
               
               DO K=0,(NBWV/2)-2
                  SPECTRE(K+1)=0.5*SPECTRE(2*K+1)+ SPECTRE(2*K+2)
     +                 + 0.5*SPECTRE(2*K+3)
               ENDDO
               
               SPECTRE(NBWV/2)=0.5*SPECTRE(NBWV-1)
               
            ELSE
***   ON REGROUPE 1/2&1--> SPECTRE(1), 11/2&2--> SPECTRE(2),
***   LE DERNIER ELEMENT SE RETROUVE SEUL.
               
               DO K=0,(NBWV/2)-2
                  SPECTRE(K+1)=SPECTRE(2*K+1)+ SPECTRE(2*K+2)
               ENDDO
               
               SPECTRE(NBWV/2)=SPECTRE(NBWV-1)
               
            ENDIF
         ENDIF
      ENDIF
      
      RETURN
*-----------------------------------------------------------------------
      
      END

      SUBROUTINE TRANS1D(AR,U1,U2,U3,MAXSIZE,N,NLOT,ICAS,IAXE,IWAY)

      IMPLICIT NONE

C     * THIS ROUTINE DOES THE SPECTRAL TRANSFORMS (FORWARD AND INVERSE)
C     * IT WORKS ON MULTIPLE ROW OR COLUMN AT A TIME.

C     * IT CAN DO TWO KINDS OF TRANSFORM:

C                 1 => DCT (COS SHIFTED)
C                 2 => DFT (REGULAR SIN & COS).
C
      INTEGER MAXSIZE,N,NLOT,ICAS,IAXE,IWAY
      REAL*8, dimension(*) :: U1,U2,U3     ! no longer used
      REAL*8 AR(MAXSIZE)

      REAL*8 R1(MAXSIZE)                   ! automatic array
      INTEGER NX,NY, NN, LOT, NIC

      REAL*8 CNORM, FAC

      REAL*8 DEL, ANGLE

      REAL*8 ZERO, HALF, ONE, TWO
      PARAMETER( ZERO = 0.0 )
      PARAMETER( HALF = 0.5 )
      PARAMETER( ONE  = 1.0 )
      PARAMETER( TWO  = 2.0 )

      LOGICAL FAST_FT
      LOGICAL LPRINT,  LDCT, LDFT
      INTEGER I, J, K, INC, JUMP
      INTEGER INCR, JUMPR, NFACT, LOTR
!       INTEGER NDX, NDR
      INTEGER NDR

!       NDX(I,J) =  1 + ( J - 1 ) * JUMP + I * INC
      NDR(I,J) =  1 + ( J - 1 ) * JUMPR + I * INCR

      LPRINT  = .FALSE.

      LDCT   = ICAS .EQ. 1
      LDFT    = ICAS .EQ. 2

      IF (LPRINT) THEN 
      PRINT *,'*********************'
      IF ( LDCT ) PRINT *,'*  QCFT8, QCFFT8    *'
      IF ( LDFT  ) PRINT *,'*   RFT8,  FFFT8    *'
      PRINT *,'*********************'
      ENDIF


C     * SELECT THE FAST OR THE SLOW FOURIER TRANSFORM 
C     * BASED ON FACTORIZATION.
      
      NFACT=N
      CALL NGFFT( NFACT )

      IF(NFACT.EQ.N) THEN 
         FAST_FT=.TRUE.
         IF ( LDCT ) CALL SETSCQR( N  , 'QCOS' )
         IF ( LDFT  ) CALL SETSCQR( N  , 'REAL' )
         
         IF (IAXE.EQ.0) THEN
            WRITE (6,*) ' '
            WRITE (6,*)'    USING FAST TRANSFORM IN X'
         ELSE
            WRITE (6,*)'    USING FAST TRANSFORM IN Y'
         ENDIF
         
      ELSE
         FAST_FT=.FALSE.
         IF (IAXE.EQ.0) THEN
            WRITE (6,*) ' '
            WRITE (6,*) '    USING SLOW TRANSFORM IN X !!! '
         ELSE
            WRITE (6,*) ' '
            WRITE (6,*) '    USING SLOW TRANSFORM IN Y !!! '
         ENDIF
      ENDIF

      CNORM = SQRT( TWO/N )
      DEL   = ACOS( - ONE )/N

      NIC  = N + 1

         IF     ( IAXE .EQ. 0 ) THEN

C     * ON VEUT PAR EXEMPLE:

C     APPEL DE TRANS1D    LISTE D'ARGUMENTS  APPEL POUR
C     DANS MAIN           DE TRANS1D         ROUTINE FFT
C     =================   =================  ===========
C
C     NX=NID    -->        N      -->        N=JUMP=N=NID
C     NY=NJD    -->        NLOT   -->        LOT=NLOT=NJD 
C                                            INC=1
            LOT  = NLOT
            INC   = 1
            JUMP  = N
            INCR  = 1
            JUMPR = N + 2
         ELSE

C     APPEL DE TRANS1D    LISTE D'ARGUMENTS  APPEL POUR
C     DANS MAIN           DE TRANS1D         ROUTINE FFT
C     =================   =================  ===========
C
C     NX=NJD    -->        N      -->           N=NJD
C     NY=NID    -->        NLOT   -->         LOT=INC=NLOT=NID 
C                                            JUMP=1
            LOT   = NLOT
            LOTR  = NLOT + 2    
            INC   = NLOT  
            JUMP  = 1
            INCR  = LOTR 
            JUMPR = 1
         ENDIF

C------------------------------------------------------------------
C
C    CASE OF QCFT8, QCFFT8
C
C------------------------------------------------------------------
         IF ( LDCT ) THEN

           IF (IWAY.EQ.-1) THEN             !  GRIDPOINT -> SPECTRAL
               IF (FAST_FT) THEN
                  CALL QCFFT8(AR, INC,JUMP,LOT,-1)
               ELSE
                  CALL QCFT8 (R1,AR,INC,JUMP,LOT,-1,N )
                  DO K=1,N*LOT
                     AR(K)=R1(K)
                  ENDDO
               ENDIF

            ELSE                             ! SPECTRAL -> GRIDPOINT

               IF (FAST_FT) THEN
                  CALL QCFFT8(AR, INC,JUMP,LOT,+1)
               ELSE
                  CALL QCFT8 (R1,AR,INC,JUMP,LOT,+1,N )
                  DO K=1,N*LOT
                     AR(K)=R1(K)
                  ENDDO
               ENDIF
            ENDIF

         ENDIF
C------------------------------------------------------------------
C
C CASE OF RFT8, FFFT8
C
C------------------------------------------------------------------
         IF ( LDFT ) THEN

            IF ( LPRINT ) THEN
               PRINT *,' ON IMPRIME AR AVANT TRANSFORME'
               DO J=1,LOT
                  PRINT *,'J = ',J,(AR(NDR(I,J)),I=0,2*(N/2)+1)
               ENDDO
            ENDIF
         
            IF (IWAY.EQ.-1) THEN
               IF (FAST_FT) THEN                   ! GRIDPOINT -> SPECTRAL
                  CALL FFFT8(AR, INCR,JUMPR,LOT,-1) 
               ELSE
                  CALL RFT8 ( R1,AR,INCR,JUMPR,LOT,-1,N)
                  DO K=1,(N+2)*LOT
                     AR(K)=R1(K)
                  ENDDO
               ENDIF
               
            ELSE                                    ! SPECTRAL -> GRIDPOINT
               IF (FAST_FT) THEN
                  CALL FFFT8(AR, INCR,JUMPR,LOT,+1) 
               ELSE
                  CALL RFT8 (R1,AR,INCR,JUMPR,LOT,+1,N)
                  DO K=1,(N+2)*LOT
                     AR(K)=R1(K)
                  ENDDO
               ENDIF
            ENDIF
            
            IF ( LPRINT ) THEN
               PRINT *,' ON IMPRIME AR APRES TRANSFORME'
               DO J=1,LOT
                  PRINT *,'J = ',J,(AR(NDR(I,J)),I=0,2*(N/2)+1)
               ENDDO
            ENDIF

         ENDIF

      RETURN 
      END

      module COMFFT8_MOD
      INTEGER, save :: N, M, NSTORE
      REAL*8, dimension(:), POINTER, save ::  SSIN, CCOS, QSIN
      logical, save :: is_allocated = .false.
      end module COMFFT8_MOD

C     $LOG: TRANS1D.FTN,V $
C     REVISION 3.1  2003/10/24 21:05:48  DUGAS
C     IMPLEMENTER DU CODE COMPATIBLE RS6000
C
C     REVISION 3.0  2000/07/24 20:39:14  ARMNRBD
C     VERSION INITIALE (DE JEAN COTE)
C     SUBROUTINE 'SETSCQR' - SETS UP COMMON 'COMFFT8X' REQUIRED BY
C                            CFFT8, SFFT8, QCFFT8, QSFFT8
C                                       AND
C                            CALLS SETFFT8 TO SET UP 'COMFFT8'
C                            REQUIRED BY FFFT8 WHICH IS USED BY 
C                            THE 4 PREVIOUS TRANSFORMS
C
CAUTHOR JEAN COTE SEPT 99
C
CARGUMENTS
C   I      - NF        - NUMBER OF GRID POINTS (LENGTH OF TRANSFORM)
C   I      - CASE      - CODED NAME OF DESIRED TRANSFORM
C                            CASE = 'COS'  FOR REAL COSINE         (CFFT8)
C                            CASE = 'SIN'  FOR REAL SINE           (SFFT8)
C                            CASE = 'QCOS' FOR REAL SHIFTED COSINE (QCFFT8)
C                            CASE = 'QSIN' FOR REAL SHIFTED SINE   (QSFFT8)
C                            CASE = 'REAL' FOR REAL PERIODIC       (FFFT8)
C
CC
C
      SUBROUTINE SETSCQR( NF, CASE )

      use COMFFT8_MOD
      IMPLICIT NONE

      INTEGER NF
      CHARACTER*(*) CASE

!       INTEGER N, M, NSTORE
!       REAL*8  SSIN, CCOS, QSIN
!       POINTER ( PTSS,SSIN(N-M-1) ), ( PTCC,CCOS(N-M-1) )
!       POINTER ( PTQS,QSIN(0:M-1) )
! 
!       COMMON    / COMFFT8X / PTSS, PTCC, PTQS, N, M, NSTORE
C
      INTEGER    I, IER
      REAL *8    DEL, ANGLE
      CHARACTER  ALLOUE*17
C
      REAL*8 ZERO, HALF, ONE, TWO, FOUR
      PARAMETER( ZERO  = 0.0 )
      PARAMETER( HALF  = 0.5 )
      PARAMETER( ONE   = 1.0 )
      PARAMETER( TWO   = 2.0 )
      PARAMETER( FOUR  = 4.0 )
C
!       DATA ALLOUE/'PAS ENCORE ALLOUE'/

!       IF (ALLOUE.NE.'PAS ENCORE ALLOUE') CALL HPDEALLC( PTSS,IER,0 )
!       IF (ALLOUE.NE.'PAS ENCORE ALLOUE') CALL HPDEALLC( PTCC,IER,0 )
!       IF (ALLOUE.NE.'PAS ENCORE ALLOUE') CALL HPDEALLC( PTQS,IER,0 )
!       IF (ALLOUE.NE.'PAS ENCORE ALLOUE') ALLOUE = 'DEJA ALLOUE'
      if(is_allocated) then
        deallocate(SSIN,QSIN,CCOS)
        is_allocated = .false.
      endif

C     N  =  LENGTH OF AUXILIARY REAL PERIODIC FOURIER TRANSFORM (FFFT8)

      IF     ( CASE .EQ. 'SIN' ) THEN
         N = NF + 1
      ELSEIF ( CASE .EQ. 'COS' ) THEN
         N = NF - 1
      ELSEIF ( CASE .EQ. 'REAL'   .OR.
     %         CASE .EQ. 'QSIN'   .OR.
     %         CASE .EQ. 'QCOS' ) THEN
         N = NF
      ELSE
         PRINT *,'ERROR IN SETSCQR -> CASE = ',CASE
         RETURN
      ENDIF
      
      M      = N/2
      NSTORE = N + 2

      allocate(SSIN(N-M-1),QSIN(0:M-1),CCOS(N-M-1))
      is_allocated = .true.

      DEL = ACOS( - ONE )/N

      DO I=1,N-M-1

         ANGLE = I * DEL
         CCOS( I ) = COS( ANGLE )
         SSIN( I ) = SIN( ANGLE )

      ENDDO

      DO I=0,M-1

         QSIN( I ) = SIN( ( I + HALF ) * DEL )

      ENDDO

      CALL SETFFT8( N )

      RETURN
      END
C     SUBROUTINE 'CFFT8' - MULTIPLE FAST REAL COSINE TRANSFORM

C     REAL COSINE TRANSFORM OF LENGTH N+1
C     THE TRANSFORM IS ITS OWN INVERSE
C     SELF INVERTING IMPLEMENTATION OF NUMERICAL RECIPES PP. 508-512

C     CREATED: SEPT 99 BY J.COTE, RPN

C     A     IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     INC   IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C          (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP  IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     LOT   IS THE NUMBER OF DATA VECTORS

C     DEFINITION OF TRANSFORM:
C     -------------------------

C     R(K) = SQRT(2/N)*SUM(I=1,...,N-1)(A(I)*COS(I*K*PI/N))
C           +SQRT(1/(2N))*( A(0) + (-1)**K * A(N) )
C
C     NOTE FOR 'A' STORED AS A(N1,N2) THEN
C
C        FOR A TRANSFORM ALONG THE FIRST DIMENSION
C
C           INC   = 1
C           JUMP  = N1
C
C        FOR A TRANSFORM ALONG THE SECOND DIMENSION
C
C           INC   = N1
C           JUMP  = 1
C
C     THE SUBROUTINE SETSCQR MUST HAVE BEEN CALLED TO SET-UP
C     THE COMMONS COMFFT8 AND COMFFT8X
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CFFT8( A, INC, JUMP, LOT )

      use COMFFT8_MOD
      IMPLICIT NONE

      INTEGER INC, JUMP, LOT
      REAL*8  A(*)

      REAL*8  AI, AS, YA, YS, XNOR
      INTEGER I, J, K, IS

      INTEGER J0, JLOT

      integer, parameter :: MAXLOT=511
      REAL*8         W(MAXLOT*NSTORE)

      REAL*8 PDWRD(2)
      REAL   PRWRD(2)

      REAL*8 ZERO, HALF, ONE, TWO, FOUR
      PARAMETER( ZERO = 0.0 )
      PARAMETER( HALF = 0.5 )
      PARAMETER( ONE  = 1.0 )
      PARAMETER( TWO  = 2.0 )
      PARAMETER( FOUR = 4.0 )

!     statement functions used for indexing
      INTEGER IJA, IJW
      IJA(I,J) = 1 + (J0+J-1)*JUMP + I*INC 
      IJW(I,J) = J + I*MAXLOT 

      XNOR = SQRT( HALF * N )

      DO 100 J0=0,LOT-1,MAXLOT
      JLOT  = MIN( MAXLOT, LOT - J0 )

      I  = 1
      IS = N - I
      DO J=1,JLOT
         AI= A( IJA(I ,J) )
         AS= A( IJA(IS,J) )
         YA = ( AS - AI )
         A( IJA( 1,J) ) = - YA * CCOS( I )
         YA = TWO *  SSIN( I ) * YA * XNOR
         YS = ( AS + AI ) * XNOR 
         W( IJW(I ,J) ) = YS + YA
         W( IJW(IS,J) ) = YS - YA
         YS = ( A( IJA(N,J) ) + A( IJA(0,J) ) ) * XNOR
         A( IJA(1,J) ) = A( IJA(1,J) ) -
     %                   HALF * ( A( IJA(N,J) ) - A( IJA(0,J) ) )
         W( IJW(0,J) ) = YS
         W( IJW(N,J) ) = YS
      ENDDO

      DO I=2,N-M-1

         IS = N - I
         DO J=1,JLOT
            AI= A( IJA(I ,J) )
            AS= A( IJA(IS,J) )
            YA = ( AS - AI )
            A( IJA( 1,J) ) = A( IJA(1,J) ) - YA * CCOS( I )
            YA = TWO *  SSIN( I ) * YA * XNOR
            YS = ( AS + AI ) * XNOR 
            W( IJW(I ,J) ) = YS + YA
            W( IJW(IS,J) ) = YS - YA
         ENDDO

      ENDDO

      DO J=1,JLOT
         A( IJA(1,J) ) = A( IJA(1,J) )/XNOR
      ENDDO

      IF ( N .EQ. 2 * M ) THEN
         DO J=1,JLOT
            W( IJW( M,J) ) = TWO * A( IJA( M,J) ) * XNOR
         ENDDO
      ENDIF

      CALL FFFT8( W, MAXLOT, 1, JLOT, -1 )

      DO J=1,JLOT
         A( IJA(0,J) ) = W( IJW(0,J) )
      ENDDO

      DO K=1,M-1
         DO J=1,JLOT
            A( IJA(2*K  ,J) ) = W( IJW(2*K  ,J) )
            A( IJA(2*K+1,J) ) = A( IJA(2*K-1,J) ) - W( IJW(2*K+1,J) )
         ENDDO
      ENDDO

      DO J=1,JLOT
         A( IJA(2*M,J) ) = W( IJW(2*M,J) )
      ENDDO

      IF ( N .NE. 2 * M ) THEN
         DO J=1,JLOT
            A( IJA(N,J) ) = A( IJA(N-2,J) ) - W( IJW(N,J) )
         ENDDO
      ENDIF

  100 CONTINUE
      RETURN
      END
C     SUBROUTINE 'QCFFT8' - MULTIPLE FAST REAL SHIFTED COSINE TRANSFORM

C     REAL SHIFTED COSINE TRANSFORM OF LENGTH N
C     IMPLEMENTATION INSPIRED BY NUMERICAL RECIPES PP. 513
C     BUT WITH A DIFFERENT NORMALIZATION

C     CREATED: SEPT 99 BY J.COTE, RPN

C     A     IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     INC   IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C          (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP  IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     LOT   IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM FOURIER TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO FOURIER

C     DEFINITION OF TRANSFORM:
C     -------------------------

C     ISIGN=+1: R(I) = SUM(K=0,...,N-1)(A(K)*COS((I+1/2)*K*PI/N))

C     ISIGN=-1: R(K) = SUM(I=0,...,N-1)(A(I)*COS((I+1/2)*K*PI/N))
C                      * ((2-DELTA(K,0))/N)

C     NOTE FOR 'A' STORED AS A(N1,N2) THEN
C
C        FOR A TRANSFORM ALONG THE FIRST DIMENSION
C
C           INC   = 1
C           JUMP  = N1
C
C        FOR A TRANSFORM ALONG THE SECOND DIMENSION
C
C           INC   = N1
C           JUMP  = 1
C
C     THE SUBROUTINE SETSCQR MUST HAVE BEEN CALLED TO SET-UP
C     THE COMMONS COMFFT8 AND COMFFT8X
C
C-----------------------------------------------------------------------
C
      SUBROUTINE QCFFT8( A, INC, JUMP, LOT, ISIGN )

      use COMFFT8_MOD
      IMPLICIT NONE

      INTEGER INC, JUMP, LOT, ISIGN
      REAL*8  A(*)

      REAL*8  AI, AS, YA, YS, C, S, RR
      INTEGER I, J, K, IS, K1, KK

      INTEGER J0, JLOT

      integer, parameter :: MAXLOT=511
      REAL*8         W(MAXLOT*NSTORE)

      REAL*8 PDWRD(2)
      REAL   PRWRD(2)

      REAL*8 ZERO, HALF, ONE, TWO, FOUR
      PARAMETER( ZERO = 0.0 )
      PARAMETER( HALF = 0.5 )
      PARAMETER( ONE  = 1.0 )
      PARAMETER( TWO  = 2.0 )
      PARAMETER( FOUR = 4.0 )

!     statement functions used for indexing
      INTEGER IJA, IJW
      IJA(I,J) = 1 + (J0+J-1)*JUMP + I*INC
      IJW(I,J) = J + I*MAXLOT

      DO 100 J0=0,LOT-1,MAXLOT
      JLOT  = MIN( MAXLOT, LOT - J0 )

      IF ( ISIGN .EQ. -1 ) THEN
C
C     TRANSFORM FROM GRIDPOINT TO FOURIER
C
         DO I = 0 , M-1

            IS = N - I - 1

            DO J=1,JLOT
               AI = A( IJA(I ,J) )
               AS = A( IJA(IS,J) )
               YS = AI + AS
               YA =  TWO *  QSIN( I ) * ( AS - AI )
               W( IJW(I ,J) ) = YS + YA
               W( IJW(IS,J) ) = YS - YA
            ENDDO

         ENDDO

         IF ( N .NE. 2 * M )  THEN
            DO J=1,JLOT
               W( IJW(M,J) ) = TWO * A( IJA(M,J) )
            ENDDO
         ENDIF

         CALL FFFT8( W, MAXLOT, 1, JLOT, -1 )

         DO J=1,JLOT
            A( IJA(0,J) ) = W( IJW(0,J) ) * HALF
         ENDDO

         DO K = 1 , M

            KK = 2*K
            K1 = KK + 1

            IF ( K .LT. M .OR. N .NE. 2 * M ) THEN
               C = CCOS( K )
               S = SSIN( K )
               DO J=1,JLOT
                  A( IJA(KK-1,J) ) = -S*W( IJW(KK,J) )+C*W( IJW(K1,J) )
                  A( IJA(KK  ,J) ) =  C*W( IJW(KK,J) )+S*W( IJW(K1,J) )
               ENDDO
            ELSE
               DO J=1,JLOT
                  A( IJA(KK-1,J) ) = -W( IJW(KK,J) )
               ENDDO
            ENDIF

         ENDDO
         IF ( N .EQ. 2 * M )  THEN
            DO J=1,JLOT
               A( IJA(N-1,J) ) = A( IJA(N-1,J) ) * HALF
            ENDDO
         ENDIF
         DO K = 2*M-3 , 1 , -2
            DO J=1,JLOT
               A( IJA(K,J) ) = A( IJA(K,J) ) + A( IJA(K+2,J) )
            ENDDO
         ENDDO


      ELSEIF ( ISIGN .EQ. +1 ) THEN
C
C     TRANSFORM FROM FOURIER TO GRIDPOINT
C
         DO J=1,JLOT
            W( IJW(0,J) ) = A( IJA(0,J) )
            W( IJW(1,J) ) = ZERO
         ENDDO

         DO K = 2 , N-1 , 2
            DO J=1,JLOT
               W( IJW(K,J) ) = A( IJA(K,J) ) * HALF
            ENDDO
         ENDDO
         DO K = 3 , 2*M-1 , 2
            DO J=1,JLOT
              W( IJW(K,J) ) = ( A( IJA(K-2,J) ) - A( IJA(K,J) ) ) * HALF
            ENDDO
         ENDDO
         IF ( N .EQ. 2 * M )  THEN
            C = ONE
         ELSE
            C = HALF
         ENDIF
         DO J=1,JLOT
            W( IJW(2*M+1,J) ) = A( IJA(2*M-1,J) ) * C
         ENDDO

         DO K = 1 , M

            IF ( K .LT. M .OR. N .NE. 2 * M ) THEN
               C = CCOS( K )
               S = SSIN( K )
            ELSE
               C = ZERO
               S = ONE
            ENDIF

            KK = 2*K
            K1 = KK + 1
            DO J=1,JLOT
               RR = W( IJW(KK,J) )
               W( IJW(KK,J) ) =  C * RR - S * W( IJW(K1,J) )
               W( IJW(K1,J) ) =  S * RR + C * W( IJW(K1,J) )
            ENDDO

         ENDDO

         CALL FFFT8( W, MAXLOT, 1, JLOT, +1 )

         DO I = 0 , M-1

            IS = N - I - 1

            DO J=1,JLOT
               YS = ( W( IJW(I ,J) ) + W( IJW(IS,J) ) ) * HALF
               YA = ( W( IJW(IS,J) ) - W( IJW(I ,J) ) ) /
     %              ( FOUR * QSIN( I ) )
               A( IJA(I ,J) ) = YS + YA
               A( IJA(IS,J) ) = YS - YA

            ENDDO

         ENDDO
         IF ( N .NE. 2 * M )  THEN
            DO J=1,JLOT
               A( IJA(M,J) ) = W( IJW(M,J) )
            ENDDO
         ENDIF

      ENDIF

  100 CONTINUE
      RETURN
      END
CCCS/R NGFFT - CALCUL DU PROCHAIN ENTIER DANS LA SUITE 4, 6, 8
C              OU QUI SE FACTORISE EN 2, 3, ET 5 QUAND > QUE 8
C              POUR FFT771 ET FFFT8, S(Q)FFT8, C(Q)FFT8
C
      SUBROUTINE NGFFT( N )

      IMPLICIT NONE
      INTEGER N
C
CAUTEUR JEAN COTE - 1990
C
CREVISION JEAN COTE - SEPT 1999, AJOUTE 4, 6, 8
C
CARGUMENTS
C   IO     - N       - EN SORTIE LE PLUS PETIT ENTIER <= N QUI FACTORISE
C
CPARAMETRES
      INTEGER L
      PARAMETER ( L = 3 )
      INTEGER K( L )
      DATA K / 2 , 3 , 5 /
C
CC
      INTEGER I,J
C

      IF ( N .LE. 8 ) THEN
         IF ( N .LE. 4 ) N = 4
         IF ( N .EQ. 5 ) N = 6
         IF ( N .EQ. 7 ) N = 8
         RETURN
      ENDIF

C<<BD MODIF: WE LOOK AT  THE NEAREST FACTORIZABLE NUMBER SMALLER THAN N INSTEAD OF GREATER
CBD==     N = N - 1
CBD==   1 N = N + 1  
      N = N + 1
    1 N = N - 1  
CBD>>
      I = N
    2 DO 3 J=1,L
         IF( MOD(I,K(J)) .EQ. 0 ) GO TO 4
    3 CONTINUE
      GO TO 1
    4 I = I/K(J)
      IF( I .NE. 1 ) GO TO 2
      RETURN
      END

C     SUBROUTINE 'QCFT8' - MULTIPLE SLOW REAL SHIFTED COSINE TRANSFORM

C     REAL SHIFTED COSINE TRANSFORM OF LENGTH N

C     CREATED: SEPT 1/99 BY J.COTE, RPN

C     R     IS THE ARRAY CONTAINING OUTPUT DATA
C     A     IS THE ARRAY CONTAINING INPUT DATA
C     INC   IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C          (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP  IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     LOT   IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM FOURIER TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO FOURIER

C     DEFINITION OF TRANSFORM:
C     -------------------------

C     ISIGN=+1: R(I) = SUM(K=0,...,N-1)(A(K)*COS((I+1/2)*K*PI/N))

C     ISIGN=-1: R(K) = SUM(I=0,...,N-1)(A(I)*COS((I+1/2)*K*PI/N))
C                      * ((2-DELTA(K,0))/N)

C     NOTE FOR 'A' STORED AS A(N1,N2) THEN
C
C        FOR A TRANSFORM ALONG THE FIRST DIMENSION
C
C           INC   = 1
C           JUMP  = N1
C
C        FOR A TRANSFORM ALONG THE SECOND DIMENSION
C
C           INC   = N1
C           JUMP  = 1
C
C-----------------------------------------------------------------------
C
      SUBROUTINE QCFT8( R, A, INC, JUMP, LOT, ISIGN, N )

C
      IMPLICIT NONE
C
      INTEGER INC, JUMP, LOT, ISIGN, N
      REAL*8  R(*), A(*)
C
      INTEGER I, J, K
C
      REAL*8 ZERO, HALF, ONE, TWO
      PARAMETER( ZERO = 0.0 )
      PARAMETER( HALF = 0.5 )
      PARAMETER( ONE  = 1.0 )
      PARAMETER( TWO  = 2.0 )

      REAL*8 DEL, ANGLE, CC, FAC
      INTEGER IJ
      IJ(I,J) = 1 + ( J - 1 ) * JUMP + I * INC

      DEL   = ACOS( - ONE )/N
C
C     THE TRANSFORM  IS ALONG THE I-TH DIRECTION
C     AND THE INDEXING IS 1+(J-1)*JUMP+ I*INC (I=0,N-1,J=1,LOT)
C                      OR 1+(J-1)*JUMP+ K*INC (K=0,N-1,J=1,LOT)
C     IN GRIDPOINT AND FOURIER SPACE RESPECTIVELY
C
      IF ( ISIGN .EQ. -1 ) THEN
C
C     TRANSFORM FROM GRIDPOINT TO FOURIER
C
         FAC = ONE/N

         DO K=0,N-1

            DO J=1,LOT
               R( IJ(K,J) ) = ZERO
            ENDDO

            ANGLE =  K * DEL

            DO I=0,N-1

               CC =  COS( ( I + HALF ) * ANGLE )

               DO J=1,LOT
                  R( IJ(K,J) ) = R( IJ(K,J) ) + CC * A( IJ(I,J) )
               ENDDO

            ENDDO

            DO J=1,LOT
               R( IJ(K,J) ) = FAC * R( IJ(K,J) )
            ENDDO

            FAC = TWO/N

         ENDDO
         
      ELSEIF( ISIGN .EQ. +1 ) THEN
C
C     TRANSFORM FROM FOURIER TO GRIDPOINT
C
         DO I=0,N-1

            DO J=1,LOT
               R( IJ(I,J) ) = A( IJ(0,J) )
            ENDDO
 
            ANGLE = ( I + HALF ) * DEL

            DO K=1,N-1

               CC =  COS( K * ANGLE )

               DO J=1,LOT
                  R( IJ(I,J) ) = R( IJ(I,J) ) + CC * A( IJ(K,J) )
               ENDDO

            ENDDO

         ENDDO

      ENDIF

      RETURN
      END
C     SUBROUTINE 'RFT8' - MULTIPLE SLOW REAL TRANSFORM

C     REAL TRANSFORM OF LENGTH N

C     CREATED: MAR 12/99 BY J.COTE, RPN

C     R     IS THE ARRAY CONTAINING OUTPUT DATA
C     A     IS THE ARRAY CONTAINING INPUT DATA
C     INC   IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C          (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP  IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     LOT   IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM FOURIER TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO FOURIER

C     ORDERING OF FOURIER COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=0, POSSIBLY B(N/2)=0 ; (2*[N/2]+2) LOCATIONS REQUIRED

C     ORDERING OF GRIDPOINT DATA:
C         X(0),X(1),X(2),...,X(N-1), ? , ? ; (2*[N/2]+2) LOCATIONS REQUIRED

C     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN

C     DEFINITION OF TRANSFORMS:
C     -------------------------

C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C               USING THE SYMMETRY THE SUM RESTRICTED TO K = 0, [N/2]

C     ISIGN=-1: A(K)= (1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C               K = 0, [N/2]

C
C     NOTE FOR 'A' STORED AS A(N1,N2) THEN
C
C        FOR A TRANSFORM ALONG THE FIRST DIMENSION
C
C           INC   = 1
C           JUMP  = N1
C
C        FOR A TRANSFORM ALONG THE SECOND DIMENSION
C
C           INC   = N1
C           JUMP  = 1
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RFT8( R, A, INC, JUMP, LOT, ISIGN, N )

C
      IMPLICIT NONE
C
      INTEGER INC, JUMP, LOT, ISIGN, N
      REAL*8  R(*), A(*)
C
      INTEGER I, J, K, KK, K1
C
      REAL*8 ZERO, HALF, ONE, TWO
      PARAMETER( ZERO = 0.0 )
      PARAMETER( HALF = 0.5 )
      PARAMETER( ONE  = 1.0 )
      PARAMETER( TWO  = 2.0 )

      REAL*8 ONEON, DEL, ANGLE, CC, SS, FAZ
      INTEGER IJ
      IJ(I,J) = 1 + ( J - 1 ) * JUMP + I * INC

      ONEON = ONE/N
      DEL   = ACOS( - ONE ) * TWO/N
C
C     THE TRANSFORM  IS ALONG THE I-TH DIRECTION
C     AND THE INDEXING IS 1+(J-1)*JUMP+I*INC (I=0,2*[N/2]+1,J=1,LOT)
C     IN THE INPUT & OUTPUT FIELDS
C
      IF ( ISIGN .EQ. -1 ) THEN
C
C     TRANSFORM FROM GRIDPOINT TO FOURIER
C
      DO K=0,2*(N/2)+1
         DO J=1,LOT
            R( IJ(K,J) ) = ZERO
         ENDDO
      ENDDO

      DO K=0,N/2

         KK = 2 * K
         K1 = KK + 1

         DO I=0,N-1

            ANGLE = K * I * DEL
            CC =   ONEON * COS( ANGLE )
            SS = - ONEON * SIN( ANGLE )

            DO J=1,LOT
               R( IJ(KK,J) ) = R( IJ(KK,J) ) + CC * A( IJ(I,J) )
               R( IJ(K1,J) ) = R( IJ(K1,J) ) + SS * A( IJ(I,J) )
            ENDDO

         ENDDO

      ENDDO

      ELSEIF ( ISIGN .EQ. +1 ) THEN
C
C    TRANSFORM FROM FOURIER TO GRIDPOINT
C
      FAZ = - ONE * MOD( N + 1, 2 )
      DO I=0,N-1

         FAZ = - FAZ
         DO J=1,LOT
            R( IJ(I,J) ) = A( IJ(0,J) ) + FAZ * A( IJ(2*(N/2),J) )
         ENDDO

         DO K=1,(N-1)/2

            KK = 2 * K
            K1 = KK + 1
            ANGLE = K * I * DEL
            CC = TWO * COS( ANGLE )
            SS = TWO * SIN( ANGLE )

            DO J=1,LOT
               R( IJ(I,J) ) = R( IJ(I,J) ) + CC * A( IJ(KK,J) )
     %                                     - SS * A( IJ(K1,J) )
            ENDDO

         ENDDO

      ENDDO

      ENDIF
         
      RETURN
      END

