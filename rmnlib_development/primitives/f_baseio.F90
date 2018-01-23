!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2018  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!
#if defined(SELF_TEST)
  program self_test   ! bare bones test
    use ISO_C_BINDING
    implicit none
    integer :: status, iun1, iun2, iun3, errors, nw
    integer*8 :: nblks, fsize
    integer :: i
    integer*8 :: wa8
    integer, external :: fnom, fclos, existe, wasize, numblks, waread2, wawrit2
    integer*8, external :: numblks64, wasize64
    integer, dimension(4) :: buf
    integer, dimension(512) :: buf512

    iun1 = 0
    status = fnom(iun1,'./fortran_formatted','SEQ+FMT+FTN+APPEND',0)
    print *,'iun1, fnom status =',iun1,status
    write(iun1,*)'ceci est un test'

    iun2 = 0
    status = fnom(iun2,'./fortran_unformatted','SEQ+UNF+FTN+APPEND',0)
    print *,'iun2, fnom status =',iun2,status
    write(iun2)(i,i=0,3)

    iun3 = 0
    status = fnom(iun3,'fortran_d77','UNF+FTN+D77+SCRATCH',5)
    print *,'iun3, fnom status =',iun3,status
    buf = 1
    write(iun3,rec=1)buf
    buf = 2
    write(iun3,rec=2)buf
    call system("ls -l $TMPDIR")
    
    status = fclos(iun1)
    print *,'iun, fclos(iun1), status =',iun1,status
    status = fclos(iun2)
    print *,'iun, fclos(iun2), status =',iun2,status
    status = fclos(iun3)
    print *,'iun, fclos(iun3), status =',iun3,status

    iun1 = 0
    status = fnom(iun1,'./fortran_wa','RND+SPARSE',0)
    print *,'iun1, fnom status(wa) =',iun1,status
    call waopen(iun1)
    buf = 1
    call wawrit(iun1,buf,1,4)
    wa8 = 1024*1024
    wa8 = 1024*wa8*4
    call wawrit64(iun1,buf,wa8,4,0)  ! sparse write at very high address
    buf = 2
    nw = wawrit2(iun1,buf,6,4)
    if(nw == 4) print *,'wawrit2 OK at 4'
    call waclos(iun1)
    call waopen(iun1)
    nblks = numblks64(iun1)
    fsize = wasize64(iun1)
    print *,'iun1(fortran_wa), numblks, fsize =',iun1,  nblks, fsize
    call waclos(iun1)

    iun2 = 0
    status = fnom(iun2,'./fortran_da','RND',0)
    print *,'iun2, fnom status(da) =',iun2,status
    call openda(iun2)
    buf512 = 1
    call writda(iun2,buf512,1,1)
    buf512 = 2
    call writda(iun2,buf512,1,2)
    call closda(iun2)

    iun3 = 0
    status = fnom(iun3,'./fortran_da','RND',0)
    print *,'iun3, fnom status(da) =',iun3,status
    call waopen(iun3)
    errors = 0
    buf512 = 0
    call waread(iun3,buf512,1,512)
    do i=1,512
      if(buf512(i) .ne. 1) errors = errors + 1
    enddo
    buf512 = 0
    nw = waread2(iun3,buf512,513,512)
    if(nw == 512) print *,'waread2 OK at 512'
    do i=1,512
      if(buf512(i) .ne. 2) errors = errors + 1
    enddo
    print *,'iun3, errors(waread) =',iun3, errors
    call waclos(iun3)
    status = fclos(iun3)
    print *,'iun3, fclos(iun3), status =',iun3,status

    iun3 = 0
    status = fnom(iun3,'./fortran_da','RND+OLD+R/O',0)
    print *,'iun3, fnom status(da) =',iun3,status
    call openda(iun3)
    errors = 0
    buf512 = 0
    call readda(iun3,buf512,1,1)
    do i=1,512
      if(buf512(i) .ne. 1) errors = errors + 1
    enddo
    buf512 = 0
    call readda(iun3,buf512,1,2)
    do i=1,512
      if(buf512(i) .ne. 2) errors = errors + 1
    enddo
    nblks = numblks(iun3)
    fsize = wasize(iun3)
    print *,'iun3, errors(readda), numblks, fsize =',iun3, errors, nblks, fsize
    call closda(iun3)

    iun3 = 0
    status = fnom(iun3,'./fortran_da','RND',0)
    print *,'iun3, fnom status(da) =',iun3,status
    call waopen(iun3)
    errors = 0
    buf512 = 0
    call waread(iun3,buf512,1,512)
    do i=1,512
      if(buf512(i) .ne. 1) errors = errors + 1
    enddo
    buf512 = 0
    call waread(iun3,buf512,513,512)
    do i=1,512
      if(buf512(i) .ne. 2) errors = errors + 1
    enddo
    print *,'iun3, errors =',iun3, errors
    call waclos(iun3)

    status = fclos(iun1)
    print *,'iun1, fclos(iun1), status =',iun1,status
    status = fclos(iun2)
    print *,'iun2, fclos(iun2), status =',iun2,status
    status = fclos(iun3)
    print *,'iun3, fclos(iun3), status =',iun3,status

    status = existe('./fortran_formatted')
    print *,'existe(./fortran_formatted) =',status
    status = existe('./fortran_unformatted')
    print *,'existe(./fortran_unformatted) =',status
    status = existe('./fortran_da')
    print *,'existe(./fortran_da) =',status
    status = existe('./fortran_wa')
    print *,'existe(./fortran_wa) =',status
    status = existe('./iepala')
    print *,'existe(./iepala) =',status

  end program
#endif
!
! iun        : if zero upon entry, fnom will find an appropriate unit number
! path       : file name (character string)
!              some_name@file_path refers to file some_name inside CMCARC archive file_path
! options    : list of + separated options (upper case or lower case)
!              STD         RPN "standard" file (implies WA+RND)
!              FTN         Fortran file (UNF, D77 may be used as sub attributes)
!              D77         Fortran direct access file (lrec must be non zero)
!              UNF         Fortran sequential unformatted file (default is formatted)
!              RND         random access file (normally used with STD)
!              WA          Word Addressable file (Big Endian) (implies RND)
!              STREAM      stream file (non Fortran, no record markers, Big Endian
!              BURP        Meteorological reports file
!              OLD         file must exist (applies to all files)
!              R/O         file is Read Only (default is Read/Write) (applies to all files) (implies OLD)
!              R/W         file is Read Write (default) (applies to all files)
!              SCRATCH     File will be removed when closed (applies to all files)
!              SPARSE      Unix WA sparse file (may be written into far beyond end of file)
!              PAGED       special type of WA file (not implemented yet)
!              REMOTE      file is on another system and accessed with a ssh (applies to all files)
!              ex.   STD+RND+OLD+R/W open existing random standard file for reading and writing
!                    FTN+UNF         open Fortran sequential file for reading and writing , create it if it does not exist
! lrec       : record length in 4 byte integers for Fortran D77 file records (should be zero otherwise)
!
  function fnom(iun, path, options, lrec) result(status)  ! moved from C to Fortran
    use ISO_C_BINDING
    implicit none
    integer, intent(INOUT) :: iun          ! iun == 0 : return usable iun upon return
    integer, intent(IN) :: lrec
    character(len=*), intent(IN) :: path, options
    integer :: status
!
    interface   ! int c_fnom(int *iun,char *nom,char *type,int lrec)
      function true_fnom(iun, path, options, lrec) result(status) bind(C,name='c_fnom')
      import :: C_INT, C_CHAR
      character(C_CHAR), dimension(*), intent(IN) :: path, options
      integer(C_INT), intent(INOUT) :: iun
      integer(C_INT), intent(IN), value :: lrec
      integer(C_INT) :: status
      end function true_fnom
    end interface
!
    character (C_CHAR), dimension(len(path)+1) :: my_path         ! add 1 position for the terminating null character
    character (C_CHAR), dimension(len(options)+1) :: my_options
!
    my_path = transfer(trim(path)//achar(0),my_path)              ! make null terminated C string from Fortran strings
    my_options = transfer(trim(options)//achar(0),my_options)
!
    status = true_fnom(iun, my_path, my_options, lrec)     ! call the C code after Fortran style to C style string conversion
!
  end function

  function get_ftn_free_unit_number() result(iun) bind(C,name='FtnFreeUnitNumber')
    implicit none
    integer :: iun
    integer, external :: get_free_unit_number

    iun = get_free_unit_number()
  end function get_ftn_free_unit_number

  function get_free_unit_number() result(iun)
    implicit none
    integer :: iun
    integer :: i
    character (len=16) :: access_mode
      iun = -1
      do i = 99, 1, -1  ! find an available Fortran unit number
	inquire(UNIT=i,ACCESS=access_mode)
	if(trim(access_mode) == 'UNDEFINED')then ! found
	  iun = i
	  exit
	endif
      enddo
    return
  end function get_free_unit_number

! this is called by the C code to set some Fortran I/O setup options (NOT USER CALLABLE)
! iun      : Fortran unit number
! path     : file name
! options  : SCRATCH APPEND OLD R/O
! lrec     : see fnom
! rndflag  : 0 sequential file, 1 random file
! unfflag  : 0 formatted file, 1 unformattted file
! lmult    : multiplier for lrec (compiler dependent)
! lng_in   : length of path
  FUNCTION qqqf7op_from_c(iun,path,options,lrec,rndflag,unfflag,lmult,lng_in) result(status) BIND(C,name='qqq_f90_options')
    use ISO_C_BINDING
    implicit none
    integer(C_INT), intent(IN), value :: iun, lrec, lmult, rndflag ,unfflag, lng_in
    character (C_CHAR), dimension(*), intent(IN) :: path, options
    integer(C_INT) :: status
  
    character(len=4096) :: name
    character(len=16) :: acc, form, action, position, fstat
    character(len=128) :: optns
    integer :: i, lng
!
    status=0
    lng = lng_in
    name(1:lng) = transfer(path(1:lng),name)    ! transfer filename to Fortran string
    i = 1
    optns = ' '
    do while(options(i) .ne. achar(0) .and. i < 128)
      optns(i:i) = options(i)
      i = i + 1
    enddo
#if defined(SELF_TEST)
print *,'qqqf7op_from_c, iun,path,lrec,rndflag,unfflag,lmult,lng',iun,"'"//name(1:lng)//"'",lrec,rndflag,unfflag,lmult,lng
print *,"      options ='"//trim(optns)//"'"
#endif
!
    if ((name(1:lng).EQ.'input') .OR. (name(1:lng).EQ.'$input')            &
     &  .OR. (name(1:lng).EQ.'output') .OR. (name(1:lng).EQ.'$output')     &
     &  .OR. (name(1:lng).EQ.'$in') .OR. (name(1:lng).EQ.'$out'))          &
     &  then
!            print *,' STDIN or STDOUT'
      return
    endif
!
    form     = 'FORMATTED'     ! default values
    acc      = 'SEQUENTIAL'
    position = 'REWIND'
    action   = 'READWRITE'
    fstat    = 'UNKNOWN'
!
    if (rndflag == 1) acc  = 'DIRECT'
    if (unfflag == 1) form = 'UNFORMATTED'
    if( index(optns,'SCRATCH',.false.) > 0 .or. index(optns,'scratch',.false.) > 0 ) then
      fstat    = 'SCRATCH'
      name     = 'NoNe'      ! name is irrelevant, make sure to give an acceptable one
      lng = 4
    endif
    if( index(optns,'APPEND',.false.) > 0 .or. index(optns,'append',.false.) > 0 ) then
      position = 'APPEND'
    endif
    if( index(optns,'OLD',.false.) > 0 .or. index(optns,'old',.false.) > 0 ) then
      fstat    = 'OLD'
    endif
    if( index(optns,'R/O',.false.) > 0 .or. index(optns,'r/o',.false.) > 0 ) then
      fstat    = 'OLD'       ! better exist if file is to be r/o
      action   = 'READ'
    endif
!
#if defined(SELF_TEST)
print *,trim(acc)//'+'//trim(form)//'+'//trim(fstat)//'+'//trim(position)//'+'//trim(action)
print *,'recl =',lrec*lmult
#endif
    if(rndflag == 1) then  ! no position nor form if random 77
      if(trim(fstat) == 'SCRATCH') then
        OPEN(iun,ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='SCRATCH',ACTION=action,RECL=lrec*lmult,ERR=77)
      else
        OPEN(iun,FILE=name(1:lng),ACCESS='DIRECT',FORM='UNFORMATTED',STATUS=fstat,ACTION=action,RECL=lrec*lmult,ERR=77)
      endif
    else                   ! no recl if not random 77
      if(trim(fstat) == 'SCRATCH') then
        OPEN(iun,ACCESS='SEQUENTIAL',FORM=form,STATUS='SCRATCH',POSITION=position,ACTION=action,ERR=77)
      else
        OPEN(iun,FILE=name(1:lng),ACCESS='SEQUENTIAL',FORM=form,STATUS=fstat,POSITION=position,ACTION=action,ERR=77)
      endif
    endif
!
    return
77  continue
print *,'error in qqqf7op_from_c'
    status = -1
    return
  end
!
  function fclos(iun) result(status)    ! "close" file opened with fnom, makes iun usable again
    use ISO_C_BINDING
    implicit none
    integer, intent(IN) :: iun
    integer :: status
    interface
      function fclose(iun) result(status) BIND(C,name='c_fclos')
        import:: C_INT
        integer(C_INT), intent(IN), value :: iun
        integer(C_INT) :: status
      end function fclose
    end interface
#if defined(SELF_TEST)
print *,'closing unit =',iun
#endif
    status = fclose(iun)
  end function
!
  FUNCTION ftnclos(iun) result(status) BIND(C,name='ftn_clos')   ! actual Fortran close (called from C)
    use ISO_C_BINDING
    implicit none
    integer(C_INT), value :: iun
    integer(C_INT) :: status
    
    status = 0
    CLOSE(iun)
    return
  end
!
  INTEGER FUNCTION LONGUEUR(NOM) ! one should use lentrim instead
    implicit none
    CHARACTER(len=*), intent(IN) :: NOM
!
    INTEGER :: LNG, I
!
    LNG = LEN(NOM)
    DO I = LNG,1,-1
      IF (NOM(I:I) .EQ. ' ') THEN
        LNG = LNG - 1
      ELSE
        exit
      ENDIF
    enddo
    LONGUEUR = LNG
    RETURN
  END
!
integer function fretour(iun)
! ARGUMENTS: in iun   unit number, ignored
! RETURNS: zero.
! Kept only for backward compatibility. NO-OP
  fretour = 0
  return
end
!
function getfdsc(iun) result(i)  ! Get file descriptor associated to unit iun
  use ISO_C_BINDING
  implicit none
  interface
    function cgetfdsc(iun) result(i) bind(C,name='c_getfdsc')
      import
      integer(C_INT), intent(IN), value :: iun
      integer(C_INT) :: i
    end function cgetfdsc
  end interface
  integer, intent(IN) :: iun
  integer :: i

  i = cgetfdsc(iun)
end function getfdsc
!
subroutine sqopen(iun)   ! Open a stream on unit iun
  use ISO_C_BINDING
  implicit none
  interface
    subroutine csqopen(iun) bind(C,name='c_sqopen')
      import
      integer(C_INT), intent(IN), value :: iun
    end subroutine csqopen
  end interface
  integer, intent(IN) :: iun

  call csqopen(iun)
end subroutine sqopen
!
subroutine sqclos(iun)   ! Close stream associated to unit iun
  use ISO_C_BINDING
  implicit none
  interface
    subroutine csqclos(iun) bind(C,name='c_sqclos')
      import
      integer(C_INT), intent(IN), value :: iun
    end subroutine csqclos
  end interface
  integer, intent(IN) :: iun

  call csqclos(iun)
end subroutine sqclos
!
subroutine sqrew(iun)   ! Rewind stream associated to unit iun
  use ISO_C_BINDING
  implicit none
  interface
    subroutine csqrew(iun) bind(C,name='c_sqrew')
      import
      integer(C_INT), intent(IN), value :: iun
    end subroutine csqrew
  end interface
  integer, intent(IN) :: iun

  call csqrew(iun)
end subroutine sqrew
!
subroutine sqeoi(iun)   ! Go to the end of stream associated to unit iun (append position)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine csqeoi(iun) bind(C,name='c_sqeoi')
      import
      integer(C_INT), intent(IN), value :: iun
    end subroutine csqeoi
  end interface
  integer, intent(IN) :: iun

  call csqeoi(iun)
end subroutine sqeoi
!
subroutine sqgetw(iun,buf,nmots)   ! get nmots words (4 bytes) from stream associated to unit iun
  use ISO_C_BINDING
  implicit none
  interface
    subroutine csqgetw(iun,buf,nmots) bind(C,name='c_sqgetw')
      import
      integer(C_INT), intent(IN), value :: iun,nmots
      integer(C_INT), intent(OUT), dimension(nmots) :: buf
    end subroutine csqgetw
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(OUT), dimension(nmots) :: buf

  call csqgetw(iun,buf,nmots)
end subroutine sqgetw
!
subroutine sqputw(iun,buf,nmots)   ! write nmots words (4 bytes) into stream associated to unit iun
  use ISO_C_BINDING
  implicit none
  interface
    subroutine csqputw(iun,buf,nmots) bind(C,name='c_sqputw')
      import
      integer(C_INT), intent(IN), value :: iun,nmots
      integer(C_INT), intent(IN), dimension(nmots) :: buf
    end subroutine csqputw
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(IN), dimension(nmots) :: buf

  call csqputw(iun,buf,nmots)
end subroutine sqputw
!
subroutine traceback_from_c() bind(C,name='f_tracebck') ! C code calls f_tracebck()
  call tracebck  ! unfortunately does nothing with most compilers
end subroutine traceback_from_c
!
! ====================================================
!     openda/closda readda/writda/checda
!     "asynchronous" random access by block routines
!     (the current implementation is SYNCHRONOUS)
! IUN(IN)     : fortran unit number
! BUF(IN/OUT) : array to write from or read into
! NS          : number of "sectors" (sector = 512 bytes)
! IS          : address of first "sector" for transfer
!               file starts at sector #1
! ====================================================
subroutine openda(iun)
  implicit none
  integer, intent(IN) :: iun

  call waopen(iun)
end subroutine openda

subroutine readda(iun,buf,ns,is)
  implicit none
  integer, intent(IN) :: iun, ns, is
  integer, intent(OUT), dimension(512,ns) :: buf
  call waread(iun,buf,(is-1)*512+1,ns*512)
end subroutine readda

subroutine writda(iun,buf,ns,is)
  implicit none
  integer, intent(IN) :: iun, ns, is
  integer, intent(IN), dimension(512,ns) :: buf

  call wawrit(iun,buf,(is-1)*512+1,ns*512)
end subroutine writda

subroutine checda(iun)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine cchecda(iun) bind(C,name='c_checda')
      import
      integer(C_INT), intent(IN), value :: iun
    end subroutine cchecda
  end interface
  integer, intent(IN) :: iun

  call cchecda(iun)
end subroutine checda

subroutine closda(iun)
  implicit none
  integer, intent(IN) :: iun

  call waclos(iun)
end subroutine closda
! ====================================================
!  waopen/waclos waread/waread64/wawrit/wawrit64
!   random access by word (4 bytes) routines
!   these routines take care of endian conversion
!   the file contents are always BIG-ENDIAN (4 bytes)
! IUN(IN)     : fortran unit number
! BUF(IN/OUT) : array to write from or read into
! NMOTS(IN)   : number of "words" to read (word = 4 bytes)
! ADR(IN)     : address of first word for transfer
!               file starts at word #1
! ====================================================
subroutine waopen(iun)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine cwaopen(iun) bind(C,name='c_waopen')
      import
      integer(C_INT), intent(IN), value :: iun
    end subroutine cwaopen
  end interface
  integer, intent(IN) :: iun

  call cwaopen(iun)
end subroutine waopen

subroutine waclos(iun)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine cwaclos(iun) bind(C,name='c_waclos')
      import
      integer(C_INT), intent(IN), value :: iun
    end subroutine cwaclos
  end interface
  integer, intent(IN) :: iun

  call cwaclos(iun)
end subroutine waclos

function waopen2(iun) result(status)
  use ISO_C_BINDING
  implicit none
  interface
    function cwaopen2(iun) result(status) bind(C,name='c_waopen2')
      import
      integer(C_INT), intent(IN), value :: iun
      integer(C_INT) :: status
    end function cwaopen2
  end interface
  integer, intent(IN) :: iun
  integer :: status

  status = cwaopen2(iun)
end function waopen2

function waclos2(iun) result(status)
  use ISO_C_BINDING
  implicit none
  interface
    function cwaclos2(iun) result(status) bind(C,name='c_waclos2')
      import
      integer(C_INT), intent(IN), value :: iun
      integer(C_INT) :: status
    end function cwaclos2
  end interface
  integer, intent(IN) :: iun
  integer :: status

  status = cwaclos2(iun)
end function waclos2

subroutine waread(iun,buf,adr,nmots)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine cwaread(iun,buf,adr,nmots) bind(C,name='c_waread')
      import :: C_INT
      integer(C_INT), intent(IN), value :: iun, nmots
      integer(C_INT), intent(IN), value :: adr
      integer(C_INT), intent(OUT), dimension(nmots) :: buf
    end subroutine cwaread
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(IN) :: adr
  integer, intent(OUT), dimension(nmots) :: buf
#if defined(SELF_TEST)
print *,'waread, adr, nmots',adr,nmots
#endif
  call cwaread(iun,buf,adr,nmots)
end subroutine waread

function waread2(iun,buf,adr,nmots) result(nread)
  use ISO_C_BINDING
  implicit none
  interface
    function cwaread2(iun,buf,adr,nmots) result(nread) bind(C,name='c_waread2')
      import :: C_INT
      integer(C_INT), intent(IN), value :: iun, nmots
      integer(C_INT), intent(IN), value :: adr
      integer(C_INT), intent(OUT), dimension(nmots) :: buf
      integer(C_INT) :: nread   ! number of "words" read
    end function cwaread2
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(IN) :: adr
  integer, intent(OUT), dimension(nmots) :: buf
  integer :: nread   ! number of "words" read
  nread = cwaread2(iun,buf,adr,nmots)
#if defined(SELF_TEST)
print *,'waread, adr, nmots, nread',adr,nmots,nread
#endif
end function waread2

subroutine waread64(iun,buf,adr,nmots,mode)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine cwaread64(iun,buf,adr,nmots,mode) bind(C,name='c_waread64')
      import
      integer(C_INT), intent(IN), value :: iun, nmots, mode
      integer(C_LONG_LONG), intent(IN), value :: adr
      integer(C_INT), intent(OUT), dimension(nmots) :: buf
    end subroutine cwaread64
  end interface
  integer, intent(IN) :: iun, nmots, mode
  integer*8, intent(IN) :: adr
  integer, intent(OUT), dimension(nmots) :: buf

  call cwaread64(iun,buf,adr,nmots,mode)
end subroutine waread64

subroutine wawrit(iun,buf,adr,nmots)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine cwawrit(iun,buf,adr,nmots) bind(C,name='c_wawrit')
      import :: C_INT
      integer(C_INT), intent(IN), value :: iun, nmots
      integer(C_INT), intent(IN), value :: adr
      integer(C_INT), intent(IN), dimension(nmots) :: buf
    end subroutine cwawrit
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(IN) :: adr
  integer, intent(IN), dimension(nmots) :: buf

#if defined(SELF_TEST)
print *,'wawrit, adr, nmots',adr,nmots
#endif
  call cwawrit(iun,buf,adr,nmots)
end subroutine wawrit

function wawrit2(iun,buf,adr,nmots) result(nwrt)
  use ISO_C_BINDING
  implicit none
  interface
    function cwawrit2(iun,buf,adr,nmots) result(nwrt) bind(C,name='c_wawrit2')
      import :: C_INT
      integer(C_INT), intent(IN), value :: iun, nmots
      integer(C_INT), intent(IN), value :: adr
      integer(C_INT), intent(IN), dimension(nmots) :: buf
      integer(C_INT) :: nwrt   ! number of "words" written
    end function cwawrit2
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(IN) :: adr
  integer, intent(IN), dimension(nmots) :: buf
  integer :: nwrt   ! number of "words" written

#if defined(SELF_TEST)
print *,'wawrit, adr, nmots',adr,nmots
#endif
  nwrt = cwawrit2(iun,buf,adr,nmots)
end function wawrit2

function numblks(iun) result(i)   ! get size of file associated with iun in units of 512 bytes
  use ISO_C_BINDING
  implicit none
  interface
    function cnumblks(iun) result(i) bind(C,name='c_numblks')
      import
      integer(C_INT), intent(IN), value :: iun
      integer(C_INT) :: i
    end function cnumblks
  end interface
  integer, intent(IN) :: iun
  integer :: i

  i = cnumblks(iun)
end function numblks

function numblks64(iun) result(i8)   ! get size of file associated with iun in units of 512 bytes
  use ISO_C_BINDING
  implicit none
  interface
    function cnumblks64(iun) result(i8) bind(C,name='c_numblks64')
      import
      integer(C_INT), intent(IN), value :: iun
      integer(C_LONG) :: i8
    end function cnumblks64
  end interface
  integer, intent(IN) :: iun
  integer*8 :: i8

  i8 = cnumblks64(iun)
end function numblks64

function wasize(iun) result(i)   ! get size of file associated with iun in words (4 bytes)
  use ISO_C_BINDING
  implicit none
  interface
    function cwasize(iun) result(i) bind(C,name='c_wasize')
      import
      integer(C_INT), intent(IN), value :: iun
      integer(C_INT) :: i
    end function cwasize
  end interface
  integer, intent(IN) :: iun
  integer :: i

  i = cwasize(iun)
end function wasize

function wasize64(iun) result(i8)   ! get size of file associated with iun in words (4 bytes)
  use ISO_C_BINDING
  implicit none
  interface
    function cwasize64(iun) result(i8) bind(C,name='c_wasize64')
      import
      integer(C_INT), intent(IN), value :: iun
      integer(C_LONG) :: i8
    end function cwasize64
  end interface
  integer, intent(IN) :: iun
  integer*8 :: i8

  i8 = cwasize64(iun)
end function wasize64

function existe(name) result(status)  ! return 1 if file 'name' exits, 0 otherwise
  use ISO_C_BINDING
  implicit none
  interface
    function cexiste(name) result(status) bind(C,name='c_existe')
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: name
      integer(C_INT) :: status
    end function cexiste
  end interface
  character(len=*), intent(IN) :: name    ! file name
  integer :: status
  character(C_CHAR), dimension(len(trim(name))+1), target :: name1

  name1 = transfer(trim(name)//achar(0),name1)
  status = cexiste(name1)
  return
end function existe

subroutine wawrit64(iun,buf,adr,nmots,mode)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine cwawrit64(iun,buf,adr,nmots,mode) bind(C,name='c_wawrit64')
      import
      integer(C_INT), intent(IN), value :: iun, nmots, mode
      integer(C_LONG_LONG), intent(IN), value :: adr
      integer(C_INT), intent(IN), dimension(nmots) :: buf
    end subroutine cwawrit64
  end interface
  integer, intent(IN) :: iun, nmots, mode
  integer*8, intent(IN) :: adr
  integer, intent(IN), dimension(nmots) :: buf

  call cwawrit64(iun,buf,adr,nmots,mode)
end subroutine wawrit64

