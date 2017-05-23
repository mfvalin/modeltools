!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
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
  program self_test
    use ISO_C_BINDING
    implicit none
    integer :: status, iun1, iun2, iun3, errors, nblks, fsize, nw
    integer :: i
    integer, external :: fnom, fclos, existe, wasize, numblks, waread2, wawrit2
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
    status = fnom(iun1,'./fortran_wa','RND',0)
    print *,'iun1, fnom status(wa) =',iun1,status
    call waopen(iun1)
    buf = 1
    call wawrit(iun1,buf,1,4)
    buf = 2
    nw = wawrit2(iun1,buf,6,4)
    if(nw == 4) print *,'wawrit2 OK at 4'
    call waclos(iun1)
    call waopen(iun1)
    nblks = numblks(iun1)
    fsize = wasize(iun1)
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
  function fnom(iun, path, options, lrec) result(status)  ! moved from C to Fortran
    use ISO_C_BINDING
    implicit none
    integer, intent(IN) :: iun, lrec
    character(len=*), intent(IN) :: path, options
    integer :: status
!
    interface   ! int c_fnom(int *iun,char *nom,char *type,int lrec)
      function true_fnom(iun, path, options, lrec) result(status) bind(C,name='c_fnom')
      import :: C_INT, C_CHAR
      character(C_CHAR), dimension(*), intent(IN) :: path, options
      integer(C_INT) :: iun
      integer(C_INT), intent(IN), value :: lrec
      integer(C_INT) :: status
      end function true_fnom
    end interface
!
    character (C_CHAR), dimension(len(path)+1) :: my_path
    character (C_CHAR), dimension(len(options)+1) :: my_options
!
    my_path = transfer(trim(path)//achar(0),my_path)
    my_options = transfer(trim(options)//achar(0),my_options)
!
    status = true_fnom(iun, my_path, my_options, lrec)
!
  end function
!
! TODO:L ajouter flags OLD R/O append scratch (passer les flags recus ?)
!
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
      OPEN(iun,FILE=name(1:lng),ACCESS=acc,FORM='UNFORMATTED',STATUS=fstat,ACTION=action,RECL=lrec*lmult,ERR=77)
    else                   ! no recl if not random 77
      OPEN(iun,FILE=name(1:lng),ACCESS=acc,FORM=form,STATUS=fstat,POSITION=position,ACTION=action,ERR=77)
    endif
!
    return
77  continue
    status = -1
    return
  end
!
  function fclos(iun) result(status)
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
  FUNCTION ftnclos(iun) result(status) BIND(C,name='ftn_clos')
    use ISO_C_BINDING
    implicit none
    integer(C_INT), value :: iun
    integer(C_INT) :: status
    
    status = 0
    CLOSE(iun)
    return
  end
!
  INTEGER FUNCTION LONGUEUR(NOM)
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
function getfdsc(iun) result(i)
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
subroutine sqopen(iun)
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
subroutine sqclos(iun)
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
subroutine sqrew(iun)
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
subroutine sqeoi(iun)
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
subroutine sqgetw(iun,buf,nmots)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine csqgetw(iun,buf,nmots) bind(C,name='c_sqgetw')
      import
      integer(C_INT), intent(IN), value :: iun,nmots
      integer(C_INT), intent(OUT) :: buf
    end subroutine csqgetw
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(OUT) :: buf

  call csqgetw(iun,buf,nmots)
end subroutine sqgetw
!
subroutine sqputw(iun,buf,nmots)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine csqputw(iun,buf,nmots) bind(C,name='c_sqputw')
      import
      integer(C_INT), intent(IN), value :: iun,nmots
      integer(C_INT), intent(IN) :: buf
    end subroutine csqputw
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(IN) :: buf

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

! subroutine waread64(iun,buf,adr,nmots,partition)
!   use ISO_C_BINDING
!   implicit none
!   interface
!     subroutine cwaread64(iun,buf,adr,nmots,partition) bind(C,name='c_waread64')
!       import
!       integer(C_INT), intent(IN), value :: iun, nmots, partition
!       integer(C_LONG_LONG), intent(IN), value :: adr
!       integer(C_INT), intent(OUT), dimension(nmots) :: buf
!     end subroutine cwaread64
!   end interface
!   integer, intent(IN) :: iun, nmots, partition
!   integer*8, intent(IN) :: adr
!   integer, intent(OUT), dimension(nmots) :: buf
! 
!   call cwaread64(iun,buf,adr,nmots,partition)
! end subroutine waread64

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

function numblks(iun) result(i)
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

function wasize(iun) result(i)
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

function existe(name) result(status)
  use ISO_C_BINDING
  implicit none
  interface
    function cexiste(name) result(status) bind(C,name='c_existe')
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: name
      integer(C_INT) :: status
    end function cexiste
  end interface
  character(len=*), intent(IN) :: name
  integer :: status
  character(C_CHAR), dimension(len(trim(name))+1), target :: name1

  name1 = transfer(trim(name)//achar(0),name1)
  status = cexiste(name1)
  return
end function existe

! subroutine wawrit64(iun,buf,adr,nmots,partition)
!   use ISO_C_BINDING
!   implicit none
!   interface
!     subroutine cwawrit64(iun,buf,adr,nmots,partition) bind(C,name='c_wawrit64')
!       import
!       integer(C_INT), intent(IN), value :: iun, nmots, partition
!       integer(C_LONG_LONG), intent(IN), value :: adr
!       integer(C_INT), intent(IN), dimension(nmots) :: buf
!     end subroutine cwawrit64
!   end interface
!   integer, intent(IN) :: iun, nmots, partition
!   integer*8, intent(IN) :: adr
!   integer, intent(IN), dimension(nmots) :: buf
! 
!   call cwawrit64(iun,buf,adr,nmots,partition)
! end subroutine wawrit64

