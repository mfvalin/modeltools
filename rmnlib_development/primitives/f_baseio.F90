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
  FUNCTION qqqf7op_from_c(iun,path,lrec,rndflag,unfflag,lmult,lng) result(status) BIND(C,name='qqq_f7op')
    use ISO_C_BINDING
    implicit none
    integer(C_INT), intent(IN), value :: iun, lrec, lmult, rndflag ,unfflag, lng
    character (C_CHAR), dimension(*), intent(IN) :: path
    integer(C_INT) :: status
  
    character(len=4096) :: name
    character(len=16) :: acc, form
!
    status=0
    name(1:lng) = transfer(path(1:lng),name)    ! transfer filename to Fortran string
!
    if ((name(1:lng).EQ.'input') .OR. (name(1:lng).EQ.'$input')            &
     &  .OR. (name(1:lng).EQ.'output') .OR. (name(1:lng).EQ.'$output')     &
     &  .OR. (name(1:lng).EQ.'$in') .OR. (name(1:lng).EQ.'$out'))          &
     &  then
!            print *,' STDIN or STDOUT'
      return
    endif
!
    acc = 'SEQUENTIAL'
    if (rndflag.eq.1) acc = 'DIRECT'
!
    form = 'FORMATTED'
    if (unfflag.eq.1) form = 'UNFORMATTED'
    OPEN(iun,FILE=name(1:lng),ACCESS=acc,FORM=form,RECL=lrec*lmult,ERR=77)
!
    return
77  continue
    status = -1
    return
  end
!
  FUNCTION ftnclos(iun) result(status) BIND(C,name='ftn_clos')
    use ISO_C_BINDING
    implicit none
    integer(C_INT) :: iun
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

subroutine waread(iun,buf,adr,nmots)
  use ISO_C_BINDING
  implicit none
  interface
    subroutine cwaread(iun,buf,adr,nmots) bind(C,name='c_waread')
      import
      integer(C_INT), intent(IN), value :: iun, nmots
      integer(C_INT), intent(IN), value :: adr
      integer(C_INT), intent(OUT), dimension(nmots) :: buf
    end subroutine cwaread
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(IN) :: adr
  integer, intent(OUT), dimension(nmots) :: buf
print *,'waread, adr, nmots',adr,nmots
  call cwaread(iun,buf,adr,nmots)
end subroutine waread

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
      import
      integer(C_INT), intent(IN), value :: iun, nmots
      integer(C_INT), intent(IN), value :: adr
      integer(C_INT), intent(IN), dimension(nmots) :: buf
    end subroutine cwawrit
  end interface
  integer, intent(IN) :: iun, nmots
  integer, intent(IN) :: adr
  integer, intent(IN), dimension(nmots) :: buf

  call cwawrit(iun,buf,adr,nmots)
end subroutine wawrit

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

