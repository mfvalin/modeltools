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

