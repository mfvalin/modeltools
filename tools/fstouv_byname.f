*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
      integer function fstouv_byname(names,iun,mode)
      implicit none
      character (len=*), intent(IN) :: names
      integer, intent(INOUT) :: iun
      character (len=*) :: mode

      integer fnom, fstouv, name_is_a_file, name_is_a_dir
      integer f_set_glob, f_next_glob
      external fnom, fstouv, name_is_a_file, name_is_a_dir
      external f_set_glob, f_next_glob

      integer, parameter :: MAXUNITS=1024
      integer, dimension(MAXUNITS) :: units
      integer :: nunits, status, i, list, istart, iend, ilen
      character (len=1024) :: globname1, globname2
      character (len=1) :: NULL

      status = -1
      nunits = 0
      units = 0
      units(1) = iun
      NULL = achar(0)

      if(0 == name_is_a_file(trim(names)//NULL)) then  ! names is a regular file
        print *,'FILE: ',trim(names)
        status = fnom(iun,names,'STD+RND',0)               ! call fnom and fstouv
        if(status /= 0) goto 111
        status = fstouv(iun,mode)
        goto 111                                           ! and return
      endif

      globname1 = trim(names)
      if(0 == name_is_a_dir(trim(globname1)//NULL)) then   ! names is a directory
        print *,'DIR: ',trim(names)
        globname1 = trim(names)//'/*.fst'                   ! open names/*.fst if names/.dir not found
        if(0 == name_is_a_file(trim(names)//'/.dir'//NULL)) then
          print *,trim(names)//'/.dir'//' FOUND'
          list = get_free_fortran_unit(99)
          open(unit=list,file=trim(names)//'/.dir',form='FORMATTED')   ! open list of names
1         read(list,*,end=2)globname2                            ! for each name
          nunits = nunits + 1                                  ! bump file count
          status = fnom(units(nunits),globname2,'STD+RND',0)   ! assign name to file
          if(status /= 0) goto 222                             ! error in fnom
          status = fstouv(units(nunits),mode)                  ! open standard file
          if(status /= 0) goto 222                             ! error in fstouv
          goto 1
2         continue                                             ! end of .dir file
          close(list)                                            ! close it
          call fstlnk(units,nunits)                            ! "link" standard files
          goto 111
        endif
      endif

      status=f_set_glob(trim(globname1)//NULL)              ! prepare to resolve globbed names
      print *,'GLOB: ',trim(globname1),status
      if(status == 0) then                                      ! matches have been found
        do while(f_next_glob(globname2,1024) == 0)
          nunits = nunits + 1
          status = fnom(units(nunits),globname2,'STD+RND',0)   ! assign name to file
          if(status /= 0) goto 333                             ! error in fnom
          status = fstouv(units(nunits),mode)                  ! open standard file
          if(status /= 0) goto 333                             ! error in fstouv
        enddo
        call fstlnk(units,nunits)                            ! "link" standard files
        call f_free_glob()
        status = 0
        goto 111                                             ! and return
      endif

!     if we get here, we have a blank or semicolon (;) separated list of names
      istart = 1
      iend = 1
      ilen = len(names)
      do while(istart < ilen)
        do while(names(iend:iend) /= '+' .and. iend <= ilen)
          iend = iend + 1
        enddo
        globname2=names(istart:iend-1)
        if(0 == name_is_a_file(trim(globname2)//NULL)) then    ! file exists ?
          nunits = nunits + 1
          status = fnom(units(nunits),globname2,'STD+RND',0)   ! assign name to file
          if(status /= 0) goto 333                             ! error in fnom
          status = fstouv(units(nunits),mode)                  ! open standard file
          if(status /= 0) goto 333                             ! error in fstouv
        endif
        istart = iend + 1
      enddo
      call fstlnk(units,nunits)
      status = 0

111   fstouv_byname = status
      return

222   close(list)               ! error, close .dir file
333   do i = 1 , nunits-1       ! errors detected, close all the standard files that we opened
        call fstfrm(units(i))
      enddo
      call f_free_glob()        ! harmless if not needed
      goto 111

      CONTAINS

      integer function get_free_fortran_unit(from)
      implicit none
      integer, intent(IN) :: from

      character (len=32) :: access_mode
      integer :: i
      get_free_fortran_unit = -1
      do i = from, 1, -1  ! find an available unit number
        inquire(UNIT=i,ACCESS=access_mode)
        if(trim(access_mode) == 'UNDEFINED')then ! found
          get_free_fortran_unit = i
          exit
        endif
      enddo
      return
      end function get_free_fortran_unit

      end function fstouv_byname
