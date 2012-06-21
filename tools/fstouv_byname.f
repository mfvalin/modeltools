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
      module ouvfrm_by_name
      integer, save :: nlists=-1
      type link_element
        integer, dimension(:), pointer :: units
        integer :: nunits
        type(link_element), pointer :: next
      end type link_element
      type(link_element), target, save :: link_list
      integer, dimension(2048), save :: iun
      integer, save :: niun=0

      contains

      integer function get_fortran_unit(from)
      implicit none
      integer, intent(IN) :: from

      character (len=32) :: access_mode
      integer :: i
      get_fortran_unit = -1
      do i = from, 1, -1  ! find an available unit number
        inquire(UNIT=i,ACCESS=access_mode)
        if(trim(access_mode) == 'UNDEFINED')then ! found
          get_fortran_unit = i
          exit
        endif
      enddo
      return
      end function get_fortran_unit

      end module ouvfrm_by_name

!     close a list of RANDOM STANDARD file(s) and "unlink" them
!     undo the work of fstouv_by_name
      integer function fstfrm_by_name(item)
      use ouvfrm_by_name
      implicit none
      integer item

      integer i
      type(link_element), pointer :: element, previous
      
      fstfrm_by_name = -1
      if(item <= 0) return

      element => link_list%next
      previous => link_list
      
      do while( associated(element) )       ! go through list of linked standard files
         if(element%units(1) == item) then  ! units(1) is the unit number for file operations
           call fstunl(element%units,element%nunits)  ! unlink the files
           do i = 1 , element%nunits
             call fstfrm(element%units(i))            ! and close them
           enddo
           deallocate(element%units)
           fstfrm_by_name = 0
           return
         endif
         previous = element
         element = element%next
      enddo

      do i = 1 , niun  ! not a linked item, scan single file list
        if(item == iun(i))then
          call fstfrm(item)
          iun(i) = 0
          fstfrm_by_name = 0
          return
        endif
      enddo
      return  ! if we get here, item was not found, ERROR
      end function fstfrm_by_name

!     open one or more RANDOM STANDARD file(s) and "link" them
      integer function fstouv_by_name(name,options)
      use ouvfrm_by_name
      implicit none
      character (len=*), intent(IN) :: name
      integer options

      include 'fstouv_by_name.inc'
      integer, parameter :: maxunits=2048
      integer, dimension(maxunits) :: units
      integer :: nunits
      integer status, i, istart, iend, ilen
      character (len=1) :: NULL
      character (len=32) :: file_mode
      character (len=*), parameter :: defdir='.dir'  ! default name of file list in directory
      character (len=2048) :: filelist, dirname, globname2
      character (len=*), parameter :: defext='*.fst'   ! default extension to look for in directory
      character (len=32) :: extension
      type(link_element), pointer :: next_list

      integer fnom, fstouv, name_is_a_file, name_is_a_dir
      integer name_is_txt_file
      integer f_set_glob, f_next_glob
      external fnom, fstouv, name_is_a_file, name_is_a_dir
      external name_is_txt_file
      external f_set_glob, f_next_glob

      if(nlists == -1) then   ! initialize internal chain of linked files units to empty state
        nullify(link_list%next)
        nlists = 0
      endif

      fstouv_by_name = -1   !  precondition fail
      units=0
      NULL = achar(0)
      extension=''
      filelist=''
      dirname=''
      if(iand(options,FSTOUV_WRITE) == FSTOUV_WRITE) then
        file_mode = 'STD+RND'
      else
        file_mode = 'STD+RND+OLD+R/O'
      endif

      if(0 == name_is_a_file(trim(name)//NULL)) then
        if(0 == name_is_txt_file(trim(name)//NULL)) then   ! text file, this is the list of file names
          filelist=name
          dirname=name
          i=len(trim(name))
          do while( i>1 .and. dirname(i:i) /= '/' )  ! keep only directory name
            dirname(i:i)=' '
            i=i-1
          enddo
        else  ! not a text file, try to open as a standard file
          print *,'FILE: ',trim(name)
          status = fnom(units(1),name,trim(file_mode),0)               ! call fnom and fstouv
          if(status /= 0) goto 111
          status = fstouv(units(1),'RND')
          if(status == 0) fstouv_by_name = units(1)
          return
        endif

      elseif(0 == name_is_a_dir(trim(name)//NULL)) then
          dirname=name//'/'
          filelist=trim(dirname)//defdir
      endif

      if(dirname /= '')then  ! we have a directory name and a file list
      endif

!     if we get here, we have a plus sign (+) separated list of names
      istart = 0
      iend = 0
      ilen = len(trim(name))
      nunits = 0
      do while(istart < ilen)
        istart = iend + 1
        iend = iend + 1
        do while(name(iend:iend) /= '+' .and. iend <= ilen)
          iend = iend + 1
        enddo
        if(name(istart:iend) == '+') cycle ! null file name
        globname2=name(istart:iend-1)
        if(0 == name_is_a_file(trim(globname2)//NULL)) then    ! file exists ?
!         check that file is a standard file with wkoffit
!         if not a standard file and LAZY flag absent, ERROR
          nunits = nunits + 1
          if(nunits > maxunits) goto 333                       ! unit table full
          status = fnom(units(nunits),globname2,file_mode,0)   ! assign name to file
          if(status /= 0) goto 333                             ! error in fnom
          status = fstouv(units(nunits),'RND')                 ! open standard file
          if(status /= 0) goto 333                             ! error in fstouv
        else                                                   ! bad file name, error
          print *,"BAD FILE NAME: '"//trim(globname2)//"'"
          status = -1
          goto 333
        endif
      enddo
      call fstlnk(units,nunits)
      status = 0

111   fstouv_by_name = status
      return

222   continue ! close(filelist)               ! error, close .dir file
333   do i = 1 , nunits         ! errors detected, close all the standard files that we opened
        call fstfrm(units(i))
      enddo
      call f_free_glob()        ! harmless if not needed
      goto 111

      end function fstouv_by_name

!     open a group of RANDOM STANDARD files and "link" them
! IN   names     : 4 possibilities
!                  a file name
!                  a directory name
!                     if directory/"dir" is found, it is assumed to contain the names of the files to open
!                     otherwise, directory/*."ext" will be used as a glob pattern
!                  a file pattern (shell style)
!                  a list of file names separated by a plus sign(+)
! IN    dir      : name of the file containing the file names (if blank, .dir)
! IN    ext      : extension of files ot open (if blank, *.fst)
! INOUT units    : units(1) in IN unless it is equal to 0
!                  units(2:nunits) is OUT
!                  will contain the list of units the files are connected to
! IN    maxunits : size of array units
! OUT   nunits   : number of files openned and linked
!
!     in case of error, the function will return a non zero error code, no file
!     will be left open (the already open ones will be closed), and no "link" operation will be performed
!
!     in case of success, the function will return zero
!     and units(1) contains the unit number to be used to address the "linked" file set
!     units(1) must be set by the caller, 
!     units(2:maxunits) will be initialized to zero before opening files
! notes:
!     this function may not be called with standard sequential files
!
!     it is the caller's responsibility to unlink and close the files when done with them
!     call fstunl(units,nunits)
!     do i = 1 , nunits
!       call fstfrm(units(i))
!     enddo
!
      integer function fstouv_byname
     %        (names,units,maxunits,nunits,dir,ext)
      use ouvfrm_by_name
      implicit none
      character (len=*), intent(IN) :: names, dir, ext
      integer, intent(IN) :: maxunits
      integer, intent(OUT) :: nunits
      integer, dimension(maxunits), intent(INOUT) :: units

      integer fnom, fstouv, name_is_a_file, name_is_a_dir
      integer f_set_glob, f_next_glob
      external fnom, fstouv, name_is_a_file, name_is_a_dir
      external f_set_glob, f_next_glob

      integer :: status, i, list, istart, iend, ilen
      character (len=1024) :: globname1, globname2, listname
      character (len=1) :: NULL
      character (len=*), parameter :: defdir='.dir'  ! default name of file list in directory
      character (len=32) :: filelist
      character (len=*), parameter :: defext='*.fst'   ! default extension to look for in directory
      character (len=32) :: extension
      character (len=*), parameter :: mode='RND'

      status = -1
      nunits = 0
      units(2:maxunits) = 0
      NULL = achar(0)
      extension = ext
      if(extension == " ") extension = defext
      filelist = dir
      if(filelist == " ") filelist = defdir

      if(0 == name_is_a_file(trim(names)//NULL)) then  ! names is a regular file
        print *,'FILE: ',trim(names)
        status = fnom(units(1),names,'STD+RND+OLD+R/O',0)               ! call fnom and fstouv
        if(status /= 0) goto 111
        status = fstouv(units(1),mode)
        nunits = 1
        goto 111                                           ! and return
      endif
!     if we get here, we have either a glob pattern, directory, or list of names
      globname1 = trim(names)
      if(0 == name_is_a_dir(trim(names)//NULL)) then   ! names is a directory
        print *,'DIR: ',trim(names)
        globname1 = trim(names)//'/'//trim(extension)        ! try names/*."ext" if names/".dir" not found
        listname  = trim(names)//'/'//trim(filelist)           ! name of file list to look for
        if(0 == name_is_a_file(trim(listname)//NULL)) then
          print *,'FOUND '//trim(listname)
          list = get_fortran_unit(99)
          open(unit=list,file=trim(listname)//NULL,form='FORMATTED')   ! open list of names
1         read(list,*,end=2)globname2                          ! for each name in .dir file
          globname2=trim(names)//'/'//globname2
          nunits = nunits + 1                                  ! bump file count
          if(nunits > maxunits) goto 222                       ! unit table full
          status = fnom(units(nunits),globname2,'STD+RND+OLD+R/O',0)   ! assign name to file
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
!     if we get here, we have either a glob pattern, directory/*.fst, or list of names
      status=f_set_glob(trim(globname1)//NULL)              ! prepare to resolve globbed names
      print *,'GLOB: ',trim(globname1),status
      if(status == 0) then                                      ! matches have been found
        do while(f_next_glob(globname2,1024) == 0)
          nunits = nunits + 1
          if(nunits > maxunits) goto 333                       ! unit table full
          status = fnom(units(nunits),globname2,'STD+RND+OLD+R/O',0)   ! assign name to file
          if(status /= 0) goto 333                             ! error in fnom
          status = fstouv(units(nunits),mode)                  ! open standard file
          if(status /= 0) goto 333                             ! error in fstouv
        enddo
        call fstlnk(units,nunits)                            ! "link" standard files
        call f_free_glob()
        status = 0
        goto 111                                             ! and return
      endif
!     if we get here, we have a plus sign (+) separated list of names
      istart = 0
      iend = 0
      ilen = len(trim(names))
      do while(istart < ilen)
        istart = iend + 1
        iend = iend + 1
        do while(names(iend:iend) /= '+' .and. iend <= ilen)
          iend = iend + 1
        enddo
        if(names(istart:iend) == '+') cycle ! null file name
        globname2=names(istart:iend-1)
        if(0 == name_is_a_file(trim(globname2)//NULL)) then    ! file exists ?
          nunits = nunits + 1
          if(nunits > maxunits) goto 333                       ! unit table full
          status = fnom(units(nunits),globname2,'STD+RND+OLD+R/O',0)   ! assign name to file
          if(status /= 0) goto 333                             ! error in fnom
          status = fstouv(units(nunits),mode)                  ! open standard file
          if(status /= 0) goto 333                             ! error in fstouv
        else                                                   ! bad file name, error
          print *,"BAD FILE NAME: '"//trim(globname2)//"'"
          status = -1
          goto 333
        endif
      enddo
      call fstlnk(units,nunits)
      status = 0

111   fstouv_byname = status
      return

222   close(list)               ! error, close .dir file
333   do i = 1 , nunits         ! errors detected, close all the standard files that we opened
        call fstfrm(units(i))
      enddo
      call f_free_glob()        ! harmless if not needed
      goto 111

      end function fstouv_byname
