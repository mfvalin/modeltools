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
!==============================================================================
      module ouvfrm_by_name   ! only used by fstfrm_by_name, fstouv_by_name, fstouv_byname, dump_ouvfrm_by_name

      integer, save :: nlists=-1
      type link_element
        integer, dimension(:), pointer :: units  ! array containing pseudo unit numbers
        integer :: nunits                        ! number of valid units in above array
        type(link_element), pointer :: next      ! pointer to next list
      end type link_element
      type(link_element), pointer, save :: link_list => NULL()

      contains
!------------------------------------------------------------------------------
      integer function get_fortran_unit() ! get the number of an unused fortran unit
      implicit none

      character (len=32) :: access_mode
      integer :: i, ios
      logical :: inuse

      get_fortran_unit = -1  ! precondition to fail
      ios = 0
      inuse = .true.
      do i = 99, 1, -1  ! find an available unit number
!        inquire(UNIT=i,ACCESS=access_mode)
!        if(trim(access_mode) == 'UNDEFINED')then ! found
        inquire(unit=i,opened=inuse,iostat=ios)
        if(.not. inuse .and. ios == 0)then ! found unit not in use, return its number
          get_fortran_unit = i
          return
        endif
      enddo
      return
      end function get_fortran_unit

      end module ouvfrm_by_name

!==============================================================================
      subroutine dump_ouvfrm_by_name()
      use ouvfrm_by_name
      implicit none
      type(link_element), pointer :: element
      integer :: i

      print *,"=======  sets of units ======="
      i = 1
      element => link_list
      do while( associated(element) )
         print *,"  set no ",i
         if(associated(element%units)) then
           print *,element%units(1:element%nunits)
         else
           print *,"element%units not associated"
         endif
         print *,"==================================="
         element => element%next
         i = i + 1
      enddo
      return
      end subroutine dump_ouvfrm_by_name

!==============================================================================
!     close a set of RANDOM STANDARD file(s) and "unlink" them
!     undo the work of fstouv_by_name
      integer function fstfrm_by_name(item)
      use ouvfrm_by_name
      implicit none
      integer item

      integer i
      type(link_element), pointer :: element, previous

      fstfrm_by_name = -1
      if(item <= 0) return

      element => link_list
      previous => null()

      do while( associated(element) )       ! go through list of linked standard files
         if(element%units(1) == item) then  ! units(1) is the unit number for file operations
           if(element%nunits > 1) then
             call fstunl(element%units,element%nunits)  ! unlink the files
           endif
           do i = 1 , element%nunits
             call fstfrm(element%units(i))            ! and close them
           enddo

           deallocate(element%units)         ! deallocate list of units
           nullify(element%units)
           if(associated(previous)) then     ! not at head of linked list
             previous%next => element%next
           else 
             link_list => element%next
           endif
           deallocate(element)               ! deallocate structure
           nullify(element)

           fstfrm_by_name = 0
           exit
         endif
         previous => element
         element => element%next
      enddo
      return 
      end function fstfrm_by_name

!==============================================================================
!     open a set of one or more RANDOM STANDARD file(s) and "link" them
! IN  name     : 5 cases
!              1) path to a standard file
!              2) path to a file containing a list of standard file names
!              3) shell style file pattern ( e.g. /a/b/c/*.fst  ./*.abc  x/y/*.Fst   *.xyz)
!              4) + separated list of standard file paths (e.g. /a/b/c/aa.fst+f/g/dd.fst+./abc.def+qwe.rty)
!              5) path to a directory
! IN options   : (user needs to include file fstouv_by_name.inc) (options may be combined by adding flags)
!                0 = no option selected
!                +FSTOUV_WRITE   open with write permission (case 1 only, ignored otherwise)
!                +FSTOUV_LIST    case 2 only, list of file names, the names of the standard files
!                                to open will be obtained by prefixing the names found in the file with 
!                                the path to the directory containning the list of file names
!                                name MUST contain a path with at least one directory element
!                                e.g. ./mylist.txt , /a/b/c/other.txt ./x/y/bla.ttx 
!                                ( something like list.txt will NOT work )
!                +FSTOUV_LAZY    if a pattern is supplied, just ignore files that are not standard files
!                                this flag is ignored unless a pattern has been used (case 3 only)
!
! case no 5:  a directory name
!     if a file called .dir is found in said directory, it is assumed to contain a list of standard file names
!     just like in case 2, that directory name is prepended to the names found in said file.
!     if such a file is not found. all files having names in *.fst in said directory will be opened
!
      integer function fstouv_by_name(name,options)
      use ouvfrm_by_name
      implicit none
      character (len=*), intent(IN) :: name
      integer options

      include 'fstouv_by_name.inc'
      integer, parameter :: maxunits=2048
      integer, dimension(maxunits) :: units
      integer :: nunits
      integer status, i, istart, iend, ilen, iunlist
      character (len=1) :: NULC
      character (len=32) :: file_mode
      character (len=*), parameter :: defdir='.dir'  ! default name of file list in directory
      character (len=2048) :: filelist, dirname, globname2
      character (len=*), parameter :: defext='*.fst'   ! default extension to look for in directory
      character (len=32) :: extension
      type(link_element), pointer :: next_list
      logical :: lazy, is_a_list

      integer fnom, fstouv, name_is_a_file, name_is_a_dir
      integer name_is_txt_file
      integer f_set_glob, f_next_glob
      external fnom, fstouv, name_is_a_file, name_is_a_dir
      external name_is_txt_file
      external f_set_glob, f_next_glob

      status = -1   !  precondition fail
      units=0
      nunits = 0
      NULC = achar(0)
      extension=''
      iunlist = -1
      file_mode = 'STD+RND+OLD+R/O'
      lazy = iand(options,FSTOUV_LAZY) == FSTOUV_LAZY
      is_a_list = iand(options,FSTOUV_LIST) == FSTOUV_LIST

      if(0 == name_is_a_file(trim(name)//NULC)) then       ! name of a list of filenames or single standard file

        if(is_a_list) then   ! text file, this is the name of a list of file names
          filelist=name
          dirname=name
          i=len(trim(name))
          do while( i>1 .and. dirname(i:i) /= '/' )  ! suppress everything after last /
            dirname(i:i)=' '
            i=i-1
          enddo
        else  ! not a list of file names, try to open as a standard file
          print *,'FILE: ',trim(name)
          if(iand(options,FSTOUV_WRITE) == FSTOUV_WRITE) 
     %                              file_mode = 'STD+RND'              ! single file to open, write mode allowed
          lazy = .false.                                               ! but lazy mode is not
          status = try_to_open(name,lazy,file_mode)
          if(status /= 0) goto 333                                     ! ERROR exit
          goto 222                ! OK, we can return
        endif

      elseif(0 == name_is_a_dir(trim(name)//NULC)) then     ! directory name
        dirname  = trim(name)//'/'
        filelist = ''
      else                                                  ! pattern or explicit + separated list of names
        filelist = name
        dirname  = ''
      endif
!
!     dirname and filelist are now set
!     filelist is blank if a directory name was given, is equal to name otherwise
!     dirname is blank if pattern or explicit + separated list of names was given,
!             contains directory name if directory was given or directory part of list of files
!
      if(dirname /= '')then  ! we have a directory name and maybe a file list
        if(filelist == '') then
          if(0 == name_is_txt_file(trim(dirname)//defdir//NULC)) then
            filelist=trim(dirname)//defdir                  ! default name for list of files
          endif
        endif
        if(filelist /= '') then    ! read file list, call fnom and fstouv for the files
          iunlist = get_fortran_unit()   ! get a fortran unit to read list of file names
        else
          filelist = trim(dirname)//defext
          dirname = ''
        endif
      endif
!
!     filelist now contains list of files, pattern, or + separated list
!     dirname has been blanked if director/default_pattern is used
!
      if(iunlist /= -1) then  ! we have a list of names from a file, lazy mode not allowed
        lazy = .false.
        open(unit=iunlist,file=trim(filelist)//NULC,form='FORMATTED')   ! open list of names
1       read(iunlist,*,end=2,err=3)globname2                             ! for each name in file
        globname2=trim(dirname)//globname2
        status = try_to_open(globname2,lazy,file_mode)
        if(status /= 0) goto 333
        goto 1
3       close(iunlist)
        status = -1
        goto 333  ! ERROR exit
2       continue  ! OK if we get here
        close(iunlist)
        goto 222

      elseif(0 == f_set_glob(trim(filelist)//NULC)) then  ! pattern with matches
          do while(f_next_glob(globname2,len(globname2)) == 0)
            status = try_to_open(globname2,lazy,file_mode)
            if(status /= 0) goto 333
          enddo
          goto 222
          call f_free_glob()
      else                            ! plus sign (+) separated list of names, lazy mode not allowed
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
          lazy = .false.
          if(0 == name_is_a_file(trim(globname2)//NULC)) then    ! file exists ?
            status = try_to_open(globname2,lazy,file_mode)
            if(status /= 0) goto 333
          else                                                   ! bad file name, error
            print *,"BAD FILE NAME: '"//trim(globname2)//"'"
            status = -1
            goto 333
          endif
        enddo
        goto 222
      endif

      return

222   if(nunits > 1) call fstlnk(units,nunits)   ! everything OK, link open files if more than one
      status = add_a_set(units,nunits)          ! and add to our list

111   fstouv_by_name = status
      call f_free_glob()        ! harmless if not needed
      return

333   status = -1
      if(iunlist /= -1) close(iunlist)
      do i = 1 , nunits         ! errors detected, close all the standard files that we opened
        call fstfrm(units(i))
        call fclos(units(i))
      enddo
      goto 111

      contains
!------------------------------------------------------------------------------
! try to open a standard file (fnom + fstouv)
! IN file   : file name
! IN lazy   : if true, do not complain if file is not standard, just ignore it
! IN mode   : opening mode for fnom
!
! a unit table full is an error
! uses units, nunits, and maxunits from the outer routine
!------------------------------------------------------------------------------
      integer function try_to_open(file,lazy,mode)
      implicit none
      logical, intent(IN) :: lazy
      character (len=*) :: file, mode

      integer :: status, iun
      integer :: fnom, fstouv, fstfrm, name_is_std_file
      external :: fnom, fstouv, fstfrm, name_is_std_file

      status = -1
      if(nunits >= maxunits) then
        print *,"ERROR: trying to link more than",maxunits," files"
        goto 444    !  table if full, ERROR
      endif

      iun = 0
      status = fnom(iun,file,mode,0)
      if(status /=0) goto 333            ! fnom failed, ERROR

      if(lazy) then
        if(0 /= name_is_std_file(trim(file)//NULC))then  ! lazy mode, not a standard file
          status = 0                                     ! no error, but call fclos
          goto 333
        else
          print *,"INFO: file '",trim(file),"' not a standard file ?"
        endif
      endif
      status = fstouv(iun,'RND')
      if(status /=0) then
        print *,"ERROR: file '",trim(file),"' not a standard file ?"
        goto 333            ! fstouv failed, ERROR
      endif

      nunits = nunits + 1
      units(nunits) = iun
      goto 444

333   if(iun /= 0) call fclos(iun)        ! call fclos if necessary
444   try_to_open = status
      return
      end function try_to_open
!------------------------------------------------------------------------------
! add a set of (linked if more than one) unit(s) to tje list
! IN units  : array of pseudo fortran unit numbers
! IN nunits : number of valid items in units
!
! the set is added at the beginning of the global table (link_list)
!------------------------------------------------------------------------------
      integer function add_a_set(units,nunits)  ! insert at head of list
      implicit none
      integer, intent(IN) :: nunits
      integer, dimension(nunits), intent(IN) :: units

      type(link_element), pointer :: element
      integer istat

      add_a_set = -1  ! precondition status to failure
      allocate(element,stat=istat)
      if(istat /= 0) return
      allocate(element%units(nunits),stat=istat)
      if(istat /= 0) then
        deallocate(element)
        return
      endif
      element%next => link_list
      element%units = units(1:nunits)
      element%nunits = nunits
      link_list => element

      add_a_set = units(1)   ! return inits(1) if insertion successful
      end function add_a_set

      end function fstouv_by_name
!==============================================================================
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
      character (len=1) :: NULC
      character (len=*), parameter :: defdir='.dir'  ! default name of file list in directory
      character (len=32) :: filelist
      character (len=*), parameter :: defext='*.fst'   ! default extension to look for in directory
      character (len=32) :: extension
      character (len=*), parameter :: mode='RND'

      status = -1
      nunits = 0
      units(2:maxunits) = 0
      NULC = achar(0)
      extension = ext
      if(extension == " ") extension = defext
      filelist = dir
      if(filelist == " ") filelist = defdir

      if(0 == name_is_a_file(trim(names)//NULC)) then  ! names is a regular file
        print *,'FILE: ',trim(names)
        status = fnom(units(1),names,'STD+RND+OLD+R/O',0)               ! call fnom and fstouv
        if(status /= 0) goto 111
        status = fstouv(units(1),mode)
        nunits = 1
        goto 111                                           ! and return
      endif
!     if we get here, we have either a glob pattern, directory, or list of names
      globname1 = trim(names)
      if(0 == name_is_a_dir(trim(names)//NULC)) then   ! names is a directory
        print *,'DIR: ',trim(names)
        globname1 = trim(names)//'/'//trim(extension)        ! try names/*."ext" if names/".dir" not found
        listname  = trim(names)//'/'//trim(filelist)           ! name of file list to look for
        if(0 == name_is_a_file(trim(listname)//NULC)) then
          print *,'FOUND '//trim(listname)
          list = get_fortran_unit()
          open(unit=list,file=trim(listname)//NULC,form='FORMATTED')   ! open list of names
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
      status=f_set_glob(trim(globname1)//NULC)              ! prepare to resolve globbed names
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
        if(0 == name_is_a_file(trim(globname2)//NULC)) then    ! file exists ?
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
