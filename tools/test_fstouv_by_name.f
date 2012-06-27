      program test_fstouv_by_name
      implicit none
      include 'fstouv_by_name.inc'
      integer status
      character *132 name
      integer iun
      integer fstouv_by_name
      integer, parameter :: MAXUNITS=1024
      integer nunits,units(MAXUNITS)

      iun=0
      print *,'=============== file1.tfs =============='
      name='file1.tfs'
      units=0
      status = fstouv_by_name(name,FSTOUV_WRITE)
      call dump_ouvfrm_by_name
      call fstfrm_by_name(status)
      call dump_ouvfrm_by_name
      print *,'===============   *.tfs   =============='
      name='*.tfs'
      units=0
      status = fstouv_by_name(name,FSTOUV_WRITE)  ! test that FSTOUV_WRITE is ignored
      call dump_ouvfrm_by_name
      call fstfrm_by_name(status)
      call dump_ouvfrm_by_name
      print *,'=============== dir/*.fst =============='
      name='mydir'
      units=0
      status = fstouv_by_name(name,0)
      call dump_ouvfrm_by_name
      stop
      print *,'=============== dir/*.Fst =============='
      name='mydir/*.Fst'
      units=0
      status = fstouv_by_name(name,0)
      call dump_ouvfrm_by_name
      print *,'=============== file+file2+.. =========='
      name='+file1.tfs++mydir/file1.Fst+file2.tfs+mydir2/file1.bbb++++'
      units=0
      status = fstouv_by_name(name,0)
      call dump_ouvfrm_by_name
      print *,'================ dir/.dir (.aaa)========'
      name='mydir2'
      units=0
      status = fstouv_by_name(name,0)
      call dump_ouvfrm_by_name
      print *,'================ dir/.Dir (.bbb)========'
      name='mydir2/.Dir'
      units=0
      status = fstouv_by_name(name,FSTOUV_LIST)
      call dump_ouvfrm_by_name
      print *,'=============== file+bad+.. =========='
      name='+file1.tfs++mydir/file1.Fst+file2.tfs+bad++++'
      units=0
      status = fstouv_by_name(name,0)
      call dump_ouvfrm_by_name
      print *,'========================================'
      stop
      end
      integer function fnom(iun,name,attr,lng)
      implicit none
      integer :: iun, lng
      character (len=*) :: name, attr
      integer, save :: baseunit=9999
      if(iun == 0) then
        iun = baseunit
        baseunit = baseunit - 1
      endif
      print *,'FNOM: iun,name,attr',iun,trim(name),' ',trim(attr)
      fnom = 0
      return
      end
      integer function fstouv(iun,mode)
      implicit none
      integer :: iun
      character (len=*) :: mode
      print *,'FSTOUV: iun,mode',iun,' ',trim(mode)
      fstouv = 0
      return
      end
      subroutine  fstlnk(units,nunits)
      integer :: nunits
      integer, dimension(nunits) :: units
!      print *,'FSTLNK: nunits=',nunits
      print *,'FSTLNK: units=',units(1:nunits)
      return
      end
      subroutine  fstunl(units,nunits)
      integer :: nunits
      integer, dimension(nunits) :: units
!      print *,'FSTLNK: nunits=',nunits
      print *,'FSTUNL: units=',units(1:nunits)
      return
      end
      subroutine  fstfrm(unit)
      integer :: unit
      print *,'FSTFRM: unit=',unit
      return
      end
      subroutine  fclos(unit)
      integer :: unit
      print *,'FCLOS: unit=',unit
      return
      end
