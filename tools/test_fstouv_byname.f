      program test_fstouv_byname
      implicit none
      integer status
      character *132 name
      integer iun
      integer fstouv_byname
      integer, parameter :: MAXUNITS=1024
      integer nunits,units(MAXUNITS)

      iun=0
      print *,'=============== file1.tfs =============='
      name='file1.tfs'
      status = fstouv_byname(name,units,MAXUNITS,nunits,'','')
      print *,'===============   *.tfs   =============='
      name='*.tfs'
      status = fstouv_byname(name,units,MAXUNITS,nunits,'','')
      print *,'=============== dir/*.fst =============='
      name='mydir'
      status = fstouv_byname(name,units,MAXUNITS,nunits,'','')
      print *,'=============== dir/*.Fst =============='
      name='mydir'
      status = fstouv_byname(name,units,MAXUNITS,nunits,'','*.Fst')
      print *,'================ dir/.dir (.aaa)========'
      name='mydir2'
      status = fstouv_byname(name,units,MAXUNITS,nunits,'','')
      print *,'================ dir/.Dir (.bbb)========'
      name='mydir2'
      status = fstouv_byname(name,units,MAXUNITS,nunits,'.Dir','')
      print *,'========================================'
      stop
      end
      integer function fnom(iun,name,attr,lng)
      implicit none
      integer :: iun, lng
      character (len=*) :: name, attr
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
      print *,'FSTLNK: nunits=',nunits
      return
      end
      subroutine  fstfrm(unit)
      integer :: unit
      print *,'FSTFRM: unit=',unit
      return
      end
