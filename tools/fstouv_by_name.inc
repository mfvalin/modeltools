      integer, parameter :: FSTOUV_LIST=1    ! this file IS a list of file names
      integer, parameter :: FSTOUV_WRITE=2   ! writing is permitted (single file case only)
      integer, parameter :: FSTOUV_CREATE=4  ! writing is permitted (single file case only)
      integer, parameter :: FSTOUV_LAZY=8    ! just ignore non RPN standard files with no error (only valid for pattern)
