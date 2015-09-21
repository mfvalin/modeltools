program demo_interp
  implicit none

  type :: weight_set
    real, dimension(:,:,:), pointer :: w
    real, dimension(:,:), pointer :: c, s
    integer, dimension(:,:,:), pointer :: ij
    integer, dimension(:,:), pointer :: n
    integer :: ni, nj, nmax
  end type
  type(weight_set) :: w

  interface
    function read_interp_weights(iun,nis,w) result(status)
      import :: weight_set
      implicit none
      integer, intent(IN) :: iun          ! weight file
      integer, intent(IN) :: nis
      type(weight_set), intent(OUT) :: w
      integer :: status
    end function
  end interface

  real *4, dimension(:,:), allocatable :: sshsrc, sshdest, sshdestref
  integer :: nargs
  character(len=1024) :: old_file, new_file
  integer :: i, status
  integer :: fstdin, fstdout
  integer, external :: fnom, fstouv, fstinf, fstlir, weighted_interp
  integer :: ni, nj, nk, nis, njs, nks
  real *8 pi, piov180
  real *8, parameter :: ONE=1.0

  pi = acos(-ONE)
  piov180 = pi/180.0
  print *,"pi=",pi,piov180
  nargs = command_argument_count()
  if(nargs /= 2) call print_usage
  call GET_COMMAND_ARGUMENT(1,old_file)
  call GET_COMMAND_ARGUMENT(2,new_file)

  print *,'u.nemo_weights '//trim(old_file)//' '//trim(new_file)

  print *,'INFO: opening input file '//trim(old_file)
  fstdin = 0
  i = fnom(fstdin,trim(old_file),'STD+RND+OLD',0)
  print *,'DEBUG: fstdin=',fstdin
  i = fstouv(fstdin,'RND')
  print *,'DEBUG: status of fstouv fstdin =',i

  print *,'INFO: opening source file '//trim(new_file)
  fstdout = 0
  i = fnom(fstdout,trim(new_file),'STD+RND+OLD',0)
  print *,'DEBUG: fstdout=',fstdout
  i = fstouv(fstdout,'RND')
  print *,'DEBUG: status of fstouv fstdout =',i

  status = fstinf(fstdout,nis,njs,nks,-1,"",-1,-1,-1,"P@","SSH")
  if(status < 0) goto 777
  allocate(sshsrc(nis,njs))
  print *,"INFO: reading source SSH"
  status = fstlir( sshsrc,fstdout,nis,njs,nks,-1,"",-1,-1,-1,"P@","SSH")
  print *,"INFO: sshsrc  min,max=",minval(sshsrc),maxval(sshsrc)

  status = read_interp_weights(fstdin,nis,w)
  ni = w%ni
  nj = w%nj

  allocate(sshdest(ni,nj),sshdestref(ni,nj))
  status = fstlir( sshdest,fstdin,ni,nj,nk,-1,"",-1,-1,-1,"P@","SSH")
  print *,"INFO: sshdest  min,max=",minval(sshdest),maxval(sshdest)

  status = weighted_interp(sshdestref,ni,nj,sshsrc,nis,njs,w%w,w%ij,w%n,w%nmax)
  call fstecr(sshdestref,sshdestref,-32,fstdin,344189600,0,0,ni,nj,1,0,0,0,'P@','SSH1','DIAGNOSTIQUE','Z',278,1298,0,0,133,.true.)
  print *,"INFO: sshsrc  min,max=",minval(sshsrc),maxval(sshsrc)
  print *,"INFO: sshdest  min,max=",minval(sshdest),maxval(sshdest)
  print *,"INFO: sshdestref  min,max=",minval(sshdestref),maxval(sshdestref)
  sshdestref = sshdestref - sshdest
  print *,"INFO: sshdestdiff  min,max=",minval(sshdestref),maxval(sshdestref)
777 continue
  call fstfrm(fstdin)
  call fstfrm(fstdout)
  stop
end program
!
function read_interp_weights(iun,nis,w) result(status)
  implicit none
  type :: weight_set
    real, dimension(:,:,:), pointer :: w
    real, dimension(:,:), pointer :: c, s
    integer, dimension(:,:,:), pointer :: ij
    integer, dimension(:,:), pointer :: n
    integer :: ni, nj, nmax
  end type
  integer, intent(IN) :: iun          ! weight file
  integer, intent(IN) :: nis          ! first dimension of output grid
  type(weight_set),intent(OUT) :: w
  integer :: status

  real *8, dimension(:,:), allocatable :: w8
  real, dimension(:,:), allocatable :: w4
  integer, dimension(:,:), allocatable :: ia, ja
  integer :: ni, nj, nk, i
  integer, external :: fstinf, fstlir
  character(len=4) :: varname
  real *8 pi, piov180
  real *8, parameter :: ONE=1.0

  pi = acos(-ONE)
  piov180 = pi/180.0
  status = -1

  status = fstinf(iun,ni,nj,nk,-1,"",-1,-1,-1,"P@","NAVG")       ! use NAVG to get dimensions
  print *,"INFO: ni,nj=",ni,nj
  allocate( w%n(ni,nj), w8(ni,nj), w4(ni,nj), ia(ni,nj), ja(ni,nj) )
  status = fstlir( w4,iun,ni,nj,nk,-1,"",-1,-1,-1,"P@","NAVG")   ! number of useful points
  w%n = w4 + 0.5
  w%nmax = maxval(w%n)
  print *,"INFO: w%n  max=",w%nmax
  allocate( w%ij(ni,nj,w%nmax), w%w(ni,nj,w%nmax), w%c(ni,nj), w%s(ni,nj) )
  do i = 1, w%nmax
    write(varname,100)'W',i
    status = fstlir( w8,iun,ni,nj,nk,-1,"",-1,-1,-1,"",varname)   ! Wnnn record
    w%w(:,:,i) = w8
    write(varname,100)'I',i
    status = fstlir( ia,iun,ni,nj,nk,-1,"",-1,-1,-1,"",varname)   ! Innn record
    write(varname,100)'J',i
    status = fstlir( ja,iun,ni,nj,nk,-1,"",-1,-1,-1,"",varname)   ! Jnnn record
    print *,"INFO: ia/ja min,max",minval(ia),maxval(ia),minval(ja),maxval(ja)
    w%ij(:,:,i) = ia + (ja - 1)*nis                               ! I and J combined into "collapsed index"
    print *,"INFO: w%w  min,max=",minval(w%w(:,:,i)),maxval(w%w(:,:,i))
  enddo
  status = fstlir( w%c,iun,ni,nj,nk,-1,"",-1,-1,-1,"P@","ANG")    ! rotation angles
  print *,"INFO: ang  min,max=",minval(w%c),maxval(w%c)
  w%s = sin(w%c * piov180)                                        ! transform into sine and cosine
  w%c = cos(w%c * piov180)
  print *,"INFO: cos  min,max=",minval(w%c),maxval(w%c)
  print *,"INFO: sin  min,max=",minval(w%s),maxval(w%s)
  w%ni = ni
  w%nj = nj

  deallocate(w8,w4,ia,ja)
  status = 0
100 format(A1,I3.3)
end function
!
function weighted_interp(d,ni,nj,s,nis,njs,w,ij,np,nmax) result(status)
  implicit none
  integer, intent(IN) :: ni, nj, nmax, nis, njs
  real, intent(OUT), dimension(ni,nj) :: d            ! output grid
  real, intent(IN), dimension(nis*njs) :: s           ! input grid
  real, intent(IN),  dimension(ni,nj,nmax) :: w       ! weights
  integer, intent(IN),  dimension(ni,nj,nmax) :: ij   ! source index table (I and J records from file combined)
  integer, intent(IN),  dimension(ni,nj) :: np        ! number of useful points
  integer :: status
  integer :: i, j, i0, in, k, maxpts
  integer, parameter :: BSIZE=16

  do j = 1 , nj
  do i0 = 1 , ni , BSIZE
    in = min(ni,i0+BSIZE)
    maxpts = maxval(np(i0:in,j))      ! not used for now
    d(i0:in,j) = 0.0
    do k = 1 , nmax
    do i = i0 , in
      d(i,j) = d(i,j) + ( w(i,j,k) * s(ij(i,j,k)) )   ! weighted sum
    enddo
    enddo
  enddo
  enddo
  status = 0
end function weighted_interp
subroutine print_usage()
  implicit none
  print *,'USAGE: u.nemo_weights old_standard_file new_standard_file'
  call qqexit(1)
end subroutine