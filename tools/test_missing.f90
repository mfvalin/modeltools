subroutine test_missing_values
!
!  rm -f missing.fst ; MISSING_VALUE_FLAGS="99.0 127 255" ./a.out
!
implicit none
integer, external :: fnom, fstluk, fstinf, fstouv, get_missing_value_flags, encode_missing_value, fstlir
integer :: junk, status, i
integer, parameter :: ASIZE=16             ! ASIZE should be kept even because of complex tests
integer, dimension(ASIZE) :: ia, ia2, work
integer, dimension(ASIZE) :: uia, uia2
integer *2, dimension(ASIZE) :: sa, sa2
integer *1, dimension(ASIZE) :: ba, ba2
integer *2, dimension(ASIZE) :: usa, usa2
integer *1, dimension(ASIZE) :: uba, uba2
real *4, dimension(ASIZE) :: fa, fa2, fa3
real *4, dimension(ASIZE) :: ca, ca2, ca3   ! used for 32 bit complex data
real *8, dimension(ASIZE) :: da, da2, da3
real *8, dimension(ASIZE) :: za, za2, za3   ! used for 64 bit complex data
integer :: im, uim
integer *2 :: sm, usm
integer *1 :: bm, ubm
real *4 :: fm
real *8 :: dm
integer :: ni, nj, nk
logical :: file_exists

status= get_missing_value_flags(fm,im,uim,dm,sm,usm,bm,ubm)
print 101,status,fm,dm,im,sm,bm,uim,usm,ubm
! 101 format(I2,2G14.5,I11,5I7,2G14.5)
101 format(I2,2G14.5,I11,5I7,4G14.5)
call set_missing_value_flags(fm,im,uim,dm,sm,usm,bm,ubm)
status= get_missing_value_flags(fm,im,uim,dm,sm,usm,bm,ubm)
print 101,status,fm,dm,im,sm,bm,uim,usm,ubm

do i=1,ASIZE
  fa(i)  = i*1.0+.5
  ca(i)  = i*1.0+.6
  ia(i)  = i-13
  uia(i) = i+22
  da(i)  = fa(i)+1.2345
  za(i)  = fa(i)+1.3456
  sa(i)  = i-13
  ba(i)  = i-13
  usa(i) = i+22
  uba(i) = i+22
enddo
fa(1)=fm   ; fa(ASIZE)=fm ; fa(2)=fm
 ca(1)=fm   ; ca(ASIZE)=fm ; ca(2)=fm
da(1)=dm   ; da(ASIZE)=dm ; da(ASIZE-1)=dm
 za(1)=dm   ; za(ASIZE)=dm ; za(ASIZE-1)=dm
ia(2)=im   ; ia(ASIZE-1)=im ; ia(3)=im   ; ia(ASIZE-2)=im
sa(2)=sm   ; sa(ASIZE-1)=sm ; sa(4)=sm   ; sa(ASIZE-3)=sm
ba(2)=bm   ; ba(ASIZE-1)=bm ; ba(1)=bm
uia(3)=uim ; uia(ASIZE-2)=uim
usa(3)=usm ; usa(ASIZE-2)=usm ; usa(ASIZE)=usm
uba(3)=ubm ; uba(ASIZE-2)=ubm
ia(ASIZE/2)=125      ! set to 127 to force error message
sa(ASIZE/2)=124      ! set to 127 to force error message
ba(ASIZE/2)=123      ! set to 127 to force error message
uia(ASIZE/2)=113     ! set to 255 to force error message
usa(ASIZE/2)=112     ! set to 255 to force error message
uba(ASIZE/2)=111     ! set to 255 to force error message

print *,'======================================'
do i=1,ASIZE
 print 101,i,fa(i),da(i),ia(i),sa(i),ba(i),uia(i),usa(i),uba(i),ca(i),za(i)
enddo
print *,'======================================'
status=-1
status=encode_missing_value(fa2,fa,ASIZE,6,16,0,0,0)   ! float
print *,'encoding status=',status
status=-1
status=encode_missing_value(da2,da,ASIZE,6,16,0,0,1)   ! double
print *,'encoding status=',status
status=-1
status=encode_missing_value(ia2,ia,ASIZE,4,8,0,0,0)    ! signed integer
print *,'encoding status=',status
status=-1
status=encode_missing_value(sa2,sa,ASIZE,4,8,0,1,0)    ! signed short
print *,'encoding status=',status
status=-1
status=encode_missing_value(ba2,ba,ASIZE,4,8,1,0,0)    ! signed byte
print *,'encoding status=',status
status=-1
status=encode_missing_value(uia2,uia,ASIZE,2,8,0,0,0)  ! unsigned integer
print *,'encoding status=',status
status=-1
status=encode_missing_value(usa2,usa,ASIZE,2,8,0,1,0)  ! unsigned short
print *,'encoding status=',status
status=-1
status=encode_missing_value(uba2,uba,ASIZE,2,8,1,0,0)  ! unsigned byte
print *,'encoding status=',status
do i=1,ASIZE
 print 101,i,fa2(i),da2(i),ia2(i),sa2(i),ba2(i),uia2(i),usa2(i),uba2(i),ca(i),za(i)
enddo
print *,'======================================'
call decode_missing_value(fa2,ASIZE,6,0,0,0)   ! float
if(all(fa==fa2))   print *,'float encoding/decoding test              passed'

call decode_missing_value(da2,ASIZE,6,0,0,1)   ! double
if(all(da==da2))   print *,'double encoding/decoding test             passed'

call decode_missing_value(ia2,ASIZE,4,0,0,0)   ! signed integer
if(all(ia==ia2))   print *,'signed integer encoding/decoding test     passed'

call decode_missing_value(sa2,ASIZE,4,0,1,0)   ! signed short
if(all(sa==sa2))   print *,'signed short encoding/decoding test       passed'

call decode_missing_value(ba2,ASIZE,4,1,0,0)   ! signed byte
if(all(ba==ba2))   print *,'signed byte encoding/decoding test        passed'

call decode_missing_value(uia2,ASIZE,2,0,0,0)  ! unsigned integer
if(all(uia==uia2)) print *,'unsigned integer encoding/decoding test   passed'

call decode_missing_value(usa2,ASIZE,2,0,1,0)  ! unsigned short
if(all(usa==usa2)) print *,'unsigned short encoding/decoding test     passed'

call decode_missing_value(uba2,ASIZE,2,1,0,0)  ! unsigned byte
if(all(uba==uba2)) print *,'unsigned byte encoding/decoding test      passed'

do i=1,ASIZE
 print 101,i,fa2(i),da2(i),ia2(i),sa2(i),ba2(i),uia2(i),usa2(i),uba2(i),ca(i),za(i)
enddo
1111 continue ! write data into standard file
inquire(file='missing.fst',EXIST=file_exists)
if(file_exists) then
  print *,"============= missing.fst exists, write test will be skipped ============"
  goto 2222 ! skip write test if file exists
endif
print *,'============= writing into standard file with and without missing values ============'
status=fnom(11,'missing.fst','STD+RND',0)
print *,'status fnom=',status
status=fstouv(11,'RND')
print *,'status fstouv=',status

call fstecr(ia,work,-8,11,0,0,0,ASIZE,1,1,1,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,4,.false.)   ! signed integer
call fstecr(ia,work,-8,11,0,0,0,ASIZE,1,1,2,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,68,.false.)  ! signed integer with missing

call fstecr(fa,work,-16,11,0,0,0,ASIZE,1,1,3,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,1,.false.)  ! float
call fstecr(fa,work,-32,11,0,0,0,ASIZE,1,1,19,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,5,.false.)  ! IEEE 32
call fstecr(ca,work,-32,11,0,0,0,ASIZE/2,1,1,119,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,8+64+128,.false.)  ! complex IEEE 32 with missing and compression
call fstecr(fa,work,-16,11,0,0,0,ASIZE,1,1,4,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,65,.false.) ! float with missing
call fstecr(fa,work,-32,11,0,0,0,ASIZE,1,1,20,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,69,.false.)  ! IEEE 32 with missing

call fstecr(uia,work,-8,11,0,0,0,ASIZE,1,1,5,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,2,.false.)  ! unsigned integer
call fstecr(uia,work,-8,11,0,0,0,ASIZE,1,1,6,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,66,.false.) ! unsigned integer with missing

call fst_data_length(8)
call fstecr(da,work,-16,11,0,0,0,ASIZE,1,1,7,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,1,.false.)  ! double
call fstecr(da,work,-64,11,0,0,0,ASIZE,1,1,17,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,5,.false.)  ! IEEE64
call fstecr(za,work,-64,11,0,0,0,ASIZE/2,1,1,117,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,8+64+128,.false.)  ! complex IEEE64 with missing and compression
call fst_data_length(8)
call fstecr(da,work,-16,11,0,0,0,ASIZE,1,1,8,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,65,.false.) ! double with missing
call fstecr(da,work,-64,11,0,0,0,ASIZE,1,1,18,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,69,.false.) ! IEEE64 with missing

call fst_data_length(2)
call fstecr(usa,work,-8,11,0,0,0,ASIZE,1,1,9,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,2,.false.)  ! unsigned short
call fst_data_length(2)
call fstecr(usa,work,-8,11,0,0,0,ASIZE,1,1,10,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,66,.false.)! unsigned short with missing

call fst_data_length(2)
call fstecr(sa,work,-8,11,0,0,0,ASIZE,1,1,11,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,4,.false.)  ! signed short
call fst_data_length(2)
call fstecr(sa,work,-8,11,0,0,0,ASIZE,1,1,12,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,68,.false.) ! signed short with missing

call fst_data_length(1)
call fstecr(ba,work,-8,11,0,0,0,ASIZE,1,1,13,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,4,.false.)  ! signed byte
call fst_data_length(1)
call fstecr(ba,work,-8,11,0,0,0,ASIZE,1,1,14,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,68,.false.) ! signed byte with missing

call fst_data_length(1)
call fstecr(uba,work,-8,11,0,0,0,ASIZE,1,1,15,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,2,.false.)  ! unsigned byte
call fst_data_length(1)
call fstecr(uba,work,-8,11,0,0,0,ASIZE,1,1,16,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,66,.false.) ! unsigned byte with missing

call fstfrm(11)
2222 continue ! read data from standard file
!  set new flag values before reading so we can see the difference 'type n' vs 'type n+64'
fm=-999.99; dm=-888.88; im=-999999; sm=-32768 ; bm=-128; uim=999999 ; usm=32767; ubm=127
print *,'============================== missing values ============================='
print *,'#     float        double            int  short   byte   uint ushort  ubyte'
call set_missing_value_flags(fm,im,uim,dm,sm,usm,bm,ubm)
print 101,status,fm,dm,im,sm,bm,uim,usm,ubm
print *,'============= reading from standard file with and without missing values ============'
status=fnom(12,'missing.fst','STD+RND+OLD',0)
print *,'status fnom=',status
status=fstouv(12,'RND')
print *,'status fstouv=',status

print *,'=========== WITHOUT missing value =========='
status=fstlir(ia2,12,ni,nj,nk,-1,' ',1,-1,-1,' ',' ')     ! integer
status=fstlir(uia2,12,ni,nj,nk,-1,' ',5,-1,-1,' ',' ')    ! unsigned integer

status=fstlir(fa2,12,ni,nj,nk,-1,' ',3,-1,-1,' ',' ')     ! float
status=fstlir(fa3,12,ni,nj,nk,-1,' ',19,-1,-1,' ',' ')     ! IEEE 32
 ca3=-1.0 ; status=fstlir(ca3,12,ni,nj,nk,-1,' ',119,-1,-1,' ',' ')     ! complex IEEE 32
if(ni /= ASIZE/2) print *,'ERROR: complex32 expected ni=',ASIZE/2,' got:',ni
call fst_data_length(8)
status=fstlir(da2,12,ni,nj,nk,-1,' ',7,-1,-1,' ',' ')     ! double
status=fstlir(da3,12,ni,nj,nk,-1,' ',17,-1,-1,' ',' ')     ! IEEE 64
 za3=-1.0 ; status=fstlir(za3,12,ni,nj,nk,-1,' ',117,-1,-1,' ',' ')     ! complex IEEE 64
if(ni /= ASIZE/2) print *,'ERROR: complex64 expected ni=',ASIZE/2,' got:',ni

call fst_data_length(2)
status=fstlir(usa2,12,ni,nj,nk,-1,' ',9,-1,-1,' ',' ')    ! unsigned short
call fst_data_length(2)
status=fstlir(sa2,12,ni,nj,nk,-1,' ',11,-1,-1,' ',' ')    ! signed short

call fst_data_length(1)
status=fstlir(ba2,12,ni,nj,nk,-1,' ',13,-1,-1,' ',' ')    ! signed byte
call fst_data_length(1)
status=fstlir(uba2,12,ni,nj,nk,-1,' ',15,-1,-1,' ',' ')   ! unsigned byte

print *,'---------------------------------------------------------------------------------------------------'
print *,'#       float  <16> double          int  short   byte   uint ushort  ubyte    float  <IEEE> double     complex 32    complex 64'
print *,'---------------------------------------------------------------------------------------------------'
do i=1,ASIZE
 print 101,i,fa2(i),da2(i),ia2(i),sa2(i),ba2(i),uia2(i),usa2(i),uba2(i),fa3(i),da3(i),ca3(i),za3(i)
enddo
print *,'============= WITH possible missing value(s) ======'
status=fstlir(ia2,12,ni,nj,nk,-1,' ',2,-1,-1,' ',' ')     ! integer with missing
status=fstlir(uia2,12,ni,nj,nk,-1,' ',6,-1,-1,' ',' ')    ! unsigned integer with missing

status=fstlir(fa2,12,ni,nj,nk,-1,' ',4,-1,-1,' ',' ')     ! float with missing
status=fstlir(fa3,12,ni,nj,nk,-1,' ',20,-1,-1,' ',' ')     ! IEE 32 with missing
call fst_data_length(8)
status=fstlir(da2,12,ni,nj,nk,-1,' ',8,-1,-1,' ',' ')     ! double with missing
status=fstlir(da3,12,ni,nj,nk,-1,' ',18,-1,-1,' ',' ')     ! double with missing

call fst_data_length(2)
status=fstlir(sa2,12,ni,nj,nk,-1,' ',12,-1,-1,' ',' ')    ! signed short with missing
call fst_data_length(2)
status=fstlir(usa2,12,ni,nj,nk,-1,' ',10,-1,-1,' ',' ')   ! unsigned short with missing

call fst_data_length(1)
status=fstlir(ba2,12,ni,nj,nk,-1,' ',14,-1,-1,' ',' ')    ! signed byte with missing
call fst_data_length(1)
status=fstlir(uba2,12,ni,nj,nk,-1,' ',16,-1,-1,' ',' ')   ! unsigned byte with missing

print *,'---------------------------------------------------------------------------------------------------'
print *,'#    float  <f16> double             int  short   byte   uint ushort  ubyte     float  <IEEE> double     complex 32    complex 64'
print *,'---------------------------------------------------------------------------------------------------'
do i=1,ASIZE
 print 101,i,fa2(i),da2(i),ia2(i),sa2(i),ba2(i),uia2(i),usa2(i),uba2(i),fa3(i),da3(i),ca3(i),za3(i)
enddo

call fstfrm(12)
stop
end
