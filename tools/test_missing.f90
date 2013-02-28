program test_missing_values
!
!  rm -f missing.fst ; MISSING_VALUE_FLAGS="99.0 127 255" ./a.out
!
implicit none
integer, external :: fnom, fstluk, fstinf, fstouv, get_missing_value_flags, encode_missing_value, fstlir
integer :: junk, status, i
integer, parameter :: ASIZE=16
integer, dimension(ASIZE) :: ia, ia2, work
integer, dimension(ASIZE) :: uia, uia2
integer *2, dimension(ASIZE) :: sa, sa2
integer *1, dimension(ASIZE) :: ba, ba2
integer *2, dimension(ASIZE) :: usa, usa2
integer *1, dimension(ASIZE) :: uba, uba2
real *4, dimension(ASIZE) :: fa, fa2
real *8, dimension(ASIZE) :: da, da2
integer :: im, uim
integer *2 :: sm, usm
integer *1 :: bm, ubm
real *4 :: fm
real *8 :: dm
integer :: ni, nj, nk

status= get_missing_value_flags(fm,im,uim,dm,sm,usm,bm,ubm)
print 101,status,fm,dm,im,sm,bm,uim,usm,ubm
100 format(I2,G14.5,I12,I12,G14.5,I12,I12,I8,I8,I8,I8)
101 format(I2,2G14.5,8I9)
call set_missing_value_flags(fm,im,uim,dm,sm,usm,bm,ubm)
status= get_missing_value_flags(fm,im,uim,dm,sm,usm,bm,ubm)
print 101,status,fm,dm,im,sm,bm,uim,usm,ubm

do i=1,ASIZE
  fa(i)  = i*1.0+.5
  ia(i)  = i-13
  uia(i) = i+22
  da(i)  = fa(i)+1.2345
  sa(i)  = i-13
  ba(i)  = i-13
  usa(i) = i+22
  uba(i) = i+22
enddo
fa(1)=fm   ; fa(ASIZE)=fm ; fa(2)=fm
da(1)=dm   ; da(ASIZE)=dm ; da(ASIZE-1)=dm
ia(2)=im   ; ia(ASIZE-1)=im ; ia(3)=im   ; ia(ASIZE-2)=im
sa(2)=sm   ; sa(ASIZE-1)=sm ; sa(4)=sm   ; sa(ASIZE-3)=sm
ba(2)=bm   ; ba(ASIZE-1)=bm ; ba(1)=bm
uia(3)=uim ; uia(ASIZE-2)=uim
usa(3)=usm ; usa(ASIZE-2)=usm ; usa(ASIZE)=usm
!uba(3)=usm ; uba(ASIZE-2)=usm
ia(ASIZE/2)=125      ! set to 127 to force error message
sa(ASIZE/2)=124      ! set to 127 to force error message
ba(ASIZE/2)=123      ! set to 127 to force error message
uia(ASIZE/2)=113     ! set to 255 to force error message
usa(ASIZE/2)=112     ! set to 255 to force error message
uba(ASIZE/2)=111     ! set to 255 to force error message

print *,'======================================'
do i=1,ASIZE
 print 101,i,fa(i),da(i),ia(i),sa(i),ba(i),uia(i),usa(i),uba(i)
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
 print 101,i,fa2(i),da2(i),ia2(i),sa2(i),ba2(i),uia2(i),usa2(i),uba2(i)
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
 print 101,i,fa2(i),da2(i),ia2(i),sa2(i),ba2(i),uia2(i),usa2(i),uba2(i)
enddo
print *,'============= writing into standard file with and without missing values ============'
status=fnom(11,'missing.fst','STD+RND',0)
print *,'status fnom=',status
status=fstouv(11,'RND')
print *,'status fstouv=',status

call fstecr(ia,work,-8,11,0,0,0,ASIZE,1,1,1,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,4,.false.)   ! signed integer
call fstecr(ia,work,-8,11,0,0,0,ASIZE,1,1,2,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,68,.false.)  ! signed integer with missing

call fstecr(fa,work,-16,11,0,0,0,ASIZE,1,1,3,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,1,.false.)  ! float
call fstecr(fa,work,-16,11,0,0,0,ASIZE,1,1,4,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,65,.false.) ! float with missing

call fstecr(uia,work,-8,11,0,0,0,ASIZE,1,1,5,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,2,.false.)  ! unsigned integer
call fstecr(uia,work,-8,11,0,0,0,ASIZE,1,1,6,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,66,.false.) ! unsigned integer with missing

call fst_data_length(8)
call fstecr(da,work,-16,11,0,0,0,ASIZE,1,1,7,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,1,.false.)  ! double
call fst_data_length(8)
call fstecr(da,work,-16,11,0,0,0,ASIZE,1,1,8,0,0,'XX','YYYY','ETIKET','X',0,0,0,0,65,.false.) ! double with missing

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
!  set new flag values before reading so we can see the difference 'type n' vs 'type n+64'
fm=-999.99; dm=-888.88; im=-999999; sm=-32768 ; bm=-128; uim=999999 ; usm=32767; ubm=127
call set_missing_value_flags(fm,im,uim,dm,sm,usm,bm,ubm)
print *,'============= reading from standard file with and without missing values ============'
status=fnom(12,'missing.fst','STD+RND+OLD',0)
print *,'status fnom=',status
status=fstouv(12,'RND')
print *,'status fstouv=',status

print *,'=========== WITHOUT missing value =========='
status=fstlir(ia2,12,ni,nj,nk,-1,' ',1,-1,-1,' ',' ')     ! integer
status=fstlir(uia2,12,ni,nj,nk,-1,' ',5,-1,-1,' ',' ')    ! unsigned integer

status=fstlir(fa2,12,ni,nj,nk,-1,' ',3,-1,-1,' ',' ')     ! float
call fst_data_length(8)
status=fstlir(da2,12,ni,nj,nk,-1,' ',7,-1,-1,' ',' ')     ! double

call fst_data_length(2)
status=fstlir(usa2,12,ni,nj,nk,-1,' ',9,-1,-1,' ',' ')    ! unsigned short
call fst_data_length(2)
status=fstlir(sa2,12,ni,nj,nk,-1,' ',11,-1,-1,' ',' ')    ! signed short

call fst_data_length(1)
status=fstlir(ba2,12,ni,nj,nk,-1,' ',13,-1,-1,' ',' ')    ! signed byte
call fst_data_length(1)
status=fstlir(uba2,12,ni,nj,nk,-1,' ',15,-1,-1,' ',' ')   ! unsigned byte
do i=1,ASIZE
 print 101,i,fa2(i),da2(i),ia2(i),sa2(i),ba2(i),uia2(i),usa2(i),uba2(i)
enddo
print *,'============= WITH possible missing value(s) ======'
status=fstlir(ia2,12,ni,nj,nk,-1,' ',2,-1,-1,' ',' ')     ! integer with missing
status=fstlir(uia2,12,ni,nj,nk,-1,' ',6,-1,-1,' ',' ')    ! unsigned integer with missing

status=fstlir(fa2,12,ni,nj,nk,-1,' ',4,-1,-1,' ',' ')     ! float with missing
call fst_data_length(8)
status=fstlir(da2,12,ni,nj,nk,-1,' ',8,-1,-1,' ',' ')     ! double with missing

call fst_data_length(2)
status=fstlir(sa2,12,ni,nj,nk,-1,' ',12,-1,-1,' ',' ')    ! signed short with missing
call fst_data_length(2)
status=fstlir(usa2,12,ni,nj,nk,-1,' ',10,-1,-1,' ',' ')   ! unsigned short with missing

call fst_data_length(1)
status=fstlir(ba2,12,ni,nj,nk,-1,' ',14,-1,-1,' ',' ')    ! signed byte with missing
call fst_data_length(1)
status=fstlir(uba2,12,ni,nj,nk,-1,' ',16,-1,-1,' ',' ')   ! unsigned byte with missing

do i=1,ASIZE
 print 101,i,fa2(i),da2(i),ia2(i),sa2(i),ba2(i),uia2(i),usa2(i),uba2(i)
enddo

call fstfrm(12)
stop
end
