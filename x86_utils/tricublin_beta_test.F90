program tricublin_d_test
  use ISO_C_BINDING
  implicit none
  include "tricublib_interface.inc"
  integer, parameter :: NI = 300
  integer, parameter :: NJ = 200
  integer, parameter :: NK = 110
  integer, parameter :: NR = 10
  integer*8, external :: rdtscp
  integer*8 :: t0, t1
  real*8 :: dx, dy, dz
  real*8, dimension(8) :: px, py, pz
  real*8, dimension(8) :: pxy
  real*8, dimension(2) :: xy
  real, dimension(-2:NI+3,-1:NJ+2,NK), target :: f1, f2, f3
  real, dimension(3,-2:NI+3,-1:NJ+2,NK), target :: f123
  real*8, dimension(NK) :: zlev
  integer :: i, j, k, ii, jj, kk, rr
  real*8 fx, fy, fz, fxyz, x, y, z, ovni, ovnj, ovnk, xx, yy, zz
  real*4, dimension(3) :: r1, r2, r3, e
  type(C_PTR) :: lv
  real, dimension(3,NI,NJ,NK) :: pxpypz
  real, dimension(NI,NJ,NK) :: expected, d, dmin, dmax, dlin
  real, dimension(3,NI,NJ,NK) :: d3
  real*8, dimension(NK) :: levels
  real :: delta, error, delta1, delta2, delta3
  real*8 :: avg, avg1, avg2, avg3
  integer :: exact, minmaxerr

  fx(x) = (x+1.0)*(x+1.1)*(x+1.2)*(x+1.3)
  fy(y) = (y+1.05)*(y+1.15)*(y+1.25)*(y+1.35)
  fz(z) = (z*1.25+1.02)
!   fxyz(x,y,z) = x + y + z
!   fxyz(x,y,z) = x * y * z
!   fxyz(x,y,z) = fx(x*ovni) * fy(y*ovnj) * fz(z*ovnk)
  fxyz(x,y,z) = fx(x*ovni) + fy(y*ovnj) * fz(z*ovnk)

  px = 0
  py = 0
  pz = 0
  ovni = 1.0_8 / (NI - 1)
  ovnj = 1.0_8 / (NJ - 1)
  ovnk = 1.0_8 / (NK - 1)
  
!   dx = .375
!   dy = .25
!   dz = .125
!   call tricubic_coeffs_d(px,py,pz,dx,dy,dz)
!   print 103,'px   = ',px(1:4),'py   = ',py(1:4),'pz   = ',pz(1:4)
!   print 103,'pzl  = ',pz(5:8)
!   print 103,'sums = ', sum(px(1:4)),sum(py(1:4)),sum(pz(1:4)),sum(pz(5:8))

  pxpypz = 0.0

  do k = 1, NK
    levels(k) = k
    do j = 1, NJ
      do i = 1, NI
        dx = .375 + i * .001
        dy = .251 + j * .001
        dz = .125 + k * .003
        xx = (i)
        yy = (j)
	zz = (k)
	f1(i,j,k) = real(fxyz(xx,yy,zz))
	f2(i,j,k) = f1(i,j,k) + 10.0
	f3(i,j,k) = f1(i,j,k) + 20.0
	f123(1,i,j,k) = f1(i,j,k)
	f123(2,i,j,k) = f2(i,j,k)
	f123(3,i,j,k) = f3(i,j,k)
! 	xx = min( max(i,2), ni - 2) + dx + 3         ! offset x de 3
! 	yy = min( max(j,2), nj - 2) + dy + 2         ! offset y de 2
	xx = min( max(i,2), ni - 2) + dx         ! offset x de 3
	yy = min( max(j,2), nj - 2) + dy         ! offset y de 2
	if(k > 1) then
	  zz = zz + dz
	else
	  zz = zz - dz
	endif
	pxpypz(1,i,j,k) = xx + 3        ! offset x de 3
	pxpypz(2,i,j,k) = yy + 2        ! offset y de 2
	pxpypz(3,i,j,k) = zz
! 	expected(i,j,k) = real(fxyz(xx-3,yy-2,zz)) ! compensation d'offset
	expected(i,j,k) = real(fxyz(xx,yy,zz))
      enddo
    enddo
  enddo
  dmin = 0
  dmax = 0
  dlin = 0

  lv = vsearch_setup_plus(levels, NK, NI+6, NJ+4, -3, -2)    ! "NI" = NI + 6, "NJ" = nj + 4 offseti = -3, offsetj = -2

  print *,'======================== 1 variable ================================'

  t0 = rdtscp()
  call tricublin_zyx1_n(d,f1(1,1,1),pxpypz,lv,NI*NJ*NK)
  t1 = rdtscp()
  print *,"cycles per point =",(t1-t0)/(NI*NJ*NK)

  exact = 0
  delta = 0.0
  avg = 0.0
  do k = 1, NK
    do j = 1, NJ
      do i = 1, NI
        error = abs(expected(i,j,k) - d(i,j,k))
        if(error == 0.0) exact = exact + 1
        delta = max(delta,error / expected(i,j,k))
        avg = avg + (error / expected(i,j,k))
        if(delta > .000001 .or. i+j+k == 3) then
          print *,'i,j,k',i,j,k
          print *,'pxpypz =',pxpypz(1,i,j,k),pxpypz(2,i,j,k),pxpypz(3,i,j,k)
          print *,'expected =',expected(i,j,k)
          print *,'result   =',d(i,j,k)
!           if(delta > .000001) stop
          if(delta > .000001) goto 111
        endif
      enddo
    enddo
  enddo
  print *,'exact =',exact,' out of',NI*NJ*NK
  print *,'%      ',real(exact)/real(NI*NJ*NK)*100
  print*,'maxerr =',delta
  print*,'avgerr =',real(avg/(NI*NJ*NK))
111 continue
  print *,'======================== 3 variables ==============================='

  t0 = rdtscp()
  call tricublin_zyx3_n(d3,f123(1,1,1,1),pxpypz,lv,NI*NJ*NK)
  t1 = rdtscp()
  print *,"cycles per point =",(t1-t0)/(NI*NJ*NK*3)

  exact = 0
  delta = 0.0
  delta1 = 0.0
  delta2 = 0.0
  delta3 = 0.0
  avg = 0.0
  avg1 = 0.0
  avg2 = 0.0
  avg3 = 0.0
  do k = 1, NK
    do j = 1, NJ
      do i = 1, NI
        error = abs(expected(i,j,k)+00.0 - d3(1,i,j,k))
        if(error == 0.0) exact = exact + 1
        delta1 = max(delta1,error / expected(i,j,k))
        avg1 = avg1 + (error / expected(i,j,k))

        error = abs(expected(i,j,k)+10.0 - d3(2,i,j,k))
        if(error == 0.0) exact = exact + 1
        delta2 = max(delta2,error / (expected(i,j,k)+10.0))
        avg2 = avg2 + error / (expected(i,j,k)+10.0)

        error = abs(expected(i,j,k)+20.0 - d3(3,i,j,k))
        if(error == 0.0) exact = exact + 1
        delta3 = max(delta3,error / (expected(i,j,k)+20.0))
        avg3 = avg3 + error / (expected(i,j,k)+20.0)

        delta = max(delta1,delta2,delta3)
        if(delta > .000001 .or. i+j+k == 3) then
          print *,'i,j,k',i,j,k
          print *,'pxpypz =',pxpypz(1,i,j,k),pxpypz(2,i,j,k),pxpypz(3,i,j,k)
          print *,'expected =',expected(i,j,k) + [0.0, 10.0, 20.0]
          print *,'result   =',d3(:,i,j,k)
!           if(delta > .000001) stop
          if(delta > .000001) goto 222
        endif
      enddo
    enddo
  enddo
  print *,'exact =',exact,' out of',NI*NJ*NK*3
  print *,'%      ',real(exact)/real(NI*NJ*NK*3)*100
  print*,'maxerr =',delta1, delta2, delta3
  print*,'avgerr =',real([avg1,avg2,avg3] / (NI*NJ*NK))
222 continue
  print *,'======================== 1 mono cubic =============================='

  d = 0
  dlin = 0
  dmin = 0
  dmax = 0
  t0 = rdtscp()
  call tricublin_mono_zyx_n(d,dlin,dmin,dmax,f1(1,1,1),pxpypz,lv,NI*NJ*NK)
  t1 = rdtscp()
  print *,"cycles per point =",(t1-t0)/(NI*NJ*NK)

  exact = 0
  delta = 0.0
  avg = 0.0
  do k = 1, NK
    do j = 1, NJ
      do i = 1, NI
        error = abs(expected(i,j,k) - d(i,j,k))
        if(error == 0.0) exact = exact + 1
        delta = max(delta,error / expected(i,j,k))
        avg = avg + (error / expected(i,j,k))
        if(delta > .000001 .or. i+j+k == 3) then
          print *,'i,j,k',i,j,k
          print *,'pxpypz =',pxpypz(1,i,j,k),pxpypz(2,i,j,k),pxpypz(3,i,j,k)
          print *,'expected =',expected(i,j,k)
          print *,'result   =',d(i,j,k),dlin(i,j,k),dmin(i,j,k),dmax(i,j,k)
!           if(delta > .000001) stop
          if(delta > .000001) goto 333
        endif
      enddo
    enddo
  enddo
  print *,'exact =',exact,' out of',NI*NJ*NK
  print *,'%      ',real(exact)/real(NI*NJ*NK)*100
  print*,'maxerr =',delta
  print*,'avgerr =',real(avg/(NI*NJ*NK))
333 continue
  print *,'======================== 1 mono linear ============================='

  exact = 0
  minmaxerr = 0
  delta = 0.0
  avg = 0.0
  do k = 1, NK
    do j = 1, NJ
      do i = 1, NI
        error = abs(expected(i,j,k) - dlin(i,j,k))
        if( dlin(i,j,k) < dmin(i,j,k) .or. dlin(i,j,k) > dmax(i,j,k) ) minmaxerr = minmaxerr + 1
        if(error == 0.0) exact = exact + 1
        delta = max(delta,error / expected(i,j,k))
        avg = avg + (error / expected(i,j,k))
!         if(delta > .000001 .or. i+j+k == 3) then
!           print *,'i,j,k',i,j,k
!           print *,'pxpypz =',pxpypz(1,i,j,k),pxpypz(2,i,j,k),pxpypz(3,i,j,k)
!           print *,'expected =',expected(i,j,k)
!           print *,'result   =',dlin(i,j,k)
!           if(delta > .000001) stop
!         endif
      enddo
    enddo
  enddo
  print *,'exact =',exact,' out of',NI*NJ*NK
  print *,'%      ',real(exact)/real(NI*NJ*NK)*100
  print*,'maxerr =',delta
  print*,'points outside of min-max range =',minmaxerr
  print*,'avgerr =',real(avg/(NI*NJ*NK))

  print *,'====================================================================='

#if defined(OLD_TEST)
  ii = 2
  jj = 3
  kk = 4
  xx = ii + 1 + dx
  yy = jj + 1 + dy
  zz = kk + 1 + dz
!   print 101,'target = ',xx,yy,zz

  print *,'===== tricublin_zyxf_beta test (cubic) ====='
  call tricublin_zyxf_beta(r3(1),f1(ii,jj,kk),px,py,pz,ni,ni*nj)
  call tricublin_zyxf_beta(r3(2),f2(ii,jj,kk),px,py,pz,ni,ni*nj)
  call tricublin_zyxf_beta(r3(3),f3(ii,jj,kk),px,py,pz,ni,ni*nj)
  e(1) = real(fxyz(xx,yy,zz))
  e(2) = e(1) + 10.0
  e(3) = e(1) + 20.0
  print 101,'expected :',e,' , got :',r3
  print 102,'rel error :',(e-r3)/r3
  print *, 'FLOPS/interp = ',3*(32*4 + 32 + 8 + 64)  ! interpolation + float -> double conversions

  print *,'===== tricublin_zyx3f_beta test (cubic) ====='
  call tricublin_zyx3f_beta(r3,f123(1,ii,jj,kk),px,py,pz,ni,ni*nj)
  e(1) = real(fxyz(xx,yy,zz))
  e(2) = e(1) + 10.0
  e(3) = e(1) + 20.0
  print 101,'expected :',e,' , got :',r3
  print 102,'rel error :',(e-r3)/r3
  print *, 'FLOPS/interp = ',3*(32*4 + 32 + 8 + 64)  ! interpolation + float -> double conversions

  print *,'===== tricublin_zyxf_beta test (linear) ====='
  kk = 1
  zz = kk + 0 + dz
  call tricublin_zyxf_beta(r3(1),f1(ii,jj,kk),px,py,pz(5),ni,ni*nj)
  call tricublin_zyxf_beta(r3(2),f2(ii,jj,kk),px,py,pz(5),ni,ni*nj)
  call tricublin_zyxf_beta(r3(3),f3(ii,jj,kk),px,py,pz(5),ni,ni*nj)
  e(1) = real(fxyz(xx,yy,zz))
  e(2) = e(1) + 10.0
  e(3) = e(1) + 20.0
  print 101,'expected :',e,' , got :',r3
  print 102,'rel error :',(e-r3)/r3
  print *, 'FLOPS/interp = ',3*(32*4 + 32 + 8 + 64)  ! interpolation + float -> double conversions

  print *,'===== tricublin_zyxf_beta timing (cubic) ====='
  kk = 4
  f2 = f1 + 10
  f3 = f2 + 10
  t0 = rdtscp()
  do rr = 1 , NR
  do KK = 1, NK - 3
    do JJ = 1, NJ - 3
      do II = 1, NI - 3
        call tricubic_coeffs_d(px,py,pz,dx,dy,dz)
        call tricublin_zyxf_beta(r1(1),f1(ii,jj,kk),px,py,pz,ni,ni*nj)
        call tricublin_zyxf_beta(r2(1),f2(ii,jj,kk),px,py,pz,ni,ni*nj)
        call tricublin_zyxf_beta(r3(1),f3(ii,jj,kk),px,py,pz,ni,ni*nj)
      enddo
    enddo
  enddo
  enddo
  t1 = rdtscp() - t0
  
  print *, 'tot cycles   = ',t1
  print *, 'per interp   = ',t1/((ni-3)*(nj-3)*(nk-3))/NR
  print *, 'FLOPS/interp = ',3*(32*4 + 32 + 8 + 64)  ! interpolation + float -> double conversions

  print *,'===== tricublin_zyx3f_beta timing (cubic) ====='
  kk = 4
  f2 = f1 + 10
  f3 = f2 + 10
  t0 = rdtscp()
  do rr = 1 , NR
  do KK = 1, NK - 3
    do JJ = 1, NJ - 3
      do II = 1, NI - 3
        call tricubic_coeffs_d(px,py,pz,dx,dy,dz)
        call tricublin_zyx3f_beta(r3,f123(1,ii,jj,kk),px,py,pz,ni,ni*nj)
      enddo
    enddo
  enddo
  enddo
  t1 = rdtscp() - t0
  
  print *, 'tot cycles   = ',t1
  print *, 'per interp   = ',t1/((ni-3)*(nj-3)*(nk-3))/NR
  print *, 'FLOPS/interp = ',3*(32*4 + 32 + 8 + 64)  ! interpolation + float -> double conversions

  if(t1 .eq. 12) then
    call tricublin_zyx1_n(r1,f1,pxpypz,lv,ni)
  endif
#endif

101 format(A,3F20.10,A,3F20.10,A,3F20.10)
103 format(A,4F20.10,/A,4F20.10,/A,4F20.10)
102 format(A,3E12.5,A,3E12.5)
end program tricublin_d_test
