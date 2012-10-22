      program stoopid
      implicit none
      include 'mpif.h'
      integer, parameter :: niter=50
      integer, parameter :: nis=3000
      integer, parameter :: njs=400
      integer, parameter :: nks=120
      integer iter,ierr,irank,iterext,isize
      real *8 rr(nks)
      real *8, dimension(niter) :: timea, timeb, timec, timed, timee, timef, time1, time2, time3
      real *8, dimension(:,:), pointer :: time1all, time2all, time3all
      real *8 tjitter,t0,t1,t2,tstart,tend
      real *8 tmax, tmin, tmean, tsum, tsum2, tdev
      real *8, dimension(:,:,:), pointer :: r,a,b

external report_cpu_binding
integer report_cpu_binding
external omp_get_thread_num
integer omp_get_thread_num
external omp_get_max_threads
integer omp_get_max_threads
external bind_thread_to_cpu
integer bind_thread_to_cpu

!print *,'MAX threads=',omp_get_max_threads()
!      call omp_set_num_threads(max(1,omp_get_max_threads()-1))
!print *,'MAX threads=',omp_get_max_threads()


      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,isize,ierr)
      if(irank==0)then
        allocate(time1all(niter,isize))
        allocate(time2all(niter,isize))
        allocate(time3all(niter,isize))
      else
        allocate(time1all(niter,1))
        allocate(time2all(niter,1))
        allocate(time3all(niter,1))
      endif
!$OMP parallel
      print *,'thread',omp_get_thread_num(),' bound to cpu',  bind_thread_to_cpu()
!$OMP end parallel

!$OMP parallel
      print *,'before from thread', 1+omp_get_thread_num(),' of',omp_get_max_threads(),' bound to',report_cpu_binding()
!$OMP end parallel

      allocate(r(nis,njs,nks))
      allocate(a(nis,njs,nks))
      allocate(b(nis,njs,nks))

      do iterext=1,3
      r = 2.1+iterext
      a = 0.5+iterext
      b = 2.3+iterext

      call computefornothing_a(r,a,b,nis,njs,nks)
      
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      tstart=mpi_wtime()
      do iter=1,niter
        timea(iter) = mpi_wtime()
        if(iterext==1) then
          call computefornothing_a(r,a,b,nis,njs,nks)
        endif
        if(iterext==2) then
          call computefornothing_b(r,a,b,nis,njs,nks)
        endif
        if(iterext==3) then
          call computefornothing_c(r,a,b,nis,njs,nks)
        endif
        timeb(iter)=mpi_wtime()
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        timec(iter)=mpi_wtime()
!print *,'iter=',iter
      enddo
      tend=mpi_wtime()
      
      timed = timeb-timea  ! compute time
      timee = timec-timeb  ! barrier time
      timef = timec-timea  ! iteration time
      call mpi_reduce(timed,time1,niter,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call mpi_gather(timed,niter,MPI_DOUBLE_PRECISION,time1all,niter,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(timee,time2,niter,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call mpi_gather(timee,niter,MPI_DOUBLE_PRECISION,time2all,niter,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(timef,time3,niter,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call mpi_gather(timef,niter,MPI_DOUBLE_PRECISION,time3all,niter,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(irank==0)then
!        do iter=1,niter
!          print 99,'iteration(',iterext,'):',iter,' time=',time1(iter),time2(iter),time3(iter)
99        format(A,I1,A,i4,a,4F12.6)
!        enddo
        print 100,'Diag(',1000*iterext,'): MAX,MIN,RATIO time for computation   =',maxval(time1(1:niter)),minval(time1(1:niter)),maxval(time1(1:niter))/minval(time1(1:niter))
        print 100,'Diag(',1000*iterext,'): MAX,MIN,RATIO time for an iteration  =',maxval(time3(1:niter)),minval(time3(1:niter)),maxval(time3(1:niter))/minval(time3(1:niter))
        print 100,'Diag(',1000*iterext,'): MAX,MIN,RATIO longest wait on barrier=',maxval(time2(1:niter)),minval(time2(1:niter)),maxval(time2(1:niter))/minval(time2(1:niter))
        print 100,'Diag(',1000*iterext,'): END to END time                      =',tend-tstart
100     format(A,I4,A,5F12.6)
        do iter=1,niter
          tmax=maxval(time1all(iter,:))   ! max compute for all MPI tiles (this iteration)
          tmin=minval(time1all(iter,:))   ! min compute for all MPI tiles (this iteration)
          tsum=sum(time1all(iter,:))
          tsum2=sum(time1all(iter,:)**2)
          tmean=tsum/isize
          timeb(iter)=tmean               ! average compute (this iteration)
          timec(iter)=sqrt(  sum(  ( time1all(iter,:)-tmean )**2  )  )   ! RMS 
          timea(iter)=timea(iter)/minval(time1all(iter,:))  ! compute jitter ratio for an iteration
          print 100,'DIAG(',iter+1000*iterext,')CPU jitter(min/max/avg/rms%):',tmin,tmax,timeb(iter),100.0*timec(iter)/timeb(iter)
        enddo
        tmax=maxval(time1all)
        tmin=minval(time1all)
        tsum=sum(time1all)
        tsum2=sum(time1all**2)
        tmean=tsum/(niter*isize)
        tdev=sqrt(  sum(  ( time1all-tmean )**2  )  )   ! RMS
        print 100,'DIAG(',1000*iterext,')TOT jitter(min/max/avg/rms%):',tmin,tmax,tmean,100*tdev/tmean
      endif

!      deallocate(r,a,b)

      enddo ! iterext

      call mpi_finalize(ierr)
      end program stoopid
!
      subroutine computefornothing_a (r,a,b,nis,njs,nks)
      implicit none
      integer nis,njs,nks
      real*8 r(nis,njs,nks)
      real*8 a(nis,njs,nks)
      real*8 b(nis,njs,nks)

      integer n,i,j,k

!print *,'nis,njs,nks=',nis,njs,nks
!$omp parallel
!$omp do
      do k=1,nks
      do j=1,njs
      do i=1,nis
         r(i,j,k)=r(i,j,k)+4*i*j*k+(a(i,j,k)-1)*b(i,j,k)
      end do
      end do
      end do
!$omp enddo
!$omp end parallel

      return
      end
!
      subroutine computefornothing_b (r,a,b,nis,njs,nks)
      implicit none
      integer nis,njs,nks
      real*8 r(nis,njs,nks)
      real*8 a(nis,njs,nks)
      real*8 b(nis,njs,nks)

      integer n,i,j,k

!print *,'nis,njs,nks=',nis,njs,nks
!$omp parallel
!$omp do
      do k=1,nks
      do j=1,njs
      do i=1,nis
         r(i,j,k)=r(i,j,k)+4*i*j*k+(r(i,j,k)-1)
      end do
      end do
      end do
!$omp enddo
!$omp end parallel

      return
      end
!
      subroutine computefornothing_c (r,a,b,nis,njs,nks)
      implicit none
      integer nis,njs,nks
      real*8 r(nis,njs,nks),rr
      real*8 a(nis,njs,nks)
      real*8 b(nis,njs,nks)

      integer n,i,j,k

!print *,'nis,njs,nks=',nis,njs,nks
!$omp parallel private(rr)
!$omp do
      do k=1,nks
      rr=1.33
      do j=1,njs
      do i=1,nis
         rr=min(1000000.D0,rr+4.0*i*j*k+(rr-1)*rr)
      end do
      end do
      r(1,1,k)=rr
      end do
!$omp enddo
!$omp end parallel

      return
      end
