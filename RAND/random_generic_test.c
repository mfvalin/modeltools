#define NSAMPLES 1000000000

#include <randomfunctions.h>
#include <mpi.h>

int main(int argc, char **argv){
  unsigned int lr;
  int i, j;
  double t0, t1, rval;
  double MPI_Wtime() ;
  unsigned int ranbuf[1200000];
  double ranbuf2[1200000];
  int pos, neg, mask, postot, negtot ;
#if defined(CYCLIC_TEST)
  int ran;
  long long count, counts ;
#endif
  double dmax, dmin, avg;
  unsigned long long *idmax, *idmin ;
  unsigned int maxpos, maxneg;
  generic_state *gen = NULL;
  int gaussdist[10];
  int index;
  int mySeed;

  MPI_Init(&argc,&argv);
  for(i=0 ; i<1200000 ; i++) ranbuf[i] = 0;
  for(i=0 ; i<1200000 ; i++) ranbuf2[i] = 0.0;
  maxpos = 0x7FFFFFFF ;
  maxneg = 0x80000000 ;
  idmax = (unsigned long long *)&dmax;
  idmin = (unsigned long long *)&dmin;
  dmax = CVTDBL_32(maxpos) ;
  dmin = CVTDBL_32(maxneg) ;
  printf("maxpos, maxneg transformed with CVTDBL_32  : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  dmax = CVTDBLS_32(maxpos) ;
  dmin = CVTDBLS_32(maxneg) ;
  printf("maxpos, maxneg transformed with CVTDBLS_32 : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);

#if defined(TEST_R250)
  mySeed = 123456;
  gen = (generic_state *)  Ran_R250_new_stream(NULL, &mySeed, 1);
#endif
#if defined(TEST_MT19937)
  gen = (generic_state *)  Ran_MT19937_new_stream(NULL, &mySeed, 1);
#endif
#if defined(CYCLIC_TEST)
  ran = IRan_generic_stream(gen);
  counts = 0;
  fprintf(stdout,"ran target = %d\n",ran);
  for(j=0 ; j<1000 ; j++){
    count = 0;
    while(ran != IRan_generic_stream(gen)) count ++ ;
    counts += count;
    fprintf(stdout,"%5d-repeat after %12Ld , running average = %12Ld\n",j+1,count,counts / (j+1)) ;
    fflush(stdout);
  }
exit(0);
#endif

  for( i=0 ; i < 1000000 ; i++) lr = IRan_generic_stream(gen);  // prime the pump

  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < NSAMPLES ; i++) lr = IRan_generic_stream(gen);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random generic integer value = %6.3f , last = %d\n",t1-t0, lr);

  for( i=0 ; i < 1000000 ; i++) rval = DRan_generic_stream(gen);
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < NSAMPLES ; i++) rval = DRan_generic_stream(gen);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random generic  double value = %6.3f , last = %g\n",t1-t0,rval);

  for( i=0 ; i < 1000000 ; i++) rval = DRanS_generic_stream(gen);
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < NSAMPLES ; i++) rval = DRanS_generic_stream(gen);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random generic  signed value = %6.3f , last = %g\n",t1-t0,rval);

  for( i=0 ; i < 10 ; i++) VecIRan_generic_stream(gen,ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_generic_stream(gen, ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random generic integer values = %6.3f , last = %d ",t1-t0,ranbuf[1000000-1]);

  postot = 0 ; negtot = 0;
  RanSetSeed_generic_stream(gen,NULL,0);
  for (j=0 ; j<100 ; j++) {
    VecIRan_generic_stream(gen,ranbuf, 1000000) ;
    mask = 1 ;
    while (mask) {
      pos = 0 ; neg = 0 ; 
      for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & mask) pos++ ; else neg++  ; 
      postot += pos ; negtot += neg ;
      mask <<= 1 ;//  printf("%5d ",pos-neg) ;
    }
  }
  printf(", pos - neg = %d\n",postot-negtot);

  for( i=0 ; i < 10 ; i++) VecDRan_generic_stream(gen, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_generic_stream(gen, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random generic  double values = %6.3f,  last = %g\n",t1-t0,ranbuf2[1000000-1]);

  for( i=0 ; i < 10 ; i++) VecDRanS_generic_stream(gen, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRanS_generic_stream(gen, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random generic  signed values = %6.3f,  last = %g\n",t1-t0,ranbuf2[1000000-1]);

  for( i=0 ; i < 1000 ; i++) VecIRan_generic_stream(gen, &ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_generic_stream(gen, &ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random generic integer values = %6.3f , last = %d\n",t1-t0,ranbuf[i+1000-2]);

  for( i=0 ; i < 1000 ; i++) VecDRan_generic_stream(gen, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_generic_stream(gen, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random generic  double values = %6.3f,  last = %g\n",t1-t0,ranbuf2[i+1000-2]);

  for( i=0 ; i < 1000 ; i++) VecDRanS_generic_stream(gen, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRanS_generic_stream(gen, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random generic  signed values = %6.3f,  last = %g\n",t1-t0,ranbuf2[i+1000-2]);

  MPI_Finalize();
  return(0);
}
