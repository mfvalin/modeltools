#if defined(NEVER_TO_BE_TRUE)
#endif

// macro used to instrument code with counters
#define INSTRUMENT(A) ;

#include <randomgeneric.h>

static uint64_t FACTOR  = 809430660 ;
#define CARRY  362436
#define NSTATE  256

typedef struct{
  REFILLBUFFUN  refill;
  RANSETSEEDFUN seed;
  IRANFUN       iran;
  DRANFUN       dran;
  DRANSFUN      drans;
  IVECRANFUN    vec_iran;
  DVECRANFUN    vec_dran;
  DVECSRANFUN   vec_drans;
  unsigned int *gauss;
  int ngauss;
  unsigned int i ;
  unsigned int c ;
  unsigned int state[NSTATE];
} mwc_state ;                  // MWC8222 generator stream control structure

static mwc_state mwc = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NSTATE - 1, CARRY };

void RanSetSeed_MWC8222(void *MWC8222, unsigned int *Seed, int nSeed)  // !InTc!
{
  unsigned int seed ;
  int i;

  mwc.i = NSTATE - 1;
  mwc.c = CARRY;

  if (nSeed == NSTATE)
  {
    for (i = 0; i < NSTATE; ++i) mwc.state[i] = (unsigned int)Seed[i];
  }else{
    seed = (Seed && (nSeed > 0)) ? Seed[0] : 0 ;
    for (i = 0 ; i < NSTATE ; ) {   /* see Knuth p.106, Table 1(16) and Numerical Recipes p.284 (ranqd1)*/
      seed = 1664525 * seed + 1013904223 ;
      if (seed <= 0) continue ;
      mwc.state[i++] = seed ;
    }
  }
}

unsigned int IRan_MWC8222(void *MWC8222)  // !InTc!
{
  uint64_t t;

  mwc.i = (mwc.i + 1) & (NSTATE - 1);
  t = FACTOR * mwc.state[mwc.i] + mwc.c;
  mwc.c = (unsigned int)(t >> 32);
  mwc.state[mwc.i] = (unsigned int)t;
  return (unsigned int)t;
}

double DRan_MWC8222(void *MWC8222)         // !InTc!  /* returns a random double (0.0 , 1.0) */
{
  uint64_t t;

  mwc.i = (mwc.i + 1) & (NSTATE - 1);
  t = FACTOR * mwc.state[mwc.i] + mwc.c;
  mwc.c = (unsigned int)(t >> 32);
  mwc.state[mwc.i] = (unsigned int)t;
  return CVTDBL_32(t);   // convert from 32 bit int to (0.0 , 1.0)
}

double DRanS_MWC8222(void *MWC8222)        // !InTc!  /* returns a random double (-1.0 , 1.0) */
{
  uint64_t t;

  mwc.i = (mwc.i + 1) & (NSTATE - 1);
  t = FACTOR * mwc.state[mwc.i] + mwc.c;
  mwc.c = (unsigned int)(t >> 32);
  mwc.state[mwc.i] = (unsigned int)t;
  return CVTDBLS_32(t);   // convert from 32 bit int to (-1.0 , 1.0)
}

void VecIRan_MWC8222(void *MWC8222, unsigned int *ran, int npts)  // !InTc!
{
  uint64_t t;
  unsigned int carry = mwc.c, index = mwc.i;
  
  for (; npts > 0; --npts, ++ran)
  {
    index = (index + 1) & (NSTATE - 1);
    t = FACTOR * mwc.state[index] + carry;
    *ran = mwc.state[index] = (unsigned int)t;
    carry = (unsigned int)(t >> 32);
  }
  mwc.c = carry;
  mwc.i = index;
}

void VecDRan_MWC8222(void *MWC8222, double *ran, int npts)  // !InTc!
{
  uint64_t t;
  unsigned int carry = mwc.c, index = mwc.i;
  
  for (; npts > 0; --npts, ++ran)
  {
    index = (index + 1) & (NSTATE - 1);
    t = FACTOR * mwc.state[index] + carry;
    mwc.state[index] = (unsigned int)t;
    *ran = CVTDBL_32(t);   // convert from 32 bit int to (0.0 , 1.0)
    carry = (unsigned int)(t >> 32);
  }
  mwc.c = carry;
  mwc.i = index;
}

void VecDRanS_MWC8222(void *MWC8222, double *ran, int npts)  // !InTc!
{
  uint64_t t;
  unsigned int carry = mwc.c, index = mwc.i;
  
  for (; npts > 0; --npts, ++ran)
  {
    index = (index + 1) & (NSTATE - 1);
    t = FACTOR * mwc.state[index] + carry;
    mwc.state[index] = (unsigned int)t;
    *ran = CVTDBLS_32(t);   // convert from 32 bit int to (-1.0 , 1.0)
    carry = (unsigned int)(t >> 32);
  }
  mwc.c = carry;
  mwc.i = index;
}
/*----------------------- END George Marsaglia MWC -------------------------*/

/*------------------------- end of MWC8222 routines ---------------------------*/

#if defined(SELF_TEST)
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
  int gaussdist[10];
  int index;

  MPI_Init(&argc,&argv);
  for(i=0 ; i<1200000 ; i++) ranbuf[i] = 0;
  for(i=0 ; i<1200000 ; i++) ranbuf2[i] = 0.0;
  maxpos = 0x7FFFFFFF ;
  maxneg = 0x80000000 ;
  idmax = (unsigned long long *)&dmax;
  idmin = (unsigned long long *)&dmin;
  dmax = CVTDBL_32(maxpos) ;
  dmin = CVTDBL_32(maxneg) ;
  printf("maxpos, maxneg transformed with CVTDBL_32new  : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  dmax = CVTDBLS_32(maxpos) ;
  dmin = CVTDBLS_32(maxneg) ;
  printf("maxpos, maxneg transformed with CVTDBLS_32new : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);

  for( i=0 ; i < 1000000 ; i++) lr = IRan_MWC8222(&mwc);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRan_MWC8222(&mwc);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random MWC8222 integer value = %6.3f\n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) rval = DRan_MWC8222(&mwc);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRan_MWC8222(&mwc);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random MWC8222 double value = %6.3f \n",t1-t0);

  for( i=0 ; i < 10 ; i++) VecIRan_MWC8222(&mwc, ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_MWC8222(&mwc, ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if((int)ranbuf[i] > 0) pos++ ; else neg++ ;
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & 1024) pos++ ; else neg++ ;
  postot = 0 ; negtot = 0;
  for (j=0 ; j<100 ; j++) {
    VecIRan_MWC8222(&mwc, ranbuf, 1000000) ;
    mask = 1 ;
    while (mask) {
      pos = 0 ; neg = 0 ; 
      for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & mask) pos++ ; else neg++  ; 
      postot += pos ; negtot += neg ;
      mask <<= 1 ;//  printf("%5d ",pos-neg) ;
    }
  }
//   printf("%d\n",postot-negtot);
  printf("time for 1E+3 x 1E+6 random MWC8222 integer values = %6.3f , pos - neg = %d\n",t1-t0,postot-negtot);

  for( i=0 ; i < 10 ; i++) VecDRan_MWC8222(&mwc, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_MWC8222(&mwc, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random MWC8222 double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecIRan_MWC8222(&mwc, &ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_MWC8222(&mwc, &ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random MWC8222 integer values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecDRan_MWC8222(&mwc, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_MWC8222(&mwc, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random MWC8222 double values = %6.3f \n",t1-t0);

  MPI_Finalize();
  return(0);
}
#endif
