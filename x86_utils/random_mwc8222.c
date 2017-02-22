#if defined(NEVER_TO_BE_TRUE)
#endif

// macro used to instrument code with counters
#define INSTRUMENT(A) ;

#include <randomgeneric.h>

/*==========================================================================
 *               ORIGINAL Copyright before code modifications              *
/*==========================================================================
 *  This code is Copyright (C) 2005, Jurgen A. Doornik.
 *  Permission to use this code for non-commercial purposes
 *  is hereby given, provided proper reference is made to:
 *  Doornik, J.A. (2005), "An Improved Ziggurat Method to Generate Normal
 *      Random Samples", mimeo, Nuffield College, University of Oxford,
 *      and www.doornik.com/research.
 *   or the published version when available.
 *  This reference is still required when using modified versions of the code.
 *  This notice should be maintained in modified versions of the code.
 *  No warranty is given regarding the correctness of this code.
 * Amendment 2017-02-22
 * Jurgen A Doornik gives permission to any branch of the government of Canada 
 * to use his random number and ziggurat code and any derivatives for
 * any purpose, whether commercial or not.
 *==========================================================================*/
/*==========================================================================
 * code modifications and refactoring done by
 * M.Valin, Feb 2017
 * Recherche en Prevision Numerique
 * Environnement Canada
 *==========================================================================*/
static uint64_t MWC_A  = 809430660 ;
#define MWC_AI 809430660
#define MWC_C  362436
#define MWC_R  256

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
  unsigned int uiState ;
  unsigned int uiCarry ;
  unsigned int auiState[MWC_R];
} mwc_state ;                  // MWC8222 generator stream control structure

static mwc_state mwc = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, MWC_R - 1, MWC_C };

// static unsigned int s_uiStateMWC = MWC_R - 1;
// static unsigned int s_uiCarryMWC = MWC_C;
// static unsigned int s_auiStateMWC[MWC_R];

void RanSetSeed_MWC8222(void *MWC8222, unsigned int *piSeed, int cSeed)  // !InTc!
{
  unsigned int seed ;
  int i;

  mwc.uiState = MWC_R - 1;
  mwc.uiCarry = MWC_C;

  if (cSeed == MWC_R)
  {
    for (i = 0; i < MWC_R; ++i) mwc.auiState[i] = (unsigned int)piSeed[i];
  }else{
    seed = piSeed && (cSeed > 0) ? piSeed[0] : 0 ;
    for (i = 0 ; i < MWC_R ; ) {   /* see Knuth p.106, Table 1(16) and Numerical Recipes p.284 (ranqd1)*/
      seed = 1664525 * seed + 1013904223 ;
      if (seed <= 0) continue ;
      mwc.auiState[i++] = seed ;
    }
  }
}

unsigned int IRan_MWC8222(void *MWC8222)  // !InTc!
{
  uint64_t t;

  mwc.uiState = (mwc.uiState + 1) & (MWC_R - 1);
  t = MWC_A * mwc.auiState[mwc.uiState] + mwc.uiCarry;
  mwc.uiCarry = (unsigned int)(t >> 32);
  mwc.auiState[mwc.uiState] = (unsigned int)t;
  return (unsigned int)t;
}

double DRan_MWC8222(void *MWC8222)         // !InTc!  /* returns a random double (0.0 , 1.0) */
{
  uint64_t t;

  mwc.uiState = (mwc.uiState + 1) & (MWC_R - 1);
  t = MWC_A * mwc.auiState[mwc.uiState] + mwc.uiCarry;
  mwc.uiCarry = (unsigned int)(t >> 32);
  mwc.auiState[mwc.uiState] = (unsigned int)t;
  return RANDBL_32new(t);   // convert from 32 bit int to (0.0 , 1.0)
}

double DRanS_MWC8222(void *MWC8222)        // !InTc!  /* returns a random double (-1.0 , 1.0) */
{
  uint64_t t;

  mwc.uiState = (mwc.uiState + 1) & (MWC_R - 1);
  t = MWC_A * mwc.auiState[mwc.uiState] + mwc.uiCarry;
  mwc.uiCarry = (unsigned int)(t >> 32);
  mwc.auiState[mwc.uiState] = (unsigned int)t;
  return RANDBLS_32new(t);   // convert from 32 bit int to (-1.0 , 1.0)
}

void VecIRan_MWC8222(void *MWC8222, unsigned int *auiRan, int cRan)  // !InTc!
{
  uint64_t t;
  unsigned int carry = mwc.uiCarry, state = mwc.uiState;
  
  for (; cRan > 0; --cRan, ++auiRan)
  {
    state = (state + 1) & (MWC_R - 1);
    t = MWC_A * mwc.auiState[state] + carry;
    *auiRan = mwc.auiState[state] = (unsigned int)t;
    carry = (unsigned int)(t >> 32);
  }
  mwc.uiCarry = carry;
  mwc.uiState = state;
}

void VecDRan_MWC8222(void *MWC8222, double *adRan, int cRan)  // !InTc!
{
  uint64_t t;
  unsigned int carry = mwc.uiCarry, state = mwc.uiState;
  
  for (; cRan > 0; --cRan, ++adRan)
  {
    state = (state + 1) & (MWC_R - 1);
    t = MWC_A * mwc.auiState[state] + carry;
    mwc.auiState[state] = (unsigned int)t;
    *adRan = RANDBL_32new(t);   // convert from 32 bit int to (0.0 , 1.0)
    carry = (unsigned int)(t >> 32);
  }
  mwc.uiCarry = carry;
  mwc.uiState = state;
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
  dmax = RANDBL_32new(maxpos) ;
  dmin = RANDBL_32new(maxneg) ;
  printf("maxpos, maxneg transformed with RANDBL_32new  : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  dmax = RANDBLS_32new(maxpos) ;
  dmin = RANDBLS_32new(maxneg) ;
  printf("maxpos, maxneg transformed with RANDBLS_32new : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);

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
