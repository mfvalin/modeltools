#if defined(NEVER_TO_BE_TRUE)
#endif
// macro used to instrument code with counters
#define INSTRUMENT(A) ;

#include <randomgeneric.h>

/*==========================================================================
 * SHR3 generators, naming consistent with the MWC8222 code
 *==========================================================================*/
/*------------------------ start of SHR3 routines --------------------------*/

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
  unsigned long jsr;
} shr3_state;                  // SHR3 generator stream control structure

#if defined(SELF_TEST)
static shr3_state shr3 = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 123456789 } ;
#endif

// #define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)

void RanSetSeed_SHR3(void *SHR3, unsigned int *piSeed, int cSeed)  // !InTc!
{
  shr3_state *state = (shr3_state *) SHR3 ;
  if(piSeed == NULL || cSeed == 0) return ; // null call, nothing to do
  state->jsr = (unsigned) *piSeed ;
}

unsigned int IRan_SHR3(void *SHR3)	  // !InTc!	/* returns a random unsigned integer */
{
  shr3_state *state = (shr3_state *) SHR3 ;
  unsigned long jz;
  unsigned long jsr;

  jsr=state->jsr ;
  jz = jsr;
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  state->jsr=jsr;
  return (jz+jsr) ;
}

double DRan_SHR3(void *SHR3)	  // !InTc!	/* returns a random double (0.0 , 1.0) */
{
  shr3_state *state = (shr3_state *) SHR3 ;
  unsigned long jz;
  unsigned long jsr;

  jsr=state->jsr ;
  jz = jsr;
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  state->jsr=jsr;
  return CVTDBL_32new(jz+jsr);   // convert from 32 bit int to (0.0 , 1.0)
}

double DRanS_SHR3(void *SHR3)	  // !InTc!	/* returns a random double (-1.0 , 1.0) */
{
  shr3_state *state = (shr3_state *) SHR3 ;
  unsigned long jz;
  unsigned long jsr;

  jsr=state->jsr ;
  jz = jsr;
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  state->jsr=jsr;
  return CVTDBLS_32new(jz+jsr);   // convert from 32 bit int to (-1.0 , 1.0)
}

void VecIRan_SHR3(void *SHR3, unsigned int *ranbuf, int n)  // !InTc!
{
  int i;
  shr3_state *state = (shr3_state *) SHR3 ;
  unsigned long jz;
  unsigned long jsr;

  jsr = state->jsr ;
  for(i=0 ; i<n ; i++) {
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    ranbuf[i] = (jz+jsr);
  };
  state->jsr=jsr;
}

void VecDRan_SHR3(void *SHR3, double *ranbuf, int n)  // !InTc!
{
  int i;
  shr3_state *state = (shr3_state *) SHR3 ;
  unsigned long jz;
  unsigned long jsr;

  jsr = state->jsr ;
  for(i=0 ; i<n ; i++){
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    ranbuf[i] = CVTDBL_32new(jz+jsr);   // convert from 32 bit int to (0.0 , 1.0)
  }
  state->jsr=jsr;
}

void VecDRanS_SHR3(void *SHR3, double *ranbuf, int n)  // !InTc!
{
  int i;
  shr3_state *state = (shr3_state *) SHR3 ;
  unsigned long jz;
  unsigned long jsr;

  jsr = state->jsr ;
  for(i=0 ; i<n ; i++){
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    ranbuf[i] = CVTDBLS_32new(jz+jsr);   // convert from 32 bit int to (0.0 , 1.0)
  }
  state->jsr=jsr;
}

/*------------------------- end of SHR3 routines ---------------------------*/

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
  dmax = CVTDBL_32new(maxpos) ;
  dmin = CVTDBL_32new(maxneg) ;
  printf("maxpos, maxneg transformed with CVTDBL_32new  : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  dmax = CVTDBLS_32new(maxpos) ;
  dmin = CVTDBLS_32new(maxneg) ;
  printf("maxpos, maxneg transformed with CVTDBLS_32new : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);

  for( i=0 ; i < 1000000 ; i++) lr = IRan_SHR3(&shr3);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRan_SHR3(&shr3);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random SHR3 integer value = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) rval = DRan_SHR3(&shr3);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRan_SHR3(&shr3);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random SHR3 double value = %6.3f \n",t1-t0);

  for( i=0 ; i < 10 ; i++) VecIRan_SHR3(&shr3, ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_SHR3(&shr3, ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if((int)ranbuf[i] > 0) pos++ ; else neg++ ;
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & 1) pos++ ; else neg++ ;
  postot = 0 ; negtot = 0;
  for (j=0 ; j<100 ; j++) {
    VecIRan_SHR3(&shr3, ranbuf, 1000000) ;
    mask = 1 ;
    while (mask) {
      pos = 0 ; neg = 0 ; 
      for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & mask) pos++ ; else neg++  ; 
      postot += pos ; negtot += neg ;
      mask <<= 1 ;//  printf("%5d ",pos-neg) ;
    }
  }
//   printf("%d\n",postot-negtot);
  printf("time for 1E+3 x 1E+6 random SHR3 integer values = %6.3f , pos - neg = %d\n",t1-t0,postot-negtot);

  for( i=0 ; i < 10 ; i++) VecDRan_SHR3(&shr3, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_SHR3(&shr3, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random SHR3 double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecIRan_SHR3(&shr3, &ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_SHR3(&shr3, &ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random SHR3 integer values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecDRan_SHR3(&shr3, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_SHR3(&shr3, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random SHR3 double values = %6.3f \n",t1-t0);

  MPI_Finalize();
  return(0);
}
#endif
