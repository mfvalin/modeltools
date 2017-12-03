#include <randomgeneric.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv){
  unsigned int lr, lr2;
  int i, j, k, l;
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
  generic_state *shr3;
  unsigned int piSeed = 123456789;
  int irep = 0;
  int count[32];

  MPI_Init(&argc,&argv);

  for(irep=0 ; irep <3 ; irep++){
  if(irep == 2) {
    printf("\n\n\n ========= SHR3 test =========\n");
    piSeed = 123456789;
    shr3 = Ran_SHR3_new_stream(NULL, &piSeed, 1);
  }
  if(irep == 1) {
    printf("\n\n\n ========= R250 test =========\n");
    piSeed = 123456789;
    shr3 = Ran_R250_new_stream(NULL, &piSeed, 1);
  }
  if(irep == 0) {
    printf("\n\n\n ========= MT19937 test =========\n");
    piSeed = 1;
    shr3 = Ran_MT19937_new_stream(NULL, &piSeed, 1);
  }

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

  for( i=0 ; i < 1000000 ; i++) lr = IRan_generic_stream(shr3);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRan_generic_stream(shr3);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random integer value = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) rval = DRan_generic_stream(shr3);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRan_generic_stream(shr3);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random double value = %6.3f \n",t1-t0);

  for( i=0 ; i < 10 ; i++) VecIRan_generic_stream(shr3, ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_generic_stream(shr3, ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if((int)ranbuf[i] > 0) pos++ ; else neg++ ;
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & 1) pos++ ; else neg++ ;
  postot = 0 ; negtot = 0;
  for(i=0 ; i<32 ; i++) count[i] = 0;
  for (j=0 ; j<1000 ; j++) {
    VecIRan_generic_stream(shr3, ranbuf, 1000000) ;
    for(k=0 ; k<1000000 ; k++) {
      l = ranbuf[k];
      for(i=0 ; i<32 ; i++) { count[i] += (l & 1) ; l >>=1 ; }
    }
//     mask = 1 ;
//     while (mask) {
//       pos = 0 ; neg = 0 ; 
//       for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & mask) pos++ ; else neg++  ; 
//       postot += pos ; negtot += neg ;
//       mask <<= 1 ;//  printf("%5d ",pos-neg) ;
//     }
  }
//   printf("%d\n",postot-negtot);
  for(i=0 ; i<32 ; i++) count[i] -= (j*k/2);
  for(i=0 ; i<32 ; i++) postot += count[i] ; 
  rval = j*k ; rval = sqrt(rval); lr = rval*sqrt(32.0); lr2 = rval;
  printf("time for 1E+3 x 1E+6 random integer values = %6.3f \n  ones - zeroes = %d (32 * %d samples) expected = +-%d (+-%d by bit)\n  by bit / 100 : ",
	 t1-t0,postot-negtot,j*k, lr,lr2);
//   printf("time for 1E+3 x 1E+6 random integer values = %6.3f \n",t1-t0);
  for(i=0 ; i<32 ; i++) printf("%d ",count[i]/100); printf("\n");

  for( i=0 ; i < 10 ; i++) VecDRan_generic_stream(shr3, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_generic_stream(shr3, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecIRan_generic_stream(shr3, &ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_generic_stream(shr3, &ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random integer values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecDRan_generic_stream(shr3, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_generic_stream(shr3, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random double values = %6.3f \n",t1-t0);

  } // for irep
  MPI_Finalize();
  return(0);
}
