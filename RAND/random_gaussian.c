/* 
 * Copyright (C) 2017 Recherche en Prevision Numerique
 * 
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

// Fortran interface, ready to be extracted
#if defined(NEVER_TRUE)
! type, bind(C) :: RANDOM_STREAM                                                          !InTf!
!   type(C_PTR) :: p                                                                      !InTf!
! end type                                                                                !InTf!

! double F_DRan_NormalZig_stream(statep *s   )                                            !InTf!
 interface                                                                                !InTf!
   function DRan_Normal_stream(stream) result(ran) bind(C,name='F_DRan_NormalZig_stream') !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRan_Normal_stream                                                        !InTf!
 end interface                                                                            !InTf!

! double F_D64Ran_NormalZig_stream(statep *s   )                                          !InTf!
 interface                                                                                !InTf!
   function D64Ran_Normal_stream(stream) result(ran) bind(C,name='F_D64Ran_NormalZig_stream') !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function D64Ran_Normal_stream                                                      !InTf!
 end interface                                                                            !InTf!

! void F_RanNormalZigSetSeed(statep *s   , int *piSeed, int cSeed)                        !InTf!
 interface                                                                                !InTf!
   subroutine RanNormalZigSetSeed(stream, piSeed, cSeed) bind(C,name='F_RanNormalZigSetSeed') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine RanNormalZigSetSeed                                                     !InTf!
 end interface                                                                            !InTf!

#endif

// to instrument the code with counters, add -DPROFILE
#if defined(PROFILE)
#define INSTRUMENT(A) A
#else
#define INSTRUMENT(A)
#endif
INSTRUMENT(static unsigned int funcalls = 0;)
INSTRUMENT(static unsigned int funquick = 0;)
INSTRUMENT(static unsigned int funtails = 0;)
INSTRUMENT(static unsigned int funwedge = 0;)
INSTRUMENT(static unsigned int funloops = 0;)
INSTRUMENT(static unsigned int funused  = 0;)

#include <randomgeneric.h>
#include <stdint.h>
/*==========================================================================
 * 
 *  slightly trimmed version of the improved Ziggurat algorithm described in
 * 
 *  Doornik, J.A. (2005), 
 *      "An Improved Ziggurat Method to Generate Normal Random Samples",
 *      mimeo, Nuffield College, University of Oxford,
 *      www.doornik.com/research.
 * 
 *  BOXAREA and TAIL values for 128 and 256 boxes from 
 * 
 *  George Marsaglia, Wai Wan Tsang (2000)
 *      "The Ziggurat Method for Generating Random Variables".
 *
 *==========================================================================*/
/*----------------------- START of gaussian generators ---------------------*/

// internal buffer size
#define ZIGBUFSIZE 256
#define ZIGBUFSIZE2 (ZIGBUFSIZE + ZIGBUFSIZE/4)

// internal buffer structure
typedef struct{
  uint32_t ir[ZIGBUFSIZE];
  uint8_t  box[ZIGBUFSIZE];
} zigbuf;

// check that internal buffer contains at least N items. if not, refill it using stream vector integer function
#define CHECKBUF(N) { if (navail-- < N) { (*zig->vec_iran)(stream,(void *)buffer,ZIGBUFSIZE2 ) ; navail = ZIGBUFSIZE-1 ; } }
#define CHECKBUF1   { if (--navail < 0) { (*zig->vec_iran)(stream,(void *)buffer,ZIGBUFSIZE2 ) ; navail = ZIGBUFSIZE-1 ; } }
#define CHECKB64(N) { if (navail-- < N) { (*zig->vec_iran)(stream,(void *)buffer,ZIGBUFSIZE2 ) ; navail = ZIGBUFSIZE2/2-1 ; } }
#define CHECKB641   { if (--navail < 0) { (*zig->vec_iran)(stream,(void *)buffer,ZIGBUFSIZE2 ) ; navail = ZIGBUFSIZE2/2-1 ; } }

// initialize ziggurat algorithm tables
static void InitZigguratMethodTables(double *redge, double *gauss, int nboxes, double tail, double boxarea){
  int i;
  double p;
  p             = exp(-.5 * tail * tail) ;   // f(x)
  redge[0]      = boxarea / p ;
  gauss[0]      = exp(-.5 * redge[0] * redge[0]);
  redge[1]      = tail ;
  gauss[1]      = exp(-.5 * redge[1] * redge[1]);
  redge[nboxes] = 0.0 ;          // top box has no acceptable rectangle
  gauss[nboxes] = 1.0 ;
  for ( i=2 ; i < nboxes ; i++) {
    redge[i] = sqrt(-2 * log(boxarea / redge[i-1] + p));  // right edge of box i
    p = exp(-.5 * redge[i] * redge[i]);                   // f(x)
    gauss[i] = p;
  }
}

// default version uses 256 boxes (~ 15% faster)
// to generate 128 box version (2 x smaller tables), add -DUSE128 when compiling
#if ! defined(USE128)
#define USE256
#endif

#if defined(USE256)

#define NBOXES1  256
static const double BOXAREA1 = 4.92867323399e-3;
static const double TAIL1    = 3.6541528853610088;
static const double TAILINV1 = 0.2736612373297583;

// ziggurat algorithm tables (256 boxes)
static double redge1[NBOXES1+1] ;  // coordinates of rectangular box right edges (box 0 is bottom box)
static double gauss1[NBOXES1+1] ;  // exp ( -.5 * (redge[i] * redge[i]) )
static int init1 = 0;

// initialize ziggurat algorithm tables (256 boxes)
static void InitZigguratMethodTables256(void)
{
  int i;
  double p;
  InitZigguratMethodTables(redge1,gauss1,NBOXES1,TAIL1,BOXAREA1);
}

// this function MUST be called at least once BEFORE calling DRan_NormalZig_stream
void RanNormalZigSetSeed(void *stream, void *values, int nvalues)  // !InTc!
{
  generic_state *zig = stream ;

  InitZigguratMethodTables256() ;
  if(zig->gauss == NULL) {   // allocate stream buffer for uniform numbers if not already done
    zig->gauss = (unsigned int *) memalign(64,ZIGBUFSIZE2*sizeof(unsigned int));
    zig->ngauss = 0;
  }
  if(zig->gauss != NULL) init1 = 1;
}

// get a gaussian distributed random number (only 32 significant bits in mantissa)
// this function will FAIL if RanNormalZigSetSeed has not been called to initialize the 
// stream and the base tables
double DRan_NormalZig_stream(void *stream)  // !InTc!
{
  generic_state *zig = stream ;
  double g, x, y, f0, f1, f2;
  int navail;
  zigbuf *buffer; 
  unsigned int i;
  INSTRUMENT(int direct ;)

  INSTRUMENT(funcalls++; ; direct = 1;)
  if(init1 == 0) { exit(1) ; return(g = -1.0e+38) ; } // miserable failure, initialization not done properly

  navail = zig->ngauss ;
  buffer = (zigbuf *) zig->gauss ;
  for(;;){
//     CHECKBUF(1);                               // need 1 uniform value
//     g = CVTDBLS_32(buffer->ir[--navail]);      // convert from 32 bit int to (-1.0 , 1.0) on the fly
    CHECKBUF1;                               // need 1 uniform value
    g = CVTDBLS_32(buffer->ir[navail]);      // convert from 32 bit int to (-1.0 , 1.0) on the fly
    i = (buffer->box[navail]) ;
    INSTRUMENT(funused += 2;)
    g = g * redge1[i];
    if (fabs(g) < redge1[i+1] ) {              // first try the rectangular boxes
      INSTRUMENT(funquick += direct;)
      break ;                                  // done
    }
    if (i == 0) {                              // bottom box: sample from the TAIL
      INSTRUMENT(funtails++ ; )
      do {
	CHECKBUF(2) ;                          // need 2 uniform values
	x = CVTDBL_32(buffer->ir[navail]);
	navail--;
	x = log(x) * TAILINV1 ;                //    / TAIL;
	y = CVTDBL_32(buffer->ir[navail]);
	y = log(y);
	INSTRUMENT(funused += 2;)
      } while (x * x + y + y > 0);
      g = (g < 0) ? x - TAIL1 : TAIL1 - x;
      break ;                                  // done
    } else {                                   // is this from the wedges?
      INSTRUMENT(funwedge ++;)
      f0 = exp(0.5 * g * g) ;
      f1 = f0 * gauss1[i] ;     // f1 = exp(-0.5 * (redge1[i] * redge1[i] - x * x) );
      f2 = f0 * gauss1[i+1] ;   // f2 = exp(-0.5 * (redge1[i + 1] * redge1[i + 1] - x * x) );
      INSTRUMENT(funused++;)
//       CHECKBUF(1);                               // need 1 uniform value
//       y = CVTDBL_32(buffer->ir[--navail]);
      CHECKBUF1;                               // need 1 uniform value
      y = CVTDBL_32(buffer->ir[navail]);
      if (f2 + y * (f1 - f2) < 1.0)  break;
      INSTRUMENT(funloops++;)
    }
  }   // for(;;) { }
  zig->ngauss = navail ;
  return(g) ;
}

// get a gaussian distributed random number (52 significant bits in mantissa)
// this function will FAIL if RanNormalZigSetSeed has not been called to initialize the 
// stream and the base tables
double D64Ran_NormalZig_stream(void *stream)  // !InTc!
{
  generic_state *zig = stream ;
  double g, x, y, f0, f1, f2;
  int navail;
  uint64_t *buffer; 
  int i;
  INSTRUMENT(int direct ;)

  INSTRUMENT(funcalls++; ; direct = 1;)
  if(init1 == 0) { exit(1) ; return(g = -1.0e+38) ; } // miserable failure, initialization not done properly

  navail = zig->ngauss ;
  buffer = (uint64_t *) zig->gauss ;
  for(;;){
    CHECKB641;                               // need 1 uniform value
    g = CVTDBLS_64(buffer[navail]);      // convert from 32 bit int to (-1.0 , 1.0) on the fly
    i = buffer[navail] & 0xFF ;
    INSTRUMENT(funused += 1;)
    g = g * redge1[i];
    if (fabs(g) < redge1[i+1] ) {              // first try the rectangular boxes
      INSTRUMENT(funquick += direct;)
      break ;                                  // done
    }
    if (i == 0) {                              // bottom box: sample from the TAIL
      INSTRUMENT(funtails++ ; )
      do {
	CHECKB64(2) ;                          // need 2 uniform values
	x = CVTDBL_64(buffer[navail]);
	navail--;
	x = log(x) * TAILINV1 ;                //    / TAIL;
	y = CVTDBL_64(buffer[navail]);
	y = log(y);
	INSTRUMENT(funused += 2;)
      } while (x * x + y + y > 0);
      g = (g < 0) ? x - TAIL1 : TAIL1 - x;
      break ;                                  // done
    } else {                                   // is this from the wedges?
      INSTRUMENT(funwedge ++;)
      f0 = exp(0.5 * g * g) ;
      f1 = f0 * gauss1[i] ;     // f1 = exp(-0.5 * (redge1[i] * redge1[i] - x * x) );
      f2 = f0 * gauss1[i+1] ;   // f2 = exp(-0.5 * (redge1[i + 1] * redge1[i + 1] - x * x) );
      INSTRUMENT(funused++;)
      CHECKB641;                               // need 1 uniform value
      y = CVTDBL_64(buffer[navail]);
      if (f2 + y * (f1 - f2) < 1.0)  break;
      INSTRUMENT(funloops++;)
    }
  }   // for(;;) { }
  zig->ngauss = navail ;
  return(g) ;
}

#endif

#if defined(USE128)

#define NBOXES0  128
static const double BOXAREA0 = 9.91256303526217e-3;
static const double TAIL0    = 3.442619855899;
static const double TAILINV0 = 0.29047645161474317;

// ziggurat algorithm tables (128 boxes)
static double redge0[NBOXES0+1] ;  // coordinates of rectangular box right edges (box 0 is bottom box)
static double gauss0[NBOXES0+1] ;  // exp ( -.5 * (redge[i] * redge[i]) )
static int init0 = 0;

// initialize ziggurat algorithm tables (128 boxes)
static void InitZigguratMethodTables128(void)
{
  int i;
  double p;
  InitZigguratMethodTables(redge0,gauss0,NBOXES0,TAIL0,BOXAREA0);
}

// this function MUST be called at least once BEFORE calling DRan_NormalZig_stream
void RanNormalZigSetSeed(void *stream, void *values, int nvalues)
{
  generic_state *zig = stream ;
  InitZigguratMethodTables128() ;
  if(zig->gauss == NULL) {   // allocate stream buffer for uniform numbers if not already done
    zig->gauss = (unsigned int *) memalign(64,ZIGBUFSIZE2*sizeof(unsigned int));
    zig->ngauss = 0;
  }
  if(zig->gauss != NULL) init0 = 1;
}

// get a gaussian distributed random number (only 32 significant bits in mantissa)
// this function will FAIL if RanNormalZigSetSeed has not been called to initialize the 
// stream and the base tables
double DRan_NormalZig_stream(void *stream)
{
  generic_state *zig = stream ;
  double g, x, y, f0, f1, f2;
  int navail;
  zigbuf *buffer; 
  unsigned int i;
  INSTRUMENT(int direct ;)

  INSTRUMENT(funcalls++; ; direct = 1;)
  if(init0 == 0) { exit(1) ; return(g = -1.0e+38) ; } // miserable failure, initialization not done properly

  navail = zig->ngauss ;
  buffer = (zigbuf *) zig->gauss ;
  for(;;){
    CHECKBUF1;                               // need 1 uniform value
    g = CVTDBLS_32(buffer->ir[navail]);      // convert from 32 bit int to (-1.0 , 1.0) on the fly
    i = (buffer->box[navail]) & 0x7F ;
    INSTRUMENT(funused += 2;)
    g = g * redge0[i];
    if (fabs(g) < redge0[i+1] ) {              // first try the rectangular boxes
      INSTRUMENT(funquick += direct;)
      break ;                                  // done
    }
    if (i == 0) {                              // bottom box: sample from the TAIL
      INSTRUMENT(funtails++ ; )
      do {
	CHECKBUF(2) ;                          // need 2 uniform values
	x = CVTDBL_32(buffer->ir[navail]);
	navail--;
	x = log(x) * TAILINV0 ;                //    / TAIL;
	y = CVTDBL_32(buffer->ir[navail]);
	y = log(y);
	INSTRUMENT(funused += 2;)
      } while (x * x + y + y > 0);
      g = (g < 0) ? x - TAIL0 : TAIL0 - x;
      break ;                                  // done
    }
    INSTRUMENT(funwedge ++;)
    f0 = exp(0.5 * g * g) ;   // is this from the wedges?
    f1 = f0 * gauss0[i] ;     // f1 = exp(-0.5 * (redge1[i] * redge1[i] - x * x) );
    f2 = f0 * gauss0[i+1] ;   // f2 = exp(-0.5 * (redge1[i + 1] * redge1[i + 1] - x * x) );
    INSTRUMENT(funused++;)
    CHECKBUF1;                               // need 1 uniform value
    y = CVTDBL_32(buffer->ir[navail]);
    if (f2 + y * (f1 - f2) < 1.0)  break;
    INSTRUMENT(funloops++;)
  }   // for(;;) { }
  zig->ngauss = navail ;
  return(g) ;
}

// get a gaussian distributed random number (52 significant bits in mantissa)
// this function will FAIL if RanNormalZigSetSeed has not been called to initialize the 
// stream and the base tables
double D64Ran_NormalZig_stream(void *stream)
{
  generic_state *zig = stream ;
  double g, x, y, f0, f1, f2;
  int navail;
  uint64_t *buffer; 
  unsigned int i;
  INSTRUMENT(int direct ;)

  INSTRUMENT(funcalls++; ; direct = 1;)
  if(init0 == 0) { exit(1) ; return(g = -1.0e+38) ; } // miserable failure, initialization not done properly

  navail = zig->ngauss ;
  buffer = (uint64_t *) zig->gauss ;
  for(;;){
    CHECKB641;                               // need 1 uniform value
    g = CVTDBLS_64(buffer[navail]);      // convert from 32 bit int to (-1.0 , 1.0) on the fly
    i = buffer[navail] & 0x7F ;
    INSTRUMENT(funused += 1;)
    g = g * redge0[i];
    if (fabs(g) < redge0[i+1] ) {              // first try the rectangular boxes
      INSTRUMENT(funquick += direct;)
      break ;                                  // done
    }
    if (i == 0) {                              // bottom box: sample from the TAIL
      INSTRUMENT(funtails++ ; )
      do {
	CHECKB64(2) ;                          // need 2 uniform values
	x = CVTDBL_64(buffer[navail]);
	navail--;
	x = log(x) * TAILINV0 ;                //    / TAIL;
	y = CVTDBL_64(buffer[navail]);
	y = log(y);
	INSTRUMENT(funused += 2;)
      } while (x * x + y + y > 0);
      g = (g < 0) ? x - TAIL0 : TAIL0 - x;
      break ;                                  // done
    } else {                                   // is this from the wedges?
      INSTRUMENT(funwedge ++;)
      f0 = exp(0.5 * g * g) ;
      f1 = f0 * gauss0[i] ;     // f1 = exp(-0.5 * (redge1[i] * redge1[i] - x * x) );
      f2 = f0 * gauss0[i+1] ;   // f2 = exp(-0.5 * (redge1[i + 1] * redge1[i + 1] - x * x) );
      INSTRUMENT(funused++;)
      CHECKB641;                               // need 1 uniform value
      y = CVTDBL_64(buffer[navail]);
      if (f2 + y * (f1 - f2) < 1.0)  break;
      INSTRUMENT(funloops++;)
    }
  }   // for(;;) { }
  zig->ngauss = navail ;
  return(g) ;
}

#endif
// get a gaussian distributed random number (full 52 bit mantissa) (deferred inplementation)
double D64RanNormalFun(void *stream)  // !InTc!
{
  return(D64Ran_NormalZig_stream(stream));
}

// Fortran entry point for initialization
void F_RanNormalZigSetSeed(void **stream, void *values, int nvalues){
  RanNormalZigSetSeed(*stream, values, nvalues);
}

// Fortran entry point using the Fortran derived type, passed by reference
double F_DRan_NormalZig_stream(void **stream){
  return(DRan_NormalZig_stream(*stream));
}

// Fortran entry point using the Fortran derived type, passed by reference
double F_D64Ran_NormalZig_stream(void **stream){
  return(D64Ran_NormalZig_stream(*stream));
}

/*------------------------ END of gaussian generators ----------------------*/
#if defined(SELF_TEST)
#if defined(FULL_TEST)
#include <mpi.h>
#endif
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
  int biggaussdist[2001];
  int index;
  generic_state *stream;
  int myseed = 123456;
  double x, p1, p2, prob, ptot;
  int est[2001];
  union{
    long l;
    double d;
  }v;
// check constants
  v.d = INVM31 ; printf("%5i %16.16lx %24.20g\n",31,v.l,v.d);
  v.d = INVM32 ; printf("%5i %16.16lx %24.20g\n",32,v.l,v.d);
  v.d = INVM33 ; printf("%5i %16.16lx %24.20g\n",33,v.l,v.d);
  v.d = INVM63 ; printf("%5i %16.16lx %24.20g\n",63,v.l,v.d);
  v.d = INVM64 ; printf("%5i %16.16lx %24.20g\n",64,v.l,v.d);
  v.d = INVM65 ; printf("%5i %16.16lx %24.20g\n",65,v.l,v.d);
//   void  RanNormalZigSetSeed(void *stream, int *piSeed, int cSeed)  ;
  void *Ran_R250_new_stream(void *clone_in, int *piSeed, int cSeed)   ;
#if defined(FULL_TEST)
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

  stream = Ran_R250_new_stream(NULL, &myseed , 1);
  RanNormalZigSetSeed(stream, &myseed, 1);
//   RanNormalZigSetSeed128(stream, &myseed, 1);
//   RanNormalZigSetSeed256(stream, &myseed, 1);

  dmin = 0.0 ; dmax = 0.0;
  for( i=0 ; i < 10 ; i++) gaussdist[i] = 0;
  for( i=0 ; i < 2001 ; i++) biggaussdist[i] = 0;
  for(j=0; j<10 ; j++) ;
  for( i=0 ; i < 1000000000 ; i++) {
#if defined(TEST64)
    rval = D64Ran_NormalZig_stream(stream);      // use C entry point
#else
    rval = DRan_NormalZig_stream(stream);      // use C entry point
#endif
    avg = avg + rval ;
    dmin = (dmin < rval) ? dmin : rval ;
    dmax = (dmax > rval) ? dmax : rval ;
    if(rval > 10.0) rval = 10.0;
    if(rval < -10.0) rval = -10.0;
    index = 1001 + rval * 100;
    biggaussdist[index] ++;
  }
  printf("for %d samples, min = %6.3f, max = %6.3f, avg = %10.7f\n",i,dmin,dmax,avg/i);
  ptot = 0.0;
  for( i=0 ; i < 2000 ; i++){
    x = (i-1001) / 100.0 ;
    p1 = exp(-0.5 * x * x)  * .01 ;
    x += .01;
    p2 =  exp(-0.5 * x * x) * .01 ;
    prob = .5 * (p1 + p2);
    ptot += prob ;
    p1 = 1000000000 * prob / 2.5066283;  // sqrt(2 * pi)
    est[i] = p1 * 1;
  }
  printf("ptot = %g, expecting 2.5066283\n",ptot);
//   for( i=0 ; i < 10 ; i++) printf("%9d ",gaussdist[i]);
//   for( i=0 ; i < 2000 ; i++) printf("%9d %9d\n",i,biggaussdist[i]);
  printf("%9s %9s %10s\n","slot","population","deviation(ppm) (center at slot 1000)");
  for( i=991 ; i <= 1010 ; i++) printf("%9d %9d %10.0f\n",i,biggaussdist[i],1000000.0*(biggaussdist[i]-est[i])*1.0/est[i]);
  printf("\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) 
#if defined(TEST64)
    rval = F_D64Ran_NormalZig_stream((void **) &stream);  // tiem Fortran entry point (costlier)
#else
    rval = F_DRan_NormalZig_stream((void **) &stream);  // tiem Fortran entry point (costlier)
#endif
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random DRan_NormalZig_stream/R250 double value = %6.3f \n",t1-t0);  // DRan_NormalZig_stream256

  t1 = 0 ; t0 = 1 ; 
  INSTRUMENT(t1 = funquick ; t0 = funcalls+1 ; )
  INSTRUMENT(printf("quick calls in gaussian generator = %7.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = funtails ;)
  INSTRUMENT(printf("tail calls in gaussian generator  = %7.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = funwedge ;)
  INSTRUMENT(printf("wedge calls in gaussian generator = %7.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = funloops ;)
  INSTRUMENT(printf("extra loops in gaussian generator = %7.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = funused ;)
  INSTRUMENT(printf("uniform random values used        = %7.3f%\n",t1 / t0 * 100.0);)
  MPI_Finalize();
#endif
  return(0);
}
#endif
