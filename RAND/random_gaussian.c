/*==========================================================================
 *               ORIGINAL Copyright before modifications                    */
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
 *==========================================================================*/

/*==========================================================================
Amendment 2017-02-22:
Jurgen A Doornik gives permission to any branch of the government of Canada 
to use his random number and ziggurat code and any derivatives 
for any purpose, whether commercial or not.
 *==========================================================================*/


// gaussian generator using the "ziggurat" method (uses one of the linear integer generators)   // !InTc!
// 128 or 256 boxes (use internal buffer for 320 32 bit integers)   // !InTc!

// Fortran interface, ready to be extracted
#if defined(NEVER_TO_BE_TRUE)

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
INSTRUMENT(unsigned int funcalls = 0;)
INSTRUMENT(unsigned int funquick = 0;)
INSTRUMENT(unsigned int funtails = 0;)
INSTRUMENT(unsigned int funwedge = 0;)
INSTRUMENT(unsigned int funloops = 0;)
INSTRUMENT(unsigned int funused  = 0;)

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
