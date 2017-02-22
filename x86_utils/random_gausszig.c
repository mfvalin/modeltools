#if defined(COMPILING_FORTRAN)
subroutine fortran_test   ! test or Fortran interfaces really
  use ISO_C_BINDING
  implicit none

  type, bind(C) :: RANDOM_STREAM                                                          !InTf!
    type(C_PTR) :: p                                                                      !InTf!
  end type                                                                                !InTf!

! double F_DRanNormalZigVec(statep *s   )                                                 !InTf!
 interface                                                                                !InTf!
   function DRanNormalZigVec(stream) result(ran) bind(C,name='F_DRanNormalZigVec')          !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRanNormalZigVec                                                          !InTf!
 end interface                                                                            !InTf!

! void F_RanNormalSetSeedZig(statep *s   , int *piSeed, int cSeed)                        !InTf!
 interface                                                                                !InTf!
   subroutine RanNormalSetSeedZig(stream, piSeed, cSeed) bind(C,name='F_RanNormalSetSeedZig') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine RanNormalSetSeedZig                                                     !InTf!
 end interface                                                                            !InTf!

! void F_RanNormalSetSeedZigVec(statep *s   , int *piSeed, int cSeed)                     !InTf!
 interface                                                                                !InTf!
   subroutine RanNormalSetSeedZigVec(stream, piSeed, cSeed) bind(C,name='F_RanNormalSetSeedZigVec')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                         !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine RanNormalSetSeedZigVec                                                  !InTf!
 end interface                                                                            !InTf!

 type(RANDOM_STREAM) :: ran250
 integer, parameter :: BUFSZ = 1000
 integer, dimension(BUFSZ) :: ibuf
 real*8, dimension(BUFSZ) :: rbuf
 return
end
#endif
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

// macro used to instrument code with counters
#define INSTRUMENT(A) ;

#include <randomgeneric.h>

static int s_cNormalInStore = 0;		     /* > 0 if a normal is in store */

static void RanVecIntU(void *STREAM, unsigned int *auiRan, int cRan)
{
// (*s_fnVecIRanu)(STREAM, auiRan, cRan);
  VecIRan_R250_stream(STREAM, auiRan, cRan);
}
static void RanSetSeed(void *STREAM, unsigned int *piSeed, int cSeed)
{
  s_cNormalInStore = 0;
// (*s_fnRanSetSeed)(STREAM, piSeed, cSeed);
  RanSetSeed_R250_stream(STREAM, piSeed, cSeed);
}

/*------------------------------ General Ziggurat --------------------------*/
#define ZIGNOR_C 128			       /* number of blocks */
#define ZIGNOR_R 3.442619855899	/* start of the right tail */
				   /* (R * phi(R) + Pr(X>=R)) * sqrt(2\pi) */
#define ZIGNOR_V 9.91256303526217e-3
#define ZIGNOR_STORE 256 * 4

/* s_adZigX holds coordinates, such that each rectangle has*/
/* same area; s_adZigR holds s_adZigX[i + 1] / s_adZigX[i] */
static double s_adZigX[ZIGNOR_C + 1], s_adZigR[ZIGNOR_C];
static int Zig_initialized = 0;

INSTRUMENT(static unsigned int zigcalls = 0;)
INSTRUMENT(static unsigned int zigquick = 0;)
INSTRUMENT(static unsigned int zigtails = 0;)
INSTRUMENT(static unsigned int zigwedge = 0;)
INSTRUMENT(static unsigned int zigloops = 0;)
INSTRUMENT(static unsigned int zigused  = 0;)
INSTRUMENT(static unsigned int zigcalls2 = 0;)
INSTRUMENT(static unsigned int zigquick2 = 0;)
INSTRUMENT(static unsigned int zigtails2 = 0;)
INSTRUMENT(static unsigned int zigwedge2 = 0;)
INSTRUMENT(static unsigned int zigloops2 = 0;)
INSTRUMENT(static unsigned int zigused2 = 0;)

static void zigNorInit(int iC, double dR, double dV)
{
  int i;	double f;

  if(Zig_initialized) return ;
  f = exp(-0.5 * dR * dR);
  s_adZigX[0] = dV / f; /* [0] is bottom block: V / f(R) */
  s_adZigX[1] = dR;
  s_adZigX[iC] = 0;

  for (i = 2; i < iC; ++i)
  {
    s_adZigX[i] = sqrt(-2 * log(dV / s_adZigX[i - 1] + f));
    f = exp(-0.5 * s_adZigX[i] * s_adZigX[i]);
  }
  for (i = 0; i < iC; ++i)
    s_adZigR[i] = s_adZigX[i + 1] / s_adZigX[i];
  Zig_initialized = 1;
}

double  DRanNormalZigVec(void *stream)  // !InTc!
{
  unsigned int i;
  double x, u, y, f0, f1;
  generic_state *zig = stream ;
  unsigned int *s_auiZigTmp, *s_auiZigRan;
  int s_cZigStored;
  unsigned char *s_auiZigBox;

  if(zig->gauss == NULL) {
    zig->gauss = (unsigned int *) memalign(64,(ZIGNOR_STORE + ZIGNOR_STORE / 4)*sizeof(unsigned int));
    zig->ngauss = 0;
//     printf("allocating buffer in gaussian stream, size=%d uints\n",ZIGNOR_STORE + ZIGNOR_STORE / 4);
  }
  if(Zig_initialized == 0) {  // initialize ziggurat internal tables if not already done
    zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V);
  }
  s_cZigStored = zig->ngauss ;
  s_auiZigTmp = zig->gauss ;
  s_auiZigRan = s_auiZigTmp ;
  s_auiZigBox = (unsigned char *) &s_auiZigTmp[ZIGNOR_STORE];

  INSTRUMENT(zigcalls++; ; direct = 1;)
  for (;;)
  {
    if (s_cZigStored <= 0) { 
      RanVecIntU(stream, s_auiZigTmp, ZIGNOR_STORE + ZIGNOR_STORE / 4); s_cZigStored = ZIGNOR_STORE ; 
    } // FillBufferZig(stream);
    --s_cZigStored;

    u = RANDBLS_32new(s_auiZigRan[s_cZigStored]);   // convert from 32 bit int to (-1.0 , 1.0) on the fly
    i = s_auiZigBox[s_cZigStored] & 0x7F;
    INSTRUMENT(zigused += 2;)
    
    if (fabs(u) < s_adZigR[i]){		 /* first try the rectangular boxes */
      INSTRUMENT(zigquick += direct;)
      zig->ngauss = s_cZigStored;
      return u * s_adZigX[i];
    }
    INSTRUMENT(direct = 0;)

    if (i == 0){						/* bottom box: sample from the tail */
      INSTRUMENT(zigtails++ ; )
      do                   //   return DRanNormalTail(ZIGNOR_R, u < 0);
      {
	if(s_cZigStored <= 1) { RanVecIntU(stream, s_auiZigTmp, ZIGNOR_STORE + ZIGNOR_STORE / 4); s_cZigStored = ZIGNOR_STORE ; }
	x = RANDBL_32new(s_auiZigRan[--s_cZigStored]);
	x = log(x) / ZIGNOR_R;
	y = RANDBL_32new(s_auiZigRan[--s_cZigStored]);
	y = log(y);
	INSTRUMENT(zigused += 2;)
      } while (-2 * y < x * x);
      zig->ngauss = s_cZigStored;
      return (u < 0) ? x - ZIGNOR_R : ZIGNOR_R - x;
    }

    INSTRUMENT(zigwedge ++;)
    x = u * s_adZigX[i];		   /* is this a sample from the wedges? */
    f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x) );
    f1 = exp(-0.5 * (s_adZigX[i + 1] * s_adZigX[i + 1] - x * x) );
    INSTRUMENT(zigused++;)
    if(s_cZigStored <= 0) { RanVecIntU(stream, s_auiZigTmp, ZIGNOR_STORE + ZIGNOR_STORE / 4); s_cZigStored = ZIGNOR_STORE ; }
    y = RANDBL_32new(s_auiZigRan[--s_cZigStored]);
    zig->ngauss = s_cZigStored;
    if (f1 + y * (f0 - f1) < 1.0)  return x;
    INSTRUMENT(zigloops++;)
  }
}
double  F_DRanNormalZigVec(statep *s)
{
  return (DRanNormalZigVec(s->p)) ;
}

void  RanNormalSetSeedZig(void *stream, unsigned int *piSeed, int cSeed)  // !InTc!
{
  generic_state *zig = stream ;

  zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V);
  RanSetSeed(stream, piSeed, cSeed);
  if(zig->gauss == NULL) {
    zig->gauss = (unsigned int *) memalign(64,(ZIGNOR_STORE + ZIGNOR_STORE / 4)*sizeof(unsigned int));
    zig->ngauss = 0;
  }
  DRanNormalZigVec(stream);
}

void F_RanNormalSetSeedZig(statep *s, unsigned int *piSeed, int cSeed)  // !InTc!
{
  RanNormalSetSeedZig(s->p, piSeed, cSeed);
}

void  RanNormalSetSeedZigVec(void *stream, unsigned int *piSeed, int cSeed)  // !InTc!
{
  RanNormalSetSeedZig(stream, piSeed, cSeed);
}

void  F_RanNormalSetSeedZigVec(statep *s, unsigned int *piSeed, int cSeed)
{
  RanNormalSetSeedZigVec(s->p, piSeed, cSeed);
}
/*--------------------------- END General Ziggurat -------------------------*/
/*---------------  Uniform to gaussian random number conversion ------------*/
//
// attempt to produce ngauss normally distributed random 64 bit floats from
// nuni 32 bit integer uniformly distributed random numberd
//
// gaussian [OUT], ngauss[IN/OUT], uniform[IN], nuni [IN]
// gaussian : output array of normal distribution random 64 bit floats
// ngauss   : [IN] number of values desired
//            [OUT] number of values yet to produce if not enough uniform values to do so
// uniform  : input array of uniform distribution 32 bit random integers
// nuni     : number of available input numbers
// the return value is the number of still usable input integers 
// (<0 if unable to produce the requested number of output numbers)
// NOTE: nuni should be at least ~ 2.1 times larger than ngauss to ensure success
int NormalFromUniform(double gaussian[], int *ngauss, int uniform[], int nuni)  // !InTc!
{
  int nout=*ngauss;
  int iz, iout;
  double u, x, y, f0, f1;

  iout = 0;
  do {
    for (;;) {                  // produce one normal distribution random number
      if(nuni < 2) goto ouch ;                 // not enough uniform randoms left

      u = RANDBLS_32new( uniform[--nuni]) ;
      iz = uniform[--nuni] & 0x7F ;

      if (fabs(u) < s_adZigR[iz]) {            // success
	gaussian[iout] = u * s_adZigX[iz] ;    // store normal random in output
	break ;
      }

      if(iz == 0){                             // bottom box: sample from the tail
	do{
	  if(nuni < 2) goto ouch ;             // not enough uniform randoms left
	  x = RANDBL_32new(uniform[--nuni]);
	  x = log(x) / ZIGNOR_R;
	  y = RANDBL_32new(uniform[--nuni]);
	  y = log(y) ;
	} while (-2 * y < x * x);
	gaussian[iout] = (u < 0) ? x - ZIGNOR_R : ZIGNOR_R - x;
	break ;
      }

      x = u * s_adZigX[iz];		       // is this a sample from the wedges? */
      f0 = exp(-0.5 * (s_adZigX[iz    ] * s_adZigX[iz    ] - x * x) );
      f1 = exp(-0.5 * (s_adZigX[iz + 1] * s_adZigX[iz + 1] - x * x) );
      if(nuni < 1) goto ouch ;                 // not enough uniform randoms left
      y = RANDBL_32new(uniform[--nuni]);
      if (f1 + y * (f0 - f1) < 1.0) {          // inside wedge
	gaussian[iout] = x ;
	break ;
      }
    }  // for (;;), point is done
    iout++ ;
  }while(iout < nout);

  *ngauss -= iout;   // should be zero at this point
  return (nuni) ;

ouch:
  *ngauss -= iout;
  return(-1);
}
/*------------------------ END of gaussian conversion ----------------------*/

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
  dmax = RANDBL_32new(maxpos) ;
  dmin = RANDBL_32new(maxneg) ;
  printf("maxpos, maxneg transformed with RANDBL_32new  : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  dmax = RANDBLS_32new(maxpos) ;
  dmin = RANDBLS_32new(maxneg) ;
  printf("maxpos, maxneg transformed with RANDBLS_32new : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);

   RanNormalSetSeedZig(R250,r250.buffer, 250);   // initialize to values already there :-)

  dmin = 0.0 ; dmax = 0.0;
  for( i=0 ; i < 10 ; i++) gaussdist[i] = 0;
  for( i=0 ; i < 100000000 ; i++) {
    rval = DRanNormalZigVec(R250);
    avg = avg + rval ;
    dmin = (dmin < rval) ? dmin : rval ;
    dmax = (dmax > rval) ? dmax : rval ;
    rval = rval > 0 ? rval : -rval ;
    index = rval;
    gaussdist[index] ++;
  }
  printf("dmin = %6.3f, dmax = %6.3f, avg = %10.7f\n",dmin,dmax,avg/i);
  for( i=0 ; i < 10 ; i++) printf("%9d ",gaussdist[i]);
  printf("\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRanNormalZigVec(R250);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random DRanNormalZigVec/R250 double value = %6.3f \n",t1-t0);  // DRanNormalZigVec

  t1 = 0 ; t0 = 1 ; 
  INSTRUMENT(t1 = zigquick ; t0 = zigcalls+1 ; )
  INSTRUMENT(printf("quick calls in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = zigtails ;)
 INSTRUMENT( printf("tail  calls in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = zigwedge ;)
  INSTRUMENT(printf("wedge calls in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = zigloops ;)
  INSTRUMENT(printf("extra loops in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = zigused ;)
  INSTRUMENT(printf("random values used in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)

  INSTRUMENT(t1 = zigquick2 ; t0 = zigcalls2+1 ; )
  INSTRUMENT(printf("quick calls in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = zigtails2 ;)
  INSTRUMENT(printf("tail  calls in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = zigwedge2 ;)
  INSTRUMENT(printf("wedge calls in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = zigloops2 ;)
  INSTRUMENT(printf("extra loops in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = zigused2 ;)
  INSTRUMENT(printf("random values used in ziggurat = %6.3f%\n",t1 / t0 * 100.0);)

  MPI_Finalize();
  return(0);
}
#endif
