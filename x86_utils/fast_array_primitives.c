#include <stdio.h>
#undef __SSExx__
#if defined(__SSE__)
#include <x86intrin.h>
#endif

#define VLEN 4
#define VSHIFT 3
#define VMASK 7

// IEEE 754 components
#define IEEE32_EXP(a) ( ((unsigned int) (a) >> 23) & 0xFF )
#define IEEE32_SGN(a) ( ((unsigned int) (a) >> 31) )
#define IEEE32(sign,exp,mantissa) ( (sign << 31) | ( (exp+127) << 23) | ((mantissa) & 0x7FFFFF) )
// minmax
#define MAX(a,b) ( (a > b) ? (a) : (b) )
#define MIN(a,b) ( (a < b) ? (a) : (b) )
// shuffle control
#define SELECT(a,b,c,d) ( a + (b<<2) + (c<<4) + (d<<6) )
// endian swap macros
#define ESWAP8x4_128(xmm)                    // 8 in 32
#define ESWAP8x4_256(ymm)
#define ESWAP16x2_128(xmm)                    // 16 in 32
#define ESWAP16x2_256(ymm)
#define ESWAP8x2_128(xmm)                    // 8 in 16
#define ESWAP8x2_256(ymm)


// get min value, max value, very rough average of 32 bit float array
// the average is only ~OK if n is a multiple of 2 * VLEN
// otherwise some points are added twice (we sum 2* np2 points)
void FloatFastQuantizeLinear(float *z_in, short *quant, int n, int nbits, float *Max, float *Min, float *Avg, float *Rescl)
{
  int i, j, np, np2, offset, maxexp, izero, mask;
  float scale, fzero;
  float *z_in0 = z_in;
#if defined(__SSE__)
  __m128  x0, x1, x2, xxmin, xxmax, xxsum, xxrng, xxsca;
  __m256  y0, yymin, yysca, point5;
  __m128i ix0, ix1;
  __m256i iy0, yymsk;
#else
  float tmax, tmin, tsum, vsum[VLEN], vmax[VLEN], vmin[VLEN];
#endif
  union {
    unsigned int i;
    float f;
  } xmax, xmin, xrng, xscl;
// initialize min, max, sum
#if defined(__SSE__)
  xxmin  = _mm_loadu_ps(&z_in[ 0]) ;  // min = max = 8 first elements of array
  xxmax  = xxmin;
  xxsum  = _mm_xor_ps(xxsum,xxsum) ;  // xxsum = 0
#else
  for(j=0 ; j<VLEN ; j++){
    vmax[j] = z_in[j] ;         // min = max = VLEN first elements of array
    vmin[j] = vmax[j] ;
    vsum[j] = 0.0 ;             // xxsum = 0
  }
#endif
  np     = ((n+VMASK) >> VSHIFT ) << VSHIFT ; // next higher multiple of 2 * VLEN
  np2    = np >> 1;                           // np/2
  offset =  n - np2 ;                         // offset to "upper" part with length a multiple of VLEN
  
  for (i=0 ; i<np2 ; i+=VLEN){
#if defined(__SSE__)
    x0    = _mm_loadu_ps(&z_in[i]) ;          // stream from beginning
    x1    = _mm_loadu_ps(&z_in[i+offset]) ;   // stream from "midpoint"
    xxsum = _mm_add_ps(xxsum,x0);
    xxmax = _mm_max_ps(xxmax,x0);
    xxmin = _mm_min_ps(xxmin,x0);
    xxsum = _mm_add_ps(xxsum,x1);
    xxmax = _mm_max_ps(xxmax,x1);
    xxmin = _mm_min_ps(xxmin,x1);
#else
    for(j=0 ; j<VLEN ; j++){
      vmax[j] = (vmax[j] > z_in[j]) ? vmax[j] : z_in[j] ;
      vmin[j] = (vmin[j] < z_in[j]) ? vmin[j] : z_in[j] ;
      vsum[j] = vsum[j] + z_in[j] ;
      vmax[j] = (vmax[j] > z_in[j+offset]) ? vmax[j] : z_in[j+offset] ;
      vmin[j] = (vmin[j] < z_in[j+offset]) ? vmin[j] : z_in[j+offset] ;
      vsum[j] = vsum[j] + z_in[j+offset] ;
    }
    z_in += VLEN;
#endif
  }
#if defined(__SSE__)
  x0 = _mm_shuffle_ps(xxmin,xxmin,SELECT(2,3,2,3)) ;  //  put elements 2 and 3 on top of 0 and 1
  x1 = _mm_shuffle_ps(xxmax,xxmax,SELECT(2,3,2,3)) ;
  x2 = _mm_shuffle_ps(xxsum,xxsum,SELECT(2,3,2,3)) ;
  xxmin = _mm_min_ps(xxmin,x0);   // min(0,2) , min(1,3), ?? , ??
  xxmax = _mm_max_ps(xxmax,x1);
  xxsum = _mm_add_ps(xxsum,x2);
  x0 = _mm_shuffle_ps(xxmin,xxmin,SELECT(1,1,1,1)) ;  // put f(1,3) on top of f(0,2)
  x1 = _mm_shuffle_ps(xxmax,xxmax,SELECT(1,1,1,1)) ;
  x2 = _mm_shuffle_ps(xxsum,xxsum,SELECT(1,1,1,1)) ;
  xxmin = _mm_min_ps(xxmin,x0);   // f( f(0,2) , f(1,3) ) , ?? , ??, ??
  xxmax = _mm_max_ps(xxmax,x1);
  xxsum = _mm_add_ps(xxsum,x2);
  xxmin = _mm_shuffle_ps(xxmin,xxmin,SELECT(0,0,0,0)) ;  // broadcast
  xxmax = _mm_shuffle_ps(xxmax,xxmax,SELECT(0,0,0,0)) ;
  xxrng = _mm_sub_ps(xxmax,xxmin) ;
  xxsum = _mm_shuffle_ps(xxsum,xxsum,SELECT(0,0,0,0)) ;
  _mm_store_ss(Min,xxmin);
  _mm_store_ss(Max,xxmax);
  _mm_store_ss(Avg,xxsum);
  *Avg /= np;
#else
// final wrapup, scalar loop over VLEN terms, only in non SSE case
  tmax = vmax[0] ;
  tmin = vmin[0] ;
  tsum = vsum[0] ;
  for(j=1 ; j<VLEN ; j++){
    tmax = (tmax > vmax[j]) ? tmax : vmax[j] ;
    tmin = (tmin < vmin[j]) ? tmin : vmin[j] ;
    tsum = tsum + vsum[j] ;    
  }
  *Max  = tmax;
  *Min  = tmin;
  *Avg  = tsum / np ;
#endif
// code to quantize floating poing into 16 bit tokens (nbits significant bits)
  mask = ~ (-1 << nbits) ;
  xmax.f = *Max;
  xmin.f = *Min;
  xrng.f = xmax.f - xmin.f;
  maxexp = MAX( IEEE32_EXP(xmax.i) , IEEE32_EXP(xmin.i) );
  maxexp = MAX( IEEE32_EXP(xrng.i) , maxexp);
  xscl.i = IEEE32(0,(nbits-1) - (maxexp-127),0);
  scale = xscl.f;   // TO BE FIXED, used to convert floats to range  0 <= X < 2**nbits -1 
  xscl.i = IEEE32(0,(maxexp-127) - (nbits-1),0); // inverse scaling factor
  *Rescl = xscl.f;
  if(IEEE32_SGN(xmax.i) + IEEE32_SGN(xmin.i) == 1 ) {   // min < 0 and max >0
    izero = (0 - *Min) * scale + .5;
    fzero = *Min + izero * *Rescl;
    *Min -= fzero ;
//     izero = (0 - *Min) * scale + .5;
//     printf("rescaled fzero = %f, from %f, new min = %f\n",*Min + izero / scale,fzero,*Min);
  }
//   printf("max min range maxexp scale mask is %f %f %f %i %f %8.8x\n",xmax.f,xmin.f,xrng.f,maxexp-127,scale,mask);
#if defined(__AVX2__) && defined(__FMA__)
  yysca = _mm256_set1_ps(scale) ;
  point5 = _mm256_set1_ps(0.5) ;
  yymin = _mm256_set1_ps(xmin.f) ;
  ix0 = (__m128i) _mm_broadcast_ss( (void *) &mask ) ;
  yymsk = _mm256_inserti128_si256(yymsk,ix0,0) ;
  yymsk = _mm256_inserti128_si256(yymsk,ix0,1) ;
  for (i=0 ; i<np2 ; i+=VLEN){
    x0    = _mm_loadu_ps(&z_in0[i]) ;          // stream from beginning
    x1    = _mm_loadu_ps(&z_in0[i+offset]) ;   // stream from "midpoint"
    y0    = _mm256_insertf128_ps(y0,x0,0) ;
    y0    = _mm256_insertf128_ps(y0,x1,1) ;
    y0    = _mm256_sub_ps(y0,yymin) ;          // x - min
    y0    = _mm256_fmadd_ps(y0,yysca,point5) ; // * scale + .5
    iy0   = _mm256_cvttps_epi32(y0) ;          // to integer (lower 16 bits used only)
    iy0   = _mm256_min_epi32(iy0,yymsk) ;
    ix0   = _mm256_extracti128_si256(iy0,0) ;  // lower 128 bits
    ix1   = _mm256_extracti128_si256(iy0,1) ;  // upper 128 bits
    ix0   = _mm_packus_epi32(ix0,ix1) ;        // pack 32 -> 16 ix0 and ix1 , then stuff both back into ix0 (ix1 in upper part)
    _mm_storel_epi64((void *)&quant[i],ix0) ;              // store lower 64 bits
    _mm_storeh_pi((void *)&quant[i+offset],(__m128)ix0) ;  // store upper 64 bits
#else
  for (i=0 ; i<np2 ; i+=VLEN){
    for(j=0 ; j<VLEN ; j++){
      quant[j] = (z_in0[j] - *Min) * scale;
      quant[j+offset] = (z_in0[j+offset] - *Min) * scale + .5;
    }
    z_in0 += VLEN;
    quant += VLEN;
#endif
  }  // for (i=0 ; i<np2 ; i+=VLEN)
}

#if defined(SELF_TEST)
#define ASIZE 1026
#define RND 0.5
#define NBITS 16
main()
{
  int i;
  float a[ASIZE+10];
  short int ia[ASIZE+10];
  int nbits = 16;
  float mi, ma, av, range, scal, mi0, fac, rescl;
  int ima, imi, exp;
  union {
    unsigned int i;
    float f;
  } m;
  int p2m16 = 0x37800000 ;  // 2.0 ** -16
  int p2p16 = 0x47800000 ;  // 2.0 ** +16
  int one   = 0x3f800000 ;  // 1.0
  int onep  = 0x3f80000f ;  // 1.0 + epsilon
  int round = 0x000080;
  int mask  = 0x00FFFF00 ;

  mi = 1.0;
  for (i=0;i<23;i++) mi = mi * 2.0;
  m.f = mi;
  printf("2 ** 23   = %8.8x , exp = %d\n\n",m.i,m.i>>23);

  m.f = 1.000001;
//   printf("%15f = %8.8x \n",m.f,((m.i) & 0x7FFFFF) | 0x800000);
  imi = m.f * mi;
//   printf(" imi            = %8.8x \n",imi);

  mi  = 1048574 ;
//   mi = mi * 8 + 7;
  m.f = mi;
  printf("%15g = %8.8x \n",mi,((m.i) & 0x7FFFFF) | 0x800000);
  exp = m.i >> 23;
  m.i = (127 + (150 - exp)) << 23;
  fac = m.f;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  if (imi > mask) imi = mask;
  printf(" fac = %g\n",fac);
  printf(" imi            = %8.8x, %12d \n\n",imi,imi);

  mi = 65535;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  if (imi > mask) imi = mask;
  m.f = mi;
  printf("%15f,  %8.8x\n",mi,(m.i & 0x00FFFFFF ) | 0x800000);
  printf(" imi            = %8.8x, %12d \n\n",imi, imi);

  mi = 4094.5;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  m.f = mi;
  printf("%15f,  %8.8x\n",mi,(m.i & 0x00FFFFFF ) | 0x800000);
  printf(" imi            = %8.8x, %12d \n\n",imi, imi);

  mi = 254.9;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  if (imi > mask) imi = mask;
  m.f = mi;
  printf("%15f,  %8.8x\n",mi,(m.i & 0x00FFFFFF ) | 0x800000);
  printf(" imi            = %8.8x, %12d \n\n",imi, imi);

  mi = 255;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  if (imi > mask) imi = mask;
  m.f = mi;
  printf("%15f,  %8.8x\n",mi,(m.i & 0x00FFFFFF ) | 0x800000);
  printf(" imi            = %8.8x, %12d \n\n",imi, imi);
//   exit (0);
  
  a[0] = 1.0;
  for (i=1 ; i<ASIZE ; i++){
    a[i] = a[i-1]*1.01;
  }
  a[4] = -8.7 ; a[7] = 65535.99 + a[4] ;
  printf("a[1,4,7,ASIZE-1] %8f %8f %8f %8f\n",a[1],a[4],a[7],a[ASIZE-2]);
  FloatFastQuantizeLinear(&a[1],&ia[1],ASIZE-2,nbits,&ma,&mi,&av,&rescl);
  printf("mi=%f, av=%f, ma=%f, rescl=%f\n",mi,av,ma,rescl);
  av = 0;
  for (i=1 ; i<ASIZE-1 ; i++) { av = av + a[i] ; } ;
  av = av / (ASIZE-2) ;
  printf("expected mi=%f, expected av=%f, expected ma=%f\n",-1.2,av,a[7]);
  printf("a[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",a[1],a[4],a[7],a[ASIZE-2]);
  printf("ia[1,4,7,ASIZE-2] %8hu %8hu %8hu %8hu\n",ia[1],ia[4],ia[7],ia[ASIZE-2]);

  m.f = ma;
  ima = m.i >> 23;
  m.f = mi;
  imi = m.i >> 23;
  imi &= 0xFF;
  if(ima > imi) m.i = (m.i >> (ima-imi)) << (ima-imi) ;
  mi0 = m.f;

  range = (ma - mi0) ;
  m.f = range ;                     // range
  m.i = m.i >> 23 ;                 // exponent of range
  m.i = m.i + 1 ;                   // range = next power of 2 > range
  m.i = (m.i - NBITS) -127 ;        // range / 2**NBITS, remove exponent bias
  m.i = 127 - m.i;                  // 1 / range , add exponent bias
  m.i = m.i << 23;                  // shift exponent into position
  scal =  m.f;                      // 1.0 / discretization interval
  printf("m.f = %f, m.i=%8.8x, range = %f, discr = %f, mi = %f, mi0 = %f\n",m.f,m.i,range,1.0/scal,mi,mi0);

//   m.i = p2m16;
//   scal = range * m.f;
  imi = (mi - mi0)*scal + RND;
  ima = (ma - mi0)*scal + RND;
  printf("range = %f, %f = %d, %f = %d\n",range,mi,imi,ma,ima);
}
#endif
