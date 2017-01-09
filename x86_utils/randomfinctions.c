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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
// #include <unistd.h>
#include <math.h>

#ifdef __LP64__
	typedef unsigned long UINT64;
	typedef          long INT64;
	#define LIT_UINT64(c) (c##ul)
	#define LIT_INT64(c)  (c##l)
#elif defined(_MSC_VER)
	typedef unsigned __int64 UINT64;
	typedef          __int64 INT64;
	#define LIT_UINT64(c) (c##ui64)
	#define LIT_INT64(c)  (c##i64)
#else 
	typedef unsigned long long UINT64;
	typedef          long long INT64;
	#define LIT_UINT64(c) (c##ull)
	#define LIT_INT64(c)  (c##ll)
#endif

#define M_RAN_INVM32	2.32830643653869628906e-010			  /* 1.0 / 2^32 */
#define M_RAN_INVM31	4.65661287307739257812e-010			  /* 1.0 / 2^31 */
#define M_RAN_INVM48	3.55271367880050092936e-015			  /* 1.0 / 2^48 */
#define M_RAN_INVM52	2.22044604925031308085e-016			  /* 1.0 / 2^52 */

#define RANDBL_32new(iRan1)          ((int)(iRan1) * M_RAN_INVM32 + (0.5 + M_RAN_INVM32 / 2))
#define RANDBLS_32new(iRan1)         ((int)(iRan1) * M_RAN_INVM31 + (M_RAN_INVM31 / 2))
#define RANDBL_48new(iRan1, iRan2)   ((int)(iRan1) * M_RAN_INVM32 + (0.5 + M_RAN_INVM48 / 2) + (int)((iRan2) & 0x0000FFFF) * M_RAN_INVM48)
#define RANDBL_52new(iRan1, iRan2)   ((int)(iRan1) * M_RAN_INVM32 + (0.5 + M_RAN_INVM52 / 2) + (int)((iRan2) & 0x000FFFFF) * M_RAN_INVM52)

/* plug-in RNG */
typedef double 		( * DRANFUN)(void);
typedef unsigned int	( * IRANFUN)(void);
typedef void   		( * IVECRANFUN)(unsigned int *, int);
typedef void   		( * DVECRANFUN)(double *, int);
typedef void   		( * RANSETSEEDFUN)(int *, int);

/*---------------------------- GetInitialSeeds -----------------------------*/
void GetInitialSeeds(unsigned int auiSeed[], int cSeed, unsigned int uiSeed, unsigned int uiMin)
{
  int i;
  unsigned int s = uiSeed;     /* may be 0 */

  for (i = 0; i < cSeed; ) {   /* see Knuth p.106, Table 1(16) and Numerical Recipes p.284 (ranqd1)*/
    s = 1664525 * s + 1013904223;
    if (s <= uiMin) continue;
    auiSeed[i++] = s;
  }
}
/*-------------------------- END GetInitialSeeds ---------------------------*/
/*==========================================================================
 * R250 ans SHR3 generators added, naming is consistent with the MWC8222 code
 *==========================================================================*/
/*------------------------ start of SHR3 addition --------------------------*/
typedef struct{
  unsigned long jsr;
} shr3_state;

static shr3_state shr3 = { 123456789 } ;

static unsigned long jsr=123456789;
// #define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)

void RanSetSeed_SHR3(int *piSeed, int cSeed)
{
  if(piSeed == NULL || cSeed == 0) return ; // null call, nothing to do
  jsr = (unsigned) *piSeed ;
}

unsigned int IRan_SHR3(void)		/* returns a random unsigned integer */
{
  unsigned long jz=jsr ; 
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  return (jz+jsr) ;
}

double DRan_SHR3(void)		/* returns a random double (0.0 , 1.0) */
{
  unsigned long jz=jsr ; 
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  return RANDBL_32new(jz+jsr);   // convert from 32 bit int to (0.0 , 1.0)
}

double DRanS_SHR3(void)		/* returns a random double (-1.0 , 1.0) */
{
  unsigned long jz=jsr ; 
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  return RANDBLS_32new(jz+jsr);   // convert from 32 bit int to (-1.0 , 1.0)
}

void VecIRan_SHR3(unsigned int *ranbuf, int n){
  int i;
  unsigned long jz;

//   for(i=0 ; i<n ; i++) ranbuf[k] = SHR3;
  for(i=0 ; i<n ; i++) {
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    ranbuf[i] = (jz+jsr);
  };
}

void VecDRan_SHR3(double *ranbuf, int n){
  int i;
  unsigned long jz;

//   for(i=0 ; i<n ; i++) ranbuf[k] = RANDBL_32new(SHR3);
  for(i=0 ; i<n ; i++){
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    ranbuf[i] = RANDBL_32new(jz+jsr);   // convert from 32 bit int to (0.0 , 1.0)
  }
}

/*------------------------- end of SHR3 addition ---------------------------*/
/*------------------------ start of R250 addition --------------------------*/
typedef struct{
  int index;
  unsigned int buffer[250];
}r250_state ;

static r250_state r250 = {
  0 ,
  {
    0x17617168 ,0x17a0d192 ,0x7cba449a ,0x86d91b38 ,0x7455bfb5 ,0x3bb194f2 ,0xd4cf89e2 ,0xee85f453 ,0x916c5ad3 ,0x8e5c32bb ,
    0x0d035201 ,0xc45d76fd ,0x04d799be ,0x409aa63c ,0x4c56a663 ,0xe0404838 ,0x88141c55 ,0xeb963899 ,0x83f34bc2 ,0xcf784dfe ,
    0x303c4f3d ,0xe2e580de ,0x9213c2a7 ,0x13a77f75 ,0xd0b5ca05 ,0xfc21d156 ,0x1e4c4c31 ,0xe4d70461 ,0xc8d3eb8a ,0x19d5a61d ,
    0x40d8f714 ,0x051457d9 ,0x60dd589d ,0xadcae8f4 ,0xae50f8e4 ,0xc0b1623d ,0x019840e2 ,0x9fd91e54 ,0x205fd793 ,0x9e17be66 ,
    0xc55b8dfa ,0x2fc0cfc5 ,0x7f292738 ,0x826e9e1f ,0x50e160e4 ,0xb3a2c76f ,0x5e03e011 ,0xe05601d4 ,0x5ef981a0 ,0x8676f0f1 ,
    0x6d8690cb ,0x8af13789 ,0xab00a515 ,0x2bcff371 ,0xd4adf100 ,0xdca5d50d ,0xc4ebad34 ,0x0d75d8f3 ,0xf77b9800 ,0x5593a16e ,
    0x1504475c ,0x453356ee ,0x1302dbd3 ,0xad1d21af ,0xcf668d00 ,0x76d68c1c ,0x91b61c7e ,0x831fbcb8 ,0x96258a86 ,0xbd1f4791 ,
    0xdf26ff85 ,0x93f5aa8d ,0xf6f6d072 ,0x55622d0f ,0x64905fb1 ,0x2d4b22c3 ,0x1ca372c3 ,0xda494a98 ,0xcfa8f513 ,0xf4737e2e ,
    0xe77f1a9c ,0x238ff343 ,0xafa9f948 ,0x20823fdf ,0x010c0e27 ,0x0020f353 ,0x31654bfd ,0x9637394d ,0x5ffd52b6 ,0x06db7185 ,
    0x76fadbb5 ,0xc90c3dd1 ,0x71364d2a ,0x49b411a8 ,0xda7556c9 ,0x86a61957 ,0xba798498 ,0x442d2c73 ,0xc89268b0 ,0x77d1be52 ,
    0x326e6b29 ,0x727a3092 ,0xba7b0780 ,0xf751def0 ,0xc1c05141 ,0x774f8101 ,0xde7495b4 ,0x2250658f ,0xd896e6e2 ,0x04a1649d ,
    0xf0dbdc2f ,0xab8a806e ,0xcccaef53 ,0x24a9daf8 ,0x9c472081 ,0xe88534a7 ,0xc6400afc ,0x5db66a6c ,0xf1bcdb68 ,0xe283ac5c ,
    0x71f6bb93 ,0x29bc3b72 ,0x6cc0d398 ,0x7ad7b430 ,0x07f8067d ,0x9bc21534 ,0xcad1a0be ,0x2eb81bf9 ,0x866e596e ,0xe6150f60 ,
    0x4f0cdac2 ,0xf14feeb3 ,0x63559552 ,0xee360d3e ,0xe950794d ,0x5674b5b7 ,0x8866636f ,0xb3bb5604 ,0xf894278a ,0xc4a13b2b ,
    0xc07f38ca ,0x1742de68 ,0xda9d902f ,0x2cb57fbc ,0x050fabca ,0x6b471777 ,0x4bdee5be ,0x60d3d78a ,0x6b31001b ,0x237156a5 ,
    0x617f2b07 ,0x4cbe6264 ,0x05610cd5 ,0x3d222b24 ,0x097ffff2 ,0x38975335 ,0x8682e0dc ,0x7e993479 ,0xe1ecf67b ,0xc619babd ,
    0x0f1ac989 ,0x1eab4e4b ,0x8c3cb3fc ,0x5787c983 ,0x74f19f89 ,0x968f257d ,0x95cc62b2 ,0x11e6bd09 ,0xd1a57e05 ,0x67358a7f ,
    0x95e23779 ,0x30efec41 ,0xe46c4803 ,0x2d2414e1 ,0x352c0fda ,0x1da1a740 ,0x28aea00b ,0xfe1dec28 ,0xae7b6c47 ,0xcfd1de31 ,
    0xa468360d ,0x544fc9e5 ,0xbcd04aa4 ,0xd2cfc38b ,0xb8a48f82 ,0xa8718902 ,0x5bd8a509 ,0x9c40dd86 ,0x6a3dadd0 ,0xd0a0d65f ,
    0xc62298c6 ,0x46393aca ,0x0b7436f2 ,0x99ddd69c ,0x839b79a7 ,0xa155be69 ,0x2e4f0458 ,0x474bd538 ,0x73d65578 ,0xa49ab70f ,
    0xbe2a3c0b ,0x69e550db ,0x9e38abcb ,0x9e483578 ,0xdabc5814 ,0x2e73f8ef ,0x4ed45df8 ,0x05f8d621 ,0x0259c01e ,0xf3927074 ,
    0xfda21b64 ,0x3476f241 ,0x9aa5d95a ,0xef86ea14 ,0x8f3fce06 ,0x8bff6bfa ,0x706ab0a2 ,0x7322f175 ,0x4e8acb27 ,0x336889cc ,
    0x373ea2e0 ,0x0cc5f5ce ,0x35a5cc68 ,0x93169549 ,0xea31a7b1 ,0x6a6569bc ,0xa776f509 ,0x5b0f310e ,0x96322244 ,0x64568c56 ,
    0x08aa6767 ,0x491799f1 ,0x17735c88 ,0x71c32f7e ,0xed0a2ec6 ,0xebd94777 ,0x9b1e1086 ,0xdc740f7a ,0x03c48151 ,0xafcb9f88 ,
    0xd835a40a ,0x21308fc2 ,0x0f459e5e ,0x0358b165 ,0x6422fa89 ,0xdd9cf11b ,0x03daccf5 ,0xec9e2bd9 ,0xe300013e ,0xa97d54e4 
  }
};

static unsigned int r250buffer[250] = {
0x17617168 ,0x17a0d192 ,0x7cba449a ,0x86d91b38 ,0x7455bfb5 ,0x3bb194f2 ,0xd4cf89e2 ,0xee85f453 ,0x916c5ad3 ,0x8e5c32bb ,
0x0d035201 ,0xc45d76fd ,0x04d799be ,0x409aa63c ,0x4c56a663 ,0xe0404838 ,0x88141c55 ,0xeb963899 ,0x83f34bc2 ,0xcf784dfe ,
0x303c4f3d ,0xe2e580de ,0x9213c2a7 ,0x13a77f75 ,0xd0b5ca05 ,0xfc21d156 ,0x1e4c4c31 ,0xe4d70461 ,0xc8d3eb8a ,0x19d5a61d ,
0x40d8f714 ,0x051457d9 ,0x60dd589d ,0xadcae8f4 ,0xae50f8e4 ,0xc0b1623d ,0x019840e2 ,0x9fd91e54 ,0x205fd793 ,0x9e17be66 ,
0xc55b8dfa ,0x2fc0cfc5 ,0x7f292738 ,0x826e9e1f ,0x50e160e4 ,0xb3a2c76f ,0x5e03e011 ,0xe05601d4 ,0x5ef981a0 ,0x8676f0f1 ,
0x6d8690cb ,0x8af13789 ,0xab00a515 ,0x2bcff371 ,0xd4adf100 ,0xdca5d50d ,0xc4ebad34 ,0x0d75d8f3 ,0xf77b9800 ,0x5593a16e ,
0x1504475c ,0x453356ee ,0x1302dbd3 ,0xad1d21af ,0xcf668d00 ,0x76d68c1c ,0x91b61c7e ,0x831fbcb8 ,0x96258a86 ,0xbd1f4791 ,
0xdf26ff85 ,0x93f5aa8d ,0xf6f6d072 ,0x55622d0f ,0x64905fb1 ,0x2d4b22c3 ,0x1ca372c3 ,0xda494a98 ,0xcfa8f513 ,0xf4737e2e ,
0xe77f1a9c ,0x238ff343 ,0xafa9f948 ,0x20823fdf ,0x010c0e27 ,0x0020f353 ,0x31654bfd ,0x9637394d ,0x5ffd52b6 ,0x06db7185 ,
0x76fadbb5 ,0xc90c3dd1 ,0x71364d2a ,0x49b411a8 ,0xda7556c9 ,0x86a61957 ,0xba798498 ,0x442d2c73 ,0xc89268b0 ,0x77d1be52 ,
0x326e6b29 ,0x727a3092 ,0xba7b0780 ,0xf751def0 ,0xc1c05141 ,0x774f8101 ,0xde7495b4 ,0x2250658f ,0xd896e6e2 ,0x04a1649d ,
0xf0dbdc2f ,0xab8a806e ,0xcccaef53 ,0x24a9daf8 ,0x9c472081 ,0xe88534a7 ,0xc6400afc ,0x5db66a6c ,0xf1bcdb68 ,0xe283ac5c ,
0x71f6bb93 ,0x29bc3b72 ,0x6cc0d398 ,0x7ad7b430 ,0x07f8067d ,0x9bc21534 ,0xcad1a0be ,0x2eb81bf9 ,0x866e596e ,0xe6150f60 ,
0x4f0cdac2 ,0xf14feeb3 ,0x63559552 ,0xee360d3e ,0xe950794d ,0x5674b5b7 ,0x8866636f ,0xb3bb5604 ,0xf894278a ,0xc4a13b2b ,
0xc07f38ca ,0x1742de68 ,0xda9d902f ,0x2cb57fbc ,0x050fabca ,0x6b471777 ,0x4bdee5be ,0x60d3d78a ,0x6b31001b ,0x237156a5 ,
0x617f2b07 ,0x4cbe6264 ,0x05610cd5 ,0x3d222b24 ,0x097ffff2 ,0x38975335 ,0x8682e0dc ,0x7e993479 ,0xe1ecf67b ,0xc619babd ,
0x0f1ac989 ,0x1eab4e4b ,0x8c3cb3fc ,0x5787c983 ,0x74f19f89 ,0x968f257d ,0x95cc62b2 ,0x11e6bd09 ,0xd1a57e05 ,0x67358a7f ,
0x95e23779 ,0x30efec41 ,0xe46c4803 ,0x2d2414e1 ,0x352c0fda ,0x1da1a740 ,0x28aea00b ,0xfe1dec28 ,0xae7b6c47 ,0xcfd1de31 ,
0xa468360d ,0x544fc9e5 ,0xbcd04aa4 ,0xd2cfc38b ,0xb8a48f82 ,0xa8718902 ,0x5bd8a509 ,0x9c40dd86 ,0x6a3dadd0 ,0xd0a0d65f ,
0xc62298c6 ,0x46393aca ,0x0b7436f2 ,0x99ddd69c ,0x839b79a7 ,0xa155be69 ,0x2e4f0458 ,0x474bd538 ,0x73d65578 ,0xa49ab70f ,
0xbe2a3c0b ,0x69e550db ,0x9e38abcb ,0x9e483578 ,0xdabc5814 ,0x2e73f8ef ,0x4ed45df8 ,0x05f8d621 ,0x0259c01e ,0xf3927074 ,
0xfda21b64 ,0x3476f241 ,0x9aa5d95a ,0xef86ea14 ,0x8f3fce06 ,0x8bff6bfa ,0x706ab0a2 ,0x7322f175 ,0x4e8acb27 ,0x336889cc ,
0x373ea2e0 ,0x0cc5f5ce ,0x35a5cc68 ,0x93169549 ,0xea31a7b1 ,0x6a6569bc ,0xa776f509 ,0x5b0f310e ,0x96322244 ,0x64568c56 ,
0x08aa6767 ,0x491799f1 ,0x17735c88 ,0x71c32f7e ,0xed0a2ec6 ,0xebd94777 ,0x9b1e1086 ,0xdc740f7a ,0x03c48151 ,0xafcb9f88 ,
0xd835a40a ,0x21308fc2 ,0x0f459e5e ,0x0358b165 ,0x6422fa89 ,0xdd9cf11b ,0x03daccf5 ,0xec9e2bd9 ,0xe300013e ,0xa97d54e4 };
static int r250_index = 250 ;

static void FillBuffer_R250_stream(r250_state *r250){
  int i;
  unsigned int *r250_buffer = r250->buffer ;
  int r250_index = r250->index;

  while(r250_index > 249) r250_index -= 250;
  for (i=0 ; i< 147 ; i++) {
    r250_buffer[ i ] = r250_buffer[ i ] ^ r250_buffer[ i + 103 ];
  }
  for (i=147 ; i<250 ; i++) {
    r250_buffer[ i ] = r250_buffer[ i ] ^ r250_buffer[ i - 147 ];
  }
  r250->index = r250_index;
}

static void FillBuffer_R250(void){
  int i;
  unsigned int *r250_buffer = r250.buffer ;

  while(r250.index > 249) r250.index -= 250;
  for (i=0 ; i< 147 ; i++) {
    r250_buffer[ i ] = r250_buffer[ i ] ^ r250_buffer[ i + 103 ];
  }
  for (i=147 ; i<250 ; i++) {
    r250_buffer[ i ] = r250_buffer[ i ] ^ r250_buffer[ i - 147 ];
  }
}

void RanSetSeed_R250(int *piSeed, int cSeed)
{
  int i;

  if (cSeed == 250 && piSeed != NULL) {
    for(i=0 ; i<250; i++) r250.buffer[ i ] = piSeed[ i ] ;
  }else{
    GetInitialSeeds(r250.buffer, 250, piSeed && cSeed ? piSeed[0] : 0, 0);
  }
}

double DRan_R250(void)		/* returns a random double (0.0 , 1.0) */
{
  register int	i, j;
  register unsigned int new_rand;
  if ( r250.index > 249 ) FillBuffer_R250();
  new_rand = r250.buffer[r250.index++] ;
  return RANDBL_32new(new_rand);   // convert from 32 bit int to (0.0 , 1.0)
}

double DRanS_R250(void)		/* returns a random double (-1.0 , 1.0) */
{
  register int	i, j;
  register unsigned int new_rand;
  if ( r250.index > 249 ) FillBuffer_R250();
  new_rand = r250.buffer[r250.index++] ;
  return RANDBLS_32new(new_rand);   // convert from 32 bit int to (0.0 , 1.0)
}

unsigned int IRan_R250(void)		/* returns a random unsigned integer */
{
  register int	i, j;
  register unsigned int new_rand;
  if ( r250.index > 249 ) FillBuffer_R250();
  new_rand = r250.buffer[r250.index++];
  return new_rand;
}

void VecIRan_R250(unsigned int *ranbuf, int n){
  int k = 0;
  int i;
  unsigned int *r250_buffer = r250.buffer ;

  while( r250.index < 250 && n > 0 ){
    ranbuf[k++] = r250_buffer[r250.index++] ;
    n-- ;
  }
  if ( n == 0 ) return;
  FillBuffer_R250_stream(&r250);     // we get here if buffer is empty before n is satisfied
  while(n >= 250){        // chunks of 250 values
    for(i=0 ; i<250 ; i++) ranbuf[k+i] = r250_buffer[i] ;
    n -= 250 ;
    k += 250 ;
    FillBuffer_R250_stream(&r250) ;
  }
  while( n > 0 ){  // n < 250 at this point
    ranbuf[k++] = r250_buffer[r250.index++] ;
    n-- ;
  }
}

void VecIRan_R250_stream(r250_state *r250, unsigned int *ranbuf, int n){
  int k = 0;
  int i;
  unsigned int *r250_buffer = r250->buffer ;
  int r250_index = r250->index ;

  while( r250_index < 250 && n > 0 ){
    ranbuf[k++] = r250_buffer[r250_index++] ;
    n-- ;
  }
  if ( n == 0 ) return;
  FillBuffer_R250_stream(r250);     // we get here if buffer is empty before n is satisfied
  while(n >= 250){        // chunks of 250 values
    for(i=0 ; i<250 ; i++) ranbuf[k+i] = r250_buffer[i] ;
    n -= 250 ;
    k += 250 ;
    FillBuffer_R250_stream(r250) ;
  }
  while( n > 0 ){  // n < 250 at this point
    ranbuf[k++] = r250_buffer[r250_index++] ;
    n-- ;
  }
  r250->index = r250_index;
}

void VecDRan_R250(double *ranbuf, int n){
  int k = 0;
  int i;
  while( r250.index < 250 && n > 0 ){
    ranbuf[k++] = RANDBL_32new(r250.buffer[r250.index++]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n-- ;
  }
  if ( n == 0 ) return;
  FillBuffer_R250();     // we get here if buffer is empty
  while(n >= 250){        // chunks of 250 values
    for(i=0 ; i<250 ; i++) ranbuf[k+i] = RANDBL_32new(r250.buffer[i]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n -= 250 ;
    k += 250 ;
    FillBuffer_R250() ;
  }
  while( n > 0 ){  // n < 250 at this point
    ranbuf[k++] = RANDBL_32new(r250.buffer[r250.index++]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n-- ;
  }
}
/*------------------------- end of R250 addition ---------------------------*/
/*------------------------ George Marsaglia MWC ----------------------------*/
#define MWC_R  256
#define MWC_A  LIT_UINT64(809430660)
#define MWC_AI 809430660
#define MWC_C  362436

typedef struct{
   unsigned int uiState ;
   unsigned int uiCarry ;
   unsigned int auiState[MWC_R];
} mwc_state ;

static mwc_state mwc;

static unsigned int s_uiStateMWC = MWC_R - 1;
static unsigned int s_uiCarryMWC = MWC_C;
static unsigned int s_auiStateMWC[MWC_R];

void RanSetSeed_MWC8222(int *piSeed, int cSeed)
{
	mwc.uiState = MWC_R - 1;
	mwc.uiCarry = MWC_C;
	
	if (cSeed == MWC_R)
	{
		int i;
		for (i = 0; i < MWC_R; ++i)
		{
			mwc.auiState[i] = (unsigned int)piSeed[i];
		}
	}
	else
	{
		GetInitialSeeds(mwc.auiState, MWC_R, piSeed && cSeed ? piSeed[0] : 0, 0);
	}
}
unsigned int IRan_MWC8222(void)
{
	UINT64 t;

	mwc.uiState = (mwc.uiState + 1) & (MWC_R - 1);
	t = MWC_A * mwc.auiState[mwc.uiState] + mwc.uiCarry;
	mwc.uiCarry = (unsigned int)(t >> 32);
	mwc.auiState[mwc.uiState] = (unsigned int)t;
    return (unsigned int)t;
}
double DRan_MWC8222(void)         /* returns a random double (0.0 , 1.0) */
{
	UINT64 t;

	mwc.uiState = (mwc.uiState + 1) & (MWC_R - 1);
	t = MWC_A * mwc.auiState[mwc.uiState] + mwc.uiCarry;
	mwc.uiCarry = (unsigned int)(t >> 32);
	mwc.auiState[mwc.uiState] = (unsigned int)t;
	return RANDBL_32new(t);   // convert from 32 bit int to (0.0 , 1.0)
}
double DRanS_MWC8222(void)        /* returns a random double (-1.0 , 1.0) */
{
	UINT64 t;

	mwc.uiState = (mwc.uiState + 1) & (MWC_R - 1);
	t = MWC_A * mwc.auiState[mwc.uiState] + mwc.uiCarry;
	mwc.uiCarry = (unsigned int)(t >> 32);
	mwc.auiState[mwc.uiState] = (unsigned int)t;
	return RANDBLS_32new(t);   // convert from 32 bit int to (-1.0 , 1.0)
}
void VecIRan_MWC8222(unsigned int *auiRan, int cRan)
{
	UINT64 t;
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
void VecDRan_MWC8222(double *adRan, int cRan)
{
	UINT64 t;
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

/*------------------- normal random number generators ----------------------*/
static int s_cNormalInStore = 0;		     /* > 0 if a normal is in store */

static DRANFUN s_fnDRanu = DRan_R250;
static DRANFUN s_fnDRanus = DRanS_R250;
static IRANFUN s_fnIRanu = IRan_R250;
static IVECRANFUN s_fnVecIRanu = VecIRan_R250;
static DVECRANFUN s_fnVecDRanu = VecDRan_R250;
static RANSETSEEDFUN s_fnRanSetSeed = RanSetSeed_R250;

double  DRanU(void)     /* returns a random double (0.0 , 1.0) */
{
    return (*s_fnDRanu)();
}
double  DRanUS(void)    /* returns a random double (-1.0 , 1.0) */
{
    return (*s_fnDRanus)();
}
unsigned int IRanU(void)
{
    return (*s_fnIRanu)();
}
void RanVecIntU(unsigned int *auiRan, int cRan)
{
    (*s_fnVecIRanu)(auiRan, cRan);
}
void RanVecU(double *adRan, int cRan)
{
    (*s_fnVecDRanu)(adRan, cRan);
}
void    RanSetSeed(int *piSeed, int cSeed)
{
   	s_cNormalInStore = 0;
	(*s_fnRanSetSeed)(piSeed, cSeed);
}
void    RanSetRan(const char *sRan)
{
   	s_cNormalInStore = 0;
	if (strcmp(sRan, "MWC8222") == 0)
	{
		s_fnDRanu = DRan_MWC8222;
		s_fnDRanus = DRanS_MWC8222;
		s_fnIRanu = IRan_MWC8222;
		s_fnVecIRanu = VecIRan_MWC8222;
		s_fnVecDRanu = VecDRan_MWC8222;
		s_fnRanSetSeed = RanSetSeed_MWC8222;
	}
	else if (strcmp(sRan, "R250") == 0)
	{
		s_fnDRanu = DRan_R250;
		s_fnDRanus = DRanS_R250;
		s_fnIRanu = IRan_R250;
		s_fnVecIRanu = VecIRan_R250;
		s_fnVecDRanu = VecDRan_R250;
		s_fnRanSetSeed = RanSetSeed_R250;
	}
	else if (strcmp(sRan, "SHR3") == 0)
	{
		s_fnDRanu = DRan_SHR3;
		s_fnDRanus = DRanS_SHR3;
		s_fnIRanu = IRan_SHR3;
		s_fnVecIRanu = VecIRan_SHR3;
		s_fnVecDRanu = VecDRan_SHR3;
		s_fnRanSetSeed = RanSetSeed_SHR3;
	}
	else
	{
		s_fnDRanu = NULL;
		s_fnDRanus = NULL;
		s_fnIRanu = NULL;
		s_fnVecIRanu = NULL;
		s_fnVecDRanu = NULL;
		s_fnRanSetSeed = NULL;
	}
}
static unsigned int IRanUfromDRanU(void)
{
    return (unsigned int)(UINT_MAX * (*s_fnDRanu)());
}
static double DRanUfromIRanU(void)
{
    return RANDBL_32new( (*s_fnIRanu)() );
}
/*----------------- end normal random number generators --------------------*/

/*------------------------------ General Ziggurat --------------------------*/
#define ZIGNOR_C 128			       /* number of blocks */
#define ZIGNOR_R 3.442619855899	/* start of the right tail */
				   /* (R * phi(R) + Pr(X>=R)) * sqrt(2\pi) */
#define ZIGNOR_V 9.91256303526217e-3

/* s_adZigX holds coordinates, such that each rectangle has*/
/* same area; s_adZigR holds s_adZigX[i + 1] / s_adZigX[i] */
static double s_adZigX[ZIGNOR_C + 1], s_adZigR[ZIGNOR_C];

static unsigned int kn[128];
static double wn[128],fn[128];

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

#define ZIGNOR_STORE 256 * 4
static unsigned int s_auiZigTmp[ZIGNOR_STORE + ZIGNOR_STORE / 4];
static unsigned char s_auiZigBox[ZIGNOR_STORE];
// static double s_adZigRan[ZIGNOR_STORE + ZIGNOR_STORE / 4];
static unsigned int *i_adZigRan = &s_auiZigTmp[ZIGNOR_STORE / 4];
static int s_cZigStored = 0;

static int FillBufferZig(){
  int i, j, k;
  RanVecIntU(s_auiZigTmp, ZIGNOR_STORE + ZIGNOR_STORE / 4);
  for (j = k = 0; j < ZIGNOR_STORE; j += 4, ++k) {
    i = s_auiZigTmp[ZIGNOR_STORE + k];
    s_auiZigBox[j + 0] = i & 0x7F;
    s_auiZigBox[j + 1] = (i>>8) & 0x7F;
    s_auiZigBox[j + 2] = (i>>16) & 0x7F;
    s_auiZigBox[j + 3] = (i>>24) & 0x7F;
  }
  s_cZigStored = j ;
}

static double DRanNormalTail(double dMin, int iNegative)
{
  double x, y;
  do
  {
    if(s_cZigStored <= 1) FillBufferZig();
//     x = DRanU();
    x = RANDBL_32new(i_adZigRan[--s_cZigStored]);
    x = log(x) / dMin;
//     y = DRanU();
    y = RANDBL_32new(i_adZigRan[--s_cZigStored]);
    y = log(y);
    INSTRUMENT(zigused += 2;)
  } while (-2 * y < x * x);
  return iNegative ? x - dMin : dMin - x;
}

static void zigNorFastInit(int iC, double dR, double dV)  // faster, but with some deficiencies
{  const double m1 = 2147483648.0;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
   int i;

/* Set up tables for RNOR */
   q=vn/exp(-.5*dn*dn);
   kn[0]=(dn/q)*m1;
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(i=126;i>=1;i--)
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(dn/tn)*m1;
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }
}

static void zigNorInit(int iC, double dR, double dV)
{
	int i;	double f;
	
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
}

double  DRanNormalZig(void)
{
	unsigned int i;
	double x, u, f0, f1;
	
	for (;;)
	{
// 		u = 2 * DRanU() - 1;
		u = RANDBLS_32new(IRanU()) ;   // convert 32 bit int to (-1.0,1.0)
		i = IRanU() & 0x7F;
		/* first try the rectangular boxes */
		if (fabs(u) < s_adZigR[i])		 
			return u * s_adZigX[i];
		/* bottom box: sample from the tail */
		if (i == 0)						
			return DRanNormalTail(ZIGNOR_R, u < 0);
		/* is this a sample from the wedges? */
		x = u * s_adZigX[i];		   
		f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x) );
		f1 = exp(-0.5 * (s_adZigX[i+1] * s_adZigX[i+1] - x * x) );
      		if (f1 + DRanU() * (f0 - f1) < 1.0)
			return x;
	}
}

double  DRanNormalZigFast(void)  // faster, but with some deficiencies
{
	int iz, hz;
	double x, y;
	const double r = ZIGNOR_R;     /* The start of the right tail */

	for (;;)
	{
		hz = (int) IRanU() ;
		iz = hz & 0x7F ;
		if (abs(hz)<kn[iz]) { return(hz*wn[iz]) ; }
		
		x=hz*wn[iz];      /* iz==0, handles the base strip */
		if(iz==0) {
		  do{ x=-log(DRanU())*0.2904764; y=-log(DRanU());} while(y+y<x*x);	/* .2904764 is 1/r */ 
		  return (hz>0)? r+x : -r-x;
		}
		/* iz>0, handle the wedges of other strips */
		if( fn[iz]+DRanU()*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;
	}
}

double  DRanNormalZigVec(void)
{
	unsigned int i, j, k, direct;
	double x, u, y, f0, f1;
	INSTRUMENT(zigcalls++; ; direct = 1;)
	for (;;)
	{
		if (s_cZigStored <= 0) FillBufferZig();
		--s_cZigStored;

		u = RANDBLS_32new(i_adZigRan[s_cZigStored]);   // convert from 32 bit int to (-1.0 , 1.0) on the fly
		i = s_auiZigBox[s_cZigStored];
		INSTRUMENT(zigused += 2;)
		
		if (fabs(u) < s_adZigR[i]){		 /* first try the rectangular boxes */
		        INSTRUMENT(zigquick += direct;)
			return u * s_adZigX[i];
		}
		INSTRUMENT(direct = 0;)

		if (i == 0){						/* bottom box: sample from the tail */
			INSTRUMENT(zigtails++ ; )
// 			return DRanNormalTail(ZIGNOR_R, u < 0);
			do
			{
				if(s_cZigStored <= 1) FillBufferZig();
				x = RANDBL_32new(i_adZigRan[--s_cZigStored]);
				x = log(x) / ZIGNOR_R;
				y = RANDBL_32new(i_adZigRan[--s_cZigStored]);
				y = log(y);
				INSTRUMENT(zigused += 2;)
			} while (-2 * y < x * x);
			return (u < 0) ? x - ZIGNOR_R : ZIGNOR_R - x;
		}

		INSTRUMENT(zigwedge ++;)
		x = u * s_adZigX[i];		   /* is this a sample from the wedges? */
		f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x) );
		f1 = exp(-0.5 * (s_adZigX[i + 1] * s_adZigX[i + 1] - x * x) );
		INSTRUMENT(zigused++;)
		if(s_cZigStored <= 0) FillBufferZig();
		y = RANDBL_32new(i_adZigRan[--s_cZigStored]);
		if (f1 + y * (f0 - f1) < 1.0)  return x;
//       		if (f1 + DRanU() * (f0 - f1) < 1.0)  return x;
		INSTRUMENT(zigloops++;)
	}
}

static unsigned int u_auiZigTmp[ZIGNOR_STORE];
static int u_cZigStored = 0;
double  DRanNormalZigFastVec(void)  // faster, but with some deficiencies
{
	int iz, hz;
	double x, y;
	const double r = ZIGNOR_R;     /* The start of the right tail */
	INSTRUMENT(int direct = 1;)
	INSTRUMENT(zigcalls2++;)
	for (;;)
	{
		if (u_cZigStored <= 0)
		{
			RanVecIntU(u_auiZigTmp, ZIGNOR_STORE);
			u_cZigStored = ZIGNOR_STORE;
		}
		--u_cZigStored;
		hz = (int) u_auiZigTmp[u_cZigStored] ;
		INSTRUMENT(zigused2++;)
		iz = hz & 0x7F ;
		if (abs(hz)<kn[iz]) { 
		  INSTRUMENT(zigquick2+=direct ; )
		  return (hz*wn[iz]) ; 
		}
		INSTRUMENT(direct = 0;)
		
		x=hz*wn[iz];      /* iz==0, handles the base strip */
		if(iz==0) {
		  INSTRUMENT(zigtails2++;)
		  do{ x=-log(DRanU())*0.2904764; y=-log(DRanU()); INSTRUMENT(zigused2 += 2;) } while(y+y<x*x);	/* .2904764 is 1/r */ 
		  return (hz>0)? r+x : -r-x;
		}
		/* iz>0, handle the wedges of other strips */
		INSTRUMENT(zigwedge2++;)
		INSTRUMENT(zigused2++;)
		if( fn[iz]+DRanU()*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;
		INSTRUMENT(zigloops2++;)
	}
}

void  RanNormalSetSeedZig(int *piSeed, int cSeed)
{
	zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V);
	RanSetSeed(piSeed, cSeed);
}
void  RanNormalSetSeedZigFast(int *piSeed, int cSeed)
{
	zigNorFastInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V);
	RanSetSeed(piSeed, cSeed);
}
void  RanNormalSetSeedZigVec(int *piSeed, int cSeed)
{
	s_cZigStored = 0;
	RanNormalSetSeedZig(piSeed, cSeed);
}
void  RanNormalSetSeedZigFastVec(int *piSeed, int cSeed)
{
	s_cZigStored = 0;
	RanNormalSetSeedZigFast(piSeed, cSeed);
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
int NormalFromUniform(double gaussian[], int *ngauss, int uniform[], int nuni)
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
main(int argc, char **argv){
  unsigned int lr;
  int i, j;
  double t0, t1, rval;
  double MPI_Wtime() ;
  unsigned int ranbuf[1200000];
  double ranbuf2[1200000];
  int pos, neg, mask, postot, negtot;
  double dmax, dmin, avg;
  unsigned long long *idmax, *idmin ;
  unsigned int maxpos, maxneg;

  MPI_Init(&argc,&argv);
  for(i=0 ; i<1200000 ; i++) ranbuf[i] = 0;
  for(i=0 ; i<1200000 ; i++) ranbuf2[i] = 0.0;
  lr = (unsigned int)mrand48();
  lr = (unsigned int)mrand48();
  lr = (unsigned int)mrand48();
  lr = (unsigned int)mrand48();
  lr = (unsigned int)mrand48();
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

   RanSetSeed_MWC8222(&lr, 1);
   RanNormalSetSeedZig(r250.buffer, 250);   // initialize to values already there :-)
   RanNormalSetSeedZigFast(r250.buffer, 250);   // initialize to values already there :-)

//   printf("static unsigned int r250.buffer[250] = {\n");
  for (i=0 ; i<25 ; i++){
    for( j=0 ; j<10 ; j++){
      lr = (unsigned int) mrand48();
//       printf("0x%8.8x ",lr);
//       if(i==24 && j==9) printf("};") ; else  printf(",") ;
      if(r250.buffer[j+10*i] != lr) exit(1);
    }
//     printf("\n");
  }
  fprintf(stderr,"TEST 1 successful\n");

  FillBuffer_R250(); FillBuffer_R250(); FillBuffer_R250(); FillBuffer_R250(); FillBuffer_R250();

  for( i=0 ; i < 1000000 ; i++) lr = IRan_MWC8222();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRan_MWC8222();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random MWC8222 integer value = %6.3f\n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) lr = IRan_R250();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRan_R250();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random R250 integer value = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) lr = IRanU();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRanU();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random IRanU/R250 integer value = %6.3f \n",t1-t0);  // IRanU

  for( i=0 ; i < 1000000 ; i++) lr = IRan_SHR3();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRan_SHR3();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random SHR3 integer value = %6.3f \n",t1-t0);
  printf("\n");

  for( i=0 ; i < 1000000 ; i++) rval = DRan_MWC8222();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRan_MWC8222();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random MWC8222 double value = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) rval = DRan_R250();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRan_R250();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random R250 double value = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) rval = DRanU();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRanU();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random DRanU/R250 double value = %6.3f \n",t1-t0);  //  DRanU

  dmin = 0.0 ; dmax = 0.0;
  for( i=0 ; i < 100000000 ; i++) {
    rval = DRanNormalZigVec();
    avg = avg + rval ;
    dmin = (dmin < rval) ? dmin : rval ;
    dmax = (dmax > rval) ? dmax : rval ;
  }
  printf("dmin = %6.3f, dmax = %6.3f, avg = %10.7f\n",dmin,dmax,avg/i);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRanNormalZigVec();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random DRanNormalZigVec/R250 double value = %6.3f \n",t1-t0);  // DRanNormalZigVec

  dmin = 0.0 ; dmax = 0.0;
  for( i=0 ; i < 100000000 ; i++) {
    rval = DRanNormalZigFastVec();
    avg = avg + rval ;
    dmin = (dmin < rval) ? dmin : rval ;
    dmax = (dmax > rval) ? dmax : rval ;
  }
  printf("dmin = %6.3f, dmax = %6.3f, avg = %10.7f\n",dmin,dmax,avg/i);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRanNormalZigFastVec();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random DRanNormalZigFastVec/R250 double value = %6.3f \n",t1-t0);  // DRanNormalZigFastVec

  for( i=0 ; i < 1000000 ; i++) rval = DRan_SHR3();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRan_SHR3();
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random SHR3 double value = %6.3f \n",t1-t0);
  printf("\n");

  for( i=0 ; i < 10 ; i++) VecIRan_MWC8222(ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_MWC8222(ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if((int)ranbuf[i] > 0) pos++ ; else neg++ ;
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & 1024) pos++ ; else neg++ ;
  postot = 0 ; negtot = 0;
  for (j=0 ; j<100 ; j++) {
    VecIRan_MWC8222(ranbuf, 1000000) ;
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
  printf("\n");

  for( i=0 ; i < 10 ; i++) VecIRan_R250(ranbuf, 1000000) ;
//   for( i=0 ; i < 10 ; i++) VecIRan_R250_stream(&r250,ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_R250(ranbuf, 1000000) ;
//     VecIRan_R250_stream(&r250,ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if((int)ranbuf[i] > 0) pos++ ; else neg++ ;
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & 1) pos++ ; else neg++ ;
  postot = 0 ; negtot = 0;
  for (j=0 ; j<100 ; j++) {
    VecIRan_R250(ranbuf, 1000000) ;
//     VecIRan_R250_stream(&r250,ranbuf, 1000000) ;
    mask = 1 ;
    while (mask) {
      pos = 0 ; neg = 0 ; 
      for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & mask) pos++ ; else neg++  ; 
      postot += pos ; negtot += neg ;
      mask <<= 1 ;//  printf("%5d ",pos-neg) ;
    }
  }
//   printf("%d\n",postot-negtot);
  printf("time for 1E+3 x 1E+6 random R250 integer values = %6.3f , pos - neg = %d\n",t1-t0,postot-negtot);
  printf("\n");

  for( i=0 ; i < 10 ; i++) VecIRan_SHR3(ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_SHR3(ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if((int)ranbuf[i] > 0) pos++ ; else neg++ ;
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & 1) pos++ ; else neg++ ;
  postot = 0 ; negtot = 0;
  for (j=0 ; j<100 ; j++) {
    VecIRan_SHR3(ranbuf, 1000000) ;
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
  printf("\n");

  for( i=0 ; i < 10 ; i++) VecDRan_MWC8222(ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_MWC8222(ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random MWC8222 double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 10 ; i++) VecDRan_R250(ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_R250(ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random R250 double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 10 ; i++) VecDRan_SHR3(ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_SHR3(ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random SHR3 double values = %6.3f \n",t1-t0);
  printf("\n");

  for( i=0 ; i < 1000 ; i++) VecIRan_MWC8222(&ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_MWC8222(&ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random MWC8222 integer values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecIRan_R250(&ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_R250(&ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random R250 integer values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecIRan_SHR3(&ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_SHR3(&ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random SHR3 integer values = %6.3f \n",t1-t0);
  printf("\n");

  for( i=0 ; i < 1000 ; i++) VecDRan_MWC8222(&ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_MWC8222(&ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random MWC8222 double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecDRan_R250(&ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_R250(&ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random R250 double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecDRan_SHR3(&ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_SHR3(&ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random SHR3 double values = %6.3f \n",t1-t0);

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
}
#endif