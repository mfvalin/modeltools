
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <limits.h>
// #include <unistd.h>
#include <math.h>

/*------------------------ multiplier constants --------------------------*/
#define M_RAN_INVM32    2.32830643653869628906e-010         /* 1.0 / 2^32 */
#define M_RAN_INVM31    4.65661287307739257812e-010         /* 1.0 / 2^31 */
#define M_RAN_INVM48    3.55271367880050092936e-015         /* 1.0 / 2^48 */
#define M_RAN_INVM52    2.22044604925031308085e-016         /* 1.0 / 2^52 */

/*------------------------    scaling macros    --------------------------*/
#define RANDBL_32new(iRan1)          ((int)(iRan1) * M_RAN_INVM32 + (0.5 + M_RAN_INVM32 / 2))
#define RANDBLS_32new(iRan1)         ((int)(iRan1) * M_RAN_INVM31 + (M_RAN_INVM31 / 2))
#define RANDBL_48new(iRan1, iRan2)   ((int)(iRan1) * M_RAN_INVM32 + (0.5 + M_RAN_INVM48 / 2) + (int)((iRan2) & 0x0000FFFF) * M_RAN_INVM48)
#define RANDBL_52new(iRan1, iRan2)   ((int)(iRan1) * M_RAN_INVM32 + (0.5 + M_RAN_INVM52 / 2) + (int)((iRan2) & 0x000FFFFF) * M_RAN_INVM52)

/* plug-in RNG */
typedef double          ( * DRANFUN)(void *);
typedef double          ( * DRANSFUN)(void *);
typedef unsigned int	( * IRANFUN)(void *);
typedef void   		( * IVECRANFUN)(void *, unsigned int *, int);
typedef void            ( * DVECRANFUN)(void *, double *, int);
typedef void            ( * DVECSRANFUN)(void *, double *, int);
typedef void            ( * RANSETSEEDFUN)(void *, unsigned int *, int);
typedef void            ( * REFILLBUFFUN)(void *);

/*------------------------ stream control structures --------------------------*/
typedef struct{                // mimic Fortran derived type (wrapped pointer to type it)
  void *p;
} statep;

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
} generic_state;               // generic part, identical at start of all stream control structures

