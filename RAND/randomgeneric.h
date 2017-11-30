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

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <limits.h>
// #include <unistd.h>
#include <math.h>

/*------------------------ multiplier constants --------------------------*/
static const double INVM31   =  4.65661287307739257812e-010 ;        /* 1.0 / 2^31 */
static const double INVM32   =  2.32830643653869628906e-010 ;        /* 1.0 / 2^32 */
static const double INVM33   =  1.1641532182693481445e-10   ;        /* 1.0 / 2^33 */
static const double INVM48   =  3.55271367880050092936e-015 ;        /* 1.0 / 2^48 */
static const double INVM52   =  2.22044604925031308085e-016 ;        /* 1.0 / 2^52 */
static const double INVM63   =  1.084202172485504434e-19    ;        /* 1.0 / 2^63 */
static const double INVM64   =  5.42101086242752217e-20     ;        /* 1.0 / 2^64 */
static const double INVM65   =  2.710505431213761085e-20    ;        /* 1.0 / 2^65 */

/*------------------------    scaling macros    --------------------------*/
#define CVTDBL_32(i1)       ((int)(i1)  * INVM32 + (0.5 + INVM33))
#define CVTDBL_64(i1)       ((long)(i1) * INVM64 + (0.5 + INVM65))
#define CVTDBLS_32(i1)      ((int)(i1)  * INVM31 + (INVM32))
#define CVTDBLS_64(i1)      ((long)(i1) * INVM63 + (INVM64))
#define CVTDBL_48(i1, i2)   ((int)(i1)  * INVM32 + (0.5 + INVM48 / 2) + (int)((i2) & 0xFFFF)  * INVM48)
#define CVTDBL_52(i1, i2)   ((int)(i1)  * INVM32 + (0.5 + INVM52 / 2) + (int)((i2) & 0xFFFFF) * INVM52)

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
#include <randomfunctions.h>
