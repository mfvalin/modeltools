#include <stdio.h>

/* default values for missing values */
static float float_missing_val=-1.0E+38;          /* very large float */
static double double_missing_val=-1.0E+38;       /* very large double */
static int int_missing_val=0x80000000;              /* largest negative 32 bit integer */
static short short_missing_val=0x8000;              /* largest negative 16 bit integer */
static signed char byte_missing_val=0x80;           /* largest negative  8 bit integer */
static unsigned int uint_missing_val=0xFFFFFFFF;    /* largest 32 bit unsigned integer */
static unsigned short ushort_missing_val=0xFFFF;    /* largest 16 bit unsigned integer */
static unsigned char ubyte_missing_val=0xFF;       /* largest  8 bit unsigned integer */

void set_plugin_missing_value_flags(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
                                      signed char *b, unsigned char *ub)
{
  float_missing_val  = *f ;
  int_missing_val    = *i ;
  uint_missing_val   = *ui ;
  double_missing_val = *d ;
  short_missing_val  = *s ;
  ushort_missing_val = *us ;
  byte_missing_val   = *b ;
  ubyte_missing_val  = *ub ;
  fprintf(stderr,"INFO: set_plugin_missing_value_flags\n");
  fprintf(stderr,"++float=%g, int=%d, uint=%u, double=%lg, short=%hd, ushort=%hu, byte=%hhd, ubyte=%hhu\n",*f,*i,*ui,*d,*s,*us,*b,*ub);
}

/* demo plugin encoders are simple copy routines */
int float_encode(float *z, float *z2, int n, int nbits) { 
  while(n-- > 0) *z++ = *z2++;
  return(0) ;
} ;
int double_encode(double *z, double *z2, int n, int nbits) { 
  while(n-- > 0) *z++ = *z2++;
  return(0) ; 
} ;
int int_encode(int *z, int *z2, int n, int nbits) { 
  while(n-- > 0) *z++ = *z2++;
  return(0) ; 
} ;
int short_encode(short *z, short *z2, int n, int nbits) { 
  while(n-- > 0) *z++ = *z2++;
  return(0) ; 
} ;
int byte_encode(signed char *z, signed char *z2, int n, int nbits) { 
  while(n-- > 0) *z++ = *z2++;
  return(0) ; 
} ;
int uint_encode(unsigned int *z, unsigned int *z2, int n, int nbits) { 
  while(n-- > 0) *z++ = *z2++;
  return(0) ; 
} ;
int ushort_encode(unsigned short *z, unsigned short *z2, int n, int nbits) { 
  while(n-- > 0) *z++ = *z2++;
  return(0) ; 
} ;
int ubyte_encode(unsigned char *z, unsigned char *z2, int n, int nbits) { 
  while(n-- > 0) *z++ = *z2++;
  return(0) ;
} ;

/* demo plugin decoders are do nothing routines */
void float_decode(float *z, int n) { return ; } ;
void double_decode(double *z, int n) { return ; } ;
void int_decode(int *z, int n) { return ; } ;
void short_decode(short *z, int n) { return ; } ;
void byte_decode(signed char *z, int n) { return ; } ;
void uint_decode(unsigned int *z, int n) { return ; } ;
void ushort_decode(unsigned short *z, int n) { return ; } ;
void ubyte_decode(unsigned char *z, int n) { return ; } ;
