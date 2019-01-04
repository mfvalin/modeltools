// generic functions
void RanSetSeed_generic_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  ;
unsigned int IRan_generic_stream(generic_state *stream)       ;
double DRan_generic_stream(generic_state *stream)       ;
double DRanS_generic_stream(generic_state *stream)       ;
void VecIRan_generic_stream(generic_state *stream, unsigned int *ranbuf, int n)       ;
void VecDRan_generic_stream(generic_state *stream, double *ranbuf, int n)       ;
void VecDRanS_generic_stream(generic_state *stream, double *ranbuf, int n)       ;

// R250 random number generator, from S.Kirkpatrick and E.  Stoll (state = 250 32 bit integers)
void RanSetSeed_R250_static(unsigned int *piSeed, int cSeed)  ;
unsigned int IRan_R250_static()	  ;
void VecIRan_R250_static(unsigned int *ranbuf, int n)  ;
void RanSetSeed_R250_stream(void *stream, unsigned int *piSeed, int cSeed)  ;
void *Ran_R250_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   ;
unsigned int IRan_R250_stream(void *stream)       ;
double DRan_R250_stream(void *stream)	  ;
double DRanS_R250_stream(void *stream)	  ;
void VecIRan_R250_stream(void *stream, unsigned int *ranbuf, int n)  ;
void VecDRan_R250_stream(void *stream, double *ranbuf, int n)  ;
void VecDRanS_R250_stream(void *stream, double *ranbuf, int n)  ;

// SHR3 3-shift-register generator (state = 32 bit integer)
void *Ran_SHR3_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   ;
void RanSetSeed_SHR3_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  ;
unsigned int IRan_SHR3_stream(generic_state *stream)	  ;
double DRan_SHR3_stream(generic_state *stream)	  ;
double DRanS_SHR3_stream(generic_state *stream)	  ;
void VecIRan_SHR3_stream(generic_state *stream, unsigned int *ranbuf, int n)  ;
void VecDRan_SHR3_stream(generic_state *stream, double *ranbuf, int n)  ;
void VecDRanS_SHR3_stream(generic_state *stream, double *ranbuf, int n)  ;

// MT19937 (Mersenne twister) generator (state = 1248 32 bit integers)
generic_state *Ran_MT19937_new_stream(generic_state *clone_in, unsigned int *piSeed, int cSeed)   ;
void RanSetSeed_MT19937_stream(generic_state *stream, unsigned int *piSeed, int cSeed)    ;
unsigned int IRan_MT19937_stream(generic_state *stream)       ;
double DRan_MT19937_stream(generic_state *stream)       ;
double DRanS_MT19937_stream(generic_state *stream)     ;
void VecIRan_MT19937_stream(generic_state *stream, uint32_t *ranbuf, int n)  ;
void VecDRan_MT19937_stream(generic_state *stream, double *ranbuf, int n)  ;
void VecDRanS_MT19937_stream(generic_state *stream, double *ranbuf, int n)  ;

// xorshift generator with 128 bit state
void *Ran_XSR128_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   ;
void RanSetSeed_XSR128_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  ;
unsigned int IRan_XSR128_stream(generic_state *stream)	  ;
double DRan_XSR128_stream(generic_state *stream)	  ;
double DRanS_XSR128_stream(generic_state *stream)	  ;
void VecIRan_XSR128_stream(generic_state *stream, unsigned int *ranbuf, int n)  ;
void VecDRan_XSR128_stream(generic_state *stream, double *ranbuf, int n)  ;
void VecDRanS_XSR128_stream(generic_state *stream, double *ranbuf, int n)  ;

// xorshiftrotate generator with 128 bit state
void *Ran_XSR128R_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   ;
void RanSetSeed_XSR128R_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  ;
unsigned int IRan_XSR128R_stream(generic_state *stream)	  ;
double DRan_XSR128R_stream(generic_state *stream)	  ;
double DRanS_XSR128R_stream(generic_state *stream)	  ;
void VecIRan_XSR128R_stream(generic_state *stream, unsigned int *ranbuf, int n)  ;
void VecDRan_XSR128R_stream(generic_state *stream, double *ranbuf, int n)  ;
void VecDRanS_XSR128R_stream(generic_state *stream, double *ranbuf, int n)  ;

// gaussian generator using the "ziggurat" method (uses one of the previous integer generators)
// 128 or 256 boxes (use internal buffer for 320 32 bit integers)
void RanNormalZigSetSeed(void *stream, void *values, int nvalues)  ;
double DRan_NormalZig_stream(void *stream)  ;
double D64Ran_NormalZig_stream(void *stream)  ;
double D64RanNormalFun(void *stream)  ;
