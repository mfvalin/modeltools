
#include <xmmintrin.h>

#pragma weak set_ftz_daz__=set_ftz_daz
void set_ftz_daz__();
#pragma weak set_ftz_daz_=set_ftz_daz
void set_ftz_daz_();
void set_ftz_daz()
{
  int mxcsr = _mm_getcsr ();
  mxcsr |= (1<<15) | (1<<11);
  mxcsr |= (1<<6);
  _mm_setcsr (mxcsr);
}
