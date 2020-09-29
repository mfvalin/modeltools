
#define QUANT_LINEAR 15
#define QUANT_EXPON0  0
#define QUANT_EXPON1  1
#define QUANT_EXPON2  2
#define QUANT_EXPON3  3
#define QUANT_EXPON4  4
#define QUANT_EXPON5  5

#if defined(FORTRAN_CODE)
  type, bind(C) :: pack_info
    type(C_PTR) :: f
  end type
#else
  typedef struct{
    void *f;
  } pack_info
#endif

#if defined(FORTRAN_CODE)
#else

// IEEE 754 components (float)
#define IEEE32_EXP(a) (( ((uint32_t) (a) >> 23) & 0xFF ) - 127 )
#define IEEE32_SGN(a) ( ((uint32_t) (a) >> 31) )
#define IEEE32_MANT(a) ( ((uint32_t) (a) & 0x7FFFFF) )
// rebuild IEEE 754 float from components
#define IEEE32(sign,exp,mantissa) ( (sign << 31) | ( ((exp+127) & 0xFF) << 23) | ((mantissa) & 0x7FFFFF) )

#endif
