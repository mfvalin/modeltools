/* capability flags used by cpu_type.c */
#define FLAG_SSE    1
#define FLAG_SSE2   2
#define FLAG_SSE3   4
#define FLAG_SSE4   8
#define FLAG_AVX   16
#define FLAG_AVX2  32
#define FLAG_FMA   64
#define FLAG_BMI  128

#if defined(IN_FORTRAN_CODE)
  interface
    function processor_has_feature(feature) result(status) bind(C,name='Cpu_has_feature')
      import :: C_INT
      integer(C_INT), intent(IN), value :: feature
      integer :: status
    end function processor_has_feature
    function get_processor_apicid() result(id)  bind(C,name='Get_cpu_apicid')
      import :: C_INT
      integer(C_INT) :: id
    end function get_processor_apicid
  end interface
#else
  int Cpu_has_feature(int feature);
  int Get_cpu_apicid();
#endif