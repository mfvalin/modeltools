#include <stdio.h>

void myfunction()
{
  printf("#define f77name(a) a\n");
}
void myfunction_()
{
  printf("#define f77name(a) a##_\n");
}
void myfunction__()
{
  printf("#define f77name(a) a##__\n");
}
void my_function()
{
  printf("#define f77_name(a) a\n");
}
void my_function_()
{
  printf("#define f77_name(a) a##_\n");
}
void my_function__()
{
  printf("#define f77_name(a) a##__\n");
}
