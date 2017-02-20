#include <stdio.h>

atomic_add(int *what,int n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   addl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

atomic_and(int *what,int n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   andl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

atomic_or(int *what,int n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   orl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

atomic_xor(int *what,int n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   xorl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}
#if defined(SELF_TEST)
main(){
 int what = 0;
 int r;
 printf("what = %d\n",what);
#pragma omp parallel for
 for(r=0 ; r<1200 ; r++) {
 atomic_add(&what,r);
//  atomic_add(&what,2);
//  atomic_xor(&what,3);
//  atomic_xor(&what,3);
//  atomic_add(&what,-4);
//  atomic_add(&what,3);
// printf("%d %3d what = %3d\n",omp_get_thread_num(),r,what);
 }
 printf("what = %d\n",what);
}
#endif
