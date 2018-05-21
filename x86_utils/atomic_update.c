#include <stdint.h>

#if defined(__amd64__)

#define INLINE static inline
#include <atomic_update.h>

uint32_t atomic_compare_exchange(uint32_t *what, uint32_t expected, uint32_t desired)
{
//    uint32_t previous;
//     asm volatile("lock cmpxchgl %2, %1"
//                  : "=a"(previous), "+m"(*what)
//                  : "q"(desired), "0"(expected));
//     return previous;
     return atomic_compare_and_swap_32(what, expected, desired);
}

int atomic_fetch_and_add(int* what, int n)
{
//   __asm__ volatile(
//       "   lock       ;\n"
//       "   xaddl %0, %1"
//       : "+r" (n), "+m" (*what) // input+output
//       : // No input-only
//       : "memory"
//     );
//     return n;
     return atomic_fetch_and_add_32(what, n);
}

uint32_t release_lock_32(uint32_t *lock, uint32_t id){  // release lock if i am the owner (id)
  if(*lock == ~id) *lock = 0;  // only if i am the owner of the lock
  return *lock ;               // return lock value, 0 means i was the owner of the lock
}

uint32_t acquire_lock_32(uint32_t *lock, uint32_t id){ // id MUST NOT BE 0xFFFFFFFF (all ones)
  while( ! atomic_compare_and_swap_32(lock, 0, ~id ) ) ; // spin loop until lock is free (value == 0)
  return *lock;
}

uint32_t test_lock_32(int32_t *lock){
  return *lock ; // return lock value, 0 means free
}

#endif

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>

static uint64_t time0=0;

static inline uint64_t rdtsc(void) {   // version rapide "out of order"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  return time0++;
#endif
}

static int32_t lockvar = 0;
main(){
  int what = 0;
  int r, n;
  int i = 1;
  uint64_t t1, t2; 
  printf("what = %d\n",what);
  t1 = rdtsc();
#pragma omp parallel for
  for(r=0 ; r<1000000 ; r++) {
//     if( acquire_lock_32(&lockvar,1) == 0) exit(1) ;
    n=atomic_fetch_and_add_32(&what,1);
//     what++;
//     atomic_xor_32(&what,3);
//     atomic_xor_32(&what,3);
//    if( release_lock_32(&lockvar,1) != 0) exit(1);
//  atomic_add_32(&what,2);
//  atomic_add_32(&what,-1);
//  atomic_add(&what,-4);
//  atomic_add(&what,3);
// printf("%d %3d what = %3d\n",omp_get_thread_num(),r,what);
  }
  t2 = rdtsc();
  printf("what = %d, expected = %d, n = %d, time = %ld cycles\n",what,1000000,n,(t2-t1)/1000000);

  what = 0;
  t1 = rdtsc();
  for(i=1 ; i < 2000000 ; i++) {n=atomic_compare_exchange( &what , i-1, i); }
  t2 = rdtsc();
  printf("what = %d, expected = %d, time = %10ld cycles\n",what,i-1,(t2-t1)/2000000);

  what=-1;
  t1 = rdtsc();
  for(i=1 ; i < 2000000 ; i++) {n=atomic_compare_exchange( &what , i-1, i); }
  t2 = rdtsc();
  printf("what = %d, expected = %d, time = %10ld cycles\n",what,-1,(t2-t1)/2000000);
}
#endif
