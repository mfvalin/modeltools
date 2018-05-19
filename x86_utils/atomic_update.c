#include <stdio.h>
#include <stdint.h>

#if defined(__amd64__)

atomic_add_32(int32_t *what,int32_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   addl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

atomic_and_32(uint32_t *what,uint32_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   andl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

atomic_or_32(uint32_t *what,uint32_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   orl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

atomic_xor_32(uint32_t *what,uint32_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   xorl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

// function cas(p : pointer to int, old : int, new : int) returns bool {
//     if *p != old {
//         return false
//     }
//     *p = new
//     return true
// }

int atomic_compare_and_swap_64(uint64_t *thevalue, uint64_t expected, uint64_t newvalue) {
	uint64_t prev;
	asm volatile("lock;\n"
		"\tcmpxchgq %1, %2;"
		: "=a"(prev)
		: "q"(newvalue), "m"(*thevalue), "a"(expected)
		: "memory");
	return prev == expected;
}

int64_t atomic_fetch_and_add_64(int64_t* ptr, int64_t value) {
	int64_t previous;
	// xadd r0, r1   - Exchange r0 and r1 loading sum into r1
	asm volatile ("lock;\n"
		"\txaddq %1, %2;"
		: "=r"(previous)
		: "0"(value), "m"(*ptr)
		: "memory");
	return previous;
}

int atomic_compare_and_swap_32(uint32_t *thevalue, uint32_t expected, uint32_t newvalue) {
	uint32_t prev;
	asm volatile("lock;\n"
		"\tcmpxchgl %1, %2;"
		: "=a"(prev)
		: "q"(newvalue), "m"(*thevalue), "a"(expected)
		: "memory");
	return prev == expected;
}

int32_t atomic_fetch_and_add_32(int32_t* ptr, int32_t value) {
	int32_t previous;
	// xadd r0, r1   - Exchange r0 and r1 loading sum into r1
	asm volatile ("lock;\n"
		"\txaddl %1, %2;"
		: "=r"(previous)
		: "0"(value), "m"(*ptr)
		: "memory");
	return previous;
}

uint32_t atomic_compare_exchange(uint32_t *what, uint32_t expected, uint32_t desired)
{
   uint32_t previous;
    asm volatile("lock cmpxchgl %2, %1"
                 : "=a"(previous), "+m"(*what)
                 : "q"(desired), "0"(expected));
    return previous;
}

int atomic_fetch_and_add(int* what, int n)
{
  __asm__ volatile(
      "   lock       ;\n"
      "   xaddl %0, %1"
      : "+r" (n), "+m" (*what) // input+output
      : // No input-only
      : "memory"
    );
    return n;
}

uint32_t release_lock_32(uint32_t *lock, uint32_t id){  // release lock if i am the owner (id)
  if(*lock == ~id) *lock = 0;  // only if i am the owner of the lock
  return *lock ;               // return lock value, 0 means i was the owner of the lock
}

uint32_t acquire_lock_32(uint32_t *lock, uint32_t id){ // id MUST NOT BE 0xFFFFFFFF (all ones)
  while( ! atomic_compare_and_swap_32(lock, 0, ~id ) ) ; // spin loop until lock is free (value == 0)
  return *lock;
}

uint32_t testlock_32(int32_t *lock){
  return *lock ; // return lock value, 0 means free
}

#endif

#if defined(SELF_TEST)

uint64_t rdtsc(void) {   // version rapide "out of order"
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

main(){
 int what = 0;
 int r, n;
 int i = 1;
 uint64_t t1, t2; 
 printf("what = %d\n",what);
 t1 = rdtsc();
#pragma omp parallel for
 for(r=0 ; r<1000000 ; r++) {
   n=atomic_fetch_and_add(&what,i);
//  atomic_add(&what,2);
//  atomic_xor(&what,3);
//  atomic_xor(&what,3);
//  atomic_add(&what,-4);
//  atomic_add(&what,3);
// printf("%d %3d what = %3d\n",omp_get_thread_num(),r,what);
 }
 t2 = rdtsc();
 printf("what = %d, n= %d, time = %ld\n",what,n,t2-t1);

  what = 0;
  t1 = rdtsc();
  for(i=1 ; i < 2000000 ; i++) {n=atomic_compare_exchange( &what , i-1, i); }
  t2 = rdtsc();
  printf("%10d %10d %10ld\n",i,what,t2-t1);

  what=-1;
  t1 = rdtsc();
  for(i=1 ; i < 2000000 ; i++) {n=atomic_compare_exchange( &what , i-1, i); }
  t2 = rdtsc();
  printf("%10d %10d %10ld\n",i,what,t2-t1);
}
#endif
