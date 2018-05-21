#if defined(__amd64__)

#if ! defined(INLINE)
  #define INLINE static inline
#endif

INLINE
atomic_add_32(int32_t *what,int32_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   addl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

INLINE
atomic_and_32(uint32_t *what,uint32_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   andl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

INLINE
atomic_or_32(uint32_t *what,uint32_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   orl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

INLINE
atomic_xor_32(uint32_t *what,uint32_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   xorl %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

INLINE
atomic_add_64(int64_t *what,int64_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   addq %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

INLINE
atomic_and_64(uint64_t *what,uint64_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   andq %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

INLINE
atomic_or_64(uint64_t *what,uint64_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   orq %1,%0 ;\n"
      : "=m" (*what)
      : "ir" (n), "m" (*what)
      :
      );
}

INLINE
atomic_xor_64(uint64_t *what,uint64_t n) {
  __asm__ __volatile__(
      "   lock       ;\n"
      "   xorq %1,%0 ;\n"
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

INLINE
int atomic_compare_and_swap_64(uint64_t *thevalue, uint64_t expected, uint64_t newvalue) {
        uint64_t prev;
        asm volatile("lock;\n"
                "\tcmpxchgq %1, %2;"
                : "=a"(prev)
                : "q"(newvalue), "m"(*thevalue), "a"(expected)
                : "memory");
        return prev == expected;
}

INLINE
int64_t atomic_fetch_and_add_64(int64_t* ptr, int64_t value) {
        int64_t previous;
        asm volatile ("lock;\n"
                "\txaddq %1, %2;"
                : "=r"(previous)
                : "0"(value), "m"(*ptr)
                : "memory");
        return previous;
}

INLINE
int atomic_compare_and_swap_32(uint32_t *thevalue, uint32_t expected, uint32_t newvalue) {
        uint32_t prev;
        asm volatile("lock;\n"
                "\tcmpxchgl %1, %2;"
                : "=a"(prev)
                : "q"(newvalue), "m"(*thevalue), "a"(expected)
                : "memory");
        return prev == expected;
}

INLINE
int32_t atomic_fetch_and_add_32(int32_t* ptr, int32_t value) {
        int32_t previous;
        asm volatile ("lock;\n"
                "\txaddl %1, %2;"
                : "=r"(previous)
                : "0"(value), "m"(*ptr)
                : "memory");
        return previous;
}

#endif