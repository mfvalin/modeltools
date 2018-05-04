#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#if defined(REPLACE_MEMALIGN)
//Convenience macro for memalign, the linux API
#define memalign(align, size) AlignedMalloc(align, size)
#endif

/**
* AlignedMalloc takes in the requested alignment and size
* We will call malloc with extra bytes for our header and the offset
* required to guarantee the desired alignment.
*/

#if defined(WITH_POINTER)

void *AlignedMalloc(size_t align, size_t size)
{
  void *mem = NULL;
  void *ptr = NULL;

  //We want it to be a power of two since align_up operates on powers of two
  assert((align & (align - 1)) == 0);

  // We require whatever user asked for PLUS space for a pointer
  // PLUS space to align pointer as per alignment requirement
  mem = malloc(size + sizeof(void*) + (align - 1));
  if(mem == NULL) return NULL ;
  
  // Location that we will return to user
  // This has space *behind* it for a pointer and is aligned
  // as per requirement
  ptr = (void**)((uintptr_t) (mem + (align - 1) + sizeof(void*)) & ~(align - 1));
  
  // Sneakily store address returned by malloc *behind* user pointer
  // void** cast is cause void* pointer cannot be decremented, cause
  // compiler has no idea "how many" bytes to decrement by
  ((void **) ptr)[-1] = mem;
#if defined(DEBUG)  
  printf("raw = %p, aligned to %ld = %p, delta = %ld\n",mem,align,ptr,ptr-mem);
#endif
  // Return user pointer
  return ptr;
}

#else

/**
* Simple macro for making sure memory addresses are aligned
* to the nearest power of two
*/
#ifndef align_up
#define align_up(num, align) \
	(((num) + ((align) - 1)) & ~((align) - 1))
#endif

//Number of bytes we're using for storing the aligned pointer offset (assumed to be < 64K)
typedef uint16_t offset_t;
#define PTR_OFFSET_SZ sizeof(offset_t)

void * AlignedMalloc(size_t align, size_t size)
{
  void * ptr = NULL;
  void *p = NULL;

  //We want it to be a power of two since align_up operates on powers of two
  assert((align & (align - 1)) == 0);

  if(align && size)
  {
    // We know we have to fit an offset value
    // We also allocate extra bytes to ensure we can meet the alignment
    uint32_t hdr_size = PTR_OFFSET_SZ + (align - 1);
    p = malloc(size + hdr_size);

    if(p)
    {
      // Add the offset size to malloc's pointer (we will always store that)
      // Then align the resulting value to the arget alignment
      ptr = (void *) align_up(((uintptr_t)p + PTR_OFFSET_SZ), align);

      //Calculate the offset and store it behind our aligned pointer
      *((offset_t *)ptr - 1) = (offset_t)((uintptr_t)ptr - (uintptr_t)p);

    } // else NULL, could not malloc
  } //else NULL, invalid arguments

#if defined(DEBUG)  
  printf("raw = %p, aligned to %ld = %p, delta = %ld\n",p,align,ptr,ptr-p);
#endif
  return ptr;
}

#endif

/**
* AlignedFree works like free(), but we work backwards from the returned
* pointer to find the correct offset and pointer location to return to free()
* Note that it is VERY BAD to call free() on an AlignedMalloc() pointer.
*/
#if defined(WITH_POINTER)

void AlignedFree(void *ptr)
{
  assert(ptr != NULL);
  // Sneak *behind* user pointer to find address returned by malloc
  // Use that address to free
#if defined(DEBUG)  
  printf("freeing raw = %p, aligned = %p\n",((void**) ptr)[-1],ptr);
#endif
  free(((void**) ptr)[-1]);
}

#else

void AlignedFree(void * ptr)
{
  assert(ptr != NULL);

  // Walk backwards from the passed-in pointer to get the pointer offset
  // We convert to an offset_t pointer and rely on pointer math to get the data
  offset_t offset = *((offset_t *)ptr - 1);

  // Once we have the offset, we can get our original pointer and call free
  void * p = (void *)((uint8_t *)ptr - offset);
#if defined(DEBUG)  
  printf("freeing raw = %p, aligned = %p\n",p,ptr);
#endif
  free(p);
}

#endif

/**
* Example Usage
*/

#ifdef SELF_TEST
int main(void)
{
    void * p = malloc(103);
    void * q = malloc(1000);
    void * r = malloc(7);

    void * x = AlignedMalloc(128, 100);
    void * y = AlignedMalloc(32, 1035);
    void * z = AlignedMalloc(64, 8);

    printf("Raw malloc pointers, no alignment enforced:\n");
    printf("\t%p, %p, %p\n", p, q, r);
    printf("\tNote: you may see 4-8 byte alignment on host PC\n");
    printf("aligned to 128: %p\n", x);
    printf("aligned to 32: %p\n", y);
    printf("aligned to 64: %p\n", z);

    AlignedFree(x), x = NULL;
    AlignedFree(y), y = NULL;
    AlignedFree(z), z = NULL;

    free(p), p = NULL;
    free(q), q = NULL;
    free(r), r = NULL;

    return 0;
}
#endif
