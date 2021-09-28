//     functions for C and FORTRAN programming
//     Copyright (C) 2019  Recherche en Prevision Numerique
// 
//     This software is free software; you can redistribute it and/or
//     modify it under the terms of the GNU Lesser General Public
//     License as published by the Free Software Foundation,
//     version 2.1 of the License.
// 
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//     Lesser General Public License for more details.

// set of inline functions used to implement locks

// Authors :: M. Valin,   RPN, 2021
//            V. Magnoux, RPN, 2021

#include <stdint.h>

void *allocate_safe_shared_memory(int32_t *sid, size_t size) ;

void ReleaseLock(volatile int32_t *lock, int32_t fence) ;
void ReleaseIdLock(volatile int32_t *lock, int32_t id, int32_t fence) ;

int32_t TryReleaseLock(volatile int32_t *lock, int32_t fence) ;
int32_t TryReleaseIdLock(volatile int32_t *lock, int32_t id, int32_t fence) ;

void AcquireLock(volatile int32_t *lock, int32_t fence) ;
void AcquireIdLock(volatile int32_t *lock, int32_t id, int32_t fence) ;

int32_t TryAcquireLock(volatile int32_t *lock, int32_t fence) ;
int32_t TryAcquireIdLock(volatile int32_t *lock, int32_t id, int32_t fence) ;

void BasicNodeBarrier(volatile int32_t *flag, volatile int32_t *count, int32_t maxcount) ;
int32_t LockOwner(volatile int32_t *lock) ;

void *setup_locks_and_barriers(void *p, uint32_t size, uint32_t maxitems) ;
uint32_t simple_node_barrier(int32_t id, int32_t maxcount) ;
uint32_t node_barrier_multi(int32_t id, int32_t maxcount) ;

static inline void full_memory_fence(){
  __asm__ volatile("mfence": : :"memory") ;
}

static inline int32_t try_acquire_idlock(volatile int32_t *lock, int32_t id) {
  __asm__ volatile("": : :"memory");
  if (__sync_val_compare_and_swap(lock, 0, (id+1)) != 0) return 0;
  full_memory_fence() ;
  return 1;
}

static inline int32_t try_acquire_lock(volatile int32_t *lock) {
  return try_acquire_idlock(lock, 1) ;
}

static inline void acquire_idlock(volatile int32_t *lock, int32_t id){
  full_memory_fence() ;
  while(__sync_val_compare_and_swap(lock, 0, (id+1)) != 0) ;
  __asm__ volatile("": : :"memory") ;
}

static inline void acquire_lock(volatile int32_t *lock){   // no id, use 1
  acquire_idlock(lock, 1) ;
}

// this will not deadlock if an attempt is made to release a lock with the wrong id
static inline int32_t try_release_idlock(volatile int32_t *lock, int32_t id){
  full_memory_fence() ;
  if(__sync_val_compare_and_swap(lock, (id+1), 0) != (id+1)) return 0 ;  // failure, not lock owner
  __asm__ volatile("": : :"memory") ;
  return 1 ;                                                             // success
}

// this will deadlock if attempt is made to release a lock with the wrong id
static inline void release_idlock(volatile int32_t *lock, int32_t id){
  if( *lock <= 0) return ;                                             // NOOP, was not even locked
  full_memory_fence() ;
  while(__sync_val_compare_and_swap(lock, (id+1), 0) != (id+1)) ;
  __asm__ volatile("": : :"memory") ;
}

// this will deadlock if attempt is made to release a lock with id other than 1
static inline void release_lock(volatile int32_t *lock){   // no id, use 1
  release_idlock(lock, 1) ;
}

// this will not deadlock if an attempt is made to release a lock with id other than 1
static inline int32_t try_release_lock(volatile int32_t *lock){   // no id, use 1
  return try_release_idlock(lock, 1) ;
}

static inline int32_t lock_owner(volatile int32_t *lock){   // find ID of lock owner
  return  (*lock - 1) ;  // < 0 if lock is free
}

static inline int32_t test_idlock(volatile int32_t *lock, int32_t id){
  return (*lock == (id+1) );   // true if locked with id
}

static inline int32_t test_lock(volatile int32_t *lock){
  return (*lock != 0 );   // true if locked with any id
}

static inline void reset_lock(volatile int32_t *lock){   // forcefully reset lock
  *lock = 0 ;
}

// count MUST be equal to zero when entering this code
// flag MUST equal to be zero or one
// flag and count must be visible to ALL participants
// maxcount is the number of participants (maxcount > number of participants will end in deadlock)
static inline void trivial_barrier(volatile int32_t *flag, volatile int32_t *count, int32_t maxcount){
  int32_t lgo, cnt ;
  lgo = *flag ;                             // get current flag value
  full_memory_fence() ;
  cnt = __sync_fetch_and_add (count, 1) ;   // increment count
  if(cnt == maxcount - 1){                  // last participant arriving at barrier
    cnt = __sync_fetch_and_add (count, -maxcount) ; // reset count to zero
    *flag = 1 - *flag ;                     // invert flag
  }else{
    while(lgo == *flag) ;                   // wait until flag is inverted
  }
}

static inline int32_t try_acquire_fence_idlock(volatile int32_t *lock, int32_t id) {
  __asm__ volatile("": : :"memory");
  if (__sync_val_compare_and_swap(lock, 0, (id+1)) != 0) return 0;
  full_memory_fence() ;
  return 1;
}

static inline int32_t try_acquire_fence_lock(volatile int32_t *lock) {
  return try_acquire_fence_idlock(lock, 1) ;
}

static inline void acquire_fence_idlock(volatile int32_t *lock, int32_t id){
  full_memory_fence() ;
  while(__sync_val_compare_and_swap(lock, 0, (id+1)) != 0) ;
  __asm__ volatile("": : :"memory") ;
}

static inline void acquire_fence_lock(volatile int32_t *lock){   // no id, use 1
  acquire_fence_idlock(lock, 1) ;
}

static inline int32_t try_release_fence_idlock(volatile int32_t *lock, int32_t id){
  full_memory_fence() ;
  if(__sync_val_compare_and_swap(lock, (id+1), 0) != (id+1)) return 0 ;  // failure, not lock owner
  __asm__ volatile("": : :"memory") ;
  return 1 ;                                                             // success
}

static inline int32_t try_release_fence_lock(volatile int32_t *lock){   // no id, use 1
  return try_release_idlock(lock, 1) ;
}

static inline void release_fence_idlock(volatile int32_t *lock, int32_t id){
  if( *lock <= 0) return ;                                             // NOOP, was not even locked
  full_memory_fence() ;
  while(__sync_val_compare_and_swap(lock, (id+1), 0) != (id+1)) ;
  __asm__ volatile("": : :"memory") ;
}

static inline void release_fence_lock(volatile int32_t *lock){   // no id, use 1
  release_fence_idlock(lock, 1) ;
}
