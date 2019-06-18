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

static inline void acquire_idlock(volatile void *lock, int32_t id){
  while(__sync_val_compare_and_swap((volatile uint32_t *)lock, 0, id) != 0) ;
}

static inline void acquire_lock(volatile void *lock){   // no id, use 1
  acquire_idlock(lock, 1) ;
}

// this will deadlock if attempt is made to release a lock with the wrong id
static inline void release_idlock(volatile void *lock, int32_t id){
  while(__sync_val_compare_and_swap((volatile uint32_t *)lock, id, 0) != id) ;
}

// this will deadlock if attempt is made to release a lock with id other than 1
static inline void release_lock(volatile void *lock){   // no id, use 1
  release_idlock(lock, 1) ;
}

static inline int32_t test_idlock(volatile void *lock, int32_t id){
  return (*(volatile uint32_t *)lock != id );   // true if locked with id
}

static inline int32_t test_lock(volatile void *lock){
  return (*(volatile uint32_t *)lock != 0 );   // true if locked with any id
}


static inline void force_reset_lock(volatile void *lock){   // forcefully reset lock
  *(volatile uint32_t *)lock = 0 ;
}
