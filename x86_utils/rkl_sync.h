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

static inline void acquire_idlock(volatile void *lock, int32_t id){
  while(__sync_val_compare_and_swap((volatile uint32_t *)lock, 0, id) != 0) ;
}

static inline void acquire_lock0(volatile void *lock){
  acquire_idlock(lock, -1) ;
}

static inline void release_idlock(volatile void *lock, int32_t id){
  while(__sync_val_compare_and_swap((volatile uint32_t *)lock, id, 0) != id) ;
}

static inline void release_lock0(volatile void *lock){
  release_idlock(lock, -1) ;
}

static inline int32_t test_idlock(volatile void *lock, int32_t id){
  return (*(volatile uint32_t *)lock != 0 );
}

static inline int32_t test_lock0(volatile void *lock){
  return (*(volatile uint32_t *)lock != 0 );
}

static inline void reset_lock(volatile void *lock){
  *(volatile uint32_t *)lock = 0 ;
}
