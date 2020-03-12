/* functions for C and FORTRAN programming
 * Copyright (C) 2019  Recherche en Prevision Numerique
 *
 * This software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
//****P* librkl/RPN kernel library memory arena management
// DESCRIPTION
// set of routines to implement named block management in a memory pool
// possibly shared by multiple threads and processes
//
//          master arena layout (there must be one and only one)
// +--------------------+--------------------+---------------------+-------------------->
// | master tables      | arena header       | symbol table        | data blocks
// +--------------------+--------------------+---------------------+-------------------->
//
//          memory arena layout (multiple arenas can coexist)
// +--------------------+---------------------+-------------------->
// | arena header       | symbol table        | data blocks
// +--------------------+---------------------+-------------------->
//
// indices are used instead of addresses because the memory arena might be mapped 
// at different addresses in different processes
//
//          data block layout
//     +----------------------------------------------------------------- +
//     |                                                                  |
//     |                                                                  v
// +-------+-------+-------+-------+.....................+-------+-------+
// |  FWD  |  IX   |  NWD  | SIGNL |  user data portion  | SIGNH |  BWD  |
// +-------+-------+-------+-------+.....................+-------+-------+
//  ^                                                                 |
//  |                                                                 |
//  +-----------------------------------------------------------------+
// FWD   : index of start of next block
// IX    : index in symbol table of this block
// NWD   : size of data portion in 64 bit units
// SIGNL : low marker (used for checking data underruns)
// SIGNH : high marker (used for checking data overruns)
// BWD   : index of start of this block
// FWD of last allocated block will point to a non existent block with FWD = 0
// 
// FWD and BWD are indices into a 64 bit unsigned integer array starting at the beginning of the memory arena
// FWD, IX, NWD, SIGNL, SIGNH, BWD are 32 bit unsigned integers
//
#if defined(NEVER_EVER_TRUE)
// EXAMPLES
//****
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <immintrin.h>

#include "memory_arena.h"

// memory store fence
#define W_FENCE asm volatile("": : :"memory"); _mm_sfence();

// memory load fence
#define R_FENCE asm volatile("": : :"memory"); _mm_lfence();

// memory load+store fence
#define M_FENCE asm volatile("": : :"memory"); _mm_mfence();

static uint32_t me = 999999999;  // identifier for this process (usually MPI rank) (alternative : getpid() )

static local_arena LA;

// interface                                                                      !InTf
//****f* librkl/memory_arena_set_id ( set memory arena id )
// Synopsis
// set owner's id (usually MPI rank) for memory arenas
//
// function memory_arena_set_id(id) result(me) BIND(C,name='memory_arena_set_id') !InTf
//   import :: C_INT                                                              !InTf
//   integer(C_INT), intent(IN), value :: id                                      !InTf
//   integer(C_INT) :: me                                                         !InTf
// end function memory_arena_set_id                                               !InTf
//
// set id for memory management arena, return identifier (-1 in case of error)
// id must be a POSITIVE INTEGER
//
// ARGUMENTS
int32_t memory_arena_set_id(uint32_t id){
//****
  if(id < 0) return -1;
  me = id + 1;
  return me;
}

// translate char string (max 8 characters) into a 64 bit unsigned token
// translation will stop at first null or space character
static inline uint64_t block_name(unsigned char *name){
  int i;
  uint64_t name64 = 0;

  for(i = 0 ; i < 8 ; i++){      // build 64 bit name token
    if(name[i] == '\0' || name[i] == ' ') break;
    name64 = (name64 << 8) | (name[i] & 0x7F);
  }
  return name64;
}

//****f* librkl/memory_arena_print_status ( print description of contents of memory arena )
// Synopsis
// dump arena header and symbol table
//
// subroutine memory_arena_print_status(mem) BIND(C,name='memory_arena_print_status') !InTf
//   import :: C_PTR                                                                  !InTf
//   type(C_PTR), intent(IN), value :: mem                                            !InTf
// end subroutine memory_arena_print_status                                           !InTf
//
// mem  : pointer to previously created and initialized memory arena
//        see  memory_arena_init
//
// ARGUMENTS
void memory_arena_print_status(void *mem){
//****
  uint64_t *mem64 = (uint64_t *) mem;
  memory_arena *ma = (memory_arena *) mem;
  symtab_entry *sym = ma->t;
  int i, j, sane;
  char name[9];
  uint64_t dname, size64;
  uint64_t *dataptr64;
  block_header *bh, *bhnext;
  block_tail *bt;

  fprintf(stderr,"Arena Header, id = %d, address = %p\n", me, ma);
  fprintf(stderr,"owner       = %8.8x\n",ma->owner);
  fprintf(stderr,"max entries = %d\n",ma->max_entries);
  fprintf(stderr,"max size    = %d\n",ma->arena_size);
  fprintf(stderr,"entries     = %d\n",ma->n_entries);
  fprintf(stderr,"first free  = %d\n",ma->first_free);

  fprintf(stderr,"\nSymbol table\n==============================================\n");
  for(i = 0 ; i < ma->n_entries ; i++){
    size64 = sym[i].data_size;
    dataptr64 = sym[i].data_index + mem64; 
    bh = (block_header *) (dataptr64);
    bt = (block_tail   *) (dataptr64 + BlockHeaderSize64 + size64);
    bhnext = (block_header *) (dataptr64 + BlockHeaderSize64 + size64 + BlockTailSize64);
    dname = sym[i].data_name;
    sane = ( bh->sign == 0xBEEFF00D ) & 
           ( bt->sign == 0xDEADBEEF ) & 
           (i == bh->ix) & 
           (bt->bwd == sym[i].data_index) &
           (bh->nwd == sym[i].data_size) &
           (bh->fwd == sym[i].data_index + sym[i].data_size + BlockHeaderSize64 + BlockTailSize64) ;
    for(j = 0; j < 9 ; j++) {
      name[j] = '\0';
      if( 0 == (dname >> 56) ) dname <<= 8;
    }
    for(j = 0; j < 8 ; j++) {
      name[j] = dname >> 56;
      dname <<= 8;
    }
    fprintf(stderr,"%4d: %4d F=%8.8x I=%8d S=%8d (%8d) %s %s %s FW=%8d FWNXT=%8d BW=%8d '%s'\n",
            i,bh->ix,sym[i].flags,sym[i].data_index,sym[i].data_size,bh->nwd,
            sane ? "T" : "F", ( bh->sign == 0xBEEFF00D ) ? "t" : "f", ( bt->sign == 0xDEADBEEF ) ? "t" : "f",
            bh->fwd, bhnext->fwd, bt->bwd, name);
  }
  fprintf(stderr,"==============================================\n");
}

//****f* librkl/memory_arena_init  ( initialize a memory arena )
// Synopsis
// initialize an already allocated 'memory arena' (node shared memory usually), 
// return id of current process
//
// function memory_arena_init(mem, nsym, size) result(me) BIND(C,name='memory_arena_init')          !InTf
//   import :: C_PTR, C_INT                                                                         !InTf
//   type(C_PTR), intent(IN), value :: mem                                                          !InTf
//   integer(C_INT), intent(IN), value :: nsym, size                                                !InTf
//   integer(C_INT) :: me                                                                           !InTf
// end function memory_arena_init                                                                   !InTf
//
// mem  : pointer to memory area (see memory_arena_init)
// size : size of memory area in 32 bit units
// nsym : size of symbol table to allocate (max number of blocks expected)
//
// ARGUMENTS
uint32_t memory_arena_init(void *mem, uint32_t nsym, uint32_t size){
//****
  memory_arena *ma = (memory_arena *) mem;
  symtab_entry *sym = ma->t;
  uint32_t size64 = size >> 1;  // round size down to 64 bit element size
  int i;
fprintf(stderr,"ma init %p, owner = %d\n", ma, ma->owner);
  while(__sync_val_compare_and_swap(&(ma->lock), 0, me) != 0); // lock memory arena

  if(ma->owner != 0) return ma->owner;                           // area already initialized, return id of initializer

  ma->owner = 0;           // initialize arena header
  ma->max_entries = nsym;
  ma->first_free = ArenaHeaderSize64 + nsym * SymtabEntrySize64;
// fprintf(stderr,"ArenaHeaderSize64 = %d, SymtabEntrySize64 = %d, nsym = %d, base = %d\n",ArenaHeaderSize64,SymtabEntrySize64,nsym,ma->first_free);
  ma->n_entries = 0;
  ma->arena_size = size64;

  for(i = 0 ; i < nsym ; i++){   // initialize symbol table to null values
    sym[i].lock       = 0;
    sym[i].flags      = 0;
    sym[i].data_index = 0;
    sym[i].data_size  = 0;
    sym[i].data_name  = 0;
  }

  ma->owner = me;  // flag area as initialized by me

  return __sync_val_compare_and_swap(&(ma->lock), me, 0); // unlock memory arena and return my id
}

// function update_local_table(mem) result(status) BIND(C,name='update_local_table')                !InTf
//   import :: C_PTR, C_INT                                                                         !InTf
//   type(C_PTR), intent(IN), value :: mem                                                          !InTf
//   integer(C_INT) :: status                                                                       !InTf
// end function update_local_table                                                                  !InTf
// update local table from master arena
// mem    : address of master arena
// return : number of arenas detected
uint32_t update_local_table(void *mem){
  master_arena *MA = (master_arena *) mem;
  memory_arena *ma = (memory_arena *) &(MA->ma);
  int i;

  while(__sync_val_compare_and_swap(&(MA->lock), 0, me) != 0); // lock master arena

  LA.lock      = 0;
  LA.master_id = me;
  LA.master_sz = MA->arena_sz;
  LA.MA        = MA;

  LA.le[0].arena_sz    = ma->arena_size;                 // memory arena associated to master arena
  LA.le[0].arena_name  = block_name((unsigned char *)(unsigned char *)"MaStEr");
  LA.le[0].ma          = ma;
fprintf(stderr,"local update, arena = %d, id = %d, address = %p, size = %ld\n",0, MA->me[0].arena_id, LA.le[0].ma, LA.le[0].arena_sz);

  for(i=1 ; i<MAX_MASTER ; i++){   // zero rest of local table
    if(LA.le[i].arena_sz == 0){                     // not initialized yet
      if(MA->me[i].arena_sz > 0) {                  // there is a shared memory segment
        LA.le[i].ma  = shmat(MA->me[i].arena_id, NULL, 0);    // attach segment, get address
        LA.le[i].arena_sz    = MA->me[i].arena_sz;
        LA.le[i].arena_name  = MA->me[i].arena_name;
      }else{                                        // no more segment, break
        break;
      }
    }
fprintf(stderr,"local update, arena = %d, id = %d, address = %p, size = %ld\n",i, MA->me[i].arena_id, LA.le[i].ma, LA.le[i].arena_sz);
  }

  __sync_val_compare_and_swap(&(MA->lock), me, 0); // unlock master arena

  return i;
}

// function master_arena_init(mem, nsym, size) result(id) BIND(C,name='master_arena_init')          !InTf
//   import :: C_PTR, C_INT                                                                         !InTf
//   type(C_PTR), intent(IN), value :: mem                                                          !InTf
//   integer(C_INT), intent(IN), value :: nsym, size                                                !InTf
//   integer(C_INT) :: id                                                                           !InTf
// end function master_arena_init                                                                   !InTf

// initialize an already allocated 'master arena' (node shared memory usually), return id of current process
// mem  : pointer to memory area
// size : size of memory area in 32 bit units
// nsym : size of symbol table to allocate (max number of blocks expected)
uint32_t master_arena_init(void *mem, uint32_t nsym, uint32_t size){
  master_arena *MA = (master_arena *) mem;
  memory_arena *ma = (memory_arena *) &(MA->ma);
  int i, status;

  size = size - MasterHeaderSize64 * 2;     // space left for memory arena proper

  while(__sync_val_compare_and_swap(&(MA->lock), 0, me) != 0); // lock master arena

  for(i=0 ; i<MAX_MASTER ; i++){   // zero all entries in master table
    MA->me[i].arena_name = 0;
    MA->me[i].arena_sz   = 0;
    MA->me[i].arena_id   = -1;
    MA->me[i].owner_id   = -1;
  }
  MA->me[0].arena_name = block_name((unsigned char *)"MaStEr");  // special name for master arena
  MA->me[0].arena_sz   = size >> 1;             // fix size entry of area 0 (arena part of master arena)
  MA->me[0].owner_id   = me;                    // creator id
// printf("MA = %p, ma = %p, delta = %ld\n",MA, ma, (void *)ma - (void *)MA);

  status = memory_arena_init(ma, nsym, size);   // initialize memory arena part of master arena

  __sync_val_compare_and_swap(&(MA->lock), me, 0); // unlock master arena

  return status;
}

// find entry in symbol table, return index if found, -1 otherwise
static inline int32_t find_block(memory_arena *ma, symtab_entry *sym, uint64_t name64){
  uint32_t i;

  for(i = 0 ; i < ma->n_entries ; i++){
    if(sym[i].data_name == name64){
      return i;
    }
  }
  return -1 ; // miserable failure
}

//****f* librkl/memory_block_find ( find a memory block in a memory arena )
// Synopsis
// find memory block 'name', return data address (NULL if not found), 
//                                  size of block (0 if not found),
//                                  block flags (0 if not found),
//
// function memory_block_find(mem, size, flags, name) result(ptr) BIND(C,name='memory_block_find') !InTf
//   import :: C_PTR, C_INT, C_CHAR                                               !InTf
//   type(C_PTR), intent(IN), value :: mem                                        !InTf
//   integer(C_INT), intent(OUT) :: size, flags                                   !InTf
//   character(C_CHAR), dimension(*), intent(IN) :: name                          !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function memory_block_find                                                 !InTf
//
// mem   : address of the managed 'memory arena'
// size  : size of the memory block (in 32 bit units)
// flags : blosk flags
// name  : name of block to find (characters beyond the 8th will be ignored)
// ptr   : local address of block
//
// ARGUMENTS
void *memory_block_find(void *mem, uint32_t *size, uint32_t *flags, unsigned char *name){
//****
  uint64_t *mem64 = (uint64_t *) mem;
  memory_arena *ma = (memory_arena *) mem;
  symtab_entry *sym = ma->t;
  void *dataptr = NULL;
  uint64_t name64 = block_name(name);
  int32_t i;

  *size = 0;         // precondition to fail
  *flags = 0;
  if(ma == NULL || name == NULL) return NULL;
  name64 = block_name(name);
  if(name64 == 0) return NULL;

  i = find_block(ma, sym, name64);
  if(i < 0) return NULL;  // name not found in symbol table

  *size = sym[i].data_size * 2;            // return size in 32 bit units
  *flags = sym[i].flags;
  dataptr = &mem64[sym[i].data_index + BlockHeaderSize64];     // pointer to actual data
  return dataptr;
}

//****f* librkl/memory_block_find_wait ( find a memory block in a memory arena ) (with wait)
// Synopsis
// same as memory_block_find, but wait until block is created (or timeout in milliseconds expires)
// (timeout(milliseconds) = -1 means infinite wait for all practical purposes, 3600000 is one hour)
//
// function memory_block_find_wait(mem, size, flags, name, timeout) result(ptr) BIND(C,name='memory_block_find_wait')  !InTf
//   import :: C_PTR, C_INT, C_CHAR                                               !InTf
//   type(C_PTR), intent(IN), value :: mem                                        !InTf
//   integer(C_INT), intent(OUT) :: size, flags                                   !InTf
//   character(C_CHAR), dimension(*), intent(IN) :: name                          !InTf
//   integer(C_INT), intent(IN), value :: timeout                                 !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function memory_block_find_wait                                            !InTf
//
// mem   : address of the managed 'memory arena'
// size  : size of the memory block (in 32 bit units)
// flags : blosk flags
// name  : name of block to find (characters beyond the 8th will be ignored)
// timeout : time to wait for block creation in milliseconds (-1 means forever)
// ptr   : local address of block
//
// ARGUMENTS
void *memory_block_find_wait(void *mem, uint32_t *size, uint32_t *flags, unsigned char *name, int timeout){
//****
  void *p = NULL;
  useconds_t delay = 1000;  // 1000 microseconds = 1 millisecond

  p = memory_block_find(mem, size, flags, name);     // does the block exist ?
  while(p == NULL && timeout > 0) {                  // no, sleep a bit and retry
    usleep(delay); timeout --;                       // decrement timeout
    p = memory_block_find(mem, size, flags, name);   // does the block exist ?
  }
// fprintf(stderr,"timeout = %d\n",timeout);
  return p;
}

//****f* librkl/memory_block_mark_init ( mark a memory block as initialized )
// Synopsis
// mark memory block 'name' as initialized, return block address if found, NULL otherwise
//
// function memory_block_mark_init(mem, name) result(ptr) BIND(C,name='memory_block_mark_init') !InTf
//   import :: C_PTR, C_CHAR                                              !InTf
//   type(C_PTR), intent(IN), value :: mem                                !InTf
//   character(C_CHAR), dimension(*), intent(IN) :: name                  !InTf
//   type(C_PTR) :: ptr                                                   !InTf
// end function memory_block_mark_init                                    !InTf
//
// mem  : address of the managed 'memory arena'
// name : name of block to mark (characters beyond the 8th will be ignored)
// ptr  : local address of block
//
// ARGUMENTS
void *memory_block_mark_init(void *mem, unsigned char *name){
//****
  uint64_t *mem64 = (uint64_t *) mem;
  memory_arena *ma = (memory_arena *) mem;
  symtab_entry *sym = ma->t;
  uint64_t name64 = block_name(name);
  int32_t i;
  void *dataptr = NULL;

  i = find_block(ma, sym, name64);
  if(i < 0) return NULL;  // name not found in symbol table

  while(__sync_val_compare_and_swap(&(sym[i].lock), 0, me) != 0); // lock block
  if(sym[i].flags == 0) sym[i].flags = me;                        // mark as initialized by me
  i = __sync_val_compare_and_swap(&(sym[i].lock), me, 0);         // unlock block

  dataptr = &mem64[sym[i].data_index + BlockHeaderSize64];        // pointer to actual data
  return dataptr;
}

//****f* librkl/memory_block_create ( create a named memory block in a memory arena )
// Synopsis
// create a named block in a managed 'memory arena'
// return start address of data (NULL in case of error)
//
// function memory_block_create(mem, size, name) result(ptr) BIND(C,name='memory_block_create') !InTf
//   import :: C_PTR, C_INT, C_CHAR                                               !InTf
//   type(C_PTR), intent(IN), value :: mem                                        !InTf
//   integer(C_INT), intent(IN), value :: size                                    !InTf
//   character(C_CHAR), dimension(*), intent(IN) :: name                          !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function memory_block_create                                               !InTf
//
// mem  : address of the managed 'memory arena'
// size : desired size of block in 32 bit units
// name : name of block to create (characters beyond the 8th will be ignored)
// ptr  : local address of created block
//
// ARGUMENTS
void *memory_block_create(void *mem, uint32_t size, unsigned char *name){
//****
  uint64_t *mem64 = (uint64_t *) mem;
  memory_arena *ma = (memory_arena *) mem;
  symtab_entry *sym = ma->t;
  uint32_t i, next;
  uint32_t fail;
  uint32_t size64 = (size + 1) >> 1;  // round size up to 64 bit element size
  uint32_t block64 = size64 + BlockHeaderSize64 + BlockTailSize64;
  block_header *bh;
  block_tail *bt;
  char *dataptr;
  uint64_t name64 = block_name(name);

  while(__sync_val_compare_and_swap(&(ma->lock), 0, me) != 0); // lock memory area

  fail  = ma->first_free + block64 + 1 > ma->arena_size;    // block larger than what we have left
  fail |= ma->n_entries == ma->max_entries;                 // symbol table is full

  if(fail){
    dataptr = NULL;
  }else{
    i = ma->n_entries;

    sym[i].lock  = 0;                     // keep lock as unlocked
    sym[i].flags = 0;                     // keep flag as uninitialized
    sym[i].data_index = ma->first_free;   // start of block
    sym[i].data_size = size64;            // data size for block
    sym[i].data_name = name64;            // data block name

    next = ma->first_free + block64;
    mem64[next] = 0;                      // fwd for next block will be 0

    bh = (block_header *) (mem64 + ma->first_free);    // start of block
    dataptr = (char *) (mem64 + ma->first_free + BlockHeaderSize64) ; // start of data in block
    bh->fwd = next;                       // next block will start there
    bh->ix = i;                           // index of this block in symbol table
    bh->nwd = size64;                     // size of data portion
    bh->sign = 0xBEEFF00D;                // marker below data

    bt = (block_tail *) (mem64 + ma->first_free + BlockHeaderSize64 + size64);
    bt->sign = 0xDEADBEEF;                // marker above data
    bt->bwd = sym[i].data_index;          // back pointer, index of start of current block

    ma->first_free = next;                // bump index of next free position
    ma->n_entries++;                      // bump number of valid entries
  }

  i = __sync_val_compare_and_swap(&(ma->lock), me, 0);         // unlock memory area
  return dataptr;
}

//****f* librkl/memory_allocate_shared ( allocate a shared memory block )
// Synopsis
// allocate a shared memory segment (on his process)
//
// function memory_allocate_shared(shmid, size) result(ptr) BIND(C,name='memory_allocate_shared') !InTf
//   import :: C_PTR, C_INT                                                       !InTf
//   integer, intent(OUT) :: shmid                                                !InTf
//   integer, intent(IN), value :: size                                           !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function memory_allocate_shared                                            !InTf
//
// shmid : shared memory id of block (set by memory_allocate_shared) (see shmget)
// size  : size of block in 32 bit units
// return local address of memory block
//
// ARGUMENTS
void *memory_allocate_shared(int *shmid, uint32_t size){    
//****
  int id = -1;
  void *shmaddr = NULL;
  size_t shmsz = size * sizeof(uint32_t);  // 32 bit units to bytes
  int err;
  struct shmid_ds dummy;

  if(me == 999999999) me = getpid();   // if not initialized, set to pid

  id = shmget(IPC_PRIVATE, shmsz, 0600);    // get a memory block, only accessible by user
  *shmid = id;                              // block id returned to caller
  if(id == -1) return NULL;
  shmaddr = shmat(id, NULL, 0);             // local address of memory block
  err = shmctl(id, IPC_RMID, &dummy);       // mark block as to be deleted when no process attached
  if(err == -1) return NULL;

  return shmaddr;     // return local address of memory block
}

// function memory_arena_create_shared(shmid, nsym, size) result(ptr) BIND(C,name='memory_arena_create_shared') !InTf
//   import :: C_PTR, C_INT                                                       !InTf
//   integer, intent(OUT) :: shmid                                                !InTf
//   integer, intent(IN), value :: nsym, size                                     !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function memory_arena_create_shared                                        !InTf

void *memory_arena_create_shared(int *shmid, uint32_t nsym, uint32_t size){
  void *shmaddr = memory_allocate_shared(shmid, size);    // request shared memory block
  int err;

  if(shmaddr == NULL) return shmaddr;             // request failed

  err = memory_arena_init(shmaddr, nsym, size);   // initialize memory arena
  if(err < 0) return NULL;

  return shmaddr;
}

// function master_arena_create_shared(shmid, nsym, size) result(ptr) BIND(C,name='master_arena_create_shared') !InTf
//   import :: C_PTR, C_INT                                                       !InTf
//   integer, intent(OUT) :: shmid                                                !InTf
//   integer, intent(IN), value :: nsym, size                                     !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function master_arena_create_shared                                        !InTf

void *master_arena_create_shared(int *shmid, uint32_t nsym, uint32_t size){
  void *shmaddr = memory_allocate_shared(shmid, size);    // request shared memory block
  int err;
  master_arena *MA;

  MA = (master_arena *) shmaddr;
  MA->lock       = 0;
  MA->arena_id   = *shmid;
  MA->arena_sz   = size >> 1;             // 64 bit units
  MA->arena_name = block_name((unsigned char *)"MaStEr");  // special name
  
printf("MA = %p, id = %d\n",MA, MA->arena_id);
  err = master_arena_init(shmaddr, nsym, size);
  if(err < 0) return NULL;
  MA->me[0].arena_id   = MA->arena_id;                    // segment id

  return MA;                       // return address of master arena
}

//****f* librkl/memory_address_from_id ( get memory address associated with shared memory segment id )
// Synopsis
// get memory address associated with shared memory segment shmid
//
// function memory_address_from_id(shmid) result(ptr) BIND(C,name='memory_address_from_id') !InTf
//   import :: C_PTR, C_INT                                                          !InTf
//   integer, intent(IN), value :: shmid                                             !InTf
//   type(C_PTR) :: ptr                                                              !InTf
// end function memory_address_from_id                                               !InTf
//
// schmid    : shared memory segment id (from memory_arena_create_shared, memory_allocate_shared, master_arena_create_shared)
// ptr       : local memory addres of said segment
// ARGUMENTS
void *memory_address_from_id(int shmid){
//****
  return shmat(shmid, NULL, 0);
}

// function memory_arena_from_master(mem) result(ptr) BIND(C,name='memory_arena_from_master') !InTf
//   import :: C_PTR                                                                 !InTf
//   type(C_PTR), intent(IN), value :: mem                                           !InTf
//   type(C_PTR) :: ptr                                                              !InTf
// end function memory_arena_from_master                                             !InTf

// get memory arena address from master arena address
// mem     : local memory address of master arena
// return  : local memory addres of memory arena from master arena
void *memory_arena_from_master(void *mem){
  master_arena *MA = (master_arena *) mem;
  return &(MA->ma);
}

// function memory_arena_from_master_id(shmid) result(ptr) BIND(C,name='memory_arena_from_master_id') !InTf
//   import :: C_PTR, C_INT                                                          !InTf
//   integer, intent(IN), value :: shmid                                             !InTf
//   type(C_PTR) :: ptr                                                              !InTf
// end function memory_arena_from_master_id                                          !InTf

// get memory address associated with shared memory segment shmid
// schmid  : master arena segment id (from master_arena_create_shared)
// return  : local memory addres of memory arena from master arena
void *memory_arena_from_master_id(int shmid){
  void *shmaddr = shmat(shmid, NULL, 0);
  master_arena *MA;

  if(shmaddr == NULL) return NULL;

  MA = (master_arena *) shmaddr;
  return &(MA->ma);
}
// end interface                                                                     !InTf
#if defined(SELF_TEST)
#include <errno.h>

#include <mpi.h>

#define NSYM 128
#define DBLK 20

int main(int argc, char **argv){
  int err, rank, size, id;
  int shmid = -1;
  void *shmaddr = NULL;
  void *p = NULL;
  int shmsz = 1024 * 1024 * 4;  // 4 MBytes
  struct shmid_ds dummy;

  err = MPI_Init(&argc, &argv);
  err = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  err = MPI_Comm_size(MPI_COMM_WORLD,&size);
  id = memory_arena_set_id(rank);
  if(rank == 0) {
//     shmid = shmget(IPC_PRIVATE, shmsz, 0600);
//     shmaddr = shmat(shmid, NULL, 0);
//     err = shmctl(shmid, IPC_RMID, &dummy);
//     system("ipcs -m");
//     if(shmaddr == (void *) -1) {
//       perror("shmat");
//       fprintf(stderr,"error attaching segment %d, %p, errno = %d\n",shmid,shmaddr,errno);
//       exit(1);
//     }
    shmaddr = memory_arena_create_shared(&shmid, NSYM, shmsz);
    memory_arena_print_status(shmaddr);
    p = memory_block_create(shmaddr, DBLK*1, "BLOCK000"); p = memory_block_mark_init(shmaddr, "BLOCK000");
    p = memory_block_create(shmaddr, DBLK*2, "BLOCK001"); p = memory_block_mark_init(shmaddr, "BLOCK001");
    p = memory_block_create(shmaddr, DBLK*3, "BLOCK002"); p = memory_block_mark_init(shmaddr, "BLOCK002");
    p = memory_block_create(shmaddr, DBLK*4, "BLOCK003"); p = memory_block_mark_init(shmaddr, "BLOCK003");
    p = memory_block_create(shmaddr, DBLK*5, "BLOCK004"); p = memory_block_mark_init(shmaddr, "BLOCK004");
  }
  err = MPI_Bcast(&shmid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  if(rank != 0) shmaddr = shmat(shmid, NULL, 0);
  fprintf(stderr,"I am process %d of %d, id = %d, shmid = %d, addr = %p\n",rank+1,size,id,shmid,shmaddr);
  err = MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 1) {
    system("ipcs -m");
    memory_arena_print_status(shmaddr);
  }
  err = MPI_Finalize();
}
#endif
