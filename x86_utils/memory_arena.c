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
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>

#include "memory_arena.h"

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
// 
//          memory arena layout
// +--------------------+---------------------+-------------------->
// | arena header       | symbol table        | data blocks
// +--------------------+---------------------+-------------------->
//
// indices are used instead of addresses because the memory arena might be mapped to different addresses 
// in different processes

static uint32_t me = 999999999;  // identifier for this process (usually MPI rank)
// interface                                                                      !InTf
// function memory_arena_set_id(id) result(me) BIND(C,name='memory_arena_set_id') !InTf
//   import :: C_INT                                                              !InTf
//   integer(C_INT), intent(IN), value :: id                                      !InTf
//   integer(C_INT) :: me                                                         !InTf
// end function memory_arena_set_id                                               !InTf

// set id for memory management arena, return identifier (-1 in case of error)
// id must be a POSITIVE INTEGER
int32_t memory_arena_set_id(uint32_t id){
  if(id < 0) return -1;
  me = id + 1;
  return me;
}

// subroutine memory_arena_print_status(mem) BIND(C,name='memory_arena_print_status') !InTf
//   import :: C_PTR                                                                  !InTf
//   type(C_PTR), intent(IN), value :: mem                                            !InTf
// end subroutine memory_arena_print_status                                           !InTf

// dump arena header and symbol table
void memory_arena_print_status(void *mem){
  uint64_t *mem64 = (uint64_t *) mem;
  arena_header *ap = (arena_header *) mem;
  symtab_entry *sym = (symtab_entry *) ( mem64 + ArenaHeaderSize64 );
  int i, j;
  char name[9];
  uint64_t dname;

  fprintf(stderr,"Arena Header, id = %d\n", me);
  fprintf(stderr,"flags       = %8.8x\n",ap->flags);
  fprintf(stderr,"max entries = %d\n",ap->max_entries);
  fprintf(stderr,"max size    = %d\n",ap->arena_size);
  fprintf(stderr,"entries     = %d\n",ap->n_entries);
  fprintf(stderr,"first free  = %d\n",ap->first_free);

  fprintf(stderr,"\nSymbol table\n");
  for(i = 0 ; i < ap->n_entries ; i++){
    dname = sym[i].data_name;
    for(j = 0; j < 9 ; j++) {
      name[j] = '\0';
      if( 0 == (dname >> 56) ) dname <<= 8;
    }
    for(j = 0; j < 8 ; j++) {
      name[j] = dname >> 56;
      dname <<= 8;
    }
    fprintf(stderr,"%4d: F=%8.8x I=%8d S=%8d '%s'\n",i,sym[i].flags,sym[i].data_index,sym[i].data_size,name);
  }
}

// function memory_arena_init(mem, nsym, size) result(id) BIND(C,name='function memory_arena_init') !InTf
//   import :: C_PTR, C_INT                                                                         !InTf
//   type(C_PTR), intent(IN), value :: mem                                                          !InTf
//   integer(C_INT), intent(IN), value :: nsym, size                                                !InTf
//   integer(C_INT) :: id                                                                           !InTf
// end function memory_arena_init                                                                   !InTf

// initialize an already allocated 'memory arena' (node shared memory usually), return id of current process
// mem  : pointer to memory area
// size : size of memory area in 32 bit units
// nsym : size of symbol table to allocate (max number of blocks expected)
uint32_t memory_arena_init(void *mem, uint32_t nsym, uint32_t size){
  uint64_t *mem64 = (uint64_t *) mem;
  arena_header *ap = (arena_header *) mem;
  symtab_entry *sym = (symtab_entry *) ( mem64 + ArenaHeaderSize64 );
  uint32_t size64 = size >> 1;  // round size down to 64 bit element size
  int i;

  while(__sync_val_compare_and_swap(&(ap->lock), 0, me) != 0); // lock memory area

  if(ap->flags != 0) return ap->flags;                           // area already initialized, return id of initializer

  ap->flags = 0;           // initialize arena header
  ap->max_entries = nsym;
  ap->first_free = ArenaHeaderSize64 + nsym * SymtabEntrySize64;
  ap->n_entries = 0;
  ap->arena_size = size64;

  for(i = 0 ; i < nsym ; i++){   // initialize symbol table
    sym[i].lock = 0;
    sym[i].flags = 0;
    sym[i].data_index = 0;
    sym[i].data_size = 0;
    sym[i].data_name = 0;
  }

  ap->flags = me;  // flag area as initilized by me

  return __sync_val_compare_and_swap(&(ap->lock), me, 0); // unlock and return my id
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

// find entry in symbol table, return index if found, -1 otherwise
static inline int32_t find_block(arena_header *ap, symtab_entry *sym, uint64_t name64){
  uint32_t i;

  for(i = 0 ; i < ap->n_entries ; i++){
    if(sym[i].data_name == name64){
      return i;
    }
  }
  return -1 ; // miserable failure
}

// function memory_block_find(mem, size, flags, name) result(ptr) BIND(C,name='') !InTf
//   import :: C_PTR, C_INT, C_CHAR                                               !InTf
//   type(C_PTR), intent(IN), value :: mem                                        !InTf
//   integer(C_INT), intent(OUT) :: size, flags                                   !InTf
//   character(C_CHAR), dimension(*), intent(IN) :: name                          !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function memory_block_find                                                 !InTf

// find memory block 'name', return data address (NULL if not found), 
//                                  size of block (0 if not found),
//                                  block flags (0 if not found),
// mem  : address of the managed 'memory arena'
// name : name of block to find (characters beyond the 8th will be ignored)
void *memory_block_find(void *mem, uint32_t *size, uint32_t *flags, unsigned char *name){
  uint64_t *mem64 = (uint64_t *) mem;
  arena_header *ap = (arena_header *) mem;
  symtab_entry *sym = (symtab_entry *) ( mem64 + ArenaHeaderSize64 );
  void *dataptr = NULL;
  uint64_t name64 = block_name(name);
  int32_t i;

  *size = 0;         // precondition to fail
  *flags = 0;
  if(ap == NULL || name == NULL) return NULL;
  name64 = block_name(name);
  if(name64 == 0) return NULL;

  i = find_block(ap, sym, name64);
  if(i < 0) return NULL;  // name not found in symbol table

  *size = sym[i].data_size * 2;            // return size in 32 bit units
  *flags = sym[i].flags;
  dataptr = &mem64[sym[i].data_index];     // pointer to actual data
  return dataptr;
}

// function memory_block_mark_init(mem, name) result(ptr) BIND(C,name='memory_block_mark_init') !InTf
//   import :: C_PTR, C_CHAR                                              !InTf
//   type(C_PTR), intent(IN), value :: mem                                !InTf
//   character(C_CHAR), dimension(*), intent(IN) :: name                  !InTf
//   type(C_PTR) :: ptr                                                   !InTf
// end function memory_block_mark_init                                    !InTf

// mark memory block 'name' as initialized, return block address if found, NULL otherwise
// mem  : address of the managed 'memory arena'
// name : name of block to mark (characters beyond the 8th will be ignored)
void *memory_block_mark_init(void *mem, unsigned char *name){
  uint64_t *mem64 = (uint64_t *) mem;
  arena_header *ap = (arena_header *) mem;
  symtab_entry *sym = (symtab_entry *) ( mem64 + ArenaHeaderSize64 );
  uint64_t name64 = block_name(name);
  int32_t i;
  void *dataptr = NULL;

  i = find_block(ap, sym, name64);
  if(i < 0) return NULL;  // name not found in symbol table

  while(__sync_val_compare_and_swap(&(sym[i].lock), 0, me) != 0); // lock block
  if(sym[i].flags == 0) sym[i].flags = me;                        // mark as initialized by me
  i = __sync_val_compare_and_swap(&(sym[i].lock), me, 0);         // unlock block

  dataptr = &mem64[sym[i].data_index];                            // pointer to actual data
  return dataptr;
}

// function memory_block_create(mem, size, name) result(ptr) BIND(C,name='memory_block_create') !InTf
//   import :: C_PTR, C_INT, C_CHAR                                               !InTf
//   type(C_PTR), intent(IN), value :: mem                                        !InTf
//   integer(C_INT), intent(IN), value :: size                                    !InTf
//   character(C_CHAR), dimension(*), intent(IN) :: name                          !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function memory_block_create                                               !InTf

// create a named block in a managed 'memory arena'
// return start address of data (NULL in case of error)
// mem  : address of the managed 'memory arena'
// size : desired size of blok in 32 bit units
// name : name of block to create (characters beyond the 8th will be ignored)
void *memory_block_create(void *mem, uint32_t size, unsigned char *name){
  uint64_t *mem64 = (uint64_t *) mem;
  arena_header *ap = (arena_header *) mem;
  symtab_entry *sym = (symtab_entry *) ( mem64 + ArenaHeaderSize64 );
  uint32_t i, next;
  uint32_t fail;
  uint32_t size64 = (size + 1) >> 1;  // round size up to 64 bit element size
  uint32_t block64 = size64 + BlockHeaderSize64 + BlockTailSize64;
  block_header *bh;
  block_tail *bt;
  char *dataptr;
  uint64_t name64 = block_name(name);

//   for(i = 0 ; i < 8 ; i++){      // build 64 bit name token
//     if(name[i] == '\0' || name[i] == ' ') break;
//     name64 = (name64 << 8) | name[i];
//   }
  

  while(__sync_val_compare_and_swap(&(ap->lock), 0, me) != 0); // lock memory area

  fail  = ap->first_free + block64 + 1 > ap->arena_size;    // block larger than what we have left
  fail |= ap->n_entries == ap->max_entries;                 // symbol table is full

  if(fail){
    dataptr = NULL;
  }else{
    ap->first_free = ap->first_free + block64;  // bump index of next free position
    i = ap->n_entries;

    sym[i].lock  = 0;                     // keep lock as unlocked
    sym[i].flags = 0;                     // keep flag as uninitialized
    sym[i].data_index = ap->first_free;   // start of block
    sym[i].data_size = size64;            // data size for block
    sym[i].data_name = name64;            // data block name

    next = ap->first_free + block64;
    mem64[next] = 0;                      // fwd for next block will be 0

    bh = (block_header *) (mem64 + ap->first_free);    // start of block
    dataptr = (char *) (mem64 + ap->first_free + BlockHeaderSize64) ; // start of data in block
    bh->fwd = next;                       // next block will start there
    bh->ix = i;                           // index of this block in symbol table
    bh->nwd = size64;                     // size of data portion
    bh->sign = 0xBEEFF00D;                // marker below data

    bt = (block_tail *) (mem64 + BlockHeaderSize64 + size64);
    bt->sign = 0xDEADBEEF;                // marker above data
    bt->bwd = sym[i].data_index;          // back pointer, index of start of current block

    ap->n_entries++;                      // bump number of valid entries
  }

  i = __sync_val_compare_and_swap(&(ap->lock), me, 0);         // unlock memory area
  return dataptr;
}

// function memory_arena_create_shared(id, nsym, size) result(ptr) BIND(C,name='memory_arena_create_shared') !InTf
//   import :: C_PTR, C_INT                                                       !InTf
//   integer, intent(OUT) :: id                                                   !InTf
//   integer, intent(IN), value :: nsym, size                                     !InTf
//   type(C_PTR) :: ptr                                                           !InTf
// end function memory_arena_create_shared                                        !InTf

void *memory_arena_create_shared(int *id, uint32_t nsym, uint32_t size){
  int shmid = -1;
  void *shmaddr = NULL;
  size_t shmsz = size * sizeof(uint32_t);  // 32 bit units to bytes
  int err;
  struct shmid_ds dummy;

  shmid = shmget(IPC_PRIVATE, shmsz, 0600);
  *id = shmid;
  if(shmid == -1) return NULL;
  shmaddr = shmat(shmid, NULL, 0);
  err = shmctl(shmid, IPC_RMID, &dummy);

  err = memory_arena_init(shmaddr, nsym, size);

  return shmaddr;
}

// function memory_arena_from_id(id) result(ptr) BIND(C,name='memory_arena_from_id') !InTf
//   import :: C_PTR, C_INT                                                          !InTf
//   integer, intent(IN), value :: id                                                !InTf
//   type(C_PTR) :: ptr                                                              !InTf
// end function memory_arena_from_id                                                 !InTf
// end interface                                                                     !InTf

// get memory address associated with shared memory segment shmid
void *memory_arena_from_id(int shmid){
  return shmat(shmid, NULL, 0);
}

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
