/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2018  Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdlib.h>
#include <stdint.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/stat.h>

#if defined(SELF_TEST)
#include <stdio.h>
#include <mpi.h>
#endif

// static int volatile __attribute__ ((aligned (64))) spins[64];
// static int volatile __attribute__ ((aligned (64))) barrs[64];
// static int volatile __attribute__ ((aligned (64))) locks[64];

static int32_t volatile *spins = NULL;
static int32_t volatile *barrs = NULL;
static int32_t volatile *locks = NULL;
static int32_t shared_id   = -1;
static int32_t shared_size = -1;
static int32_t *shared_ptr = NULL;

int32_t setup_locks_and_barriers(void *ptr, uint32_t size){
  size_t siz;
  int i;

  siz = (size / 192) * 192;   // need a multiple of 48 int32_t
  if(siz <= 0) return -1;     // not enough space

  if(ptr == NULL) ptr = malloc(siz); // allocate if NULL
  if(ptr == NULL) return -1;         // allocation failed

  siz = siz / 4 / 3;          // in int32_t units (multiple of 16)
  spins = (int32_t *)ptr;
  barrs = spins + siz;
  locks = barrs + siz;
  for(i=0 ; i<siz ; i++) {
    spins[i] = 0;
    barrs[i] = 0;
    locks[i] = 0;
  }
  
  return 0;
}

int32_t shared_locks_and_barriers_id(){
  return shared_id;
}

int32_t shared_locks_and_barriers_size(){
  return shared_size;
}

void *shared_locks_and_barriers_ptr(){
  return shared_ptr;
}

// allocate a  shared memory segment of size Bytes 
// segment will be preventively marked for deletion (only works on linux)
// return address and id of shared memory segment (in sid)
void *allocate_safe_shared_memory(int32_t *sid, uint32_t size){
  size_t siz = size;
  int32_t id;
  void *ptr;
  struct shmid_ds shm_buf;

  id  = shmget(IPC_PRIVATE,siz,IPC_CREAT|S_IRUSR|S_IWUSR);  // allocate shared memory segment
  if(id < 0) return NULL;         // error if id < 0
  ptr = shmat(id,NULL,0);
  if(ptr == NULL) return NULL;    // something is very wrong
  shmctl(id,IPC_RMID,&shm_buf);   // mark segment for deletion preventively (only works on linux)
  *sid  = id;                     // return shared memory segment id
  return ptr;                     // return address of shared memory segment
}

// return address and id of shared memory segment of size Bytes (rounded up to multiple of 192)
// that has been preventively marked for deletion
void *setup_shared_locks_and_barriers(int32_t *sid, uint32_t size){
  void *ptr;
  int32_t status;

  if(shared_ptr != NULL) return NULL; // setup already done

  size = ((size + 191) / 192) * 192;   // need a multiple of 48 int32_t (192 bytes)

  ptr = allocate_safe_shared_memory(sid, size);

  shared_id = *sid;
  shared_ptr = ptr;
  shared_size = size / 3 / 4;  // number of locks/spins/barriers (3 tables of int32_t)

  status = setup_locks_and_barriers(ptr, size);
  if(status != 0) return NULL;    // error during setup

  return ptr;                     // return address of shared memory segment
}

int32_t wait_barrier(int barrier, int maxcount){
 int32_t count;

 if(spins == NULL || barrs == NULL) return -1;

 count = __sync_fetch_and_add (barrs+barrier, 1);
 if(count == maxcount-1) {
  spins[barrier] = maxcount;
 }else{
  while(spins[barrier] != maxcount);
 }
 return 0;
}

int32_t reset_barrier(int barrier){

 if(spins == NULL || barrs == NULL) return -1;

 spins[barrier] = 0;
 barrs[barrier] = 0;
 return 0;
}

int32_t test_lock(int lock){

 if(spins == NULL || barrs == NULL) return -1;

 return locks[lock]-1;
}

int32_t release_lock(int lock, int me){
  int status;

 if(spins == NULL || barrs == NULL) return -1;

//  printf("attempting release of lock %d by %d owned by %d\n",lock,me,locks[lock]-1);
 return (__sync_val_compare_and_swap(locks+lock, me+1, 0) == me+1);
}

int32_t acquire_lock(int lock, int me){
  int count = 0;

 if(spins == NULL || barrs == NULL) return -1;

//  printf("attempting acquisition of lock %d owned by %d\n",lock,locks[lock]-1);
 while(__sync_val_compare_and_swap(locks+lock, 0, me+1) != 0) {
//    if(count++ == 0) printf("process %s waiting for lock %d\n",me,lock);
 }
//  printf("lock %d acquired by %d, lock = %d\n",lock,me,locks[lock]);
 return 0;
}
#if defined(SELF_TEST)

static uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

main(int argc, char **argv){
  int ierr, i;
  int localrank, localsize;
  int globalrank, globalsize;
  int peerrank, peersize;
  MPI_Comm MY_World = MPI_COMM_NULL;
  MPI_Comm MY_Peers = MPI_COMM_NULL;
  MPI_Comm temp_comm;
  int myhost, myhost0, myhost1;
  void *ptr;
  int32_t sid, size;
  uint64_t t0 , t1, t2;
  double tmin, tmax, tavg, tmp;

  ierr = MPI_Init( &argc, &argv );
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&globalrank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&globalsize);

  myhost  = gethostid();
  myhost0 = myhost & 0x7FFFFFFF  ; // lower 31 bits
  myhost1 = (myhost >> 31) & 0x1 ; // upper bit
  ierr = MPI_Comm_split(MPI_COMM_WORLD,myhost0,globalrank,&temp_comm) ;  // split WORLD using the lower 31 bits of host id , weight=rank in base
  ierr = MPI_Comm_split(temp_comm     ,myhost1,globalrank,&MY_World) ;   // re split using the upper bit of host id , weight=rank in base
  ierr = MPI_Comm_rank(MY_World,&localrank);                             // rank of this PE on this SMP node
  ierr = MPI_Comm_size(MY_World,&localsize);                             // number of PEs on this SMP node
  ierr = MPI_Comm_split(MPI_COMM_WORLD,localrank,globalrank,&MY_Peers) ; // communicator with PES of same local rank on other SMP nodes
  ierr = MPI_Comm_rank(MY_Peers,&peerrank);
  ierr = MPI_Comm_size(MY_Peers,&peersize);
  if(localrank == 0) printf("PEs in node = %d, peers in MY_Peers = %d\n",localsize,peersize);
  size = 3*4096;

  if(localrank == 0){
    ptr = allocate_safe_shared_memory(&sid, size);
    ierr = MPI_Bcast(&sid,1,MPI_INTEGER,0,MY_World);
    ierr = setup_locks_and_barriers(ptr, size);
#if ! defined(TIMING)
    printf("setup status = %d\n",ierr);
#endif
    ierr = MPI_Barrier(MY_World);

#if ! defined(TIMING)
    ierr = acquire_lock(1,localrank);
    if(ierr == 0) printf("lock acquired by %d\n",localrank);
    sleep(1);
    ierr = release_lock(1,localrank);
    printf("lock %d release status = %d\n",1,ierr);
    if(ierr != 0) printf("lock release successful\n");
#endif
  }else{
    ierr = MPI_Bcast(&sid,1,MPI_INTEGER,0,MY_World);
    ptr = shmat(sid,NULL,0);
    ierr = setup_locks_and_barriers(ptr, size);  
#if ! defined(TIMING)
    printf("setup status = %d\n",ierr);
#endif
    ierr = MPI_Barrier(MY_World);
#if ! defined(TIMING)
    ierr = acquire_lock(1,localrank);
    if(ierr == 0) printf("lock acquired by %d\n",localrank);
    sleep(1);
    ierr = release_lock(1,localrank);
    printf("lock %d release status = %d\n",1,ierr);
    if(ierr != 0) printf("lock release successful\n");
#endif
  }
    t0 = rdtsc();
    for(i=0 ; i<100 ; i++){
      ierr = acquire_lock(2,localrank);
      ierr = release_lock(2,localrank);
    }
    t1 = rdtsc();
//     printf("lock acquire/release time = %d\n",(t1-t0)/100);
    tmp = (t1-t0)/100;
    ierr = MPI_Allreduce(&tmp,&tmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&tmp,&tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&tmp,&tavg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    tavg = tavg / globalsize;
    if(globalrank == 0) printf("barrier min, max, avg = %9f, %9f, %9f\n",tmin,tmax,tavg);
    t0 = rdtsc();
    for(i=0 ; i<100 ; i++){
//       ierr = reset_barrier(i+1);
//       ierr = wait_barrier(i, localsize);
      ierr = reset_barrier(2);
      ierr = wait_barrier(1, localsize);
      ierr = reset_barrier(3);
      ierr = wait_barrier(2, localsize);
      ierr = reset_barrier(1);
      ierr = wait_barrier(3, localsize);
    }
    t1 = rdtsc();
//     printf("barrier time = %d\n",(t1-t0)/300);
    tmp = (t1-t0)/300;
    ierr = MPI_Allreduce(&tmp,&tmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&tmp,&tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&tmp,&tavg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    tavg = tavg / globalsize;
    if(globalrank == 0) printf("barrier min, max, avg = %9f, %9f, %9f\n",tmin,tmax,tavg);

  ierr = MPI_Finalize();
}
#endif
