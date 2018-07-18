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
#include <unistd.h>
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

typedef struct{
  int32_t volatile spin;
  int32_t volatile barr;
  int32_t volatile busy;
  int32_t volatile lock;
  int32_t volatile flag;
  int32_t volatile pad[11];
} comm_area;

static comm_area *node;
static int32_t volatile *spins = NULL;
static int32_t volatile *barrs = NULL;
static int32_t volatile *locks = NULL;
static int32_t volatile *busys = NULL;
static int32_t shared_id   = -1;
static int32_t shared_size = -1;
static int32_t *shared_ptr = NULL;

int32_t setup_locks_and_barriers(void *ptr, uint32_t size){
  size_t siz;
  int i;

  siz = (size / 256) * 256;   // need a multiple of 64 int32_t
  if(siz <= 0) return -1;     // not enough space

  if(ptr == NULL) ptr = malloc(siz); // allocate if NULL
  if(ptr == NULL) return -1;         // allocation failed

  siz = siz / 4 / 4;          // in int32_t units (multiple of 16)
  spins = (int32_t *)ptr;
  barrs = spins + siz;
  locks = barrs + siz;
  busys = locks + siz;
  node  = (comm_area *) ptr;
  for(i=0 ; i<siz ; i++) {
    spins[i] = 0;
    barrs[i] = 0;
    locks[i] = 0;
    busys[i] = 0;
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

// id       : identifier for this thread/process  ( 0 <= id < maxcount )
// maxcount : number of threads/processes for this barrier
void node_barrier00(int32_t id, int32_t maxcount){
  int32_t count, i, answer;

  if(maxcount < 2) return;      // barrier with 1 thread/process is a NO-OP

  barrs[id] = 1 ^ barrs[id];
  answer = ( barrs[id] == 0 ) ? 0 : maxcount ;
//   printf("answer(%d) = %d\n",id,answer);
  while(1){
    count = barrs[0];
    for(i=1 ; i<maxcount ; i++) count += barrs[i];
    break;
//     if(count == answer) break;
  }
}

// id       : identifier for this thread/process  ( 0 <= id < maxcount )
// maxcount : number of threads/processes for this barrier
void node_barrier0(int32_t id, int32_t maxcount){
  int32_t count;

  if(maxcount < 2) return;      // barrier with 1 thread/process is a NO-OP
  if(id == 0){
    barrs[id+1] = 1;
    while(barrs[id] != 2);
    barrs[id]  = 0;
  }else if(id == maxcount-1){
    while(barrs[id] != 1);
    barrs[id]  = 0;
    barrs[id-1]  = 2;
  }else{
    while(barrs[id] != 1);
    barrs[id+1] = 1;
    while(barrs[id] != 2);
    barrs[id-1]  = 2;
    barrs[id]  = 0;
  }

}

void node_barrier_multi(int32_t id, int32_t maxcount){  // version with 2 level "flag flip"
  int lgo, lgo2, count;
  int group, base, maxpop, lastgrp;
  int maxsiz = 6;              // eventually maxsiz = table1[maxcount]

  if(maxcount < 2) return;     // trivial case, no need for fancy footwork

  maxpop = maxsiz; group = 0; base = 0;
  lastgrp = (maxcount + maxpop - 1) / maxpop ;       // eventually  lastgrp = table2[maxcount]
  while(id >= base+maxpop) { group++ ; base += maxpop ; }  // which group of maxsiz does id belong to ? 
  maxpop = maxcount - base;
  maxpop = (maxpop > maxsiz) ? maxsiz : maxpop;
//   printf("maxpop = %d, group = %d, lastgrp = %d\n", maxpop, group, lastgrp);

  lgo  = node[group].flag;
  lgo2 = node[lastgrp].flag;
  count = __sync_fetch_and_add (&(node[group].barr), 1);        // bump group count
  if(count == maxpop - 1){                   // group count has been reached
      count = __sync_fetch_and_add (&(node[lastgrp].barr), 1);  // bump lastgrp count
      if(count == lastgrp - 1){              // lastgrp count has been reached
	node[lastgrp].barr = 0;              // reset lastgrp count to zero
	node[lastgrp].flag = 1 - lgo2;       // set lastgrp done flag
      }else{
	while(lgo2 == node[lastgrp].flag) ;  // wait till lastgrp done
      }
      node[group].barr = 0;                  // reset group count to zero
      node[group].flag = 1 - lgo;            // set group done flag
  }else{
      while(lgo == node[group].flag) ;       // wait till group done
  }
}

void node_barrier_single(int32_t id, int32_t maxcount){  // version with "flag flip"
  int lgo, count;

  if(maxcount < 2) return;     // trivial case, no need for fancy footwork

  lgo = barrs[0];
  count = __sync_fetch_and_add (barrs+16, 1);
  if(count == maxcount - 1){
    barrs[16] = 0;
    barrs[0] = 1 - barrs[0];
  }else{
    while(lgo == barrs[0]) ;
  }
}

#define POST_RCV 0x40000000
// static inline get_from(int me, int partner){
//   barrs[me] = (POST_RCV+partner);     // post request for partner's data
//   while(barrs[me] != 0);              // wait for partmer's acknowledge
// }
// 
// static inline put_to(int me, int partner){
//   while(barrs[partner] != (POST_RCV+me));  // wait till partner posts request from me
//   barrs[partner] = 0;                      // send data acknowledge to partner
// }
// 
// static inline get_put(int me, int partner){
//   barrs[me] = (POST_RCV+partner);          // post request for partner's data
//   while(barrs[partner] != (POST_RCV+me));  // wait till partner posts request from me
//   barrs[partner] = 0;                      // send data acknowledge to partner
//   while(barrs[me] != 0);                   // wait for partmer's acknowledge
// }

static inline get_from(int me, int partner){
  node[me].barr = (POST_RCV+partner);     // post request for partner's data
  while(node[me].barr != 0);              // wait for partmer's acknowledge
}

static inline put_to(int me, int partner){
  while(node[partner].barr != (POST_RCV+me));  // wait till partner posts request from me
  node[partner].barr = 0;                      // send data acknowledge to partner
}

static inline get_put(int me, int partner){
  node[me].barr = (POST_RCV+partner);          // post request for partner's data
  while(node[partner].barr != (POST_RCV+me));  // wait till partner posts request from me
  node[partner].barr = 0;                      // send data acknowledge to partner
  while(node[me].barr != 0);                   // wait for partmer's acknowledge
}

static inline get_put_m(int me, int src, int dst){
//   printf("(%d) send to %d, get from %d\n",me, src, dst);
  node[me].barr = (POST_RCV+src);          // post request for partner's data
  while(node[dst].barr != (POST_RCV+me));  // wait till partner posts request from me
  node[dst].barr = 0;                      // send data acknowledge to partner
  while(node[me].barr != 0);                   // wait for partmer's acknowledge
}

void node_barrier_2(int32_t id, int32_t maxcount){  // version with "messages"
  int mask, src, dst;

  if(maxcount < 2) return;     // trivial case, no need for fancy footwork

  mask = 1;
  while(mask < maxcount){
    dst = (id + mask) % maxcount;
    src = (id - mask + maxcount) % maxcount;
    get_put_m(id, src, dst);
    mask <<= 1;
  }
//   printf("(%d) DONE\n",id);
}

void node_barrier(int32_t id, int32_t maxcount){  // version with "messages"
  int mask, src, dst;

  if(maxcount < 2) return;     // trivial case, no need for fancy footwork

  mask = 1;
  while(mask < maxcount){
    dst = (id + mask) % maxcount;
    src = (id - mask + maxcount) % maxcount;
    get_put_m(id, src, dst);
    mask <<= 1;
  }
//   printf("(%d) DONE\n",id);
}

void node_barrier_0(int32_t id, int32_t maxcount){  // version with "messages"
  int powerof2, n, rest, gap, partner;

  if(maxcount < 2) return;     // trivial case, no need for fancy footwork

  powerof2 = 1;
  n = 2;
  while(n <= maxcount) {       // get largest power of 2 <= maxcount
    powerof2 = n ;
    n <<= 1 ; 
  }
  rest = maxcount - powerof2;  // the rest above largest power of 2

  if(id < powerof2){
    if(id < rest) {            // get from upper part partner if there is one
      partner = id + powerof2;
      get_from(id, partner);
    }
    for(gap=1 ; gap < powerof2 ; gap<<=1){
      partner = id ^ gap;
      get_put(id, partner);
    }
    if(id < rest) {            // put to upper part partner if there is one
      partner = id + powerof2;
      put_to(id, partner);
    }
  }else{
    partner = id - powerof2;
    put_to(id, partner);
    get_from(id, partner);
  }
}

void node_barrier1(int32_t id, int32_t maxcount){
  int32_t count;
  int32_t group;
  int32_t halfcount;
  int32_t grouplimit;

  if(maxcount < 2) return;      // barrier with 1 thread/process is a NO-OP

  count = __sync_fetch_and_add (barrs, 1);  // increment barrier count unpon entry
  if(count == maxcount -1){
    barrs[0] = 0;
  }
  while(barrs[0] != 0) ; // printf("count = %d\n",barrs[0]);
  return;

  group      = (id & 1) ^ 1;            // 0 (id is odd) or 1 (id is even)
  halfcount  = maxcount >> 1;
  grouplimit = halfcount + (maxcount & group) - 1;   // 1 less than group population
  // the even group population is one more than the odd group if maxcount is odd

  spins[group+2] = 0;
  count = __sync_fetch_and_add (barrs+group, 1);  // increment barrier count unpon entry
  if(count == grouplimit) {           // the last member of the (even/odd) group ?
    spins[group] = 1;                 // yes, set the group spin flag to true
    barrs[group] = 0;                 // and reset the barrier count of the group
  }
  while( (spins[0] & spins[1]) == 0); // wait until every thread/process has joined
}

void node_barrier2(int32_t id, int32_t maxcount){
  int32_t count;
  int32_t group;
  int32_t halfcount;
  int32_t grouplimit;

  if(maxcount < 2) return;      // barrier with 1 thread/process is a NO-OP

  count = __sync_fetch_and_add (barrs+1, 1);  // increment barrier count unpon entry
  if(count == maxcount -1){
    barrs[1] = 0;
  }
  while(barrs[1] != 0) ; // printf("count = %d\n",barrs[0]);
  return;

  group      = (id & 1) ^ 1;            // 0 (id is odd) or 1 (id is even)
  halfcount  = maxcount >> 1;
  grouplimit = halfcount + (maxcount & group) - 1;   // 1 less than group population
  // the even group population is one more than the odd group if maxcount is odd

  spins[group+4] = 0;
  count = __sync_fetch_and_add (barrs+group, 1);  // increment barrier count unpon entry
  if(count == grouplimit) {           // the last member of the (even/odd) group ?
    spins[group+2] = 1;                 // yes, set the group spin flag to true
    barrs[group] = 0;                 // and reset the barrier count of the group
  }
  while( (spins[2] & spins[3]) == 0); // wait until every thread/process has joined
}

void node_barrier3(int32_t id, int32_t maxcount){
  int32_t count;
  int32_t group;
  int32_t halfcount;
  int32_t grouplimit;

  if(maxcount < 2) return;      // barrier with 1 thread/process is a NO-OP

  count = __sync_fetch_and_add (barrs+2, 1);  // increment barrier count unpon entry
  if(count == maxcount -1){
    barrs[2] = 0;
  }
  while(barrs[2] != 0) ; // printf("count = %d\n",barrs[0]);
  return;

  group      = (id & 1) ^ 1;            // 0 (id is odd) or 1 (id is even)
  halfcount  = maxcount >> 1;
  grouplimit = halfcount + (maxcount & group) - 1;   // 1 less than group population
  // the even group population is one more than the odd group if maxcount is odd

  spins[group+0] = 0;
  count = __sync_fetch_and_add (barrs+group, 1);  // increment barrier count unpon entry
  if(count == grouplimit) {           // the last member of the (even/odd) group ?
    spins[group+4] = 1;                 // yes, set the group spin flag to true
    barrs[group] = 0;                 // and reset the barrier count of the group
  }
  while( (spins[4] & spins[5]) == 0); // wait until every thread/process has joined
}

int32_t wait_barrier(int barrier, int id, int maxcount){
  int32_t count;

  if(maxcount < 2) return 0;      // barrier with 1 thread/process is a NO-OP

  if(spins == NULL || barrs == NULL) return -1;

  if(id ==0) spins[barrier] = 0;  // lowest id thread/process sets spin flag to false
				  // initial state of spin flags MUST BE 0
  count = __sync_fetch_and_add (barrs+barrier, 1);
  if(count == maxcount-1) {
    barrs[barrier] = 0;            // last one to arrive resets the barrier count to zero
    spins[barrier] = 1;            // and sets the global spin flag to true
  }else{
    while(spins[barrier] == 0);
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
  return (__sync_val_compare_and_swap(&(node[lock].lock), me+1, 0) == me+1);
}

int32_t acquire_lock(int lock, int me){
  int count = 0;

  if(spins == NULL || barrs == NULL) return -1;

//  printf("attempting acquisition of lock %d owned by %d\n",lock,locks[lock]-1);
  while(__sync_val_compare_and_swap(&(node[lock].lock), 0, me+1) != 0) {
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
  size =globalsize *4096;   // 4 KBytes per process

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
    if(globalrank == 0) printf("lock min, max, avg = %9f, %9f, %9f\n",tmin,tmax,tavg);
    t0 = rdtsc();
    for(i=0 ; i<100 ; i++){
      node_barrier_multi(localrank, localsize);
      node_barrier_multi(localrank, localsize);
      node_barrier_multi(localrank, localsize);
    }
    t1 = rdtsc();
//     printf("barrier time = %d\n",(t1-t0)/300);
    tmp = (t1-t0)/300;
    ierr = MPI_Allreduce(&tmp,&tmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&tmp,&tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&tmp,&tavg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    tavg = tavg / globalsize;
    if(globalrank == 0) printf("SMP barrier min, max, avg = %9f, %9f, %9f\n",tmin,tmax,tavg);

    t0 = rdtsc();
    for(i=0 ; i<100 ; i++){
      ierr = MPI_Barrier(MY_World);
      ierr = MPI_Barrier(MY_World);
      ierr = MPI_Barrier(MY_World);
    }
    t1 = rdtsc();
//     printf("barrier time = %d\n",(t1-t0)/300);
    tmp = (t1-t0)/300;
    ierr = MPI_Allreduce(&tmp,&tmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&tmp,&tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&tmp,&tavg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    tavg = tavg / globalsize;
    if(globalrank == 0) printf("MPI barrier min, max, avg = %9f, %9f, %9f\n",tmin,tmax,tavg);

  ierr = MPI_Finalize();
}
#endif
