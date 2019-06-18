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

//****P* librkl/RPN kernel library locks and barriers
// DESCRIPTION
// set of routines to implement locks and barriers that can be used between threads belonging
// to the same process or between processes running on the same shared memory node (MPI, etc...)
//
#if defined(NEVER_EVER_TRUE)
// EXAMPLES
  =============  increment under lock =============

  static volatile int my_lock = 0;
  int kount;
  ......
  AcquireLock(&my_lock);  // acquire lock
  kount++;                // increment variable under lock
  ReleaseLock(&my_lock);  // release lock

  ......
  AcquireIdLock(&my_lock, id);  // acquire lock number id (0 -> maxitems-1 , see setup_locks_and_barriers)
  kount++;                // increment variable under lock
  ReleaseIdLock(&my_lock, id);  // release lock number id

  =============  get a block of shared memory =============

  void *ptr
  int sid, localrank;
  int size = ....;
    .....
    
  ierr = MPI_Comm_rank(MY_World,&localrank);         // rank of this PE on this SMP node
  if(localrank == 0){
    ptr = allocate_safe_shared_memory(&sid, size);   // create a shared memory segment and get it's address
    ierr = MPI_Bcast(&sid,1,MPI_INTEGER,0,MY_World); // broadcast the segment's id
  }else{                                             // other PEs
    ierr = MPI_Bcast(&sid,1,MPI_INTEGER,0,MY_World); // get the segment's id
    ptr = shmat(sid,NULL,0);                         // get the shared memory segment's address
  }
  ierr = MPI_Barrier(MY_World);                      // wait for all PEs to have attached the segment

  =============  setup of a block of shared memory for barriers and locks =============

  int *next;
  next = (int *) setup_locks_and_barriers(ptr, size, 16);  // 16 locks, 16 barriers, point next to area that follows
  // only on rank 0 normally

  =============  an intra-node barrier =============

  ierr = MPI_Comm_rank(MY_World,&localsize);         // how many PEs in MY_World on this SMP node
    .....
  simple_node_barrier(0, localsize);                 // use barrier 0 (0 -> maxitems-1 , see setup_locks_and_barriers)
//****
#endif

#include <stdlib.h>
#include <stdint.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/stat.h>

#include <rkl_sync.h>

// interface   !InTf!

// default number of locks/barriers = 16
#define MAXLOCKID 15
static uint32_t maxlockid = MAXLOCKID;

// locks
static uint32_t volatile lock0[MAXLOCKID+1] __attribute__ ((aligned (64)));
static uint32_t volatile *locks;       // in case we want to allocate to another address

// barrier (partial) done flags
static uint32_t volatile barrs0[MAXLOCKID+1] __attribute__ ((aligned (64)));
static uint32_t volatile *barrs;       // in case we want to allocate to another address

// barrier (partial) counts
static uint32_t volatile count0[MAXLOCKID+1] __attribute__ ((aligned (64)));
static uint32_t volatile *count;       // in case we want to allocate to another address

//****f* librkl/allocate_safe_shared_memory
// Synopsis
// allocate a  shared memory segment of size Bytes 
// segment will be preventively marked for deletion (only works on linux)
// return address and id of shared memory segment (in sid)
// in case of error id is <0 and the returned address is NULL
//
// ARGUMENTS
void *allocate_safe_shared_memory(int32_t *sid, uint32_t size)   // !InTc!
//****
{
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

//****f* librkl/setup_locks_and_barriers
// Synopsis
// use memory pointed to by p to set up control tables for maxitems locks and barriers
// if there is not enough memory, return NULL
//
// p        : base address
// size     : amount of available memory in bytes
// maxitems : number of locks/barriers to allocate (will be rounded down to a multiple of 16)
// return the next available address (above tables)
//
// ARGUMENTS
void *setup_locks_and_barriers(void *p, uint32_t size, uint32_t maxitems)   // !InTc!
//****
{
  maxitems = (maxitems >> 4) << 4;  // multiple of 16 x 4 byte items
  if(size < (maxitems * 3 * sizeof(uint32_t)) ) return NULL;   // not enough memory for tables, OOPS
  maxlockid = maxitems - 1;         // maximum index
  locks = (uint32_t *) p;
  barrs = locks + maxitems;
  count = barrs + maxitems;
  return (void *)(count + maxitems);     // next free address (after lock/barrier tables)
}

//****f* librkl/simple_node_barrier
// Synopsis
// implement a barrier between threads or processes (on the same SMP node)
// unsophisticated version (uses "flag flip")
// id       : barrier (0 -> maxitems-1) is to be used (see setup_locks_and_barriers)
// maxcount : number of threads/processes for this barrier
//
// the function will return 0 upon success, 1 in case of error (invalid id)
//
// Fortran interface
//   function simple_node_barrier(id, maxcount) result(status)   !InTf!
//     import :: C_INT                                           !InTf!
//     integer(C_INT), intent(IN), value :: id, maxcount         !InTf!
//     integer(C_INT) :: status                                  !InTf!
//   end function simple_node_barrier                            !InTf!
// ARGUMENTS
uint32_t simple_node_barrier(int32_t id, int32_t maxcount)   // !InTc!
//****
{
  int lgo, cnt;

  if(id < 0 || id >  maxlockid) return 1;  // invalid id
  if(maxcount == 1) return 0;              // trivial case, no need for fancy footwork

  lgo = barrs[id];                         // save current flag value
  cnt = __sync_fetch_and_add (count+id, 1);   // increment barrier count
  if(cnt == maxcount - 1){                 // last participant arriving at barrier
    cnt = __sync_fetch_and_add (count+id, -maxcount);  // reset count to zero
    barrs[id] = 1 - barrs[id];             // invert flag
  }else{                                   // all other participants
    while(lgo == barrs[id]) ;   // wait while flag is not inverted
  }
  return 0;
}

#if defined(TEST_MULTI)
// more sophisticated version with 2 level "flag flip" (may be faster for large counts)
// id is this PEs identifier (0 < maxcount)
uint32_t node_barrier_multi(int32_t id, int32_t maxcount)   // !InTc!
{  
  int lgo, lgo2, cnt;
  int group, base, allgrps, grpcnt, maxid;

  if(maxcount == 1) return 0;            // trivial case, no need for fancy footwork

  maxid = maxcount -1 ;                  // 0 -> maxcount -1
  allgrps = (maxid >> 3) + 1 ;           // ordinal of global group (last group + 1)
  if(allgrps  >  maxlockid) return 1;    // OOPS, too many groups
  group = id >> 3 ;                      // which group of 8 does id belong to ? 
  base = group << 3 ;                    // group number * 8 (lowest id of group)
  grpcnt = maxid - base ;
  if(grpcnt > 7) grpcnt = 7 ;            // population of my group - 1 (max 7)
  

  lgo  = barrs[group];                   // flag to signal that my group is done
  lgo2 = barrs[allgrps];                 // flag to signal that all groups are done
  cnt = __sync_fetch_and_add (&count[group], 1);        // bump my group count
  if(cnt == grpcnt ){                    // last pe in a group
      cnt = __sync_fetch_and_add (&count[allgrps], 1);  // bump allgrps count
      if(cnt == allgrps - 1){                // all groups are done, i am the last one
	count[allgrps] = 0;              // reset allgrps count to zero
	barrs[allgrps] = 1 - lgo2;       // set allgrps done flag
      }else{
	while(lgo2 == barrs[allgrps]) ;  // wait till all groups done
      }
      count[group] = 0;                  // reset my group count to zero
      barrs[group] = 1 - lgo;            // set my group done flag
  }else{
      while(lgo == barrs[group]) ;       // wait till my group done
  }
  return 0;
}
#endif

//****f* librkl/AcquireLock
// Synopsis
// acquire a lock, using 4 byte variable at address lock
// the variable pointed to by lock must have been initialized to zero
// lock   : address of lock variable
//
// ARGUMENTS
void AcquireLock(volatile int32_t *lock)   // !InTc!
//****
{
  acquire_lock(lock) ;
}

//****f* librkl/AcquireIdLock
// Synopsis
// acquire a lock, using 4 byte variable at address lock
// the variable pointed to by lock must have been initialized to zero
// and will be set to id when lock is acquired
// lock   : address of lock variable
// id     : identifier for this thread/process
//
// ARGUMENTS
void AcquireIdLock(volatile int32_t *lock, int32_t id)   // !InTc!
//****
{
  acquire_idlock(lock, id) ;
}

//****f* librkl/ReleaseLock
// Synopsis
// release a lock acquired via AcquireLock using 4 byte variable at address lock
// attempting to release a lock that is not acquired will result in a deadlock
//
// ARGUMENTS
void ReleaseLock(volatile int32_t *lock)   // !InTc!
//****
{
  release_lock(lock) ;
}

//****f* librkl/ReleaseIdLock
// Synopsis
// release a lock acquired via AcquireIdLock using 4 byte variable at address lock
// attempting to release a lock that was not acquired with this id will result in a deadlock
// lock   : address of lock variable
// id     : identifier for this thread/process
//
// ARGUMENTS
void ReleaseIdLock(volatile int32_t *lock, int32_t id)   // !InTc!
//****
{
  release_idlock(lock, id) ;
}

// end interface   !InTf!

#if defined(SELF_TEST)

static uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv){
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
  uint64_t t0 , t1;
  double tmin, tmax, tavg, tmp;
  volatile int32_t *kount;
  double tall[2048];

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

  if(localrank == 0){                                // root PE
    ptr = allocate_safe_shared_memory(&sid, size);   // create a shared memory segment and get it's address
    ierr = MPI_Bcast(&sid,1,MPI_INTEGER,0,MY_World); // broadcast the segment's id
  }else{                                             // other PEs
    ierr = MPI_Bcast(&sid,1,MPI_INTEGER,0,MY_World); // get the segment's id
    ptr = shmat(sid,NULL,0);                         // get the shared memory segment's address
  }
  ierr = MPI_Barrier(MY_World);                      // wait for all PEs to have attached the segment

  kount = (int *) setup_locks_and_barriers(ptr, size, 16);
  kount[1] = 0;
  ierr = MPI_Barrier(MY_World);                // kount[1] is now available to all PEs

// performance and correctness test of increment under lock
  t0 = rdtsc();
  for(i=0 ; i<100 ; i++){
    AcquireLock(kount);
    kount[1]++;
    ReleaseLock(kount);  // release_idlock
  }
  t1 = rdtsc();
  ierr = MPI_Barrier(MY_World);   // wait until all PEs are done
  tmp = (t1-t0)/100;
  ierr = MPI_Allreduce(&tmp,&tmin,1,MPI_DOUBLE,MPI_MIN,MY_World);  // collect minimum value
  ierr = MPI_Allreduce(&tmp,&tmax,1,MPI_DOUBLE,MPI_MAX,MY_World);  // collect maximum value
  ierr = MPI_Allreduce(&tmp,&tavg,1,MPI_DOUBLE,MPI_SUM,MY_World);  // sum of timings to compute average
  ierr = MPI_Gather(&tmp,1,MPI_DOUBLE,&tall,1,MPI_DOUBLE,0,MY_World);
  tavg = tavg / localsize;
  if(kount[1] != 100 * localsize) {
    printf("ERROR: kount = %d, expected = %d\n",kount[1], 100 * localsize);
    exit(1);
  }
  if(globalrank == 0){
    printf("lock min, max, avg = %6.0f, %6.0f, %6.0f, kount = %d\n",tmin,tmax,tavg,kount[1]);
    for(i = 0 ; i < localsize ; i++){
      printf("%6.0f ",tall[i]);
    }
    printf("\n\n");
  }

// performance test 
  t0 = rdtsc();
  for(i=0 ; i<100 ; i++){
    simple_node_barrier(0, localsize);  // 3 barriers in a row to really stress things
    simple_node_barrier(0, localsize);
    simple_node_barrier(0, localsize);
  }
  t1 = rdtsc();
  tmp = (t1-t0)/300;
  ierr = MPI_Allreduce(&tmp,&tmin,1,MPI_DOUBLE,MPI_MIN,MY_World);
  ierr = MPI_Allreduce(&tmp,&tmax,1,MPI_DOUBLE,MPI_MAX,MY_World);
  ierr = MPI_Allreduce(&tmp,&tavg,1,MPI_DOUBLE,MPI_SUM,MY_World);
  ierr = MPI_Gather(&tmp,1,MPI_DOUBLE,&tall,1,MPI_DOUBLE,0,MY_World);
  tavg = tavg / globalsize;
  if(globalrank == 0){
    printf("simple SMP barrier min, max, avg = %6.0f, %6.0f, %6.0f\n",tmin,tmax,tavg);
    for(i = 0 ; i < localsize ; i++){
      printf("%6.0f ",tall[i]);
    }
    printf("\n\n");
  }

#if defined(TEST_MULTI)
test of more complex barrier algorithm
  t0 = rdtsc();
  for(i=0 ; i<100 ; i++){
    node_barrier_multi(localrank, localsize);
    node_barrier_multi(localrank, localsize);
    node_barrier_multi(localrank, localsize);
  }
  t1 = rdtsc();
  tmp = (t1-t0)/300;
  ierr = MPI_Allreduce(&tmp,&tmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  ierr = MPI_Allreduce(&tmp,&tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ierr = MPI_Allreduce(&tmp,&tavg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  tavg = tavg / globalsize;
  if(globalrank == 0) printf("multi SMP barrier min, max, avg = %6.0f, %6.0f, %6.0f\n",tmin,tmax,tavg);
#endif

  t0 = rdtsc();
  for(i=0 ; i<100 ; i++){
    ierr = MPI_Barrier(MY_World);
    ierr = MPI_Barrier(MY_World);
    ierr = MPI_Barrier(MY_World);
  }
  t1 = rdtsc();
  tmp = (t1-t0)/300;
  ierr = MPI_Allreduce(&tmp,&tmin,1,MPI_DOUBLE,MPI_MIN,MY_World);
  ierr = MPI_Allreduce(&tmp,&tmax,1,MPI_DOUBLE,MPI_MAX,MY_World);
  ierr = MPI_Allreduce(&tmp,&tavg,1,MPI_DOUBLE,MPI_SUM,MY_World);
  tavg = tavg / globalsize;
  if(globalrank == 0) printf("MPI barrier min, max, avg = %6.0f, %6.0f, %6.0f\n",tmin,tmax,tavg);

  ierr = MPI_Finalize();
  if(ierr != MPI_SUCCESS) return 1;
  return 0;
}
#endif
