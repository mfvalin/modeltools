/* 
 * Copyright (C) 2018  Environnement Canada
 *
 * This is free software; you can redistribute it and/or
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
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/types.h>
#include <sys/shm.h>
#include <circular_buffer.h>

#define SPACE_AVAILABLE(in,out,limit)  ((in < out) ? out-in-1 : limit-in+out-1)

#define DATA_AVAILABLE(in,out,limit)  ((in > out) ? in-out : limit-out+in-1)

// the Fortran interfaces
#if defined(THIS_BETTER_NEVER_BE_TRUE)
interface
  function circular_buffer_init(p, nbytes) result(limit) bind(C,name='circular_buffer_init')
    import :: C_PTR, C_INT
    type(C_PTR), intent(IN) :: p
    type(C_INT), intent(IN) :: nbytes
  end function circular_buffer_init
end interface
#endif

static inline void move_integers(int *dst, int*src, int n){
  memcpy(dst, src, sizeof(int)*n);
}

// initialize a circular buffer
// nbytes is the size in bytes of the memory area
// return 0 upon success, -1 otherwise
int circular_buffer_init(circular_buffer_p p, int32_t nbytes){
  if(nbytes < 4096 || p == NULL) return -1;   // area is too small
  p->m.first = 0;
  p->m.in    = 0;
  p->m.out   = 0;
  p->m.limit = (nbytes - sizeof(fiol_management)) / sizeof(int);
  return 0;
}

// create and initialize a circular buffer of size nbytes in "shared memory"
// return the "shared memory segment" address of the circular buffer upon success, NULL otherwise
// shmid will be set to the shared memory id of the "shared memory segment upon success, -1 otherwise
circular_buffer_p circular_buffer_create_shared(int32_t *shmid, int32_t nbytes){
  void *t;
  size_t sz = nbytes;
  int id;
  struct shmid_ds ds;
  int status;

  *shmid = -1;
  if(nbytes < 64*1024) return NULL;
  id = shmget(IPC_PRIVATE, sz, IPC_CREAT);   // create shared memory segment
  if(id == -1) return NULL;                  // error occurred
  t = shmat(id, NULL, 0);                    // attach shared memory segment
  if( t == (void *) -1) return NULL;         // error occurred
  status = shmctl(id, IPC_RMID, &ds);        // mark segment for deletion (ONLY SAFE ON LINUX)
  if(status != 0) return NULL;               // this should not fail
  *shmid = id;
  return (circular_buffer_init((circular_buffer_p)t, nbytes) != 0) ?NULL : (circular_buffer_p)t;
}

// detach from a "shared memory segment" circular buffer 
// return 0 upon success, -1 otherwise
int circular_buffer_dertach_shared(circular_buffer_p p){
  return shmdt(p) ;   // detach from "shared memory segment" creeated by circular_buffer_create_shared
}

// create and initialize a circular buffer of size nbytes
// return the address of the circular buffer upon success, NULL otherwise
circular_buffer_p circular_buffer_create(int32_t nbytes){
  circular_buffer_p t;
  size_t sz = nbytes;

  if(nbytes < 4096) return NULL;
  t = (circular_buffer_p ) malloc(sz);
  return (circular_buffer_init(t, nbytes) != 0) ? NULL : t;
}

// returns the current number of empty slots available
int circular_buffer_space_available(circular_buffer_p p){
  int  *inp = &(p->m.in);
  int  *outp = &(p->m.out);
  int in, out, limit;

  limit = p->m.limit;
  in = *inp;
  out = *outp;
  return SPACE_AVAILABLE(in,out,limit);
}

// wait until at least n empty slots are available for inserting data
// returns the actual number of empty slots available
int circular_buffer_wait_space_available(circular_buffer_p p, int n){
  int volatile *inp = &(p->m.in);
  int volatile *outp = &(p->m.out);
  int in, out, limit, navail;

  limit = p->m.limit;
  navail = 0;
  while(navail <n){
    in = *inp;
    out = *outp;
    navail = SPACE_AVAILABLE(in,out,limit);
  }
  return navail;
}

// returns the current number of data tokens available
int circular_buffer_data_available(circular_buffer_p p){
  int  *inp = &(p->m.in);
  int  *outp = &(p->m.out);
  int in, out, limit;

  limit = p->m.limit;
  in = *inp;
  out = *outp;
  return DATA_AVAILABLE(in,out,limit);
}

// returns a pointer to the  start of the data buffer
// useful in conjunction with circular_buffer_data_snoop
int32_t *circular_buffer_data_buffer(circular_buffer_p p){
  return  p->data;
}

// get the address of the first position in the circular data buffer
int32_t *circular_buffer_start(circular_buffer_p p){
  return p->data;  // start of data buffer
}

// returns a pointer to the "in" position, assumes that the caller knows the start of data buffer
// n1 slots available at "in", n2 slots available at "start"
int32_t *circular_buffer_advance_in(circular_buffer_p p, int32_t *n1, int32_t *n2){
  int  *inp = &(p->m.in);
  int  *outp = &(p->m.out);
  int in, out, limit;

  limit = p->m.limit;
  in = *inp;
  out = *outp;
  if(in < out) {
    *n1 = out - in - 1; // available at "in"
    *n2 = 0;            // nothing available at beginning of buffer
  }else{                // in >= out
    *n1 = limit - in;   // "in" -> "limit"
    *n2 = out;          // available at beginning of buffer (technically out - first)
  }
  return p->data+in;
}

// returns a pointer to the "out" position, assumes that the caller knows the start of data buffer
// n1 tokens available at "out", n2 tokens available at "start"
int32_t *circular_buffer_advance_out(circular_buffer_p p, int32_t *n1, int32_t *n2){
  int  *inp = &(p->m.in);
  int  *outp = &(p->m.out);
  int in, out, limit;

  limit = p->m.limit;
  in = *inp;
  out = *outp;
  if(in < out) {
    *n1 = limit - out;  // available after "out"
    *n2 = in;           // available at beginning of buffer (technically in - first)
  }else{
    *n1 = in - out;     // "out" -> "in"
    *n2 = 0;            // nothing at beginning of buffer
  }
  return p->data+out;
}

// wait until at least n data tokens are available for extracting data
// returns the actual number of data tokens available
int circular_buffer_wait_data_available(circular_buffer_p p, int n){
  int volatile *inp = &(p->m.in);
  int volatile *outp = &(p->m.out);
  int in, out, limit, navail;

  limit = p->m.limit;
  navail = 0;
  while(navail <n){
    in = *inp;
    out = *outp;
    navail = DATA_AVAILABLE(in,out,limit);
  }
  return navail;
}


// atomic extraction of n tokens into the dst array
// returns the number of data tokens available after this operation
int circular_buffer_atomic_get(circular_buffer_p p, int *dst, int n){
  int volatile *inp = &(p->m.in);
  int volatile *outp = &(p->m.out);
  int *buf = p->data;
  int in, out, limit, navail, ni;

  // wait until enough data is available
  limit = p->m.limit;
  navail = 0; in = 0 ; out = 0;
  while(navail <n){
    in = *inp;
    out = *outp;
    navail = DATA_AVAILABLE(in,out,limit);
  }

  if(out < in){         // 1 segment
    move_integers(dst, buf+out, n);
    out += n;
  }else{                // 1 or 2 segments
    ni = n > (limit-out) ? (limit-out) : n;
    move_integers(dst, buf+out, ni);
    n -= ni;
    out += ni;
    dst += ni;
    if(out >= limit) out = 0;
    move_integers(dst, buf+out, n);
    out += n;
  }
  *outp = out;
  in = *inp;
  out = *outp;
  return DATA_AVAILABLE(in,out,limit);
}

// atomic insertion of n tokens from the src array
// returns the number of empty slots available after this operation
int circular_buffer_atomic_put(circular_buffer_p p, int *src, int n){
  int volatile *inp = &(p->m.in);
  int volatile *outp = &(p->m.out);
  int *buf = p->data;
  int in, out, limit, navail, ni;

  // wait until there is enough room to insert data
  limit = p->m.limit;
  navail = 0; in = 0 ; out = 0;
  while(navail <n){
    in = *inp;
    out = *outp;
    navail = SPACE_AVAILABLE(in,out,limit);
  }

  if(in < out){         // 1 segment
    move_integers(buf+in, src, n);
    in += n;
  }else{                // 1 or 2 segments
    ni = n > (limit-in) ? (limit-in) : n;
    move_integers(buf+in, src, ni);
    n -= ni;
    in += ni;
    src += ni;
    if(in >= limit) in = 0;
    move_integers(buf+in, src, n);
    in += n;
  }
  *inp = in;
  in = *inp;
  out = *outp;
  return SPACE_AVAILABLE(in,out,limit);
}
