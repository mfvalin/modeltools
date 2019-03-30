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
#include <string.h>
#include <circular_buffer.h>

#define SPACE_AVAILABLE(in,out,limit)  ((in < out) ? out-in-1 : limit-in+out-1)

#define DATA_AVAILABLE(in,out,limit)  ((in > out) ? in-out : limit-out+in-1)

static inline void move_integers(int *dst, int*src, int n){
  memcpy(dst, src, sizeof(int)*n);
}

// initialize a circular buffer
// nbytes is the size in bytes of the memory area
// returns the  maximum number of tokens that can fit in the circular buffer
int circular_buffer_init(circular_buffer_p p, int32_t nbytes){
  if(nbytes < 4096) return -1;   // area is too small
  p->m.first = 0;
  p->m.in    = 0;
  p->m.out   = 0;
  p->m.limit = (nbytes - sizeof(fiol_management)) / sizeof(int);
  return p->m.limit - 1;   // maximum number of tokens that the circular buffer may contain
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
