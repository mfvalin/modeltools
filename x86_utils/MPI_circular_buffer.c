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
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <circular_buffer.h>
#include <mpi.h>

typedef struct{    // values from partners
  int first;
  int in;
  int out;
  int limit;
}remote_entry;

typedef struct{
  void *winbase;
  MPI_Win c_win;
  MPI_Fint f_win;
  uint32_t winsize;
  remote_entry *rtab;
} remote_window;

typedef struct{
  remote_window r;
  circular_buffer c;
} mpi_circular_buffer;

// create and initialize a collective set of circular buffers
mpi_circular_buffer *MPI_circular_buffer_init(MPI_Fint comm, int32_t nbytes){
  MPI_Aint fullsize = nbytes;
  void *baseptr;
  int i, ier ;
  MPI_Aint offset0;
  MPI_Win c_win;
  MPI_Comm c_comm = MPI_Comm_f2c(comm);
  mpi_circular_buffer *rw;
  int comm_rank = 0;
  int comm_size = 0;

  ier = MPI_Comm_rank(c_comm, &comm_rank);
  if(ier != MPI_SUCCESS) return NULL;
  ier = MPI_Comm_size(c_comm, &comm_size);
  if(ier != MPI_SUCCESS) return NULL;

  ier = MPI_Alloc_mem(fullsize, MPI_INFO_NULL, &baseptr);
  if(ier != MPI_SUCCESS) return NULL;

  rw = (mpi_circular_buffer *) baseptr;
  rw->r.winbase = baseptr;
  rw->r.winsize = nbytes;
  rw->r.rtab = (remote_entry *) malloc(sizeof(remote_entry) * comm_size);  // table to get info from remote partners
  if(rw->r.rtab == NULL) return NULL;

  offset0 = sizeof(remote_window) / sizeof(int32_t);  // offset of circular buffer itself into window
  rw->c.first = offset0 + 4;
  rw->c.in    = offset0 + 4;
  rw->c.out   = offset0 + 4;
  rw->c.limit = offset0 + 4 + (nbytes - sizeof(mpi_circular_buffer)) / sizeof(int32_t);

  ier = MPI_Win_create(baseptr, fullsize, 4, MPI_INFO_NULL, c_comm, &c_win);  // create window
  if(ier != MPI_SUCCESS) return NULL;

  rw->r.c_win = c_win;
  rw->r.f_win = MPI_Win_c2f(c_win);
  for(i=0 ; i<comm_size-1 ; i++){                   // get first, in, out, limit from partners
    if(i != comm_rank){       // not ME
      ier = MPI_Win_lock(MPI_LOCK_SHARED,i,0,c_win);
      if(ier != MPI_SUCCESS) return NULL;
      ier = MPI_Get(rw->r.rtab + i, 4, MPI_INTEGER, i, offset0, 4, MPI_INTEGER, c_win);
      if(ier != MPI_SUCCESS) return NULL;
      ier = MPI_Win_unlock(i,c_win);
      if(ier != MPI_SUCCESS) return NULL;
    }else{                    // ME
      rw->r.rtab[i].first = rw->c.first;
      rw->r.rtab[i].in    = rw->c.in;
      rw->r.rtab[i].out   = rw->c.out;
      rw->r.rtab[i].limit = rw->c.limit;
    }
  }
  ier = MPI_Barrier(c_comm);  // everybody has all the info about everybody
  if(ier != MPI_SUCCESS) return NULL;
  return rw;
}
