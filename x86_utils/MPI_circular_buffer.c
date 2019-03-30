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

#define MAX_WINDOWS 16

typedef struct{
  fiol_management   *rtab; // management table for remote circular buffers
  circular_buffer_p *ltab; // pointers to local circular buffers
  void *winbase;           // base address of local communication window
  uint32_t *winsize;       // table of partner window sizes
  MPI_Win c_win;           // C window
  MPI_Win c_comm;          // C communicator
  MPI_Fint f_win;          // Fortran window
  MPI_Fint f_comm;         // Fortran communicator
  uint32_t my_size;        // size of local window
} remote_window;

// create and initialize a remote set of circular buffers
// on each PE, a circular buffer of size nbytes will be created for each member of the communicator
// the function return a pointer ot an opaque object describing the buffer set
// sizes and offsets are in int32_t units
// the way things are set up, a PE could send a remote message to itself, buffers are allocated
// and control information is set appropriately
remote_window *MPI_circular_buffer_create(MPI_Fint comm, int32_t nbytes){
  MPI_Aint fullsize;
  void *baseptr;
  int i, ier ;
  MPI_Win c_win;
  MPI_Comm c_comm = MPI_Comm_f2c(comm);
  remote_window *rw;
  int comm_rank = 0;
  int comm_size = 0;

  ier = MPI_Comm_size(c_comm, &comm_size);   // size of communicator (number of partners)
  if(ier != MPI_SUCCESS) return NULL;
  ier = MPI_Comm_rank(c_comm, &comm_rank);   // self rank in communicator
  if(ier != MPI_SUCCESS) return NULL;

  fullsize = nbytes * comm_size;             // nbytes per partner
  ier = MPI_Alloc_mem(fullsize, MPI_INFO_NULL, &baseptr);
  if(ier != MPI_SUCCESS) return NULL;        // allocation failed
  rw->winsize = (uint32_t *) malloc(comm_size * sizeof(uint32_t));  // table of sizes (size may differ on each PE)
  if(rw->winsize == NULL) return NULL;       // allocation failed

  rw = (remote_window *) malloc(sizeof(remote_window));  // allocate remote window object
  if(rw == NULL) return NULL;                // allocation failed
  rw->winbase = baseptr;                     // base address for communication window
  rw->my_size = nbytes / sizeof(int32_t);    // my size (in int32_t units)

  rw->rtab = (fiol_management_p)   malloc(sizeof(fiol_management)   * comm_size);  // table of remote partners control info
  if(rw->rtab == NULL) return NULL;          // allocation failed
  rw->ltab = (circular_buffer_p *) malloc(sizeof(circular_buffer_p) * comm_size);  // table of pointers to local circular buffers
  if(rw->ltab == NULL) return NULL;          // allocation failed
  
  ier = MPI_Allgather(&rw->my_size, 1, MPI_INTEGER, rw->winsize, 1, MPI_INTEGER, c_comm);  // get size info from all partners
  if(ier != MPI_SUCCESS) return NULL;
  for(i=0 ; i<comm_size-1 ; i++){            // compute first, in, out, limit for partners and self using size
    rw->rtab[i].first = 0;                   // initialize control information for my remote circular buffer on partner i
    rw->rtab[i].in    = 0;
    rw->rtab[i].out   = 0;
    rw->rtab[i].limit = rw->winsize[i] - sizeof(fiol_management)/sizeof(int32_t);
    rw->ltab[i] = (circular_buffer_p) ( (int32_t *) baseptr + (i * rw->my_size) );  // pointer to local circular buffer i
    rw->ltab[i]->m.first = 0;                // initialize control information for local circular buffer i
    rw->ltab[i]->m.in    = 0;
    rw->ltab[i]->m.out   = 0;
    rw->ltab[i]->m.limit = rw->my_size - sizeof(fiol_management)/sizeof(int32_t);
    rw->winsize[i] = rw->winsize[i] * comm_rank;   // transform into offset into window for my buffer on partner i
  }
  // base offset for my buffer on partner of rank i is ( offset = rw->winsize[i] )
  // data buffer will be located at offset + 4 ,
  // first at offset + 0, in at offset + 1, out at offset + 2, limit at offset + 3
  // first is expected to be zero.
  // a later MPI_Get of 4 integers at offset + 0 wil get current values of first, in , out, limit

  ier = MPI_Win_create(baseptr, fullsize, 4, MPI_INFO_NULL, c_comm, &c_win);  // create window
  if(ier != MPI_SUCCESS) return NULL;

  rw->c_win = c_win;                 // C window 
  rw->f_win = MPI_Win_c2f(c_win);    // Fortran window

  return (remote_window *) rw;
}
