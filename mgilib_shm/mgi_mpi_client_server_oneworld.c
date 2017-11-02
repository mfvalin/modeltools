/* useful routines for C and FORTRAN programming
 * Copyright (C) 2017  Division de Recherche en Prevision Numerique, Environnement Canada
 * Copyright (C) 2017  Centre ESCER UQAM
 *
 * This code is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this code; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#include "mpi.h" 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define TAG 123456
#define MAX_DATA 128

int main( int argc, char **argv ) 
{
  int size, rank;
  MPI_Comm mpi_ichan, mpi_chan;
  MPI_Win window;
  MPI_Aint winsize;
  MPI_Status status;
  MPI_Aint TargetDisp;
  int dispunit = sizeof(int);
  int memptr[1024*1024];
  int handshake, answer, tag, i;
  int local[1024];
  

  MPI_Init( &argc, &argv ); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("I am process %d of %d\n",rank+1,size);
  winsize = 1024*1024*dispunit;
  tag = 12345;

  if(rank == 0) {  // server
    handshake = rank;
    MPI_Intercomm_create(MPI_COMM_SELF, 0, MPI_COMM_WORLD, size/2, TAG, &mpi_ichan);
    MPI_Intercomm_merge(mpi_ichan, 0, &mpi_chan);
    MPI_Comm_size(mpi_chan, &size);
    MPI_Comm_rank(mpi_chan, &rank);
    printf("server has rank %d of %d\n",rank,size);
    MPI_Barrier(mpi_chan);
    MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, mpi_chan, &window);
    printf("server window created\n");
    for (i=0 ; i < winsize/dispunit ; i++) memptr[i] = (-i) ;
    printf("server[1] = %d, server[1024] = %d\n",memptr[1],memptr[1024]);
    printf("server window initialized\n");
    MPI_Barrier(mpi_chan);
    MPI_Sendrecv(&handshake, 1, MPI_INTEGER, 1-rank, tag, &answer, 1, MPI_INTEGER, 1-rank, tag, mpi_chan, &status);
    printf("server has received answer = %d\n",answer);
    MPI_Barrier(mpi_chan);
    printf("server[1] = %d, server[1024] = %d\n",memptr[1],memptr[1024]);
    MPI_Win_free(&window);
  }
  if(rank == size/2) {  // client
    handshake = rank;
    MPI_Intercomm_create(MPI_COMM_SELF, 0, MPI_COMM_WORLD, 0, TAG, &mpi_ichan);
    MPI_Intercomm_merge(mpi_ichan, 1, &mpi_chan);
    MPI_Comm_size(mpi_chan, &size);
    MPI_Comm_rank(mpi_chan, &rank);
    printf("client has rank %d of %d\n",rank,size);
    MPI_Barrier(mpi_chan);
    MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, mpi_chan, &window);
    for (i=0 ; i<1024 ; i++) local[i] = i;
    printf("client window created\n");
    MPI_Barrier(mpi_chan);
    MPI_Sendrecv(&handshake, 1, MPI_INTEGER, 1-rank, tag, &answer, 1, MPI_INTEGER, 1-rank, tag, mpi_chan, &status);
    printf("client has received answer = %d\n",answer);
    TargetDisp = 1;
    MPI_Win_lock(MPI_LOCK_SHARED,1-rank,0,window);
    MPI_Put(local, 1024, MPI_INTEGER, 1-rank, TargetDisp, 1024, MPI_INTEGER, window);
    MPI_Win_unlock(1-rank,window);
    MPI_Barrier(mpi_chan);
    MPI_Win_free(&window);
  }

  
  MPI_Finalize();
}