#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv){
  int rank, size, nw;
  MPI_Win window;
  int buffer[4096];
  int wbuf[1024*1024];
  int i, j;
  MPI_Aint TargetDisp, winsize;
  int index;
  int dispunit = sizeof(int);
  int src, dst, errors;

  MPI_Init( &argc, &argv );

  winsize = sizeof(wbuf) / dispunit;
  MPI_Win_create(&wbuf[0], winsize, dispunit, MPI_INFO_NULL, MPI_COMM_WORLD, &window);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank == 0){
    wbuf[0] = 1;
    for(i=1 ; i < sizeof(wbuf) / sizeof(int) ; i++) wbuf[i] = 9999;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    index = 1;
    for(i = 1 ; i< size ; i++){
      rank = wbuf[index];
      nw = wbuf[index+1];
      errors = 0;
      for(j=0 ; j< nw ; j++) if(wbuf[index+2+j] != ((rank << 16) + j)) errors++;
      printf("from rank %3d, n = %3d, at %4d, first = %8.8x, last = %8.8x, errors = %3d\n", rank, nw, index, wbuf[index+1+1], wbuf[index+1+nw], errors);
      index = index + nw + 2;
    }
    MPI_Win_free(&window);
  }else{
    buffer[0] = rank;
    buffer[1] = rank;
    for(i=0 ; i<rank ; i++) buffer[i+2] = (rank << 16) + i;
    MPI_Barrier(MPI_COMM_WORLD);

    TargetDisp = 0;
    src = rank + 2;  // will use rank + 2 integers in remote buffer
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,window);
    MPI_Get(&dst, 1, MPI_INTEGER, 0, TargetDisp, 1, MPI_INTEGER, window);
    MPI_Accumulate(&src, 1, MPI_INTEGER,  0, TargetDisp, 1, MPI_INTEGER, MPI_SUM, window);
//     MPI_Get_accumulate(&src, 1, MPI_INTEGER, &dst, 1, MPI_INTEGER, 0, TargetDisp, 1, MPI_INTEGER, MPI_SUM, window);
//     MPI_Fetch_and_op(&src, &dst, MPI_INTEGER, 0, TargetDisp, MPI_SUM, window);
    MPI_Win_unlock(0,window);
    printf("rank %3d, adding %3d, storing at %5d\n", rank, src, dst);

    TargetDisp = dst;
    nw = rank + 2;
    MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
    MPI_Put(&buffer[0], nw, MPI_INTEGER, 0, TargetDisp, nw, MPI_INTEGER, window);
    MPI_Win_unlock(0,window);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_free(&window);
  }

  MPI_Finalize();
}

