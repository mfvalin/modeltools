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
#define MAX_DATA 128
#define MAXSIZE 512*1024
#define FATAL 1

typedef struct{
  int control;    // flags
  int first;      // start of buffer
  int in;         // insertion index
  int out;        // extraction index
  int limit;      // end of buffer + 1
  int data[1];    // place holder, start of data buffer
} mpi_channel;

static mpi_channel *memptr;

int MPI_Create_named_port(char *publish_name);
int MPI_Unpublish_named_port( char *service_name);
int MPI_Close_named_port(char *publish_name);
int MPI_Connect_to_named_port(char *publish_name, MPI_Comm *server, MPI_Comm *local);
int MPI_Accept_on_named_port(char *publish_name, MPI_Comm *client, MPI_Comm *local);

int main( int argc, char **argv ) 
{ 
  MPI_Comm server, client1, client2, local, local1, local2; 
  MPI_Win window, window1, window2;
  MPI_Aint winsize;
  int dispunit = sizeof(int);
  MPI_Aint TargetDisp;
  double t1, t2;
//   void *memptr;
  MPI_Status status; 
  char port_name[MPI_MAX_PORT_NAME]; 
  char port_name_plus[MPI_MAX_PORT_NAME+128]; 
  double buf[MAX_DATA]; 
  int data[MAXSIZE+10];
  int    size, rank, localrank; 
  char *servicename="demo_mgi";
  int i, tag;
  int n=1;
  int lockloop;
  int world_size = 0;

  MPI_Init( &argc, &argv ); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
#if defined(SERVER)
  world_size++;
#endif
#if defined(CLIENT)
  world_size += 2;
#endif
  if (size != world_size) {
    printf("ERROR: MPI WORLD size must be %d for this test\n",world_size);
    MPI_Finalize();  
    exit(0);
  }else{
    printf("INFO: WORLD size is %d as expected\n",world_size);
  }
#if defined(SERVER)
  if(rank == 0) {   // server
    winsize = 1024*1024*dispunit + sizeof(mpi_channel);
    MPI_Alloc_mem(winsize,MPI_INFO_NULL,&memptr);
    memptr->control = 0;
    memptr->first = 1;
    memptr->in = 1;
    memptr->out = 1;
    memptr->limit = 1024*1024;

    MPI_Create_named_port(servicename);

    MPI_Accept_on_named_port(servicename, &client1, &local1);
    MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local1, &window1);
    MPI_Comm_rank(local1, &localrank);
    printf("server rank in local1 = %d\n",localrank);

    MPI_Accept_on_named_port(servicename, &client2, &local2);
    MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local2, &window2);
    MPI_Comm_rank(local2, &localrank);
    printf("server rank in local2 = %d\n",localrank);

    printf("INFO: all clients are now connected to server\n"); 

    for(i=0 ; i<MAXSIZE ; i++) memptr->data[i] = -1;
    // receive a single two sided message from clients 1 and 2
    MPI_Recv( buf, MAX_DATA, MPI_DOUBLE, 1, MPI_ANY_TAG, local1, &status ); 
    printf("message received from client 1\n");
    MPI_Recv( buf, MAX_DATA, MPI_DOUBLE, 1, MPI_ANY_TAG, local2, &status ); 
    printf("message received from client 2\n");

    printf("before barrier(local1)\n");
    MPI_Barrier(local1);
    printf("before barrier(local2)\n");
    MPI_Barrier(local2);
    MPI_Barrier(MPI_COMM_WORLD);
// one sided action starts here
#if defined(LOCKLOOP)
    lockloop = 100000 ; // 100 ms
    while(lockloop-- >0){
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window1);
      MPI_Win_unlock(0,window1);
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window2);
      MPI_Win_unlock(0,window2);
//      usleep(1);
    } 
#endif
    printf("before barrier(client1)\n");
    printf("server data[0] = %d, data[MAXSIZE-1] = %d\n",memptr->data[0],memptr->data[MAXSIZE-1]);

// getting ready to terminate
    MPI_Barrier(client1);
    MPI_Comm_disconnect( &client1 ); 
    MPI_Win_free(&window1);
//     MPI_Comm_free( &client1 );

    printf("before barrier(client2)\n");
    printf("server data[0] = %d, data[MAXSIZE-1] = %d\n",memptr->data[0],memptr->data[MAXSIZE-1]);
    MPI_Barrier(client2);

    printf("server data[0] = %d, data[MAXSIZE-1] = %d\n",memptr->data[0],memptr->data[MAXSIZE-1]);

    MPI_Comm_disconnect( &client2 ); 
    MPI_Win_free(&window2);
//     MPI_Comm_free( &client2 );
    MPI_Close_named_port(servicename);
    MPI_Unpublish_named_port(servicename);
    MPI_Unpublish_named_port(servicename);
    MPI_Free_mem(memptr);
    printf("INFO: server shutting down\n");
  }
#endif
#if defined(CLIENT)
  if(rank > world_size-3) {  // clients 1 and 2
    MPI_Connect_to_named_port(servicename, &server, &local);
    winsize = 1024*1024*dispunit;
    MPI_Alloc_mem(winsize,MPI_INFO_NULL,&memptr);
    MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local, &window);
    MPI_Comm_rank(local, &localrank);
    printf("client rank in local = %d\n",localrank);
 
    tag = 2; /* Action to perform */ 
//         MPI_Send( buf, n, MPI_DOUBLE, 0, tag, server ); 
    MPI_Send( buf, n, MPI_DOUBLE, 0, tag, local ); 
    printf("message sent to server\n");

    printf("before barrier(local)\n");
    MPI_Barrier(local);

    if(rank == world_size-2){   // one sided put to server, send receive with client 2
      for(i=0 ; i<MAXSIZE ; i++) data[i] = i;
      MPI_Barrier(MPI_COMM_WORLD);
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 5;
      MPI_Put(data, MAXSIZE, MPI_INTEGER, 0, TargetDisp, MAXSIZE, MPI_INTEGER, window);  // put to server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after put, time = %G\n",1000*(t2-t1));
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 0;
      MPI_Put(data, 100000, MPI_INTEGER, 0, TargetDisp, 100000, MPI_INTEGER, window);  // put to server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after put, time = %G\n",1000*(t2-t1));
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 0;
      MPI_Put(data, 10000, MPI_INTEGER, 0, TargetDisp, 10000, MPI_INTEGER, window);  // put to server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after put, time = %G\n",1000*(t2-t1));
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 0;
      MPI_Put(data, 1000, MPI_INTEGER, 0, TargetDisp, 1000, MPI_INTEGER, window);  // put to server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after put, time = %G\n",1000*(t2-t1));
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 0;
      MPI_Put(data, 1, MPI_INTEGER, 0, TargetDisp, 1, MPI_INTEGER, window);  // put to server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after put, time = %G\n",1000*(t2-t1));
      MPI_Sendrecv(buf, 1, MPI_INTEGER,world_size-1-rank, tag, buf+1, 1, MPI_INTEGER, world_size-1-rank, tag, MPI_COMM_WORLD, &status);  // tell partner it is done
    }
    if(rank == world_size-1){   // send receive with client 1, one sided get from server, check that data has been received
      for(i=0 ; i<MAXSIZE ; i++) data[i] = -1;
      MPI_Barrier(MPI_COMM_WORLD);
      printf("before sendrecv\n");
      MPI_Sendrecv(buf, 1, MPI_INTEGER,world_size-1-rank, tag, buf+1, 1, MPI_INTEGER, world_size-1-rank, tag, MPI_COMM_WORLD, &status);  // wait until partner is done
      printf("after sendrecv\n");
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 5;
      MPI_Get(data, MAXSIZE, MPI_INTEGER, 0, TargetDisp, MAXSIZE, MPI_INTEGER, window);  // get from server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after get time = %G\n",1000*(t2-t1));
      printf("read data[0] = %d, data[MAXSIZE-1] = %d\n",data[0],data[MAXSIZE-1]);
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 0;
      MPI_Get(data, 100000, MPI_INTEGER, 0, TargetDisp, 100000, MPI_INTEGER, window);  // get from server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after get time = %G\n",1000*(t2-t1));
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 0;
      MPI_Get(data, 10000, MPI_INTEGER, 0, TargetDisp, 10000, MPI_INTEGER, window);  // get from server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after get time = %G\n",1000*(t2-t1));
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 0;
      MPI_Get(data, 1000, MPI_INTEGER, 0, TargetDisp, 1000, MPI_INTEGER, window);  // get from server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after get time = %G\n",1000*(t2-t1));
      t1 = MPI_Wtime();
      MPI_Win_lock(MPI_LOCK_SHARED,0,0,window);
      TargetDisp = 0;
      MPI_Get(data, 1, MPI_INTEGER, 0, TargetDisp, 1, MPI_INTEGER, window);  // get from server memory
      MPI_Win_unlock(0,window);
      t2 = MPI_Wtime();
      printf("after get time = %G\n",1000*(t2-t1));
    }
    printf("before barrier(server)\n");
    MPI_Barrier(server);
    MPI_Comm_disconnect( &server ); 
    MPI_Win_free(&window);
    MPI_Free_mem(memptr);
//     MPI_Comm_free(&server);
    printf("INFO: client shutting down\n");
  }
#endif
  MPI_Finalize(); 
} 


