#include "mpi.h" 
#define MAX_DATA 128
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

int main( int argc, char **argv ) 
{ 
  MPI_Comm client1, client2, local1, local2; 
  MPI_Win window1, window2;
  MPI_Aint winsize;
  int dispunit = sizeof(int);
//   void *memptr;
  MPI_Status status; 
  char port_name[MPI_MAX_PORT_NAME]; 
  double buf[MAX_DATA]; 
  int    size, rank; 
  char *servicename="demo_mgi";

  MPI_Init( &argc, &argv ); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  if (size != 1) error(FATAL, "Server too big");

  winsize = 1024*1024*dispunit + sizeof(mpi_channel);
  MPI_Alloc_mem(winsize,MPI_INFO_NULL,&memptr);
  memptr->control = 0;
  memptr->first = 1;
  memptr->in = 1;
  memptr->out = 1;
  memptr->limit = 1024*1024;

  MPI_Open_port(MPI_INFO_NULL, port_name); 
  printf("INFO: server available at %s\n",port_name); 
  MPI_Publish_name( servicename, MPI_INFO_NULL, &port_name[0]);

  MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client1 );
  MPI_Intercomm_merge(client1, 0, &local1);
  MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local1, &window1);

  MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client2 ); 
  MPI_Intercomm_merge(client2, 0, &local2);
  MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local2, &window2);

  printf("INFO: all clients connected to server\n"); 
  MPI_Unpublish_name( servicename, MPI_INFO_NULL, &port_name[0]);

            MPI_Recv( buf, MAX_DATA, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, client1, &status ); 
            printf("something done\n");
            MPI_Recv( buf, MAX_DATA, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, client2, &status ); 
            printf("something done\n");
    printf("before barrier(client1)\n");
    MPI_Barrier(client1);
    MPI_Win_free(&window1);
    MPI_Comm_disconnect( &client1 ); 
//     MPI_Comm_free( &client1 );
    printf("before barrier(client2)\n");
    MPI_Barrier(client2);
    MPI_Win_free(&window2);
    MPI_Comm_disconnect( &client2 ); 
//     MPI_Comm_free( &client2 );
    MPI_Close_port(port_name);
//     MPI_Free_mem(memptr);
    printf("INFO: server shutting down\n");
    MPI_Finalize(); 
} 


