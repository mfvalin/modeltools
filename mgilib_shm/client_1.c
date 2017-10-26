#include <mpi.h> 
#include <stdio.h>
#include <stdlib.h>
// client for separate demo of server + client
// expects 2 processes
// mpirun -n 2 client 'port string frtom previously started server'
// aprun  -n 2 client 'port string frtom previously started server'
int main( int argc, char **argv ) 
{ 
  MPI_Comm server, client1, client2, local, local1, local2;
  int size, rank;
  char *port_name;

  MPI_Init( &argc, &argv ); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  port_name = argv[1];

  if (size != 2) {
    printf("ERROR: MPI WORLD size must be 2 for this test\n");
    MPI_Finalize();  
    exit(0);
  }

    printf("client %d: port name = '%s'\n",rank,port_name);
    sleep(rank+1);  // necessary delay before connect to ensure that accept has been posted before connect
    MPI_Comm_connect( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &server );

    printf("before barrier(server)\n");
    MPI_Barrier(server);
    MPI_Comm_disconnect( &server ); 
    printf("INFO: client shutting down\n");
  MPI_Finalize(); 
} 

