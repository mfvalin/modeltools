#include <mpi.h> 
#include <stdio.h>
#include <stdlib.h>
// combined demo server+client
// expects 3 processes
// mpirun -n 3 .....
// aprun  -n 3 .....
int main( int argc, char **argv ) 
{ 
  MPI_Comm server, client1, client2, local, local1, local2;
  int size, rank;
  char port_name[MPI_MAX_PORT_NAME];

  MPI_Init( &argc, &argv ); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

  if (size != 3) {
    printf("ERROR: MPI WORLD size must be 3 for this test\n");
    MPI_Finalize();  
    exit(0);
  }

  if(rank == 0) {   // server

    MPI_Open_port(MPI_INFO_NULL, port_name);
    MPI_Bcast(port_name,MPI_MAX_PORT_NAME,MPI_CHARACTER,0,MPI_COMM_WORLD);
    printf("server: port name = '%s'\n",port_name);

    MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client1 );
    MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client2 );

    printf("before barrier(client1)\n");
    MPI_Barrier(client1);
    MPI_Comm_disconnect( &client1 );

    printf("before barrier(client2)\n");
    MPI_Barrier(client2);
    MPI_Comm_disconnect( &client2 );
 
    MPI_Close_port(port_name);
    printf("INFO: server shutting down\n");

  }

  if(rank > 0) {  // clients 1 and 2

    MPI_Bcast(port_name,MPI_MAX_PORT_NAME,MPI_CHARACTER,0,MPI_COMM_WORLD);
    printf("client %d: port name = '%s'\n",rank,port_name);
    sleep(rank);  // necessary delay before connect to ensure that accept has been posted before connect
    MPI_Comm_connect( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &server );

    printf("before barrier(server)\n");
    MPI_Barrier(server);
    MPI_Comm_disconnect( &server ); 
    printf("INFO: client shutting down\n");
  }
  MPI_Finalize(); 
} 

