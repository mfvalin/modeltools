#include <mpi.h> 
#include <stdio.h>
#include <stdlib.h>
// client for separate demo of server + client
// expects 1 process
// mpirun -n 1 server
// aprun  -n 1 server
// the port string will be printed and will be given to the client
int main( int argc, char **argv ) 
{ 
  MPI_Comm server, client1, client2, local, local1, local2;
  int size, rank;
  char port_name[MPI_MAX_PORT_NAME];

  MPI_Init( &argc, &argv ); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

  if (size != 1) {
    printf("ERROR: MPI WORLD size must be 1 for this test\n");
    MPI_Finalize();  
    exit(0);
  }


    MPI_Open_port(MPI_INFO_NULL, port_name);
    printf("server: port name = '%s'\n",port_name);

    MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client1 );  // accept connection from client1
    MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client2 );  // accept connection from client2

    printf("before barrier(client1)\n");
    MPI_Barrier(client1);
    MPI_Comm_disconnect( &client1 );

    printf("before barrier(client2)\n");
    MPI_Barrier(client2);
    MPI_Comm_disconnect( &client2 );
 
    MPI_Close_port(port_name);
    printf("INFO: server shutting down\n");

  MPI_Finalize(); 
} 

