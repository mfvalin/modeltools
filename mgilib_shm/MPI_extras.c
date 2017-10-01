/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2017  Division de Recherche en Prevision Numerique, Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

int MPI_Unpublish_name( char *service_name, MPI_Info info,  char *port_name)
{
  char filename[4096];
  char *home = getenv("HOME");
  char *gossipdir = getenv("GOSSIP_MPI_DIR");

  if(gossipdir == NULL) {
    snprintf(filename,4096,"%s/%s/%s",getenv("HOME"),".gossip/MPI",service_name);
  }else{
    snprintf(filename,4096,"%s/%s",gossipdir,service_name);
  }
  unlink(filename);
  printf("INFO: unpublished '%s' as '%s'\n",service_name,filename);
  return MPI_SUCCESS;
}

int MPI_Publish_name( char *service_name, MPI_Info info,  char *port_name)
{
  FILE *gossip;
  char filename[4096];
  char filenew[4096];
  char *home = getenv("HOME");
  char *gossipdir = getenv("GOSSIP_MPI_DIR");

  if(gossipdir == NULL) {
    snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip");
    mkdir(filename,0755);
    snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip/MPI");
    mkdir(filename,0755);
    snprintf(filename,4096,"%s/%s/%s.new",getenv("HOME"),".gossip/MPI",service_name);
    snprintf(filenew,4096,"%s/%s/%s"     ,getenv("HOME"),".gossip/MPI",service_name);
  }else{
    snprintf(filename,4096,"%s/%s.new",gossipdir,service_name);
    snprintf(filenew,4096,"%s/%s"     ,gossipdir,service_name);
  }
  unlink(filename);
  unlink(filenew);

  gossip = fopen(filename,"w");
  fprintf(gossip,"%s",port_name);
  fclose(gossip);

  link(filename,filenew);
  unlink(filename);
  return MPI_SUCCESS;
}

int MPI_Lookup_name( char *service_name, MPI_Info info, char *port_name)
{
  char filename[4096];
  int fd, nc;
  int wait=0;
  char *home = getenv("HOME");
  char *gossipdir = getenv("GOSSIP_MPI_DIR");

  if(gossipdir == NULL) {
    snprintf(filename,4096,"%s/%s/%s",getenv("HOME"),".gossip/MPI",service_name);
  }else{
    snprintf(filename,4096,"%s/%s",gossipdir,service_name);
  }

  while( (fd=open(filename,0)) < 0) { wait++ ; usleep(1000); }
  nc=read(fd,port_name,MPI_MAX_PORT_NAME);
  close(fd);
  port_name[nc]='\0';
  printf("MPI_Lookup_name: wait time = %d msec\n",wait);

  return MPI_SUCCESS;
}

int MPI_Create_named_port(char *publish_name)
{
  char port_name[MPI_MAX_PORT_NAME]; 
  MPI_Open_port(MPI_INFO_NULL, port_name); 
  printf("INFO: port available at %s\n",port_name); 
  return(MPI_Publish_name(publish_name, MPI_INFO_NULL, &port_name[0]));
}

MPI_Comm MPI_Connect_to_named_port(char *publish_name)  // connect to published port and verify connection
{
  char port_name[MPI_MAX_PORT_NAME];
  MPI_Comm server, local;
  int handshake = 1;  // client sends 1, expects to receive 0
  int answer = -1;
  int localrank, localsize;
  int tag = 123456;
  MPI_Status status; 

  MPI_Lookup_name(publish_name, MPI_INFO_NULL, &port_name[0]);
  printf("connecting to server at '%s'\n",port_name); 
  MPI_Comm_connect( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &server ); 
  MPI_Intercomm_merge(server, 1, &local);

  MPI_Comm_rank(local, &localrank);
  MPI_Comm_size(local, &localsize);
  if(localsize != 2) return MPI_COMM_NULL;

  MPI_Barrier(local);
  MPI_Sendrecv(&handshake, 1, MPI_INTEGER, 1-localrank, tag, &answer, 1, MPI_INTEGER, 1-localrank, tag, local, &status);
  return (answer != 0) ? MPI_COMM_NULL : local ;  // error if not corrrect answer from server
}

MPI_Comm MPI_Accept_on_named_port(char *publish_name)  // accept on published port and verify connection
{
  char port_name[MPI_MAX_PORT_NAME];
  MPI_Comm client, local; 
  int handshake = 0;  // server sends 0, expects to receive 1
  int answer = -1;
  int localrank, localsize;
  int tag = 123456;
  MPI_Status status; 

  MPI_Lookup_name(publish_name, MPI_INFO_NULL, &port_name[0]);
  printf("server at '%s' accepting connection\n",port_name); 
  MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client );
  MPI_Intercomm_merge(client, 0, &local);

  MPI_Comm_rank(local, &localrank);
  MPI_Comm_size(local, &localsize);
  if(localsize != 2) return MPI_COMM_NULL;

  MPI_Barrier(local);
  MPI_Sendrecv(&handshake, 1, MPI_INTEGER, 1-localrank, tag, &answer, 1, MPI_INTEGER, 1-localrank, tag, local, &status);
  return (answer != 1) ? MPI_COMM_NULL : local ;  // error if not corrrect answer from client
}
