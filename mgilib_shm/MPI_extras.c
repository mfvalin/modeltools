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
#define MPI_ERROR MPI_ERR_LASTCODE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#define MAX_NAMES 128
static int last_name=-1;
static char *public_names[MAX_NAMES];
static char *port_names[MAX_NAMES];
static char is_published[MAX_NAMES];

static int name_lookup(char *public_name){
  int i;

  for (i=0 ; i<=last_name ; i++) {
    if(public_names[i] == NULL) continue;
    if(strcmp(public_name,public_names[i]) == 0) {
      printf("DEBUG: %s found in slot %d\n",public_name,i);
      return(i);
    }
  }
  return -1;
}

static int name_insert(char *public_name, char *port_name){
  size_t maxlen = MPI_MAX_PORT_NAME;
  size_t stlen;
  char *str1;
  char *str2;

  if(last_name+1 >= MAX_NAMES) return -1;
// printf("public_name(%ld), port_name(%ld)\n",strnlen(public_name,maxlen),strnlen(port_name,maxlen));
  str1 = malloc(10+strnlen(public_name,maxlen));
  str2 = malloc(10+strnlen(port_name  ,maxlen));
  if(str1 == NULL || str2 == NULL) return -1;

  public_names[last_name+1] = str1;
  port_names[last_name+1]   = str2;
  is_published[last_name+1] = 1;
  while(*public_name) { *str1++ = *public_name++ ; } ; *str1 = '\0';
  while(*port_name)   { *str2++ = *port_name++   ; } ; *str2 = '\0';
// return -1;
//   strncpy(str1,public_name,maxlen);
//   strncpy(str2,port_name  ,maxlen);
// printf("DEBUG: in name_insert\n");
  last_name++;
  return last_name;
}

int MPI_Unpublish_name(char *service_name, MPI_Info info, char *port_name)
{
  char filename[4096];
  char *home = getenv("HOME");
  char *gossipdir = getenv("GOSSIP_MPI_DIR");
  int target = name_lookup(service_name);

  if(target != -1) {
    if(is_published[target] == 0) {
      printf("WARNING: '%s' already unpublished\n",service_name);
      return MPI_SUCCESS;   // already unpublished
    }
    is_published[target] = 0;
  }

  if(gossipdir == NULL) {
    snprintf(filename,4096,"%s/%s/%s",getenv("HOME"),".gossip/MPI",service_name);
  }else{
    snprintf(filename,4096,"%s/%s",gossipdir,service_name);
  }
  unlink(filename);
  printf("INFO: unpublished '%s' as '%s'\n",service_name,filename);
  return MPI_SUCCESS;
}

int MPI_Unpublish_named_port( char *service_name)
{
  return MPI_Unpublish_name(service_name,MPI_INFO_NULL,"Not Used");
}

int MPI_Publish_name( char *service_name, MPI_Info info,  char *port_name)
{
  FILE *gossip;
  char filename[4096];
  char filenew[4096];
  char *home = getenv("HOME");
  char *gossipdir = getenv("GOSSIP_MPI_DIR");
  int target = name_lookup(service_name);

  if(target != -1) return MPI_SUCCESS; // already published

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

  target = name_insert(service_name, port_name);
printf("DEBUG: after name_insert\n");
  if(target == -1) {
    printf("WARNING: cache table for named ports is full\n");
  }
  printf("INFO: published '%s' as '%s'\n",service_name,filenew);
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

int MPI_Close_named_port(char *publish_name)
{
  char port_name[MPI_MAX_PORT_NAME];
  MPI_Lookup_name(publish_name, MPI_INFO_NULL, &port_name[0]);
  MPI_Close_port(port_name);
  printf("INFO: named port '%s' closed\n",publish_name);
}

int MPI_Create_named_port(char *publish_name)
{
  char port_name[MPI_MAX_PORT_NAME]; 
  MPI_Open_port(MPI_INFO_NULL, port_name); 
  printf("INFO: named port '%s' available at %s\n",publish_name,port_name); 
  return(MPI_Publish_name(publish_name, MPI_INFO_NULL, &port_name[0]));
}

int MPI_Connect_to_named_port(char *publish_name, MPI_Comm *server, MPI_Comm *local)  // connect to published port and verify connection
{
  char port_name[MPI_MAX_PORT_NAME];
  int handshake = 1;  // client sends 1, expects to receive 0
  int answer = -1;
  int localrank, localsize;
  int tag = 123456;
  MPI_Status status; 

  MPI_Lookup_name(publish_name, MPI_INFO_NULL, &port_name[0]);
  printf("INFO: connecting to server at '%s'\n",port_name); 
  MPI_Comm_connect( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, server ); 
  printf("INFO: connected to server at '%s'\n",port_name); 
  MPI_Intercomm_merge(*server, 1, local);

  MPI_Comm_rank(*local, &localrank);
  MPI_Comm_size(*local, &localsize);
  if(localsize != 2) {
    printf("ERROR: size of port communicator must be 2\n");
    return MPI_ERROR;
  }

  MPI_Barrier(*local);
  MPI_Sendrecv(&handshake, 1, MPI_INTEGER, 1-localrank, tag, &answer, 1, MPI_INTEGER, 1-localrank, tag, *local, &status);
  if(answer != 0) {
    printf("ERROR: bad handshake answer, expected 0, got %d\n",answer);
    return MPI_ERROR;
  }
  return MPI_SUCCESS ;
}

int MPI_Accept_on_named_port(char *publish_name, MPI_Comm *client, MPI_Comm *local)  // accept on published port and verify connection
{
  char port_name[MPI_MAX_PORT_NAME];
  int handshake = 0;  // server sends 0, expects to receive 1
  int answer = -1;
  int localrank, localsize;
  int tag = 123456;
  MPI_Status status; 

  MPI_Lookup_name(publish_name, MPI_INFO_NULL, &port_name[0]);
  printf("INFO: server at '%s' accepting connection\n",port_name); 
  MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, client );
  printf("INFO: server at '%s' accepted connection\n",port_name); 
  MPI_Intercomm_merge(*client, 0, local);

  MPI_Comm_rank(*local, &localrank);
  MPI_Comm_size(*local, &localsize);
  if(localsize != 2) {
    printf("ERROR: size of port communicator must be 2\n");
    return MPI_ERROR;
  }

  MPI_Barrier(*local);
  MPI_Sendrecv(&handshake, 1, MPI_INTEGER, 1-localrank, tag, &answer, 1, MPI_INTEGER, 1-localrank, tag, *local, &status);
  if(answer != 1) {
    printf("ERROR: bad handshake answer, expected 1, got %d\n",answer);
    return MPI_ERROR;
  }
  printf("INFO: handshake successful\n");
  return MPI_SUCCESS ;
}