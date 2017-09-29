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

  snprintf(filename,4096,"%s/%s/%s",getenv("HOME"),".gossip/MPI",service_name);
  unlink(filename);
  printf("INFO: unpublished '%s' as '%s'\n",service_name,filename);
  return MPI_SUCCESS;
}

int MPI_Publish_name( char *service_name, MPI_Info info,  char *port_name)
{
  FILE *gossip;
  char filename[4096];
  char filenew[4096];

  snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip");
  mkdir(filename,0755);
  snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip/MPI");
  mkdir(filename,0755);
  snprintf(filename,4096,"%s/%s/%s.new",getenv("HOME"),".gossip/MPI",service_name);
  unlink(filename);
  snprintf(filenew,4096,"%s/%s/%s",getenv("HOME"),".gossip/MPI",service_name);
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

  snprintf(filename,4096,"%s/%s/%s",getenv("HOME"),".gossip/MPI",service_name);

  while( (fd=open(filename,0)) < 0) { wait++ ; usleep(1000); }
  nc=read(fd,port_name,MPI_MAX_PORT_NAME);
  close(fd);
  port_name[nc]='\0';
  printf("MPI_Lookup_name: wait time = %d msec\n",wait);

  return MPI_SUCCESS;
}

