/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2017  Division de Recherche en Prevision Numerique, Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
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
#if defined(ALWAYS_COMPILE_DOCUMENTATION)
  the following environment variables are used in this package 

  GOSSIP_MPI_TIMEOUT    : timeout in seconds when accepting/requesting connection
  GOSSIP_MPI_DIR        : base directory for GOSSIP files
  HOME
  MPI_WORLD_NAME        : "name" assigned to this MPI world (normally set by r.run_in_parallel or equivalent)

  extra MPI helper routines dealing with connections

#endif

#include <mpi.h>

#include "mgi.h"
#define MPI_ERROR MPI_ERR_LASTCODE

#define MAX_PORT_EXTRA 128
#define MPI_MAX_PORT_NAME_PLUS MPI_MAX_PORT_NAME+MAX_PORT_EXTRA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/shm.h>

#if defined(DEBUG)
#define DEBUG_(a)    a
#else
#define DEBUG_(a)    ;
#endif

#define MAX_NAMES 128
// typedef struct{
//   char *name;     // published name
//   char *wname;    // MPI world name "UnKnOwN" by default
//   char *port;     // MPI port name
//   void *addr;     // address in shared memory
//   int  hostid;    // host id
//   int  memid;     // shared memory id
//   int  rank;      // rank in MPI world
//   int  pub;       // published flag
// } name_entry;
// static name_entry name_table[MAX_NAMES];  // name and attributes table

static int last_name=-1;                  // index of last published name

static char *public_names[MAX_NAMES];     // published name
static char *port_names[MAX_NAMES];       // MPI port name
static void *memaddr[MAX_NAMES];
static int  shmem_id[MAX_NAMES];
static char is_published[MAX_NAMES];

static char *world_name = NULL;        // name associated with this MPI world (normally set by r.run_in_parallel or equivalent)
static char *gossip_dir = NULL;        // base gossip directory (NULL means use default value)
static char *home_dir = NULL;          // user's home directory
static long host_id = 0;               // machine host id
static int rank_in_world = -1;         // rank in MPI_COMM_WORLD
static int must_init = 1;              // if nonzero, package is not initialized and init_vars() should be called
static int maxwait = 100;              //  accept connection timeout 100 seconds by default

static void init_vars(){  // initialize some basic strings and variables
  char *tmp;
  tmp = getenv("GOSSIP_MPI_TIMEOUT");             // timeout in seconds for accept connection
  if(tmp) maxwait = atoi(tmp);
  gossip_dir = getenv("GOSSIP_MPI_DIR");
  home_dir   = getenv("HOME");                    // this environment variable should always exist
  world_name = getenv("MPI_WORLD_NAME");
  if(world_name == NULL) world_name = "?";        // if environment variable not set
  host_id = gethostid();
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_in_world);
  must_init = 0;
}

// build_gossip_name : internal service routine used to build file names for GOSSIP MPI routines
//
// char buf[bufmax] will contain the relevant file name
// prefix is only used if environment variable GOSSIP_MPI_DIR is not defined (and should start with / if not blank)
// name : basename for file name
// ext  : file name extension (should start with .)
//
static int build_gossip_name(char *buf, int bufmax, const char *prefix, const char *name, const char *ext, int mkdirs)
{
  if(must_init) init_vars() ;

  if(gossip_dir == NULL) {   // use default directory names, environment variable GOSSIP_MPI_DIR was not defined
    if(mkdirs){              // create relevant directories in case they do not already exist
      snprintf(buf,bufmax,"%s/%s",home_dir,".gossip");  // $HOME/.gossip
      mkdir(buf,0755);
      snprintf(buf,bufmax,"%s/%s%s",home_dir,".gossip",prefix);  // prefix should start with / if not ""
      mkdir(buf,0755);
    }
    snprintf(buf,bufmax,"%s/%s/%s%s%s",home_dir,".gossip",prefix,name,ext ? ext : "");  // ext should start with .
  }else{                     // use value from environment variable GOSSIP_MPI_DIR, directory MUST exist
    snprintf(buf,bufmax,"%s/%s%s",gossip_dir,name,ext ? ext : "");
  }
}

// name_lookup : internal service routine to find a name in the internal name cache table
static int name_lookup(const char *public_name){
  int i;

  for (i=0 ; i<=last_name ; i++) {
    if(public_names[i] == NULL) continue;                           // blank or suppressed entry
    if(strcmp(public_name,public_names[i]) == 0) {
      DEBUG_( printf("DEBUG: '%s' found in slot %d\n",public_name,i) ) ;
      return(i);  // return index to published name
    }
  }
  return -1;   // not found
}

// name_insert :  insert a new (name value) pair into the internal name cache table
// port_name can be a multi line collection of items
static int name_insert(const char *public_name,const char *port_name){
  size_t maxlen = MPI_MAX_PORT_NAME + 128;  // add room for host id, host name, MPI world name, memory id, and MPI rank
  char *str1;
  char *str2;

  if(last_name+1 >= MAX_NAMES) {
    printf("ERROR (name_insert): name cache table is full (max size = %d)\n",MAX_NAMES);
    return -1;
  }

  if(name_lookup(public_name) >= 0) {
    printf("ERROR (name_insert): '%s' already defined\n",public_name);
    return -1; // already in tables
  }

  str1 = malloc(10+strnlen(public_name,maxlen));   // allocate space for names
  if(str1 == NULL) return -1;
  str2 = malloc(10+strnlen(port_name  ,maxlen));
  if(str2 == NULL) {
    free(str1);
    return -1;
  }

  last_name++;    // bump table size index

  public_names[last_name] = str1;
  while(*public_name) { *str1++ = *public_name++ ; } ; *str1 = '\0';  // copy published name

  port_names[last_name] = str2;
  while(*port_name)   { *str2++ = *port_name++   ; } ; *str2 = '\0';  // copy composite 'port name'

  is_published[last_name] = 1;          // name is published
  shmem_id[last_name]     = -1;         // no shared memory tag is associated with name
  memaddr[last_name]      = NULL;       // no shared memory region is associated with name

  return last_name;  // return index to published name
}

// remove a name from the internal name cache table and "unpublish" it
// this function OVERLOADS the MPI library function MPI_Unpublish_name
// arguments info and port_name are there for compatibility with the MPI library function
// and are not used
int MPI_Unpublish_name(constchar *service_name, MPI_Info info, constchar *port_name)
{
  char filename[4096];

  int target = name_lookup(service_name);   // is service_name already published

  if(target != -1) {                        // name exists
    if(is_published[target] == 0) {         // and already unpublished
      printf("WARNING: '%s' already unpublished\n",service_name);
      return MPI_SUCCESS;   // already unpublished
    }
    is_published[target] = 0;
  }

  build_gossip_name(filename, sizeof(filename), "/MPI", service_name, NULL, 1);   // file used to publish 'service_name'
  unlink(filename);                                                               // remove file
  printf("INFO: unpublished '%s' as '%s'\n",service_name,filename);

  return MPI_SUCCESS;
}

// remove a name from the internal name cache table
// calls MPI_Unpublish_name with appropriate extra arguments
int MPI_Unpublish_named_port(const char *service_name)
{
  return MPI_Unpublish_name((constchar *) service_name,MPI_INFO_NULL,"Not Used");
}

// add a name from the internal name cache table and "publish" it
// this function OVERLOADS the MPI library function MPI_Publish_name
// argument info is there for compatibility with the MPI library function and is not used
int MPI_Publish_name(constchar *service_name, MPI_Info info, constchar *port_name)
{
  FILE *gossip;
  char filename[4096];
  char filenew[4096];

  int target = name_lookup(service_name);

  if(target != -1) return MPI_SUCCESS; // already published

  build_gossip_name(filename, sizeof(filename), "/MPI", service_name, ".new", 1);
  build_gossip_name(filenew , sizeof(filenew ), "/MPI", service_name,   NULL, 0);   // file used to publish 'service_name'
  unlink(filename);   // in case there is an old file with this name
  unlink(filenew);

  gossip = fopen(filename,"w");       // write port_name into file used to publish (temporary name)
  fprintf(gossip,"%s",port_name);
  fclose(gossip);

  link(filename,filenew);             // atomic creation of file used to publish
  unlink(filename);

  target = name_insert(service_name, port_name);    // insert into name cache table
// DEBUG_( printf("DEBUG: after name_insert\n") ) ;
  if(target == -1) {
    printf("WARNING: cache table for named ports is full\n");
  }
  printf("INFO: published '%s' as '%s'\n",service_name,filenew);
  return MPI_SUCCESS;
}

// add a name into the internal name cache table and "publish" it
// calls MPI_Publish_name with appropriate extra argument
int MPI_Publish_named_port(const char *service_name, const char *port_name)
{
  return MPI_Publish_name((constchar *) service_name, MPI_INFO_NULL, (constchar *) port_name);
}

// get port_name string associated with name service_name
// this function OVERLOADS the MPI library function MPI_Lookup_name
// argument info is there for compatibility with the MPI library function and is not used
int MPI_Lookup_name(constchar *service_name, MPI_Info info, char *port_name)
{
  char filename[4096];
  int fd, nc;
  int wait=0;
  char *t;

  build_gossip_name(filename, sizeof(filename), "/MPI", service_name, NULL, 1);

  port_name[0] = '\0' ;
  while( (fd=open(filename,0)) < 0) { wait++ ; usleep(1000); }  // microsleep 1 millisecond
  nc=read(fd,port_name,MPI_MAX_PORT_NAME_PLUS-3);
  close(fd);
  port_name[nc+0]='\0'; // make sure that there are three null bytes at end of string
  port_name[nc+1]='\0';
  port_name[nc+2]='\0';
  t = port_name;
  while( *t != '\0' ) {
    if(*t == '\n') *t = '\0' ;  // replace newlines in string with nulls (the string really ends at @@)
    t++;
  }
  printf("MPI_Lookup_name: wait time = %d msec\n",wait);

  return MPI_SUCCESS;
}

// close MPI port published under name 'publish_name' if it exists
int MPI_Close_named_port(const char *publish_name)
{
  char port_name[MPI_MAX_PORT_NAME_PLUS];
  char *t = port_name;

  port_name[0] = '\0';
  MPI_Lookup_name((constchar *) publish_name, MPI_INFO_NULL, port_name);  // get MPI port name
  if(port_name[0] != '\0') MPI_Close_port(port_name);
  printf("INFO: named port '%s' closed\n",publish_name);
}

// create and publish a connection "port" (will be retrieved using appropriate lookup function)
// if no_mpi_port != 0, create MPI port
// shmid is shared memory tag of shared memory area
// publish_name is the name under which the "port" will be published
int MPI_Create_named_port(const char *publish_name, int shmid, int no_mpi_port)
{
  char port_name[MPI_MAX_PORT_NAME+1];
  char buf[MPI_MAX_PORT_NAME_PLUS];   // enough room to add hostid, shmid, MPI world name, and rank in MPI_COMM_WORLD
  int len;

  if(must_init) init_vars() ;

  snprintf(port_name,sizeof(port_name),"/dev/null");
  if(! no_mpi_port) MPI_Open_port(MPI_INFO_NULL, port_name);
  port_name[MPI_MAX_PORT_NAME] = '\0';
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_in_world);

  snprintf(buf, sizeof(buf),"%s\n%32s %6d\n %16d %16d\n", port_name, world_name, rank_in_world, host_id, shmid);
  buf[sizeof(buf)-1] = '\0';   // three null bytes guaranteed at end
  buf[sizeof(buf)-2] = '\0';
  buf[sizeof(buf)-3] = '\0';

  printf("INFO: named port '%s' available at %s\n",publish_name,buf); 
  return( MPI_Publish_name((constchar *) publish_name, MPI_INFO_NULL, buf) );
}

// split 3 line token returned by MPI_Lookup_name() into MPI port name, shared memory tag, rank in WORLD
// if not in same MPI world, rank = -1
// if not on same host, shared memory tag = -1
static int MPI_split_port_name(char *port_name, int *shmid, int *rank){  
  char worldname[MAX_PORT_EXTRA];
  char *t = port_name;
  int nitems;
  int hid;

  if(must_init) init_vars() ;

  *shmid = -1;
  *rank = -1;

  while(*t != '\0') t++; t++;  // position after first null character
  nitems = sscanf(port_name,"%s %d",worldname,rank);  // world name and rank in world
  if(worldname[0] == '?' || world_name[0] == '?' || strncmp(worldname,world_name,MAX_PORT_EXTRA) != 0 ) *rank = -1;  // not the same world name

  while(*t != '\0') t++; t++;  // position after second null character
  nitems = sscanf(port_name,"%d %d",&hid,shmid);      // hostid and shared memory tag
  if(hid != host_id) *shmid = -1;                     // not the right host
  
}

// post a request for connection to "publish_name"
// request includes MPI world name , machine host id, and rank in MPI world
//
// returns 0 if OK, -1 if error or timeout
static int MPI_Request_accept_on_named_port(const char *publish_name)
{
  char filereq[4096];
  char filenam[4096];
  char buffer[MAX_PORT_EXTRA];
  int wait = 0;
  int fd;

  if(must_init) init_vars() ;

  build_gossip_name(filereq, sizeof(filereq), "/MPI", publish_name, ".req", 0);  // request for connection (final)
  build_gossip_name(filenam, sizeof(filenam), "/MPI", publish_name, "", 0);      // request for connection (temporary)

  snprintf(buffer,sizeof(buffer),"%32s %d %d\n",world_name,host_id,rank_in_world); // name of MPI world, host id, rank in MPI world
  
  while( (fd = open(filenam,O_WRONLY+O_CREAT+O_EXCL,0700)) == -1) {  // try to create temporary file
    if(wait*10 >= maxwait) return -1; // timeout (default 100 seconds)
    wait++;
    usleep(100000);  // wait 100 msec
  }
  write(fd,buffer,strnlen(buffer,MAX_PORT_EXTRA));  // write request into temporary file
  close(fd);

  while((fd = open(filereq,O_RDONLY)) >= -1) {  // wait until filereq ceases to exist (previous request has been accepted by server)
    close(fd);
    if(wait*10 >= maxwait) {
      unlink(filenam);    // cancel tentative request
      return -1;          // timeout (default 100 seconds)
    }
    wait++;
    usleep(100000);       // wait 100 msec
  }
  
  link(filenam,filereq);  // atomic move of temporary file to request file
  unlink(filenam);        // remove temporary file

  wait = 0;
  while((fd = open(filereq,O_RDONLY)) >= -1) {  // wait until filereq ceases to exist (request has been accepted by server)
    close(fd);
    if(wait*10 >= maxwait) {
      unlink(filereq);    // cancel tentative request
      return -1;          // timeout (default 100 seconds)
    }
    wait++;
    usleep(100000);       // wait 100 msec
  }

  return(0);   // request has been posted and accepted, MPI_Connect_to_named_port should succeed
}

// returns
// -1 if error
//  0 if timeout
//  2 if link established
//    *arena == NULL , local and client not MPI_COMM_NULL  : MPI communicator link (one or two MPI worlds)
//    *arena != NULL , local and client both MPI_COMM_NULL : shared memory link, arena points to local address of shared memory area
//    *arena == NULL , local and client both MPI_COMM_NULL : we have a problem
// 
int MPI_Connect_to_named_port(const char *publish_name, MPI_Comm *server, MPI_Comm *local, void **arena)  // connect to published port and verify connection
{
  char port_name[MPI_MAX_PORT_NAME_PLUS];
  int handshake = 1;  // client sends 1, expects to receive 0
  int answer = -1;
  int localrank, localsize;
  MPI_Status status; 
  char filereq[4096];
  int shmid, rank;
  int rstatus;

  *arena  = NULL ;   // will be set to non NULL only if useful
  *server = MPI_COMM_NULL;
  *local  = MPI_COMM_NULL;

  if( MPI_Lookup_name((constchar *) publish_name, MPI_INFO_NULL, port_name) != MPI_SUCCESS ){        // get port name
    printf("ERROR: lookup of '%s' unsuccessful\n",publish_name);
    return -1;
  }

  MPI_split_port_name(port_name, &shmid, &rank);

  rstatus = MPI_Request_accept_on_named_port(publish_name);  // request for connection accept if no other choice
  if(rstatus != 0) {
    printf("WARNING: timeout while requesting connection to '%s'\n",publish_name);
    return 0;
  }

  if(shmid > 0){                                          // same host, shared memory segment is accessible
    *arena = shmat(shmid, NULL, 0);                       // attach to shared memory
    return (*arena == NULL) ? -1 : 2 ;   // no need for handshake if using shared memory, return error if shmat returned NULL
  }else{
    if(rank >= 0){                                        // same world
      MPI_Intercomm_create(MPI_COMM_SELF, 0, MPI_COMM_WORLD, rank, INTERCOMM_TAG, server);
    }else{                                                // not the same world and not common host, use remote MPI port
      printf("INFO: connecting to server at '%s'\n",port_name); 
      MPI_Comm_connect( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, server ); 
      printf("INFO: connected to server at '%s'\n",port_name); 
    }
  }

  MPI_Intercomm_merge(*server, 1, local);  // create intracommunicator between client and server

  MPI_Comm_rank(*local, &localrank);
  MPI_Comm_size(*local, &localsize);
  if(localsize != 2) {
    printf("ERROR: size of client-server communicator must be 2\n");
    return -1;
  }

  MPI_Barrier(*local);
  MPI_Sendrecv(&handshake, 1, MPI_INTEGER, 1-localrank, INTERCOMM_TAG, &answer, 1, MPI_INTEGER, 1-localrank, INTERCOMM_TAG, *local, &status);
  if(answer != 0) {
    printf("ERROR: bad handshake answer, expected 0, got %d\n",answer);
    return -1;
  }
  return 2 ;
}

// returns number of partners (normally 2, 1 client + 1 server)
// -1 if error
//  0 if timeout (no request for connection after 100 msec)
//  2 if link established
//    *arena == NULL , local and client not MPI_COMM_NULL  : MPI communicator link (one or two MPI worlds)
//    *arena != NULL , local and client both MPI_COMM_NULL : shared memory link
//    *arena == NULL , local and client both MPI_COMM_NULL : we have a problem
//  publish_name : name under which the connection has been published (in)
//  client       : intercommunicator (out)
//  local        : intracommunicator (out) (will be used to communicate with client)
//  arena        : address in memory of communication buffer if applicable (out)
//  timeout      : timeout in milliseconds (max 1000 seconds) (in)
int MPI_Accept_on_named_port(const char *publish_name, MPI_Comm *client, MPI_Comm *local, void **arena, int timeout)  // accept on published port and verify connection
{
  char port_name[MPI_MAX_PORT_NAME_PLUS];
  int handshake = 0;  // server sends 0, expects to receive 1
  int answer = -1;
  int localrank, localsize;
  MPI_Status status; 
  char filereq[4096];
  FILE *req;
  int wait = 0;
  int nc;
  int hid, shmid, rank;
  char worldname[MAX_PORT_EXTRA];

  *arena  = NULL ;   // will be set to non NULL only if useful
  *client = MPI_COMM_NULL;
  *local  = MPI_COMM_NULL;

  build_gossip_name(filereq, sizeof(filereq), "/MPI", publish_name, ".req", 0);  // request for connection
  req = fopen(filereq,"r");
  while(req == NULL){  // loop until filereq appears, then remove it (timeout may cause a premature return)
    if(wait > 0 || timeout == 0) return 0; // timeout
    wait++;
    req = fopen(filereq,"r");
    usleep(timeout*1000);  // wait 100 msec
  }
  nc = fscanf(req,"%32s%d%d%d",worldname,&rank,&hid,&shmid);  // get world name, rank, hostid, and shared memory tag
  fclose(req);
  unlink(filereq);    // remove request for connection

  if(hid == host_id && shmid != -1) {
    *arena = shmat(shmid, NULL, 0);      // use shared memory tag to connect if same host and shared memory tag is valid
    return (*arena == NULL) ? -1 : 2 ;   // no need for handshake if using shared memory, return error if shmat returned NULL
  }else{
    if(worldname[0] != '?' && world_name[0] != '?' || strncmp(worldname,world_name,32) == 0) {  // in same MPI world
      MPI_Intercomm_create(MPI_COMM_SELF, 0, MPI_COMM_WORLD, rank, INTERCOMM_TAG, client); // use rank in world of requester to create intercommunicator
    }else{   // not in same MPI world, use inter world MPI port
      MPI_Lookup_name((constchar *) publish_name, MPI_INFO_NULL, port_name);
      printf("INFO: server at '%s' accepting connection\n",port_name); 
      MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, client );
      printf("INFO: server at '%s' accepted connection\n",port_name); 
    }
  }

  MPI_Intercomm_merge(*client, 0, local);  // create intracommunicator between client and server

  MPI_Comm_rank(*local, &localrank);
  MPI_Comm_size(*local, &localsize);
  if(localsize != 2) {
    printf("ERROR: size of client-server communicator must be 2\n");
    return -1;
  }

  MPI_Barrier(*local);
  MPI_Sendrecv(&handshake, 1, MPI_INTEGER, 1-localrank, INTERCOMM_TAG, &answer, 1, MPI_INTEGER, 1-localrank, INTERCOMM_TAG, *local, &status);
  if(answer != 1) {
    printf("ERROR: bad handshake answer, expected 1, got %d\n",answer);
    return -1;
  }
  printf("INFO: handshake successful\n");
  return 2 ;
}

// data         : where to store data
// n            : number of data words
// disp         : offset into the remote window (in words)
// rankoftarget : rank of target process in communicator associated with window
// window       : MPI window
// lock         : a 'lock' is needed (single request)
int MPI_Get_words_simple(void *data, int n, int disp, int rankoftarget, MPI_Win window, int lock){   // get 4 byte words from remote window
  MPI_Aint TargetDisp = disp;
  int value;
  if(lock) MPI_Win_lock(MPI_LOCK_SHARED, rankoftarget, 0, window);
  value = MPI_Get(data, n, MPI_INTEGER, rankoftarget, TargetDisp, n, MPI_INTEGER, window);
  if(lock) MPI_Win_unlock(rankoftarget, window);
  return value;
}

// data         : where to fetch data from
// n            : number of data words
// disp         : offset into the remote window (in words)
// rankoftarget : rank of target process in communicator associated with window
// window       : MPI window
// lock         : a 'lock' is needed (single request)
int MPI_Put_words_simple(void *data, int n, int disp, int rankoftarget, MPI_Win window, int lock){   // put 4 byte words into remote window
  MPI_Aint TargetDisp = disp;
  int value;
  if(lock) MPI_Win_lock(MPI_LOCK_SHARED, rankoftarget, 0, window);
  value = MPI_Put(data, n, MPI_INTEGER, rankoftarget, TargetDisp, n, MPI_INTEGER, window);
  if(lock) MPI_Win_unlock(rankoftarget, window);
  return value;
}

// dst          : remote value BEFORE op
// src          : value to "op" at destination
// disp         : offset into the remote window (in words)
// rankoftarget : rank of target process in communicator associated with window
// window       : MPI window
// op           : MPI operator to apply at remote destination (add, min, max, ...)
// lock         : a 'lock' is needed (single request)
int MPI_Fetch_and_op_int_simple(void *src, void *dst, int disp, int rankoftarget, MPI_Win window, MPI_Op op, int lock){  // fetch and op, 1 integer
  MPI_Aint TargetDisp = disp;
  int value;
  if(lock) MPI_Win_lock(MPI_LOCK_SHARED, rankoftarget, 0, window);
  MPI_Get(dst, 1, MPI_INTEGER, rankoftarget, TargetDisp, 1, MPI_INTEGER, window);               // fetch before op
  MPI_Accumulate(src, 1, MPI_INTEGER,  rankoftarget, TargetDisp, 1, MPI_INTEGER, op, window);   // op
//   value = MPI_Fetch_and_op(src, dst, MPI_INTEGER, rankoftarget, TargetDisp, op, window);
  if(lock) MPI_Win_unlock(rankoftarget, window);
  return value;
}
