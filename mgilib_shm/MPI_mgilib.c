#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/shm.h>
#include <fcntl.h>
#include <string.h>
#include <malloc.h>

static int DEBUG=1;

// the following 3 routines are temporarily static and will have to be moved to their own source file
static int MPI_mgi_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name);
static int MPI_mgi_Publish_name(const char *service_name, MPI_Info info, const char *port_name);
static int MPI_mgi_Lookup_name(const char *service_name, MPI_Info info, char *port_name);

// the following 2 routines and associated declarations will have to be moved to their own source file
typedef void (*fptr)(void);
#define MAX_AT_TABLE 10
static fptr table[MAX_AT_TABLE];
static int nf=-1;
int MPI_Finalize();
int at_MPI_Finalize(fptr callback);

/*
 *     data buffer layout
 * 
 *     first : points to first element of buffer            (never updated)
 *     limit : points one past the last element in buffer   (never updated)
 *     in    : points to insertion position                 (updated by process that inserts data)
 *     out   : points to extraction position                (updated by process that extracts data)
 * 
 *     empty buffer, in = out
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | +
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                                                     ^
 *      | first          in | out                                                                 | limit
 *     ( 0 data elements in buffer, room for (limit - first - 1)
 * 
 *     buffer neither empty nor full
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     | | | | | | | | | | |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d| | | | | | | | | | | | | | | | +
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                     ^                               ^
 *      | first             | out                                 | in                            | limit
 *     (in - out) data elements in buffer, room for (limit - in) + (out - in -1) more
 * 
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     |d|d|d|d|d|d|d|d|d|d| | | | | | | | | | | | | | | | | | | |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d+
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                     ^                               ^
 *      | first             | in                                  | out                           | limit
 *     (in - first) + limit - out) data elements in buffer, room for (out - in + 1)
 * 
 * 
 *     full buffer, in = (out -1) (modulo limit)
 *     +-----------------------------------------------------------------------------------------+
 *     |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d| |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d+
 *     +-----------------------------------------------------------------------------------------+
 *      ^                                                       ^ ^                               ^
 *      | first                                              in | | out                           | limit
 */

int MPI_mgi_close(int channel);
int MPI_mgi_init(const char *alias);
int MPI_mgi_create(const char *alias, char mode, int cpl);
int MPI_mgi_create_endf(int cpl, MPI_Fint f_comm1, MPI_Fint f_comm2);
int MPI_mgi_create_end(int cpl, MPI_Comm comm1, MPI_Comm comm2);
int MPI_mgi_create_beginf(int cpl, MPI_Fint f_comm1, MPI_Fint f_comm2);
int MPI_mgi_create_begin(int cpl, MPI_Comm comm1, MPI_Comm comm2);
void MPI_mgi_closeall(void);
int MPI_mgi_term(void);
int MPI_mgi_open(int mpi_channel, unsigned char *mode);
int MPI_mgi_read(int channel, unsigned char *data, int nelm, unsigned char *dtyp);
int MPI_mgi_write(int channel, unsigned char *data, int nelm, unsigned char *dtyp);
void *MPI_mgi_memptr(int channel);
static int MPI_Can_Publish_name(const char *service_name, int test);

#if defined(WITH_FORTRAN_MGI)
ftnword f77_name (mgi_init) (char *channel_name, F2Cl lname);
ftnword f77_name (mgi_open) (ftnword *f_chan, char *mode, F2Cl lmode);
ftnword f77_name (mgi_read) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
ftnword f77_name (mgi_write) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
ftnword f77_name (mgi_clos) (ftnword *f_chan);
ftnword f77_name (mgi_term) (void);
#endif

static int in_closeall = 0 ;
static int closeall_done = 0 ;

static int debug_rank = -1;          // normally only needed for self test mode

int MPI_Finalize(){                  // interceptor, there are things we may want to do at MPI_Finalize time
  while (nf >= 0) {
if(DEBUG) printf("DEBUG %d: custom MPI_Finalize, nf = %d\n",debug_rank, nf);
    (*table[nf--])() ; // call registered routines to be called at mpi_finalize in reverse order
  }
  return PMPI_Finalize();            // now call the real MPI finalize routine
}
int at_MPI_Finalize(fptr callback){
  if(nf >= MAX_AT_TABLE-1) return -1;
  table[++nf] = callback;
if(DEBUG) printf("DEBUG %d: inserted callback %d\n",debug_rank,nf);
  return(0);
}

#define CHANNEL_UNKNOWN -1
#define CHANNEL_IDLE 0
#define CHANNEL_ACTIVE 1
#define CHANNEL_STOPPED 2

// r_control ... r_limit will be the image of control ... limit from remote memory arena
typedef struct{
  int control;    // flags
  int first;      // start of buffer
  int in;         // insertion index
  int out;        // extraction index
  int limit;      // end of buffer + 1

  int r_control;  // remote flags
  int r_first;    // remote first
  int r_in;       // remote in
  int r_out;      // remote out
  int r_limit;    // remote limit

  int data[1];  // place holder, start of data buffer
} arena;

typedef struct{
  char alias[32];         // short channel name
  char *channel_name;     // full MPI channel name
  char *port_name;        // MPI port name (from MPI_Open_port)
  arena *winbuf;          // local one sided window memory address   (local memory arena)
  MPI_Win window;         // one sided MPI window communicator
  MPI_Comm global, local; // inter-communicator, intra-communicator, both expected to have same size membership (2 members) 
  int is_server;          // 0 for client, 1 for server
  int is_active;          // 1 if active, 0 if inactive
  int winsize;            // window size in KBytes
  int thispe;             // my rank in intra-communicator for this channel
  int otherpe;            // rank of remote process in intra-communicator for this channel
  char mode;              // read or write
  char allocmem;          // 1 if allocated with malloc, 0 if node shared memory
  char pe0;               // 1 if PE0 of application, 0 otherwise
  char cpl;               // 1 if coupler PE, 0 otherwise
} MPI_mgi_channel ;

#define MAXALIAS 31
#define MAX_CHANNELS 8
static MPI_mgi_channel mpi_channel_table[MAX_CHANNELS];
static int last_mpi_channel=-1;

void *MPI_mgi_memptr(int channel){
  if(channel >= MAX_CHANNELS) return NULL;
  return mpi_channel_table[last_mpi_channel].winbuf;
}

static int timeout = 1000000;  // 1000 seconds by default
int MPI_mgi_timeout(int new_timeout){
  int old = timeout;
  timeout = new_timeout;
  return old;
}

int MPI_mgi_term() // close the books (we prefer to wait for MPI_Finalize in order to avoid undue wait for close)
{
//   if i am a PE0, push arena->control = CHANNEL_STOPPED to all remotes (lazy way would be set my own to CHANNEL_STOPPED)
//   MPI_mgi_closeall() ;
}

// cpl  : -1 coupler process, 0 PE0 from application >0 PEn from application
// comm used to cummunicate between PE0 and coupler process
int MPI_mgi_create_begin(int cpl, MPI_Comm comm1, MPI_Comm comm2)  // call before any call to MPI_mgi_create
{
  return(0);
}

int MPI_mgi_create_begin_f(int cpl, MPI_Fint f_comm1, MPI_Fint f_comm2)  // call before any call to MPI_mgi_create
{
  MPI_Comm comm1 = MPI_Comm_f2c(f_comm1);
  MPI_Comm comm2 = MPI_Comm_f2c(f_comm2);

  return MPI_mgi_create_begin(cpl, comm1, comm2);
}

// cpl  : -1 coupler process, 0 PE0 from application >0 PEn from application
// comm used to cummunicate between PE0 and coupler process
int MPI_mgi_create_end(int cpl, MPI_Comm comm1, MPI_Comm comm2)  // call after last call to MPI_mgi_create
{
  int i, status, id, data_length;
  int dispunit = sizeof(int);
  MPI_Aint winsize;
  char *memptr;
  struct shmid_ds shmbuf;
  MPI_Status mpi_stat;
  char mode;
  arena *arenaptr;
  int window_size;

//   if(cpl > 0) return(0);  // nothing to do on normal nodes
  if(comm2 == MPI_COMM_NULL) return(0);  // nothing to do on normal nodes

  for(i=0 ; i<=last_mpi_channel ; i++) {

    mode = mpi_channel_table[i].mode;
    if(mode != 'r' && mode != 'R' && mode != 'w' && mode != 'W' ) continue ;  // we have nothing to do with this channel

//     if(cpl == 0){      // PE 0
    if(comm1 != MPI_COMM_NULL){      // PE 0
      if(mpi_channel_table[i].mode == 'R'){
        mpi_channel_table[i].is_active = 1;
if(DEBUG) printf("DEBUG %d: receiving shmid\n",debug_rank);
        status = MPI_Recv(&id, 1, MPI_INTEGER, 0, 0, comm2, &mpi_stat);    // get shared memory segment id from coupler process
if(DEBUG) printf("DEBUG %d: received shmid\n",debug_rank);
        mpi_channel_table[i].winbuf = shmat(id,NULL,0);                    // attach to it, mpi_channel_table[i].winbuf => shared memory segment
        mpi_channel_table[i].window = MPI_WIN_NULL;                        // no one sided window on PE0 if reading
if(DEBUG) printf("DEBUG %d: PE0 attached to id=%d, address=%p, '%s'\n",debug_rank,id,mpi_channel_table[i].winbuf,mpi_channel_table[i].alias);
      }
      if(mpi_channel_table[i].mode == 'W'){
        mpi_channel_table[i].is_active = 1;
        MPI_mgi_Lookup_name(mpi_channel_table[i].channel_name, MPI_INFO_NULL, mpi_channel_table[i].port_name);  // get port name published under name channel_name
if(DEBUG) printf("DEBUG %d: connecting to '%s'\n",debug_rank,mpi_channel_table[i].port_name);
        status = MPI_Comm_connect(mpi_channel_table[i].port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &mpi_channel_table[i].global );   // connect to port
if(DEBUG) printf("DEBUG %d: connected to '%s'\n",debug_rank,mpi_channel_table[i].port_name);
        MPI_Intercomm_merge(mpi_channel_table[i].global, 1, &mpi_channel_table[i].local);
        MPI_Comm_rank(mpi_channel_table[i].local,&mpi_channel_table[i].thispe);
        mpi_channel_table[i].otherpe = 1 - mpi_channel_table[i].thispe;
        mpi_channel_table[i].winbuf = (arena *) malloc(sizeof(arena));
        mpi_channel_table[i].allocmem = 1;
        arenaptr = mpi_channel_table[i].winbuf;            // local one-sided memory area
        arenaptr->control = CHANNEL_UNKNOWN;
        arenaptr->first   = -1;
        arenaptr->in      = -1;
        arenaptr->out     = -1;
        arenaptr->limit   = -1;
        arenaptr->r_control = CHANNEL_UNKNOWN;            // remote copy of one-sided memory area parameters
        arenaptr->r_first   = -1;
        arenaptr->r_in      = -1;
        arenaptr->r_out     = -1;
        arenaptr->r_limit   = -1;

        winsize = sizeof(arena);
        status = MPI_Win_create(mpi_channel_table[i].winbuf, winsize, dispunit, MPI_INFO_NULL, mpi_channel_table[i].local, &mpi_channel_table[i].window);  // create window 
if(DEBUG) printf("DEBUG %d: PE0 short window buffer at %p '%s'\n",debug_rank,mpi_channel_table[i].winbuf,mpi_channel_table[i].alias);
      }
    }else{             // coupler PE
      if(mpi_channel_table[i].mode == 'R'){
        mpi_channel_table[i].is_active = 1;
        id = shmget(IPC_PRIVATE,mpi_channel_table[i].winsize,IPC_CREAT|S_IRUSR|S_IWUSR);  // allocate mpi_channel_table[i].winsize sized shared memory segment
        mpi_channel_table[i].winbuf = shmat(id,NULL,0);                                   // mpi_channel_table[i].winbuf => shared memory segment
if(DEBUG) printf("DEBUG %d: CPL created id=%d, address=%p '%s'\n",debug_rank,id,mpi_channel_table[i].winbuf,mpi_channel_table[i].alias);
        shmctl(id,IPC_RMID,&shmbuf);                                                      // mark segment for deletion
        arenaptr = mpi_channel_table[i].winbuf;            // local one-sided memory area
        data_length = &(arenaptr->data[0]) - &(arenaptr->control);   // offset of first data element with respect to arena
        arenaptr->control = CHANNEL_ACTIVE;             // setup of own memory arena
        arenaptr->first   = data_length;                // will use memptr[first] to memptr[limit-1]
        arenaptr->in      = arenaptr->first;            // in = out = first, buffer starts empty
        arenaptr->out     = arenaptr->in;
        window_size       = mpi_channel_table[i].winsize / sizeof(int);
        arenaptr->limit   = window_size - data_length;
        arenaptr->r_control = CHANNEL_UNKNOWN;          // mark remote arena control block status as unknown
        arenaptr->r_first   = -1;
        arenaptr->r_in      = -1;
        arenaptr->r_out     = -1;
        arenaptr->r_limit   = -1;

if(DEBUG) printf("DEBUG %d: sending shmid\n",debug_rank);
        status = MPI_Send(&id, 1, MPI_INTEGER, 1, 0,  comm2);                             // send segment id to PE0
if(DEBUG) printf("DEBUG %d: sent shmid, accepting on '%s'\n",debug_rank,mpi_channel_table[i].port_name);
        status = MPI_Comm_accept(mpi_channel_table[i].port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &mpi_channel_table[i].global );  // accept from port
if(DEBUG) printf("DEBUG %d: accepted on '%s'\n",debug_rank,mpi_channel_table[i].port_name);
        MPI_Intercomm_merge(mpi_channel_table[i].global, 0, &mpi_channel_table[i].local);
        MPI_Comm_rank(mpi_channel_table[i].local,&mpi_channel_table[i].thispe);
        mpi_channel_table[i].otherpe = 1 - mpi_channel_table[i].thispe;
        winsize = sizeof(arena) + mpi_channel_table[i].winsize * dispunit;
        status = MPI_Win_create(mpi_channel_table[i].winbuf, winsize, dispunit, MPI_INFO_NULL, mpi_channel_table[i].local, &mpi_channel_table[i].window);
      }
      if(mpi_channel_table[i].mode == 'W'){
        mpi_channel_table[i].window = MPI_WIN_NULL;
        mpi_channel_table[i].winbuf = NULL;
        continue ; // nothing to do for coupler in write mode
      }
    }
  }
}

int MPI_mgi_create_end_f(int cpl, MPI_Fint f_comm1, MPI_Fint f_comm2)  // call after last call to MPI_mgi_create
{
  MPI_Comm comm1 = MPI_Comm_f2c(f_comm1);
  MPI_Comm comm2 = MPI_Comm_f2c(f_comm2);

  return MPI_mgi_create_end(cpl, comm1, comm2);
}
// create one channel (it is assumed that no mgi activity has started)
// this channel will be opened as "alias"
// application must call this routine once for each "alias" channel name it uses
// mode : 'r' or 'w'
// cpl  : -1 coupler process, 0 PE0 from application >0 PEn from application
//
// export MGI_MPI_CFG=" mastername n : size_1 name_1 : ... : size_n name_n "
//  mastername         : experiment name (MUST be unique at any given time)
//  name_1 ...  name_n : channel names 
//
// return value : 0 if open + publish successful
//                1 MPI channel already published via another alias
//               -1 ERROR
int MPI_mgi_create(const char *alias, char mode, int cpl) 
{
  int status, status2 ;
  char *cfg ;
  char *service_name ;
  char mastername[128] ;
  int i, ordinal ;

  if(cpl > 0) return(0);      // nothing to do on normal nodes
  if(mode != 'r' && mode != 'R' && mode != 'w' && mode != 'W' ) return(0) ; // neither read nor write, nothing to do

  if(last_mpi_channel == -1){     // initialize global channel table from environment variable MGI_MPI_CFG (do only once)
    cfg = getenv("MGI_MPI_CFG") ;
    if(cfg == NULL) return(-1) ;
    sscanf(cfg,"%64s%d",mastername,&last_mpi_channel) ;   // get mastername and number of mpi channels
// if(DEBUG) printf("DEBUG: MPI channel prefix = '%s', mpi channels = %d, alias = '%s'\n",mastername,last_mpi_channel,alias);
    last_mpi_channel--;
    for(i=0 ; i<=last_mpi_channel ; i++) {    // get size, alias1, alias2 for all channels, build channel name
      while(*cfg != ':') {                    // skip to : delimiter
        if(*cfg == '\0') return(-1);          // premature termination of configuration string
        cfg++;
      }
      cfg++;
      sscanf(cfg,"%d%31s",&mpi_channel_table[i].winsize,mpi_channel_table[i].alias);  // get channel size and short name
      mpi_channel_table[i].winsize *= 1024;  // size is in KBytes
// if(DEBUG) printf("DEBUG: cfg = '%s'\n",cfg);
      mpi_channel_table[i].channel_name = (char *) malloc(257);
      memset(mpi_channel_table[i].channel_name,0,257);
      snprintf(mpi_channel_table[i].channel_name,257,"%s_%s",mastername,mpi_channel_table[i].alias);   // build full channel name
      mpi_channel_table[i].port_name = (char *) malloc(MPI_MAX_PORT_NAME+1);
      memset(mpi_channel_table[i].port_name,0,MPI_MAX_PORT_NAME+1);

// if(DEBUG) printf("DEBUG %d: i = %d, master = '%s', channel = '%s', size = %d\n",
//        debug_rank,i,
//        mpi_channel_table[i].channel_name,
//        mpi_channel_table[i].alias,
//        mpi_channel_table[i].winsize);

      mpi_channel_table[i].winbuf    = NULL ;                        // window not created yet
      mpi_channel_table[i].window    = MPI_WIN_NULL ;
      mpi_channel_table[i].global    = MPI_COMM_NULL ;
      mpi_channel_table[i].local     = MPI_COMM_NULL ;
      mpi_channel_table[i].is_server = 0 ;
      mpi_channel_table[i].is_active = 0 ;
      mpi_channel_table[i].thispe    = -1 ;
      mpi_channel_table[i].otherpe   = -1 ;
      mpi_channel_table[i].mode      = 'U' ;
      mpi_channel_table[i].allocmem  = 0;
      mpi_channel_table[i].pe0       = 0;
      mpi_channel_table[i].cpl       = 0;
    }   // end for i
    at_MPI_Finalize(MPI_mgi_closeall) ;// setup for at_MPI_Finalize (close everything at finalize time)
  }

  ordinal = -1 ;
  for(i=0 ; i<=last_mpi_channel ; i++) {
    if( strncmp(alias,mpi_channel_table[i].alias,MAXALIAS) == 0 ) {   // alias found
      mpi_channel_table[i].mode =  toupper(mode);
      ordinal = i ;
      service_name = mpi_channel_table[i].channel_name ; // MPI channel name is 'mastername_alias'
    }
  }
// if(DEBUG) printf("DEBUG %d: channel '%s', mode '%c', ordinal = %d\n",debug_rank,alias,mode,ordinal);
  if(ordinal < 0) return(-1) ;       // alias not found , OOPS 
  if(cpl == 0) mpi_channel_table[i].pe0 = 1;    // PE0 from application
  if(cpl == -1) mpi_channel_table[i].cpl = 1;   // this is a coupler PE
  if ( mpi_channel_table[ordinal].mode == 'W' ) return(0) ;  // nothing more to do in write mode
  if ( cpl >= 0) return(0) ;  // nothing more to do for PE0 -> PEn

  status2 = MPI_Can_Publish_name(service_name,0) ;    // is service_name already published 
// if(DEBUG) printf("DEBUG %d: to be published , status2 = %d\n",debug_rank, status2);
  if(status2 < 0) return(1);                          // yes, OOPS
  
// if(DEBUG) printf("DEBUG %d: opening port\n",debug_rank);
  status = MPI_Open_port(MPI_INFO_NULL, mpi_channel_table[ordinal].port_name);    // create port, get port name
  if (status == MPI_SUCCESS) {                    // publish port under name channel_name
    status = MPI_mgi_Publish_name(service_name, MPI_INFO_NULL, mpi_channel_table[ordinal].port_name);
if(DEBUG) printf("DEBUG %d: ordinal,server,status = %d %d %d,published port = '%s'\n",
       debug_rank,ordinal,mpi_channel_table[ordinal].is_server,status,mpi_channel_table[ordinal].port_name);
  }
  if(status != MPI_SUCCESS) {              // port creation or publication failed
    free(mpi_channel_table[ordinal].channel_name) ; mpi_channel_table[ordinal].channel_name = NULL ;
//     free(mpi_channel_table[ordinal].alias1)       ; mpi_channel_table[ordinal].alias1 = NULL ;
//     free(mpi_channel_table[ordinal].alias2)       ; mpi_channel_table[ordinal].alias2 = NULL ;
    free(mpi_channel_table[ordinal].port_name)    ; mpi_channel_table[ordinal].port_name = NULL ;
    return(-1) ; 
  }
  mpi_channel_table[ordinal].is_server = 1 ;  // set status to server only if open port and publish successful

  return(0);  // everything O.K., port created and successfully published
}

int MPI_mgi_close(int channel) // close a MPI channel (it is assumed that there is no pending activity)
{
  int status;
//   int channel = mpi_channel ;  // convert pseudo channel to MPI channel

// if(DEBUG) printf("DEBUG %d: closing channel %d '%s'\n",debug_rank,channel,mpi_channel_table[channel].alias);
  if(channel >= MAX_CHANNELS || channel < 0) return -1;             // invalid channel number
  if(mpi_channel_table[channel].channel_name == NULL) return -1;    // channel not/no longer active
  if(mpi_channel_table[channel].is_active == 0) return -1;          // channel never active or already closed

  mpi_channel_table[channel].is_active = -1;                        // mark channel as inactive but not totally closed
  //
  // if i am writing (PE0), set remode control flag to CHANNEL_STOPPED with MPI_Put (should have been done by mgi_term)
  // if i am reading (CPL), it is an ERROR if arena->control if CHANNEL_ACTIVE (would actually be stuck in waiting loop)
  //
if(DEBUG) printf("DEBUG %d: channel %d marked as inactive\n",debug_rank,channel);
  if(in_closeall == 0) return 0;                                    // wait for closeall to really close the books
  mpi_channel_table[channel].is_active = 0;                         // mark channel as inactive and totally closed
  if(mpi_channel_table[channel].window == MPI_WIN_NULL) {
if(DEBUG) printf("DEBUG %d: no window on channel = %d '%s'\n",debug_rank,channel,mpi_channel_table[channel].alias);
    return 0;   // no window = no communicator
  }
// if(DEBUG) printf("DEBUG %d: barrier on local, channel = %d '%s'\n",debug_rank,channel,mpi_channel_table[channel].alias);
  MPI_Barrier(mpi_channel_table[channel].local);
  if(mpi_channel_table[channel].window != MPI_WIN_NULL){
// if(DEBUG) printf("DEBUG %d: freeing window, channel = %d, win = %p\n",debug_rank,channel,mpi_channel_table[channel].winbuf);
    status = MPI_Win_free(&mpi_channel_table[channel].window);               // free window communicator
    mpi_channel_table[channel].window = MPI_WIN_NULL;
if(DEBUG) printf("DEBUG %d: freed window, channel = %d, status = %d, win = %p\n",debug_rank,channel,status,mpi_channel_table[channel].winbuf);
     if(mpi_channel_table[channel].allocmem) {
// if(DEBUG) printf("DEBUG %d: freing window memory, channel = %d, status = %d\n",debug_rank,channel,status);
       free(mpi_channel_table[channel].winbuf);
//     status = MPI_Free_mem(mpi_channel_table[channel].winbuf);                        // free memory associated with window
if(DEBUG) printf("DEBUG %d: freed window memory, channel = %d, status = %d\n",debug_rank,channel,status);
     }
    mpi_channel_table[channel].winbuf = NULL;
  }

  MPI_Barrier(mpi_channel_table[channel].local);

// if(DEBUG) printf("DEBUG %d: disconnecting local, channel = %d\n",debug_rank,channel);
  MPI_Comm_disconnect( &mpi_channel_table[channel].local );         // disconnect intra-communicator
// if(DEBUG) printf("DEBUG %d: disconnecting global, channel = %d\n",debug_rank,channel);
  MPI_Comm_disconnect( &mpi_channel_table[channel].global );        // disconnect iner-communicator

  if(mpi_channel_table[channel].is_server) {                        // on server, unpublish and close port
    MPI_mgi_Unpublish_name(mpi_channel_table[channel].channel_name, MPI_INFO_NULL, "no_port_name");
    MPI_Close_port(mpi_channel_table[channel].port_name);
  }
  MPI_Can_Publish_name(mpi_channel_table[channel].channel_name,1) ; // remove lock file if not already done

  free(mpi_channel_table[channel].channel_name);                    // free channel name array
  mpi_channel_table[channel].channel_name = NULL;
  free(mpi_channel_table[channel].port_name);                       // free port name array
  mpi_channel_table[channel].port_name = NULL;
  mpi_channel_table[channel].is_server =  0;                        // invalidate server flag
  mpi_channel_table[channel].thispe    = -1;                        // invalidate ranks
  mpi_channel_table[channel].otherpe   = -1;
if(DEBUG) printf("DEBUG %d: closed channel %d\n",debug_rank,channel);
  return 0;
}

void MPI_mgi_closeall(void)  // close all channels if not already done
{
  int i;
// if(DEBUG) printf("DEBUG %d: entering closeall\n",debug_rank);
  if(closeall_done) return ;  // already done
  in_closeall = 1 ;
  for(i=0 ; i<=last_mpi_channel ; i++) MPI_mgi_close(i);
  in_closeall = 0 ;
  closeall_done = 1 ;
// if(DEBUG) printf("DEBUG %d: exiting closeall\n",debug_rank);
}

//
// this function MUST BE CALLED ONLY BY A RANK 0 PROCESS
//
// MGI_MPI_CFG=" prefix n : size1 aname_1 bname_1 : ... : aname_n bname_n "
//
// this function has essentially become : get pseudo channel number associated with alias
// even if alias matches alias1, odd if alias matches alias2
// -1 is returned upon error
int MPI_mgi_init(const char *alias)
{
  int i; // , mpi_channel, bump;

//   mpi_channel = -1 ;
//   bump = -1;
  for(i=0 ; i<=last_mpi_channel ; i++){
// if(DEBUG) printf("DEBUG: alias, alias1, alias2 = '%s'%s'%s'\n",alias,mpi_channel_table[i].alias1,mpi_channel_table[i].alias2);
    if( strncmp(mpi_channel_table[i].alias,alias,MAXALIAS) == 0 && mpi_channel_table[i].mode != 'U') {
      return(i) ;
    }
//     if( strncmp(mpi_channel_table[i].alias2,alias,MAXALIAS) == 0 ) bump = 1 ;
//     if( bump != -1 ){
//       mpi_channel = i ;
//       break ;
//     }
  }
//   return( (bump < 0) ? bump : mpi_channel*2 + bump) ;   // return -1 if alias not found
  return(-1);
}

//
// this function MUST BE CALLED ONLY BY A RANK 0 PROCESS
//
// mpi_channel  channel number obtained from MPI_mgi_init
// mode         'R' or 'W' (will be enforced by MPI_mgi_read and MPI_mgi_write
// result       -1 in case of failure, otherwise channel number (>-0) for read/write/close calls
//
// setup of the one sided communication window, setup of intercommunicator after accept/connect
// setup of local first/in/out/limit, get remote first/in/out/limit
//
// read(R)/write(W) mode is really for backwards compatibility, the MPI channel is really bi-directional
// each alias subchannel is treated as uni-directional
int MPI_mgi_open(int mpi_channel,unsigned char *mode)
{
  MPI_Win window;
  MPI_Comm local;
  void *remote;
  MPI_Aint TargetDisp;
  int data_length ;
  arena *arenaptr ;
  int otherpe, i ;
  int window_size ;

  if(mpi_channel < 0 || mpi_channel > last_mpi_channel ) return(-1) ;     // bad channel number
  if(toupper(*mode) != 'R' && toupper(*mode) != 'W') return(-1) ;         // bad mode
  if(mpi_channel_table[mpi_channel].is_active != 1) {
    printf("ERROR: MPI channel not initialized properly\n");
    return(-1);
  }
  arenaptr = mpi_channel_table[mpi_channel].winbuf;
if(DEBUG) printf("DEBUG %d: MPI_mgi_open : channel = %d, mode = '%c'\n",debug_rank,mpi_channel,*mode);
  if(toupper(*mode) == 'W') return mpi_channel ;
if(DEBUG) printf("DEBUG %d: local window parameters %d %d %d %d %d\n",debug_rank,arenaptr->control,arenaptr->first,arenaptr->in,arenaptr->out,arenaptr->limit);

//   MPI_Barrier(mpi_channel_table[mpi_channel].local);    // sync with remote process, so that we are sure to get correct control ... limit values from it
// if(DEBUG) printf("DEBUG : remote window parameters %d %d %d %d %d\n",arenaptr->r_control,arenaptr->r_first,arenaptr->r_in,arenaptr->r_out,arenaptr->r_limit);

// get remote memory arena configuration from "otherpe" (later we will "put" remote in and "get" remote out)
//   data_length = &(arenaptr->r_control) - &(arenaptr->control);
//   TargetDisp = 0;        //  what we get is right at the beginning of memory area, length = data_length integers
//   remote = &(arenaptr->r_control) ;
//   MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
//   MPI_Get(remote, data_length, MPI_INTEGER, otherpe, TargetDisp, data_length, MPI_INTEGER, window); // get remote control ... limit
//   MPI_Win_unlock(otherpe,window);

// if(DEBUG) printf("DEBUG : remote window parameters %d %d %d %d %d\n",arenaptr->r_control,arenaptr->r_first,arenaptr->r_in,arenaptr->r_out,arenaptr->r_limit);
//   MPI_Barrier(mpi_channel_table[mpi_channel].local);    // sync with remote process done

  return mpi_channel ;                                                  // return mpi channel number;
}
int MPI_mgi_open_old(int mpi_channel,unsigned char *mode)
{
  char *port_name;
  char *channel_name ;
  MPI_Win window;
  int rank;
  MPI_Comm global, local;
  MPI_Aint winsize;
  int dispunit = sizeof(int);
  void *memptr;
  void *remote;
  MPI_Aint TargetDisp;
  int data_length ;
  arena *arenaptr ;
  int otherpe, i ;
  int server ;
  int window_size ;
//   int mpi_channel = channel;      // get MPI channel number
  int status ;

  if(mpi_channel < 0 || mpi_channel > last_mpi_channel ) return(-1) ;     // bad channel number
  if(toupper(*mode) != 'R' && toupper(*mode) != 'W') return(-1) ;         // bad mode
  if(mpi_channel_table[mpi_channel].is_active != 0) {
    printf("ERROR: attempt to open an already open MPI channel\n");
    return(-1);
  }

  mpi_channel_table[mpi_channel].mode = toupper(*mode);    // set mode for sub channel (backwards compatibility)
  channel_name = mpi_channel_table[mpi_channel].channel_name ;   // true_name of MPI channel addresses as alias1 or alias2
  window_size  = mpi_channel_table[mpi_channel].winsize ;        // size in "int" units
  winsize      = dispunit * window_size ;                        // size in bytes

  server = mpi_channel_table[mpi_channel].is_server;             // get server flag

  if(! server) mpi_channel_table[mpi_channel].port_name=malloc(MPI_MAX_PORT_NAME+1);  // not "server", allocate space for name lookup
  port_name = mpi_channel_table[mpi_channel].port_name;

  if(server){
// if(DEBUG) printf("DEBUG : accept on '%s'\n",port_name);
    status = MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &global );  // wait for client to connect (collective over communicator)
if(DEBUG) printf("DEBUG : accepted on '%s' %d\n",port_name,status);
    MPI_Intercomm_merge(global, 0, &local);                                  // create communicator for one-sided window
  }else{
    MPI_mgi_Lookup_name(channel_name, MPI_INFO_NULL, port_name);             // get port name published under name channel_name
// if(DEBUG) printf("DEBUG : connect to '%s' = '%s'\n",channel_name,port_name);
    status = MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &global ); // wait until connected to server (collective over communicator)
if(DEBUG) printf("DEBUG : connected to '%s' = '%s' %d\n",channel_name,port_name,status);
    MPI_Intercomm_merge(global, 1, &local);                                  // create communicator for one-sided window
  }
// if(DEBUG) printf("DEBUG : intercom merged\n");
  MPI_Comm_rank(local,&mpi_channel_table[mpi_channel].thispe);               // size of communicator local should be 2
  otherpe = 1 - mpi_channel_table[mpi_channel].thispe;
// if(DEBUG) printf("DEBUG : otherpe = %d\n",otherpe);
  mpi_channel_table[mpi_channel].otherpe = otherpe;                          // otherpe is 0 if i am 1, 1 if i am 0
  mpi_channel_table[mpi_channel].global = global;
  mpi_channel_table[mpi_channel].local  = local;
  mpi_channel_table[mpi_channel].is_active = 1;

// if(DEBUG) printf("DEBUG : allocating %ld bytes\n",winsize);
  status = MPI_Alloc_mem(winsize,MPI_INFO_NULL,&memptr);                              // allocate local memory for one-sided window
// if(DEBUG) printf("DEBUG : memory at %p %d\n",memptr,status);
  mpi_channel_table[mpi_channel].winbuf = memptr;                            // local memory address of 1 sided window
  status = MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local, &window);  // create one-sided window (local buffer to receive remote writes)
if(DEBUG) printf("DEBUG : window created %d\n",status);
  mpi_channel_table[mpi_channel].window = window;

  arenaptr = (arena *)memptr ;                                               // local one-sided memory area
// setup of own memory arena
  data_length = &(arenaptr->data[0]) - &(arenaptr->control);   // offset of first data element with respect to arena
  arenaptr->control = CHANNEL_ACTIVE;
  arenaptr->first   = data_length;                // will use memptr[first] to memptr[limit-1]
  arenaptr->in      = arenaptr->first;            // in = out = first, buffer starts empty
  arenaptr->out     = arenaptr->in;
  arenaptr->limit   = window_size ;

  MPI_Barrier(local);    // sync with remote process, so that we are sure to get correct control ... limit values from it

if(DEBUG) printf("DEBUG : local window parameters %d %d %d %d %d\n",arenaptr->control,arenaptr->first,arenaptr->in,arenaptr->out,arenaptr->limit);
// if(DEBUG) printf("DEBUG : remote window parameters %d %d %d %d %d\n",arenaptr->r_control,arenaptr->r_first,arenaptr->r_in,arenaptr->r_out,arenaptr->r_limit);

// get remote memory arena configuration from "otherpe" (later we will "put" remote in and "get" remote out)
  data_length = &(arenaptr->r_control) - &(arenaptr->control);
  TargetDisp = 0;        //  what we get is right at the beginning of memory area, length = data_length integers
  remote = &(arenaptr->r_control) ;
  MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
  MPI_Get(remote, data_length, MPI_INTEGER, otherpe, TargetDisp, data_length, MPI_INTEGER, window); // get remote control ... limit
  MPI_Win_unlock(otherpe,window);

if(DEBUG) printf("DEBUG : remote window parameters %d %d %d %d %d\n",arenaptr->r_control,arenaptr->r_first,arenaptr->r_in,arenaptr->r_out,arenaptr->r_limit);
  MPI_Barrier(local);    // sync with remote process done

  return mpi_channel ;                                                  // return mpi channel number;
}

// read nelm data elements of type *dtyp (I/R/D/C) into data from pseudo_channel
int MPI_mgi_read(int channel,unsigned char *data, int nelm, unsigned char *dtyp){
  MPI_Datatype mpitype;
  int nitems = nelm;
  int ntok;
  int avail = 0;
  arena *memory;
  int *buffer;
  useconds_t wait = 1000;   // 1 millisecond
  int sleepcount = 0;
  int meta;
  MPI_mgi_channel *local;
  int first, in, out, limit, navail, temp;
  size_t nbytes, nbytes1, nbytes2;
  unsigned char rtyp = '?';
  int nmeta, otherpe;
//   int channel = pseudo_channel ;      // get MPI channel number
  int *idata = (int *) data;
  MPI_Aint TargetDisp = 0;
  MPI_Win window;
  arena *winbuf;

  if(channel > last_mpi_channel || channel < 0) return -1;          // invalid channel number
  if(mpi_channel_table[channel].channel_name == NULL) return -1;    // channel not/no longer active
  if(mpi_channel_table[channel].is_active != 1) return -1;          // channel not/no longer active
  if(mpi_channel_table[channel].mode != 'R') return -1;             // wrong direction

  switch(*dtyp){    // tokens are 32 bit integers
    case 'R':       // 32 bit floating point data
      ntok = nitems;
      nbytes = nitems << 2;   // nitems * 4
      mpitype = MPI_FLOAT;
      break;
    case 'I':      // 32 bit integer data
      ntok = nitems;
      nbytes = nitems << 2;   // nitems * 4
      mpitype = MPI_INT;
      break;
    case 'C':      // character data
      ntok = (nitems + 3) >> 2;  //  (nitems rounded up to a multiple of 4) / 4
      nbytes = nitems;
      mpitype = MPI_CHAR;
      break;
    case 'D':     // 64 bit floating point data (we might want to downgrade on write, upgrade on read)
      ntok = nitems << 1;     // nitems * 2
      nbytes = nitems << 3;   // nitems * 8
      mpitype = MPI_DOUBLE;
      break;
    default :     // anything else is an error
      fprintf(stderr,"ERROR: MPI_mgi_read, bad type '%c', valid types are I/R/D/C\n",*dtyp);
      return -1;
      break;
  }
  memory = (arena *)mpi_channel_table[channel].winbuf;  // local memory arena
  buffer = (int *) memory;

if(DEBUG) printf("DEBUG : before read, memory->first = %d, memory->in = %d, memory->out = %d, memory->limit = %d \n",
                 memory->first,memory->in,memory->out,memory->limit);
// if(DEBUG) printf("DEBUG : memory = %p, buffer = %p, channel = %d\n",memory,buffer,channel);
// we need 1 metadata element + ntok data elements
  while(1){                      // wait until there is enough data in local buffer to satisfy request
    first = memory->first;       // index of start of data buffer
    in    = memory->in;          // will be updated by remote writer process
    out   = memory->out;         // will be updated by this process
    limit = memory->limit;       // bufer[limit] is one slot past end of buffer
// if(DEBUG) printf("DEBUG : in=%d, out=%d, navail=%d\n",in,out,in-out-1);
    navail = (in >= out) ? (in - out) : ((limit - out) + (in - first));  // number of available data elements
//     if(in == out) navail = 0;
// if(DEBUG) printf("DEBUG : in=%d, out=%d, navail = %d\n",in,out,navail);
    if(navail >= ntok) break; // enough data elements to satisfy request

    winbuf = mpi_channel_table[channel].winbuf;
    window = mpi_channel_table[channel].window;
    otherpe = mpi_channel_table[channel].otherpe;
    MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
// if(DEBUG) printf("DEBUG : locked window(in) at %d\n",otherpe);
//     TargetDisp =  &(winbuf->r_in) - &(winbuf->control) ;
//     MPI_Get(&temp, 1, MPI_INTEGER, otherpe, TargetDisp, 1, MPI_INTEGER, window);   // get 'my in' from remote server
// if(DEBUG) printf("DEBUG : got remote temp=%d, in=%d\n",temp,in);
//     in = temp;
    MPI_Win_unlock(otherpe,window);
// if(DEBUG) printf("DEBUG : unlocked window\n");
    usleep(wait);                   // not enough, wait for remote writer to update "in" (wait 1 ms by default)
    if(++sleepcount > timeout) {    // max timeout exceeded (default is 1000 seconds)
if(DEBUG) printf("DEBUG : timeout\n");
      return 0;
    }
  }   // local data buffer : buffer[first] to buffer[limit-1];

  if(out >= limit) out = first;  // wrap around
  meta = buffer[out];            // get metadata
  // validate metadata ( type + 'length' << 8 )   'length' = real length & 0xFFFFFF
  nmeta = meta >> 8;
  rtyp  = meta & 0xFF;
if(DEBUG) printf("DEBUG : sleep=%d, in=%d, out=%d, rtyp = %c, nmeta = %d\n",sleepcount,in,out,rtyp,nmeta);
  if(rtyp != *dtyp || nmeta != (nelm & 0xFFFFFF)) {
    fprintf(stderr,"ERROR: MPI_mgi_read, type/length mismatch, expected %d'%c', got %d,'%c'\n",nelm,*dtyp,nmeta,rtyp);
    return -1;
  }
  out++;
  if(out >= limit) out = first;

  if(in > out){                                            // data is in 1 piece
    memcpy(data,&(buffer[out]),nbytes);                    // buffer[out] -> buffer[out + ntok -1] : ntok tokens
    out += ntok;                                           // got data, update out
  }else{                                                   // data is split in 2 pieces
    nbytes1 = (limit - out) << 2;                          // (limit - out) * 4 bytes for 1st piece
    memcpy(data,&(buffer[out]),nbytes1);                   // buffer[out] -> buffer[limit-1] : limit - out tokens
    nbytes2 = nbytes - nbytes1;                            // nbytes = 1st piece
    memcpy(&(data[nbytes1]),&(buffer[first]),nbytes2);     // buffer[first] -> buffer[out-1] : ntok - (limit - out) tokens
    out = first + (ntok - (limit - out));                  // data has been copied, it is safe to update out
  }
  memory->out = out ;                                      // update out, remote partner will update itself if need be
if(DEBUG) printf("DEBUG : after read, in=%d, out=%d\n",in,out);
  return nelm;     // return number of elements read
}

// write nelm data elements of type *dtyp (I/R/D/C) from data into pseudo_channel 
int MPI_mgi_write(int channel, unsigned char *data, int nelm, unsigned char *dtyp){
  MPI_Datatype mpitype;
  int nitems = nelm;
  int ntok;
  int avail = 0;
  int first, in, out, limit, navail;
  size_t nbytes, nbytes1, nbytes2;
  int sleepcount = 0;
  int wait = 1000;
  arena *winbuf;
  arena temparena;
  MPI_Win window;
  MPI_Aint TargetDisp = 0;
  int meta;
  int *idata = (int *)data;
  int otherpe, status;
//   int channel = pseudo_channel ;      // get MPI channel number

  if(channel > last_mpi_channel || channel < 0) return -1;          // invalid channel number
  if(mpi_channel_table[channel].channel_name == NULL) return -1;    // channel not/no longer active
  if(mpi_channel_table[channel].is_active != 1) return -1;          // channel not/no longer active
  if(mpi_channel_table[channel].mode != 'W') return -1;             // wrong direction

if(DEBUG) printf("DEBUG : writing %d items of type %c from %p\n",nelm,*dtyp,data);
  switch(*dtyp){    // tokens are 32 bit integers
    case 'R':       // 32 bit floating point data
      ntok = nitems;
      nbytes = nitems << 2;   // nitems * 4
      mpitype = MPI_FLOAT;
      break;
    case 'I':      // 32 bit integer data
      ntok = nitems;
      nbytes = nitems << 2;   // nitems * 4
      mpitype = MPI_INT;
      break;
    case 'C':      // character data
      ntok = (nitems + 3) >> 2;  //  (nitems rounded up to a multiple of 4) / 4
      nbytes = nitems;
      mpitype = MPI_CHAR;
      break;
    case 'D':     // 64 bit floating point data (we might want to downgrade on write, upgrade on read)
      ntok = nitems << 1;     // nitems * 2
      nbytes = nitems << 3;   // nitems * 8
      mpitype = MPI_DOUBLE;
      break;
    default :     // anything else is an error
      fprintf(stderr,"ERROR: MPI_mgi_write, bad type '%c', valid types are I/R/D/C\n",*dtyp);
      return -1;
      break;
  }
  winbuf = mpi_channel_table[channel].winbuf;
  window = mpi_channel_table[channel].window;
  otherpe = mpi_channel_table[channel].otherpe;
// we may have to get remote out index to update r_out if not enough space
  while(1){      // wait until there is enough space in remote buffer to satisfy request (ntok + 1) data tokens
    if(winbuf->r_control = CHANNEL_UNKNOWN){   // not initialized yet, get remote area parameters
      MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
      MPI_Get(&temparena, 5, MPI_INTEGER, otherpe, 0, 5, MPI_INTEGER, window);   // get control info from remote server
      MPI_Win_unlock(otherpe,window);
      winbuf->r_control = temparena.control ;
      winbuf->r_first   = temparena.first ;
      winbuf->r_in      = temparena.in ;
      winbuf->r_out     = temparena.out ;
      winbuf->r_limit   = temparena.limit ;
    }
    first  = winbuf->r_first;
    in     = winbuf->r_in;
    out    = winbuf->r_out;
    limit  = winbuf->r_limit;
    navail = (in < out) ? (out - in - 1) : ((limit - in) + (out - first) - 1);  // number of available spaces
// if(DEBUG) printf("DEBUG : navail = %d, ntok = %d\n",navail,ntok);
    if(navail >= ntok + 1) break;  // enough space is available, OK for remote put

    TargetDisp = &(winbuf->r_out) - &(winbuf->r_control);  // not enough space, get remote 'out' to see if there is more
    MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
// if(DEBUG) printf("DEBUG : locked window(out) at %d\n",otherpe);
    MPI_Get(&out, 1, MPI_INTEGER, otherpe, TargetDisp, 1, MPI_INTEGER, window);   // get 'out' from remote server
// if(DEBUG) printf("DEBUG : got remote out=%d\n",out);
    MPI_Win_unlock(otherpe,window);
    winbuf->r_out = out;              // update local copy of remote out
// if(DEBUG) printf("DEBUG : unlocked window\n");

    navail = (in < out) ? (out - in - 1) : ((limit - in) + (out - first - 1));  // number of available tokens
    if(navail >= ntok) break;  // enough space is available

    usleep(wait);              // not enough space, remote reader needs to further increment increment 'out'
    if(sleepcount++ > timeout)  return ntok ;    // max timeout exceeded (default is ~1000 seconds)
  }

  meta = (*dtyp) & 0xFF;
  meta = meta | (nelm & 0xFFFFFF) << 8 ; // build and send 32 bit token of metadata header
  TargetDisp = in;
  MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
// if(DEBUG) printf("DEBUG : locked window(meta) at %d\n",otherpe);
  status = MPI_Put(&meta,1,MPI_INTEGER,otherpe,TargetDisp,1,MPI_INTEGER,window);      // remote write metadata
if(DEBUG) printf("DEBUG : metadata written %d, status = %d, sleep=%d\n", meta,status,sleepcount);
  MPI_Win_unlock(otherpe,window);
// if(DEBUG) printf("DEBUG : unlocked window\n");
  if(++in >= limit) in = first;    // wrap around

  MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
// if(DEBUG) printf("DEBUG : locked window(data) at %d\n",otherpe);
  TargetDisp = in;
  if(out <= in){
    MPI_Put(idata, ntok, MPI_INTEGER, otherpe, TargetDisp, ntok, MPI_INTEGER, window);     // remote write data
    in += ntok;
  }else{
    MPI_Put(idata, limit-in, MPI_INTEGER, otherpe, TargetDisp, limit-in, MPI_INTEGER, window);   // remote write data (part 1)
    TargetDisp = first;
    MPI_Put(&(idata[limit-in]), ntok-(limit-in), MPI_INTEGER, otherpe, TargetDisp, ntok-(limit-in), MPI_INTEGER, window);  // remote write data (part 2)
    in = first + (ntok-(limit-in));
  }
  MPI_Win_unlock(otherpe,window);
// if(DEBUG) printf("DEBUG : %d data elements written\n", nelm);

  winbuf->r_in = in;
// if(DEBUG) printf("DEBUG : remote in index updated = %d\n",winbuf->r_in);

  TargetDisp = &(winbuf->r_in) - &(winbuf->r_control);                  // offset of in in arena
  MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
// if(DEBUG) printf("DEBUG : locked window(in) at %d\n",otherpe);
  MPI_Put(&in,1,MPI_INTEGER,otherpe,TargetDisp,1,MPI_INTEGER,window);   // update in on remote server
if(DEBUG) printf("DEBUG : remote in written %d\n",in);
  MPI_Win_unlock(otherpe,window);
// if(DEBUG) printf("DEBUG : unlocked window\n");
  return 0;
}

static int MPI_mgi_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name)
{
  char filename[4096];
  char *mpi_mgi_home = getenv("MGI_MPI_HOME") ;

  if(mpi_mgi_home != NULL) {        // MPI channel files directory
    if( mpi_mgi_home[0] == '/' ){
      snprintf(filename,4096,"%s/%s.channel",mpi_mgi_home,service_name);
    }else{
      snprintf(filename,4096,"%s/%s/%s.channel",getenv("HOME"),mpi_mgi_home,service_name);
    }
  }else{
    snprintf(filename,4096,"%s/%s/%s.channel",getenv("HOME"),".gossip/MPI",service_name);
  }
  unlink(filename);                                  // remove channel file
  MPI_Can_Publish_name(service_name,1) ;             // remove lock file

  return MPI_SUCCESS;
}

static int MPI_Can_Publish_name(const char *service_name, int test)
{
  FILE *gossip;
  char filename[4096];
  char filenew[4096];
  int status;
  char *mpi_mgi_home = getenv("MGI_MPI_HOME") ;

  if(mpi_mgi_home != NULL) {
    if( mpi_mgi_home[0] == '/' ){
      snprintf(filename,4096,"%s/%s.lock",mpi_mgi_home,service_name);
    }else{
      snprintf(filename,4096,"%s/%s/%s.lock",getenv("HOME"),mpi_mgi_home,service_name);
    }
  }else{
    snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip");
    mkdir(filename,0755);
    snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip/MPI");   // make sure that ~/.gossip/MPI directory exists
    mkdir(filename,0755);
    snprintf(filename,4096,"%s/%s/%s.lock",getenv("HOME"),".gossip/MPI",service_name);
  }
// if(DEBUG & (test ==0)) printf("DEBUG: would be publishing in %s\n",filename);
  if(test == 0) {
    status = open(filename,O_WRONLY+O_CREAT+O_EXCL,00700);   // test mode, try to create lock file
// if(DEBUG) printf("DEBUG %d: status from create '%s' = %d\n",debug_rank,filename,status);
  }else{
    status = unlink(filename);        // cleanup mode, remove lock file
// if(DEBUG) printf("DEBUG %d: status from unlink '%s' = %d\n",debug_rank,filename,status);
  }

  return status;
}

static int MPI_mgi_Publish_name(const char *service_name, MPI_Info info, const char *port_name)
{
  FILE *mpi_gossip;
  char filename[4096];
  char filenew[4096];
  char *mpi_mgi_home = getenv("MGI_MPI_HOME") ;

  if(mpi_mgi_home != NULL) {        // MPI channel files directory
    if( mpi_mgi_home[0] == '/' ){
      snprintf(filename,4096,"%s/%s.new",    mpi_mgi_home,service_name);
      snprintf(filenew, 4096,"%s/%s.channel",mpi_mgi_home,service_name);
    }else{
      snprintf(filename,4096,"%s/%s/%s.new",    getenv("HOME"),mpi_mgi_home,service_name);
      snprintf(filenew, 4096,"%s/%s/%s.channel",getenv("HOME"),mpi_mgi_home,service_name);
    }
if(DEBUG) printf("DEBUG %d: publishing in %s\n",debug_rank,filenew);
  }else{                            //  default directory for MPI channel files
if(DEBUG) printf("DEBUG %d: publishing in %s, MGI_MPI_HOME='%s'\n",".gossip/MPI",debug_rank,mpi_mgi_home);
    snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip");
    mkdir(filename,0755);
    snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip/MPI");
    mkdir(filename,0755);
    snprintf(filename,4096,"%s/%s/%s.new",    getenv("HOME"),".gossip/MPI",service_name);
    snprintf(filenew, 4096,"%s/%s/%s.channel",getenv("HOME"),".gossip/MPI",service_name);
  }
  unlink(filename);
  unlink(filenew);

  mpi_gossip = fopen(filename,"w");
  fprintf(mpi_gossip,"%s",port_name);
  fclose(mpi_gossip);

  link(filename,filenew);
  unlink(filename);
  return MPI_SUCCESS;
}

static int MPI_mgi_Lookup_name(const char *service_name, MPI_Info info, char *port_name)
{
  char filename[4096];
  int fd, nc;
  int wait=0;
  char *mpi_mgi_home = getenv("MGI_MPI_HOME") ;

  if(mpi_mgi_home != NULL) {        // MPI channel files directory
    if( mpi_mgi_home[0] == '/' ){
      snprintf(filename,4096,"%s/%s.channel",mpi_mgi_home,service_name);   // absolute path
    }else{
      snprintf(filename,4096,"%s/%s/%s.channel",getenv("HOME"),mpi_mgi_home,service_name);  // path relative to home
    }
  }else{
    snprintf(filename,4096,"%s/%s/%s.channel",getenv("HOME"),".gossip/MPI",service_name);
  }
if(DEBUG) printf("DEBUG %d: looking for '%s'\n",debug_rank,filename);
  while( (fd=open(filename,0)) < 0) { wait++ ; usleep(1000); }
  nc=read(fd,port_name,MPI_MAX_PORT_NAME);
  close(fd);
  port_name[nc]='\0';
  printf("MPI_Lookup_name: wait time = %d msec\n",wait);

  return MPI_SUCCESS;
}
#if defined(COUPLED)

static char *MGI_MPI_CFG="MGI_MPI_CFG=mgitest 4 : 1024 yin2opa : 256 opa2yin : 512 yan2opa : 784 opa2yan " ;
static char *MGI_MPI_HOME="MGI_MPI_HOME=.gossip/MPI_TEST" ;

int main(int argc, char **argv){
  typedef struct{
    char *name;
    char mode;
  } chdef;
  int size, rank, color, status, i, activeareas;
  int i0, i1, i2, i3;
  MPI_Comm comm_cpl, comm_grid;
  char mode;
  arena *ptr;
  int writbuffer[100];
  int readbuffer[100];
#if defined(YIN)
  chdef ch0 = {"yin2opa" , 'w'};
  chdef ch1 = {"opa2yin" , 'r'};
  chdef ch2 = {"yan2opa" , 'U'};
  chdef ch3 = {"opa2yan" , 'U'};
#endif
#if defined(YAN)
  chdef ch0 = {"yan2opa" , 'w'};
  chdef ch1 = {"opa2yan" , 'r'};
  chdef ch2 = {"yin2opa" , 'U'};
  chdef ch3 = {"opa2yin" , 'U'};
#endif
#if defined(OPA)
  chdef ch0 = {"yin2opa" , 'r'};
  chdef ch1 = {"opa2yin" , 'w'};
  chdef ch2 = {"yan2opa" , 'r'};
  chdef ch3 = {"opa2yan" , 'w'};
#endif

  MPI_Init( &argc, &argv );
  putenv(MGI_MPI_CFG);
  putenv(MGI_MPI_HOME);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(size < 2) exit(1);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); debug_rank = rank;

//   MPI_Comm_split(MPI_COMM_WORLD, color, key, MPI_Comm *newcomm)
  color = (rank<2) ? 1 : 0;     // Coupler PE + PE0
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm_cpl);
  if(color == 0) comm_cpl = MPI_COMM_NULL;
//   comm_cpl = MPI_COMM_WORLD ;

  color = (rank>0) ? 1 : 0;     // grid PEs
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm_grid);
  if(color == 0) comm_grid = MPI_COMM_NULL;

  MPI_mgi_create_begin(rank - 1, MPI_COMM_WORLD, comm_cpl);
  MPI_mgi_create(ch0.name, ch0.mode, rank - 1);
  MPI_mgi_create(ch1.name, ch1.mode, rank - 1);
  MPI_mgi_create(ch2.name, ch2.mode, rank - 1);
  MPI_mgi_create(ch3.name, ch3.mode, rank - 1);
  MPI_mgi_create_end(rank - 1, comm_grid, comm_cpl);
  printf("DEBUG %d:  channel creation done\n",debug_rank);

  MPI_Barrier(MPI_COMM_WORLD); printf("DEBUG %d: Barrier 1\n",debug_rank);
//==============================================================================
  if(rank == 1){   // PE0
    i0 = MPI_mgi_init(ch0.name); printf("DEBUG %d: name = '%s', number = %d\n",debug_rank,ch0.name,i0);
    if(i0 != -1) MPI_mgi_open(i0,&ch0.mode);
    i1 = MPI_mgi_init(ch1.name); printf("DEBUG %d: name = '%s', number = %d\n",debug_rank,ch1.name,i1);
    if(i1 != -1) MPI_mgi_open(i1,&ch1.mode);
    i2 = MPI_mgi_init(ch2.name); printf("DEBUG %d: name = '%s', number = %d\n",debug_rank,ch2.name,i2);
    if(i2 != -1) MPI_mgi_open(i2,&ch2.mode);
    i3 = MPI_mgi_init(ch3.name); printf("DEBUG %d: name = '%s', number = %d\n",debug_rank,ch3.name,i3);
    if(i3 != -1) MPI_mgi_open(i3,&ch3.mode);
  }
//==============================================================================
  MPI_Barrier(MPI_COMM_WORLD); printf("DEBUG %d: Barrier 2\n",debug_rank);

//   printf("DEBUG %d: MAX_CHANNELS = %d\n",debug_rank,MAX_CHANNELS);
  for(i=0 ; i<MAX_CHANNELS ; i++){
//     printf("DEBUG %d: channel %d\n",debug_rank,i);
    mode = mpi_channel_table[i].mode;
    if(mode == 'r' || mode == 'R' || mode == 'w' || mode == 'W' ) {  // only possible on PE0 or Coupler PE
      mpi_channel_table[i].is_active = 1 ;
      if(mode == 'r' || mode == 'R'){
        ptr = mpi_channel_table[i].winbuf;
        ptr->control = CHANNEL_ACTIVE;
        printf("DEBUG %d: read channel %d activated\n",debug_rank,i);
      }
      printf("DEBUG %d:  channel %d '%s' marked as active, mode ='%c'\n",debug_rank,i,mpi_channel_table[i].alias,mode);
    }
  }
//==============================================================================
  if(rank == 0){   // Coupler PE
//     sleep(1);
    printf("DEBUG %d:  Coupler waiiiiiiiiiting\n",debug_rank);
    while(1){
      activeareas = 0;
      for(i=0 ; i<MAX_CHANNELS ; i++){
        mode = mpi_channel_table[i].mode;
        if(mode == 'r' || mode == 'R'){
          ptr = mpi_channel_table[i].winbuf;
          if(ptr->control == CHANNEL_ACTIVE) activeareas++;
          MPI_Win_lock(MPI_LOCK_SHARED,mpi_channel_table[i].thispe,0,mpi_channel_table[i].window);  // active loop for "passive" :-) one sided
          MPI_Win_unlock(mpi_channel_table[i].thispe,mpi_channel_table[i].window);
        }
      }
      printf("DEBUG %d:  Coupler waiting on %d channel(s)\n",debug_rank,activeareas);
      if(activeareas == 0) break;
      sleep(1);
    }
  }
//==============================================================================
  if(rank == 1){   // PE0
    for(i=0 ; i<100 ; i++) writbuffer[i] = i;
    for(i=0 ; i<100 ; i++) readbuffer[i] = 9999;
    printf("DEBUG %d:  PE0 sleeeeeeeeping\n",debug_rank);
    for(i=0 ; i<MAX_CHANNELS ; i++){
      mode = mpi_channel_table[i].mode;
      if(mode == 'w' || mode == 'W'){
        printf("DEBUG %d: writing into channel %s\n",debug_rank,mpi_channel_table[i].alias);
      }
    }
    for(i=0 ; i<MAX_CHANNELS ; i++){
      mode = mpi_channel_table[i].mode;
      if(mode == 'r' || mode == 'R'){
        printf("DEBUG %d: reading from channel %s\n",debug_rank,mpi_channel_table[i].alias);
      }
    }
    sleep(10);
    for(i=0 ; i<MAX_CHANNELS ; i++){
      mode = mpi_channel_table[i].mode;
      if(mode == 'r' || mode == 'R'){
        ptr = mpi_channel_table[i].winbuf;
        ptr->control = CHANNEL_STOPPED;
        printf("DEBUG %d: read channel %d stopped\n",debug_rank,i);
      }
    }
  }
//==============================================================================
  MPI_Finalize();
  return 0;
}
#endif
#if defined(SELF_TEST)
int main(int argc, char **argv){
  arena test_arena;
  arena *arenaptr = &test_arena;
  int size, rank, localrank, status;
  int data_offset =  &(arenaptr->data[0]) - &(arenaptr->control);
  MPI_Comm global, local ;
  int test_data = 999999;
  MPI_Status mpistat;
  int channel1, channel2;
  int tbuf[100], tbuf2[100];
  int i;
  arena *zz;

  for(i=0 ; i<10 ; i++) tbuf[i] = 1000+i;
  for(i=0 ; i<10 ; i++) tbuf2[i] = 100+i;

  printf("data offset = %d\n",data_offset);
  MPI_Init( &argc, &argv );
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  debug_rank = rank;

  channel1 = MPI_mgi_create("aname");
if(DEBUG) printf("DEBUG %d: create aname = %d\n",rank,channel1);
  channel2 = MPI_mgi_create("bname");
if(DEBUG) printf("DEBUG %d: create bname = %d\n",rank,channel2);
if(DEBUG) printf("DEBUG %d: server = %d\n",rank,mpi_channel_table[0].is_server);

  channel1 = MPI_mgi_init("aname");
  channel2 = MPI_mgi_init("bname");
if(DEBUG) printf("DEBUG %d: channel1 = %d, channel2 = %d\n",rank,channel1,channel2);
  status = MPI_mgi_open(channel1,"R") ;
if(DEBUG) printf("DEBUG %d: aname status = %d\n",rank,status);
  status = MPI_mgi_open(channel2,"W") ;    // deliberate error
if(DEBUG) printf("DEBUG %d: bname status = %d\n",rank,status);

// cheating a bit because open is simulated
//   if(mpi_channel_table[0].is_server){
//     MPI_Comm_accept(mpi_channel_table[0].port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &global );  // wait for client to connect (collective over communicator)
// if(DEBUG) printf("DEBUG %d: accepted connection\n",rank);
//     MPI_Intercomm_merge(global, 0, &local);                                  // create communicator for one-sided window
// if(DEBUG) printf("DEBUG %d: merged connection\n",rank);
//   }else{
// if(DEBUG) printf("DEBUG %d: lookup '%s'\n",rank,mpi_channel_table[0].channel_name);
//     MPI_mgi_Lookup_name(mpi_channel_table[0].channel_name, MPI_INFO_NULL, mpi_channel_table[0].port_name);  // get port name published under name channel_name
// if(DEBUG) printf("DEBUG %d: server is '%s'\n",rank,mpi_channel_table[0].port_name);
//     MPI_Comm_connect(mpi_channel_table[0].port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &global ); // wait until connected to server (collective over communicator)
// if(DEBUG) printf("DEBUG %d: connected to server ;%s'\n",rank,mpi_channel_table[0].port_name);
//     MPI_Intercomm_merge(global, 1, &local);                                  // create communicator for one-sided window
// if(DEBUG) printf("DEBUG %d: merged to server ;%s'\n",rank,mpi_channel_table[0].port_name);
//   }
  local = mpi_channel_table[0].local ;
  global = mpi_channel_table[0].global ;
  MPI_Comm_rank(local, &localrank);
//   mpi_channel_table[0].is_active = 1 ; 
//   mpi_channel_table[0].local = local;
//   mpi_channel_table[0].global = global;
// end of cheating, open simulated
  
  printf("rank %d of %d\n",rank+1,size);

//   sleep(2);
//   MPI_Barrier(local);
  MPI_Sendrecv(&localrank, 1, MPI_INTEGER, 1-localrank, 123, 
	       &test_data, 1, MPI_INTEGER, 1-localrank, 123,
	       local, &mpistat);
  printf("localrank %d received %d\n",localrank,test_data);
//   MPI_Barrier(local);

  if(localrank == 0) {

printf("DEBUG %d: writing\n",rank);
    status = MPI_mgi_write(channel1,(void *)tbuf, "I", 10);
printf("DEBUG %d: write status = %d\n",rank,status);
    status = MPI_mgi_write(channel1,(void *)tbuf, "R", 5);
printf("DEBUG %d: write status = %d\n",rank,status);
//   MPI_Barrier(local);

  }else{

//   MPI_Barrier(local);
    sleep(5);
printf("DEBUG %d: reading\n",rank);
    status = MPI_mgi_read(channel1,(void *)tbuf2, "I", 10);
printf("DEBUG %d: read status = %d\n",rank,status);
printf("DEBUG : data read");for(i=0 ; i<10 ; i++)printf("%5d ",tbuf2[i]); printf("\n");

    status = MPI_mgi_read(channel1,(void *)tbuf2, "R", 5);
printf("DEBUG %d: read status = %d\n",rank,status);
printf("DEBUG : data read");for(i=0 ; i<5 ; i++)printf("%5d ",tbuf2[i]); printf("\n");

  }

  MPI_Barrier(local);
  zz = (arena *)mpi_channel_table[0].winbuf;
  printf("DEBUG :"); printf(" %c %d ",zz->data[0]&0xFF,zz->data[0]>>8); for(i=1 ; i<20 ; i++) printf("%5d ",zz->data[i]); printf("\n");
  printf("DEBUG : in = %d, out = %d\n",zz->in,zz->out);
  status = MPI_mgi_close(channel1) ;
printf("DEBUG %d: close channel 1 status = %d\n",rank,status);
  status = MPI_mgi_close(channel2) ;  // deliberate error
printf("DEBUG %d: close channel 2 status = %d\n",rank,status);
  MPI_mgi_Unpublish_name("mastername_0", MPI_INFO_NULL, "");
  MPI_Finalize();
  return 0;
}
#endif
