/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 * Copyright (C) 2014       ESCER center, UQAM
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* MGILIB.C

   Functions written for programs that run in "parallel" and need
   to exchange information simultaneously. MGI stands for Model Gossip
   Interface. Each program using these functions will have to compile
   with this library. There are 6 functions:

   MGI_INIT
   to initialize a channel using the name given and open the socket.
   It will also allocate a writing buffer dynamically and if all is
   successful, it will return a channel number (socket descriptor).

   MGI_OPEN
   to open the channel in a certain mode:
   'R' for read: returns 0 to signal that data been written.
                 or returns nblks to be read
   'W' for write: returns 1 if open is ok.
   'S' for storing a restart file:returns 1 if open is ok.

   MGI_READ
   to read from a channel that is open for READ mode.
   It accepts the following type of data:
   'C' for CHARACTER
   'I' for INTEGER
   'R' for REAL
   'D' for REAL*8
   It returns the number of blocks left to read from channel.

   MGI_WRITE
   to write to a channel that is open for WRITE mode.
   It accepts the same type of data as MGI_READ.
   It returns the number of blocks written to channel.

   MGI_CLOS
   to close the mode of a channel and check to make sure all is
   transmitted as requested. It returns the status of the data
   file after it is closed.

   MGI_TERM
   to delete the PID file that was created in the beginning and
   to release all the memory allocated dynamically. It closes all
   the filepipes therefore, breaking all the pipe connections with
   the other programs.

   ***NOTE: These functions are written to keep enough bits for the
   equivalent of an integer/float in C or integer/real*4 in FORTRAN.
   In other words, you will lose some precision with the 64-bit
   compilation unless real*8 is used.

   Revision: March/April 2014, M.Valin, UQAM

   A new capability has been added without altering the API for
   MGI_INIT / MGI_OPEN / MGI_WRITE / MGI_READ / MGI_CLOS / MGI_TERM
   shared memory (shmget / shmat ...) can now be used as a communication
   medium when both ends of the "channel" reside on the same host
   the "shared memory channel" related files will be found inside the
   $HOME/.gossip/SHM directory ($HOME/.gossip for the TCP/IP mode).
   When using the shared memory mode, the gossip library goes unused
   and could even be stubbed.
   Contrary to TCP/IP channels, there is NO LOSS OF PRECISION when sending
   real*8 / double data.
   Some code simplification has been done.
   Error messages have been slightly reworked.

   new entry point:
   -PrintMgiError (C)
   -print_mgi_error(Fortran)
           print the text associated with an error code on the standard error
           file descriptor

   new programs:
   -mgi_shmem (C) is used to create the shared memory area associated with
           a "shared memory channel".
    mgi_shmem (C) will also monitor "shared memory channel" activity if needed.

   -mgi_test (C) is used to test a "shared memory channel" using the C API
           and a scenario file.

   -f_mgi_test (Fortran) is used to test a "shared memory channel" using the
           Fortran API and a scenario file.

   Revision: September 2015, M.Valin, UQAM
   If the code is compiled with -DWITHOUT_GOSSIP the MGI library will only support
   the shared memory mode.
   A bit of code refactoring has been done to make it easier to follow.
   Fortran callable mgi_read_c function has been added because there are two(2) character
   strings instead of one(1) in the interface.

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <rpnmacros.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <string.h>
#include <ctype.h>
#include "mgi.h"

#define MEMCPY_LIMIT 10

#if ! defined(WITHOUT_GOSSIP)
#include <gossip.h>
#endif

/* error codes header file */
#include <cgossip.h>

static int perfmon = 0;  /* if 1 , collect performance stats */
static int printed = 1;
static long long time_in;

static long long bytes_read = 0;       /* bytes read */
static long long time_read = 0;        /* time spent reading */
static long long call_read = 0;        /* number of calls to mgi_read */
static long long wait_read = 0;        /* time spent waiting on a read (approximate) */

static long long bytes_write = 0;    /* bytes written */
static long long time_write = 0;     /* time spent writing */
static  long long call_write = 0;    /* number of calls to mgi_write */
static long long wait_write = 0;     /* time spent waiting on a write (approximate) */

static channel chn[MAX_CHANNELS];

static int ichan = 0;
static int init = 0;
static int SIG_ACTIVE = 1;

#ifdef NOT_USED
static char PID_file[MAX_STR];
static char *mgidir;
#endif

static int *intBuffer;

#ifdef NOT_USED
static void getmgidir ();
static int makepidfile ();
static void removepidfile ();
static int validchan (int chan);
static void strcopy (char *s, char *t, int charlen);
#endif

static char *mgierrors[]={
  "INIT_ERROR",
  "SERVER_ERROR",
  "CONNECTION_ERROR",
  "READ_ERROR",
  "WRITE_ERROR",
  "READ_TIMEOUT",
  "WRITE_TIMEOUT",
  "READ_TYPE_ERROR",
  "WRITE_TYPE_ERROR",
  "DATA_LENGTH_ERROR",
  "SEND_COMMAND_ERROR"
};

void PrintMgiError(int code)
{
  int index = -1-code;
  if(index>=0 && index <=10) {
    fprintf(stderr,"MGI ERROR: %d %s\n",code,mgierrors[index]);
  }
}
#pragma weak print_mgi_error__=print_mgi_error
#pragma weak print_mgi_error_=print_mgi_error
void print_mgi_error__(int *code);
void print_mgi_error_(int *code);
void print_mgi_error(int *code)
{
  PrintMgiError(*code);
}

ftnword f77_name (mgi_init) (char *channel_name, F2Cl lname);
ftnword f77_name (mgi_open) (ftnword *f_chan, char *mode, F2Cl lmode);
ftnword f77_name (mgi_read) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
ftnword f77_name (mgi_write) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
ftnword f77_name (mgi_clos) (ftnword *f_chan);
ftnword f77_name (mgi_term) ();
void f77_name (mgi_set_timeout) (ftnword *chan, ftnword *timeout);

#if ! defined(WITHOUT_GOSSIP)

extern int connect_to_subchannel_by_name (char *channel, char *subchannel, char *mode);
extern int GET_ack_nack (int socket, char *message);
extern int write_record (int fclient, void *buf, int longueur, int tokensize);
extern void *read_record (int fclient, void *buf, int *longueur, int maxlongueur, int tokensize);
extern char *get_gossip_dir (int display);

extern void init_client_table ();
extern void set_client_timeout (int fclient, int timeout);
extern int get_client_timeout (int fclient);
extern int close_channel (int fclient, char *channel);

extern int get_stream_timeout(int gchannel);
extern int send_command(char *buf);
extern int retry_connect(int chan);
extern int get_timeout_signal(int gchannel);
extern int signal_timeout(int gchannel);

//ftnword f77_name (mgi_read_oob) ();   /* read out of band */
//ftnword f77_name (mgi_write_oob) ();   /* write out of band */
#endif

#pragma weak mgi_perf_print_=mgi_perf_print
#pragma weak mgi_perf_print__=mgi_perf_print
void mgi_perf_print_();
void mgi_perf_print__();
void mgi_perf_print(){
  double tim, speed, temp;
  double tim2, speed2;

  if(printed) return;
  tim = 0.0;
  speed = 0.0;
  fprintf(stderr,"\n====== MGI performance report ======\n");

  tim = time_read * .000001;       /* in seconds */
  tim2 = tim - wait_read * .000001;  /* subtract time spent waiting */
  temp = bytes_read * .000001;     /* in MegaBytes */
  speed = (tim <= 0.0) ? 0.0 : temp / tim ;
  speed2 = (tim2 <= 0.0) ? 0.0 : temp / tim2 ;
  fprintf(stderr,"%10Ld bytes read,    %8Ld calls to mgi_read,  %10.6f(%10.6f) seconds : %10.6f(%10.6f) MBytes/sec\n",
          bytes_read,call_read,tim,tim2,speed,speed2);

  tim = time_write * .000001;       /* in seconds */
  tim2 = tim - wait_write * .000001;  /* subtract time spent waiting */
  temp = bytes_write * .000001;     /* in MegaBytes */
  speed = (tim <= 0.0) ? 0.0 : temp / tim ;
  speed2 = (tim2 <= 0.0) ? 0.0 : temp / tim2 ;
  fprintf(stderr,"%10Ld bytes written, %8Ld calls to mgi_write, %10.6f(%10.6f) seconds : %10.6f(%10.6f) MBytes/sec\n",
          bytes_write,call_write,tim,tim2,speed,speed2);
  fprintf(stderr,"\n====================================\n");
  printed = 1;
}

#pragma weak mgi_perf_on_=mgi_perf_on
#pragma weak mgi_perf_on__=mgi_perf_on
void mgi_perf_on_();
void mgi_perf_on__();
void mgi_perf_on(){
  if(perfmon == 0) atexit(mgi_perf_print);
  perfmon = 1;
  printed = 0;
}

#pragma weak mgi_perf_off_=mgi_perf_off
#pragma weak mgi_perf_off__=mgi_perf_off
void mgi_perf_off_();
void mgi_perf_off__();
void mgi_perf_off(){
  perfmon = 0;
}

static long long mgi_time(){
  struct timeval tv;
  long long temp;

  gettimeofday(&tv,NULL);
  temp = tv.tv_sec;
  temp *= 1000000;
  temp += tv.tv_usec;

  return(temp);
}
/***********************************************************************************************/

static int validchan( int chan )
/*
 *   validate the channel number
 *   must be greater than zero and less than or equal to ICHAN
 *   channel must be active (TCP/IP or shared memory)
 */

{
  if(chan < 0 || chan >= MAX_CHANNELS) return (CONNECTION_ERROR);   /* a priori invalid channel number */
  if ( chn[chan].buffer != NULL )      return ( 0);   /* TCP/IP channel */
  if ( chn[chan].shmbuf != NULL )      return ( 0);   /* shared memory channel */

  return(-1);
}
/***********************************************************************************************/

void f77_name (mgi_set_timeout) (ftnword *chan, ftnword *timeout)
{
  int channel = *chan;
  if(validchan(channel) != 0) return; /* invalid channel */

#if ! defined(WITHOUT_GOSSIP)
  set_client_timeout(chn[channel].gchannel, (int) *timeout);
#endif
  chn[channel].timeout = *timeout;  /* timeout value useful locally for shared memory communications */
}

static void strcopy( char *s, char *t, F2Cl charlen )
/* to copy a string given by a fortran routine and place the NULL
 *        character at the end of the true (charlen) length of the string */
{
  int i;

  i = 0;
  while ( (*s++ = *t++) != ' ' && i++ < charlen);
  if (*s-- == ' ') *s = '\0';
  else *s++ = '\0';
}
/* --------------------------------------------------------------------------- */
/* #define DEBUG */

/*********************************************************************************************/
#ifdef NOTUSED

void f77_name (mgi_nosig) ()
     /* to disable the signals between filepipes */
{
  /* SIG_ACTIVE = 0; */
  fprintf(stderr,"MGI_NOSIG: deprecated call\n");
}
/* to copy a string given by a fortran routine,
   the space character is ignored, it is taken
   here as the end of the string */
static void strcopy_( char *s1, char *s2, int lengths1, int lengths2 )
{
  int i = 0;

  while ( (*s1++ = *s2++) != ' ' && i++ < lengths2 );

}

int check_ends(unsigned char *s1,unsigned char *s2, int s1length, int s2length, int i )
{
  if(*s2 == ' ' )
  {
    return (*(s1 - 1) - *(s2 - 1));
  }
  else
  {
    if(i == s2length)
      return (*(s1 - 1) - *(s2 - 1));
    else
      return -1;
  }
}

/* to compare two strings given by a fortran routine, considering
   the space character as the end of the string*/
static int f_strcmp( unsigned char *s1, unsigned char *s2, int s1length, int s2length )
{
  int i = 0;
  int length;

  if(s1length <= s2length)
    length = s1length;
  else
    length = s2length;

  while( i < length && *s1 != ' ' && *s2 != ' ' && *s1 == *s2 && *s1 != '\0' && *s2 != '\0')
    {
      s1++;
      s2++;
      i++;
    }


  if(*s1 == ' ')
    {
      fprintf(stderr, "mgilib2::f_strcmp(), before return if(*s1 == ' '), s1 => %s\n ", s1);
      return check_ends(s1, s2, s1length, s2length, i);

    }
  else if(*s2 == ' ')
    {
      fprintf(stderr, "mgilib2::f_strcmp(), before return if(*s2 == ' '), s2 => %s\n ", s2);
      return 2 + check_ends(s2, s1, s2length, s1length, i);
    }

  return (*s1 - *s2);

}

static void getmgidir()
     /* to get the value of the environment variable "MGI_DIR" */
{
  if ( (mgidir = getenv("MGI_DIR")) == NULL)
    {
      fprintf(stderr,"Environment variable 'MGI_DIR' undefined --\n");
      /* exit(1); */
    }
}

static int makepidfile()
     /* to make the PID file */
{
  char stuff[MAX_STR];
  sprintf( PID_file, "%s/%d", mgidir, getpid() );
  sprintf( stuff, "%s/%s", mgidir, "PROCS" );
  fprintf(stderr, "linking :%s: to :%s:\n", PID_file, stuff );
  return ( link ( stuff, PID_file ) );
}

static void removepidfile()
     /* to remove the PID file */
{
  fprintf(stderr, "removing %s\n", PID_file );
  unlink( PID_file );
}
#endif

#if ! defined(WITHOUT_GOSSIP)
/********************************************************************************************
 * nb = bwrite(chan,buffer,nelem,dtype) : send data to channel "chan"
 * chan   : channel number
 * buffer : data to send
 * nelem  : number of data pieces
 * dtype  : data type 'C', 'I', 'R', 'D'
 *
 * return value : number of bytes not sent,   > 0 : error
 ********************************************************************************************/
static int bwrite ( int chan, void *buffer, int nelem, char *dtype )
{
  int nb, ier;

  fd_set wfds, rfds;
  struct timeval tv;

  FD_ZERO(&wfds);
  FD_SET(chn[chan].gchannel, &wfds);

  tv.tv_sec = get_stream_timeout(chn[chan].gchannel);
  tv.tv_usec = 0;

  if (select(chn[chan].gchannel + 1, NULL, &wfds, NULL, &tv))
    {
      ier = send_command_to_server(chn[chan].gchannel, "WRITE");

    }
  else
    {
      /* return ier = signal_timeout(chn[chan].gchannel); */
      return WRITE_TIMEOUT;
    }


  if(ier < 0 )
    {
      fprintf(stderr,"ERROR: (bwrite) unable to send write command\n");
      /* return -1; */
      return WRITE_ERROR;
    }

#ifdef DEBUG
  fprintf(stderr,"mgilib2::bwrite(), ==\n");
#endif
  nb = 0 ;
  if(*dtype == 'I' || *dtype == 'R')
    {
      nb = write_record(chn[chan].gchannel, (unsigned char *)buffer, nelem, sizeof(int));
    }

  else if(*dtype == 'D')
    {
      nb = write_record(chn[chan].gchannel, (char *)buffer, nelem, sizeof(double));
    }

  else if(*dtype == 'C')
    {
      nb = write_record(chn[chan].gchannel, buffer, nelem, 1);
    }

  /* get_ack_nack(chn[chan].gchannel); */

  FD_ZERO(&rfds);
  FD_SET(chn[chan].gchannel, &rfds);

  tv.tv_sec = get_stream_timeout(chn[chan].gchannel);
  tv.tv_usec = 0;

  if (select(chn[chan].gchannel + 1, &rfds, NULL, NULL, &tv))
    {
      get_ack_nack(chn[chan].gchannel);

    }
  else
    {
      fprintf(stderr, "ERROR: (bwrite) timeout = %d, get_ack_nack(), else\n", (int) tv.tv_sec);
      /* return ier = signal_timeout(chn[chan].gchannel); */
      return  WRITE_TIMEOUT;
    }


  return nb;
}
#endif
/********************************************************************************************
 * close a channel and signal that it can be opened in another mode
 *
 * status = mgi_clos(channel_number)
 * status == 0 : no error
 ********************************************************************************************/
ftnword f77_name (mgi_clos) (ftnword *f_chan)
{
  int ier = 0, chan;
  char buf[1024];
  chan = (int) *f_chan;
  int force=0;

  if(chan > 1000){
    chan -=1000;
    force = 1;
  }

  if(validchan(chan) != 0) return(CONNECTION_ERROR); /* invalid channel */

  if(chn[chan].shmbuf != NULL) {  /* shared memory channel  */

    if(chn[chan].mode == 'W'){   /* channel is open for write */
      if(chn[chan].shmbuf->write_status == MGI_SHM_ACTIVE) {
        chn[chan].shmbuf->write_status = MGI_SHM_IDLE;    /* mark as initialized but not open for write */
      }else{
        fprintf(stderr,"ERROR: (mgi_clos) inconsistent write status on shm channel '%s'\n",chn[chan].name);
        if(force) chn[chan].shmbuf->write_status = MGI_SHM_IDLE;
      }
    }

    if(chn[chan].mode == 'R'){   /* channel is open for read */
      if(chn[chan].shmbuf->read_status  == MGI_SHM_ACTIVE) {
        chn[chan].shmbuf->read_status  = MGI_SHM_IDLE;    /* mark as initialized but not open for read  */
      }else{
        fprintf(stderr,"ERROR: (mgi_clos) inconsistent read status on shm channel '%s'\n",chn[chan].name);
        if(force) chn[chan].shmbuf->read_status = MGI_SHM_IDLE;
      }
    }
    fprintf(stderr,"INFO: (mgi_clos) shared memory channel '%s' is closed \n", chn[chan].name);
    chn[chan].mode = ' ';      /* mark channel as no longer open */
    return 0;
  }  /* shared memory channel  */

#if defined(WITHOUT_GOSSIP)
  return(CONNECTION_ERROR); /* shared memory channel or else !! */
#else
  if(chn[chan].gchannel != 0)   /* TCP/IP channel */
    {
      snprintf(buf, 1023, "%s %s", "END", chn[chan].name);
      ier = send_command(buf);
      fprintf(stderr,"INFO: (mgi_clos) subchannel '%s' is closed \n", chn[chan].name);
    }

   if(chn[chan].buffer)
    {
      free(chn[chan].buffer);
      chn[chan].buffer = NULL;
    }
  return ier;
#endif
}
int C_mgi_clos (int c_chan){
  ftnword f_chan=c_chan;
  return (f77_name (mgi_clos) (&f_chan)) ;
}
/********************************************************************************************
 * close all channels, detach from shared memory if shared memory channel
 * return value: 0 if no error, != 0 if error
 ********************************************************************************************/
ftnword f77_name (mgi_term) ()
{
  int chan, ier = -1;

  for (chan = 0; chan <= ichan; chan++)
    {
      if(chn[chan].shmbuf != NULL) {  /* shared memory channel  */
        /* get number of attaches, if == 1 i am last and will dump rest of buffer into resatart file */
        if(chn[chan].shmbuf->write_status == MGI_SHM_IDLE && chn[chan].shmbuf->read_status == MGI_SHM_IDLE) {  /* both close operations done */
          if(chn[chan].shmbuf->in != chn[chan].shmbuf->out){                               /* buffer is not empty */
            fprintf(stderr,"WARNING: (mgi_term) shared memory channel '%s' is not empty\n", chn[chan].name);
            /* open channel file for exclusive writing and dump buffer */
          }
        }
        fprintf(stderr,"INFO: (mgi_term) shared memory channel '%s' fully closed\n", chn[chan].name);
        ier = 0;
        ier = shmdt(chn[chan].shmbuf) ; /* detach from shared memory */
        chn[chan].shmbuf = NULL;        /* make it obvious that process has detached */
        continue;
      }

#if defined(WITHOUT_GOSSIP)
      continue;
#else
      if(chn[chan].name && strcmp((char *)chn[chan].name, "") && chn[chan].gchannel > 0)    /* TCP/IP channel */
	{
	  ier = send_command("END");
          fprintf(stderr,"INFO: (mgi_term) subchannel '%s' fully closed\n", chn[chan].name);

	  if(chn[chan].buffer)
	    {
	      free(chn[chan].buffer);
	      chn[chan].buffer = NULL;
	    }
	}
#endif
    }

  return ier;
}
int C_mgi_term(){
  return( f77_name (mgi_term) () );
}
/********************************************************************************************
 * channel = mgi_init(channel_name)
 * channel_name : Fortran string
 * channel : channel number ( < 0 if error)
 *
 ********************************************************************************************/
ftnword f77_name (mgi_init) (char *channel_name, F2Cl lname)
     /* To initialize a channel given a channel_name.
	It will return a number to represent this channel (1 to MAX_CHANNELS-1 */
{
  int chan, items;
  char env_var_name[1024];
  char shm_fil_name[1024];
  char *env_var_value;
  FILE *FD;

  if (init == 0)
    {
      init = 1;
    }

#ifdef DEBUG
  fprintf(stderr,"MGI_INIT ** \n");
#endif
  ichan++;
  if (ichan >= MAX_CHANNELS)
    {
      fprintf(stderr,"ERROR: (mgi_init) Too many channels assigned; MAX = %d\n", MAX_CHANNELS);
      ichan--;
      /* return -1; */
      return INIT_ERROR;
    }
  else
    {
      chan = ichan;
      if (lname < MAX_NAME)
      {
	strcopy(chn[chan].name, channel_name, lname);
      }
      else
      {
        fprintf(stderr,"ERROR: (mgi_init) Length of channel name > %d chars.\n",MAX_NAME-1);
        return INIT_ERROR;   /* return -1; */
      }
      chn[chan].fd_data = -1;

      /* initialize channel */
      chn[chan].msgno_W = 0;
      chn[chan].msgno_R = 0;
      chn[chan].nblks = 0;
      chn[chan].mode = ' ';
      chn[chan].pos = 0;
      chn[chan].gchannel = 0;
      chn[chan].shmid = -1;   /* set shared memory id to unused */
      chn[chan].shmbuf = NULL;
      chn[chan].buffer = NULL;
      chn[chan].timeout = 1000 ; /* default timeout is 1000 seconds */

      snprintf(env_var_name,sizeof(env_var_name)-1,"SHM_%s",chn[chan].name);
      env_var_name[sizeof(env_var_name)-1]='\0';
      env_var_value = getenv(env_var_name);
      if(env_var_value != NULL) {                /* it is a shared memory channel */
        chn[chan].shmid = atoi(env_var_value);   /* get shared memory segment ID */
      }else{  /* look for a channel id file */
        snprintf(shm_fil_name,sizeof(shm_fil_name),"%s/.gossip/SHM/%s.id",getenv("HOME"),chn[chan].name);
        FD=fopen(shm_fil_name,"r");
        if(FD != NULL){
          items = fscanf(FD,"%d",&chn[chan].shmid);  /* get shared memory segment ID */
          if(items != 1) {
            fprintf(stderr,"ERROR: (mgi_init) invalid shared memory channel id \n");
            return INIT_ERROR;
          }
          fclose(FD);
        }
      }
      if(chn[chan].shmid != -1) {                /* it is a shared memory channel */
        chn[chan].shmbuf = shmat(chn[chan].shmid, NULL, 0);
        if(chn[chan].shmbuf == (void *) -1) { /* cannot attach shared memory segment */
          fprintf(stderr,"ERROR: (mgi_init)  channel '%s': Cannot attach shared memory segment %d\n",chn[chan].name,chn[chan].shmid);
          return INIT_ERROR;
        }
        /* chn[chan].shmbuf->limit must not be 0 if memory segment was properly initialized at creation  */
      }
#if defined(WITHOUT_GOSSIP)
      else {
        fprintf(stderr,"ERROR: (mgi_init) cannot find channel %s\n",chn[chan].name);
        return INIT_ERROR;
      }
#endif

      if (SIG_ACTIVE)
      {
        fprintf(stderr,"INFO: (mgi_init) Initializing %s channel: '%s' \n", (chn[chan].shmbuf == NULL) ? "TCP/IP" : "shared memory" , chn[chan].name);
      }
      if ((intBuffer = (int *) malloc(BUFSIZE * sizeof(int))) == NULL)
      {
        fprintf(stderr,"ERROR: (mgi_init)  channel '%s': Cannot allocate memory for intBuffer\n",chn[chan].name);
        return INIT_ERROR;
      }

    chn[chan].buffer = intBuffer;
    }

  return(chan);
}
int C_mgi_init (char *channel_name){
  F2Cl lname = strlen(channel_name);
  return (f77_name (mgi_init) (channel_name, lname)) ;
}
/********************************************************************************************
 * open a previously initialized channel (mgi_init) for read/write/store
 * status = mgi_open(channel,mode)
 * channel : channel number obtained from mgi_init
 * mode    : 'R', 'W', 'S'
 *
 * status  : channel number if valid channel number, error code otherwise
 ********************************************************************************************/
ftnword f77_name (mgi_open) (ftnword *f_chan, char *mode, F2Cl lmode)
     /* to open a channel in mode "mode"; where mode can be:
	'R' for reading
	'W' for writing
	'S' for storing
     */
{
  int chan;
  chan = (int) *f_chan;
  mgi_shm_buf *shm;
  int force=0;

  if(chan > 1000){
    chan -=1000;
    force = 1;
  }
  if(validchan(chan) != 0) return(CONNECTION_ERROR); /* invalid channel */
  if (*mode == 'W')
    {
      if(chn[chan].shmbuf == NULL) {   /* not a shared memory channel */
#if defined(WITHOUT_GOSSIP)
        return(CONNECTION_ERROR); /* invalid channel */
#else
        chn[chan].gchannel = connect_to_subchannel_by_name( get_gossip_dir(0), chn[chan].name, "write" );

        if( chn[chan].gchannel < 0 )
          chn[chan].gchannel = retry_connect( chan );
#endif
      } else {   /* it is  a shared memory channel, initialize it for write  */
        shm = chn[chan].shmbuf;
        if(shm->write_status != MGI_SHM_IDLE && force == 0) return CONNECTION_ERROR;  /* not initialized for write or already connected for write */
        shm->write_status = MGI_SHM_ACTIVE;  /* mark as connected for write */
        chn[chan].gchannel = 100000+chan ;   /* fake fd */
        chn[chan].mode = 'W';
#if defined(WITHOUT_GOSSIP)
        chn[chan].timeout =  1000000 ; /* default 1000 second timeout */
#else
        chn[chan].timeout =  get_client_timeout(chn[chan].gchannel);     /* get timeout value for this channel */
#endif
      }

    }
  else if (*mode == 'R')
    {
      if(chn[chan].shmbuf == NULL) {   /* not a shared memory channel */
#if defined(WITHOUT_GOSSIP)
        return(CONNECTION_ERROR); /* invalid channel */
#else
        chn[chan].gchannel = connect_to_subchannel_by_name( get_gossip_dir(0), chn[chan].name, "read" );

        if( chn[chan].gchannel < 0 )
          chn[chan].gchannel = retry_connect( chan );
#endif
      } else {   /* it is  a shared memory channel, initialize it for read */
        shm = chn[chan].shmbuf;
        if(shm->read_status != MGI_SHM_IDLE && force == 0) return CONNECTION_ERROR;  /* not initialized for write or already connected for write */
        shm->read_status = MGI_SHM_ACTIVE;  /* mark as connected for read */
        chn[chan].gchannel = 100000+chan ;   /* fake fd */
        chn[chan].mode = 'R';
#if defined(WITHOUT_GOSSIP)
        chn[chan].timeout =  1000000 ; /* default 1000 second timeout */
#else
        chn[chan].timeout =  get_client_timeout(chn[chan].gchannel);     /* get timeout value for this channel */
#endif
      }
    }
  else if (*mode == 'S')
    { /* store mode (for restart files) */
      chn[chan].mode = 'S';
      chn[chan].nblks = 0;
      chn[chan].msgno_W++;
      chn[chan].pos = 0;
    }

  if(chn[chan].gchannel < 0)
    {
      fprintf(stderr, "ERROR: (mgi_open) Connection Failed, the Server may be down \n" );
      /* exit(-1); */
      return CONNECTION_ERROR;
    }

#if ! defined(WITHOUT_GOSSIP)
  /* initialize timeout table */
  init_client_table( chn[chan].gchannel );
#endif
  if (SIG_ACTIVE)
  {
    fprintf(stderr,"INFO: (mgi_open) Opening channel: '%s' , mode = '%c'\n", chn[chan].name, *mode) ;
  }
  return chan;
}
int C_mgi_open (int c_chan, char c_mode){
  F2Cl lmode=1;
  ftnword f_chan=c_chan;
  char mode=c_mode;
  return( f77_name (mgi_open) (&f_chan, &mode, lmode) );
}

#if ! defined(WITHOUT_GOSSIP)
/********************************************************************************************
 * connection helpers
 * call mgi_set_retry_connect(number_of_retries)
 ********************************************************************************************/
/* if connection to server fails          */
/* default:  retry 10 times after a sleep */
/* else use user value                    */

int USER_TRY_CONNECT = 10;
void f77_name (mgi_set_retry_connect) (ftnword *try_nbr)
{
  printf( "INFO: (mgi_open) setting try to connect USER_TRY_CONNECT: '%d' times\n", (int) *try_nbr );
  if((int) *try_nbr > 0 && (int) *try_nbr < 10)
    USER_TRY_CONNECT = (int) *try_nbr;
}

int mgi_get_retry_connect(int chan)
{

  return USER_TRY_CONNECT;

}

/* if connection to server fails          */
/* default: retry 10 times after a sleep  */
/* interval of 10 secs                    */
int retry_connect( int chan )
{
  int PING_INTERVAL = 10;
  int ping_ord0 = mgi_get_retry_connect(chan);
  int ping_ord =  mgi_get_retry_connect(chan);

  while( chn[chan].gchannel < 0 && ping_ord > 0 )
    {
      sleep( PING_INTERVAL );
      fprintf(stderr, "INFO: (mgi_open) Connection to Server Failed,  retry: %d/%d \n", ping_ord0 - ping_ord + 1, ping_ord0 );
      chn[chan].gchannel = connect_to_subchannel_by_name( get_gossip_dir(0), chn[chan].name, "write" );
      ping_ord--;
    }
  return chn[chan].gchannel;

}
#endif
/********************************************************************************************
 * write into shared memory buffer, character version
 ********************************************************************************************/
static int shm_write_c(mgi_shm_buf *shm,void *buf,int nelem,int timeout){
  volatile int in, out, inplus;
  int limit;
  int ntok=nelem;
  unsigned char *str = (unsigned char *) buf;
  unsigned int token;
  int maxiter = timeout*1000 ; /* 1000 iterations is one second */
  int iter;

//fprintf(stderr,"DEBUG: shm_write_c, writing %d characters\n",nelem);
  in = shm->in ; out = shm->out ; limit = shm->limit;
  bytes_write += ntok;
  while(ntok > 0){
    inplus = (in+1 > limit) ? 0 : in+1 ;
    iter = maxiter;
    while(inplus == out && iter > 0) {       /* shared memory circular buffer is full */
      shm->in = in;              /* update in pointer in shared memory */
      usleep(1000);              /* sleep for 1 millisecond */
      wait_write += 1000;
      out = shm->out;            /* update out pointer from shared memory */
      iter--;                    /* decrement timeout counter */
    }
    if(iter <= 0) { time_write += (mgi_time() - time_in) ; return(WRITE_TIMEOUT); }
    token = *str++ ;
    token <<= 8 ; if(ntok >  2) token |= *str++ ;
    token <<= 8 ; if(ntok >  1) token |= *str++ ;
    token <<= 8 ; if(ntok >  0) token |= *str++ ;
    ntok = ntok -4;
    if(ntok < 0) ntok = 0;
    shm->data[in] = token;
    in = inplus;
  }
  shm->in = in;  /* update in pointer in shared memory */
  time_write += (mgi_time() - time_in) ;
  return(ntok);
}

/********************************************************************************************
 * write into shared memory buffer, integer/real/double/character version
 ********************************************************************************************/
static int shm_write(mgi_shm_buf *shm,void *buf,int nelem,int type,int timeout){
  volatile int in, out, inplus;
  int limit;
  int ntok, navail;
  unsigned int *buffer = (unsigned int *) buf ;
  int maxiter = timeout*1000 ; /* 1000 iterations is one second */
  int iter, i;

  if(shm->write_status != MGI_SHM_ACTIVE) return(WRITE_ERROR) ; /* not properly setup to write */

  time_in = mgi_time();
  if(type == 'C') return ( shm_write_c(shm,buf,nelem,timeout) );                    /* characters */
  if(type != 'I' && type != 'R' && type != 'D') return(WRITE_ERROR) ;       /* unsupported type */
  ntok = nelem;
  if(type == 'D') ntok = nelem*2 ;                 /* 8 byte tokens */
  bytes_write += ntok*4;

  in = shm->in ; out = shm->out ; limit = shm->limit;
//  fprintf(stderr,"DEBUG: Write ntok=%d\n",ntok);
  while(ntok > 0){
    inplus = (in+1 > limit) ? 0 : in+1 ;
    iter = maxiter;
    while(inplus == out && iter > 0) {       /* shared memory circular buffer is full */
      shm->in = in;              /* update in pointer in shared memory */
      usleep(1000);              /* sleep for 1 millisecond */
      wait_write += 1000;
      out = shm->out;            /* update out pointer from shared memory */
      iter--;                    /* decrement timeout counter */
      if(iter <= 0) { time_write += (mgi_time() - time_in) ; return(WRITE_TIMEOUT); }
    }
    if(in >= out){   /* can write into  in -> limit except if out == 0 (can only use in -> limit-1)*/
      navail = (limit+1) - in;
      if(out == 0) navail--;
    }else{          /* can write into  in -> out-1   */
      navail = (out-1) - in;
    }
    navail = (navail > ntok) ? ntok : navail;             /* only need to write ntok tokens */
    if(navail < MEMCPY_LIMIT) {        /* special case for writing few tokens */
      for (i=0 ; i<navail ; i++) shm->data[in+i] = buffer[i];
    }else{
      memcpy(&shm->data[in],buffer,sizeof(int)*navail);     //    shm->data[in] = *buffer;
    }
//    fprintf(stderr,"DEBUG: Write navail=%d, buffer = %d, %d, %d\n",navail,buffer[0],shm->data[in],in);
    in += navail;
    if(in > limit) in = 0;  //    in = inplus;
    ntok -= navail;         //    ntok--;
    buffer += navail;       //    buffer++;
  }
  shm->in = in;  /* update in pointer in shared memory */
  time_write += (mgi_time() - time_in) ;
  return(ntok*sizeof(unsigned int));
}
int ShmWriteBuf(mgi_shm_buf *shm,void *buf,int nelem,int type,int timeout){
  return( shm_write(shm,buf,nelem,type,timeout) );
}
/********************************************************************************************
 * read from shared memory buffer, character version
 ********************************************************************************************/
static int shm_read_c(mgi_shm_buf *shm,void *buf,int nelem, int len,int timeout){
  volatile int in, out;
  int limit;
  int ntok=nelem;
  unsigned char *str = (unsigned char *) buf;
  unsigned int token;
  int pad = len - nelem;
  int maxiter = timeout*1000 ; /* 1000 iterations is one second */
  int iter;
//fprintf(stderr,"DEBUG: shm_read_c reading %d characters\n",nelem);
  in = shm->in ; out = shm->out ; limit = shm->limit;
  bytes_read += ntok;
  while(ntok > 0){
    iter = maxiter;
    while(in == out && iter > 0) {           /* shared memory circular buffer is empty */
      shm->out = out;            /* update out pointer in shared memory */
      usleep(1000);              /* sleep for 1 millisecond */
      wait_read += 1000;
      in = shm->in;              /* update in pointer from shared memory */
      iter--;                    /* decrement timeout counter */
    }
    if(iter <= 0) { time_read += (mgi_time() - time_in) ; return(READ_TIMEOUT); }
    token = shm->data[out] ;
    if(ntok >  3) {*str++ = (token>>24)&0xFF ; *str++ = (token>>16)&0xFF ; *str++ = (token>>8)&0xFF ; *str++ = token&0xFF ; } ;
    if(ntok == 3) {*str++ = (token>>24)&0xFF ; *str++ = (token>>16)&0xFF ; *str++ = (token>>8)&0xFF ; };
    if(ntok == 2) {*str++ = (token>>24)&0xFF ; *str++ = (token>>16)&0xFF ; };
    if(ntok == 1) {*str++ = (token>>24)&0xFF ; };
    ntok = ntok -4;
    if(ntok < 0) ntok = 0;
    out = (out+1 > limit) ? 0 : out+1;   /* bump out */
  }
  shm->out = out;  /* update out pointer in shared memory */
  while(pad-- > 0) *str++ = ' ';
  time_read += (mgi_time() - time_in) ;
  return(nelem);
}

/********************************************************************************************
 * read from shared memory buffer, integer/real/double/character version
 ********************************************************************************************/
static int shm_read(mgi_shm_buf *shm,void *buf,int nelem,int type, int len,int timeout){
  volatile int in, out;
  int limit;
  int ntok, navail;
  unsigned int *buffer = (unsigned int *) buf ;
  int maxiter = timeout*1000 ; /* 1000 iterations is one second */
  int iter, i;

  if(shm->read_status != MGI_SHM_ACTIVE) return(READ_ERROR) ; /* not properly setup to read */

  time_in = mgi_time();
  if(type == 'C') return ( shm_read_c(shm,buf,nelem,len,timeout) );        /* characters */
  if(type != 'I' && type != 'R' && type != 'D') return(READ_ERROR) ;       /* unsupported type */
  ntok = nelem ;
  if(type == 'D') ntok = nelem*2 ;                 /* 8 byte tokens */
  bytes_read += ntok*4;

  in = shm->in ; out = shm->out ; limit = shm->limit;
//  fprintf(stderr,"DEBUG: Read ntok=%d\n",ntok);
  while(ntok > 0){
    iter = maxiter;
    while(in == out && iter > 0) {           /* shared memory circular buffer is empty */
      shm->out = out;            /* update out pointer in shared memory */
      usleep(1000);              /* sleep for 1 millisecond */
      wait_read += 1000;
      in = shm->in;              /* update in pointer from shared memory */
      iter--;                    /* decrement timeout counter */
      if(iter <= 0) { time_read += (mgi_time() - time_in) ; return(READ_TIMEOUT); }
    }
    /* at this point out cannot be equal to in */
    if(out < in){   /* can get out -> in-1 */
      navail = in - out;
    }else{  /* can get out -> limit */
      navail = (limit + 1) - out;
    }
    navail = (navail > ntok) ? ntok : navail;             /* only need to read ntok tokens */
    if(navail < MEMCPY_LIMIT) {   /* special case for reading few tokens */
      for (i=0 ; i<navail ; i++) buffer[i] = shm->data[out+i];
    }else{
      memcpy(buffer,&shm->data[out],sizeof(int)*navail);    //    *buffer = shm->data[out] ;
    }
//    fprintf(stderr,"DEBUG: Read navail=%d, buffer = %d, %d, %d\n",navail,buffer[0],shm->data[out],out);
    out += navail;                     /* bump out */
    out = (out > limit) ? 0 : out;     //    out = (out+1 > limit) ? 0 : out+1;
    ntok -= navail;                    //    ntok--;
    buffer += navail;                  //    buffer++;
  }
  shm->out = out;  /* update out pointer in shared memory */
  time_read += (mgi_time() - time_in) ;
  return(nelem);
}
int ShmReadBuf(mgi_shm_buf *shm,void *buf,int nelem,int type, int len,int timeout){
  return( shm_read(shm,buf,nelem,type,len,timeout) );
}
/********************************************************************************************
 * status = mgi_write(channel,buffer,nelem,dtype)
 * channel   : channel number obtained from mgi_init
 * buffer    : data to be written to channel (integer/real/double/character) intent(IN)
 * nelem     : number of data pieces
 * dtype     : dta type 'C', 'I', 'R', 'D'
 *
 * status    : number of bytes not written or -1 if invalid channel number
 ********************************************************************************************/
ftnword f77_name (mgi_write) (ftnword *f_chan, void *buffer, ftnword *f_nelem, char *dtype, F2Cl ltype)
     /* to write elements from "buffer" into the specified channel
	opened for WRITEMODE. It actually writes

	The following data types (dtype) are accepted:
	'C': character
	'I': integer
	'R': real
	'D': real*8 ; note that only the precision of a real would be kept
     */
{
  int nb, chan, nelem;
  int lnblnk_();
  int len_typ;
  int timeout;
  int status;

  chan = (int) *f_chan;
  nelem = (int) *f_nelem;
  char *tmpstr;

  if(validchan(chan) != 0) return(CONNECTION_ERROR); /* invalid channel */

  if( nelem <= 0 )
    {
      fprintf(stderr,"ERROR: (mgi_write) cannot write data with length = %d\n", nelem);
      return WRITE_ERROR;
    }

  if( chn[chan].gchannel < 0 )
    {
      fprintf(stderr,"ERROR: (mgi_write) cannot connect to server using descriptor: '%d'!!!\n", chn[chan].gchannel);
      return WRITE_ERROR;
    }

  if(chn[chan].shmbuf != NULL) {                              /* shared memory channel  */
    timeout = chn[chan].timeout;
    len_typ = nelem << 8 ;
    len_typ |= *dtype ;
    status = shm_write(chn[chan].shmbuf,&len_typ,1,'I',timeout) ;  /* record length + type */
    if(status != 0) return(status);

    call_write++;
    if(*dtype != 'C'){
      return( shm_write(chn[chan].shmbuf,buffer,nelem                    ,*dtype,timeout) );
    }else{
      return( shm_write(chn[chan].shmbuf,buffer,(nelem<ltype)?nelem:ltype,*dtype,timeout) );
    }
  }

  /* if we get here, it is a "gossip" TCP/IP channel */
#if defined(WITHOUT_GOSSIP)
  return(WRITE_ERROR);  /* MUST be a shared memory channel */
#else

  if ( *dtype == 'C' )
    {
      nelem = ( *f_nelem < ltype ) ? (int) *f_nelem:ltype;

      tmpstr = (char *)malloc(nelem + 1);

#ifdef DEBUG
      fprintf(stderr,"MGI_WRITE: data type = %c, elts Nbr = %d, strlen = %d,  subchannel = %s\n", dtype[0], nelem, ltype, chn[chan].name);
#endif

      strncpy( tmpstr, (char *)buffer, nelem);
      tmpstr[nelem] = '\0';

      if ((nb = bwrite(chan, (unsigned char *)tmpstr, nelem, dtype)) > 0)
	{
          fprintf(stderr,"ERROR: (mgi_write) channel '%s': %d bytes written (character) \n", chn[chan].name, nb);
	  free( tmpstr );
	  return WRITE_ERROR;     /* return number of bytes not sent instead ? */
	}
      free( tmpstr );
    }


  else if (*dtype == 'I' || *dtype == 'R' || *dtype == 'D' )
    {
      chn[chan].nblks++;

#ifdef DEBUG
      fprintf(stderr,"MGI_WRITE: data type = %c, elts Nbr = %d, subchannel = %s\n", dtype[0], nelem, chn[chan].name);
      fprintf(stderr,"MGI_WRITE: data type = %s\n", dtype);
      fprintf(stderr,"MGI_WRITE: elts Nbr = %d\n", nelem);
#endif

      if ((nb = bwrite(chan, (unsigned char *)buffer, nelem, dtype)) > 0)
	{
          fprintf(stderr,"ERROR: (mgi_write) channel '%s': only %d bytes written\n", chn[chan].name, nb);
	  return WRITE_ERROR;     /* return number of bytes not sent instead ? */
	}
    }

  else
    {
      fprintf(stderr,"ERROR: (mgi_write) channel '%s': Unknown data type: %c\n", chn[chan].name, *dtype);
      return WRITE_TYPE_ERROR;
    }

  if(nb == WRITE_TIMEOUT)
    {
      if(get_timeout_signal(chn[chan].gchannel))
	{
	  if (*dtype == 'C')
            fprintf(stderr, "ERROR: (mgi_write) TIMEOUT for write '%d of Character data' \n", nelem);

	  else if(*dtype == 'I')
            fprintf(stderr, "ERROR: (mgi_write) TIMEOUT for write '%d of Integer data'\n", nelem);

	  else if(*dtype == 'R')
            fprintf(stderr, "ERROR: (mgi_write) TIMEOUT for write '%d of Real data' \n", nelem);

	  else if(*dtype == 'D')
            fprintf(stderr, "ERROR: (mgi_write) TIMEOUT for write '%d of Double data' \n", nelem);

	  return signal_timeout(chn[chan].gchannel);
	}
    }

  return nb;
#endif
}
ftnword f77_name(mgi_write_c) (ftnword *f_chan, void *buffer, ftnword *f_nelem, char *dtype, F2Cl lbuf, F2Cl ltype)
{
  return( f77_name(mgi_write) (f_chan, buffer, f_nelem, dtype, lbuf) );
}
int C_mgi_write (ftnword *f_chan, void *buffer, int c_nelem, char c_dtype){
  F2Cl ltype=1;
  char dtype=c_dtype;
  ftnword f_nelem=c_nelem;
  return( f77_name (mgi_write) (f_chan, buffer, &f_nelem, &dtype, ltype) );
}
/********************************************************************************************
 * status = mgi_read(channel,buffer,nelem,dtype)
 * channel   : channel number obtained from mgi_init
 * buffer    : data to be read from channel (integer/real/double/character) intent(OUT)
 * nelem     : number of data pieces
 * dtype     : dta type 'C', 'I', 'R', 'D'
 *
 * status    : number of bytes read or negative error code
 ********************************************************************************************/
ftnword f77_name (mgi_read) (ftnword *f_chan, void *buffer, ftnword *f_nelem, char *dtype, F2Cl ltype)

     /* to read elements directly from the data file related to the
	specified channel into "buffer". The channel must be opened for
	READMODE only.
	The following data types (dtype) are accepted:
	'C': character
	'I': integer (int)
	'R': real    (float)
	'D': real*8  (double)
     */
{
  int ier, chan, nelem;
  int lt = ltype;
  int len_typ;
  char typ;
  int timeout;
  int status;

  chan = (int) *f_chan;
  nelem = (int) *f_nelem;

  if(validchan(chan) != 0) return(CONNECTION_ERROR); /* invalid channel */

  if(nelem <= 0)
    {
      fprintf(stderr,"ERROR: (mgi_read) cannot read data with length = %d\n", nelem);
      /* return -1; */
      return DATA_LENGTH_ERROR;
    }

  bzero(buffer, nelem);  /* nelem * size_of_element would me more appropriate */

  if(chn[chan].shmbuf != NULL){                                  /* this is shared memory channel   */
    timeout = chn[chan].timeout;
    status = shm_read(chn[chan].shmbuf,&len_typ,1,'I',1,timeout);
    if(status <= 0) return(status);
    if( ((len_typ >> 8) != nelem) || ((len_typ & 0xFF) != *dtype) ) {
      typ = len_typ & 0xFF ;
      fprintf(stderr,"ERROR: (mgi_read) length/type mismatch, expected: %d/%c, got: %d/%c\n",nelem,*dtype,len_typ >> 8,typ);
      return(READ_ERROR);
    }
    call_read++;
    if(*dtype != 'C'){
      return( shm_read(chn[chan].shmbuf,buffer,nelem,*dtype,nelem,timeout) ) ;
    }else{
      return( shm_read(chn[chan].shmbuf,buffer,nelem,*dtype,lt   ,timeout) ) ;
    }
  }

/* if we get here, it is a "gossip" TCP/IP channel */
#if defined(WITHOUT_GOSSIP)
  return(READ_ERROR);  /* MUST be a shared memory channel */
#else

  ier = send_command_to_server(chn[chan].gchannel, "READ");

  if(ier < 0)
    {
      fprintf(stderr,"ERROR: (mgi_read) unable to send write command for channel: '%s'\n", chn[chan].name);
      return SEND_COMMAND_ERROR;
    }

  if (*dtype == 'I')
    { /* integer */

//      fprintf(stderr, "MGI_READ: 'Integer', elts Nbr = %d, channel = '%s'\n", nelem, chn[chan].name);

      buffer = (int *)read_record( chn[chan].gchannel, (int *)buffer, &nelem, nelem, sizeof(int) );

      if(buffer != NULL)
	{
	  get_ack_nack( chn[chan].gchannel );
	  return nelem;
	}
      else
	{
	  if( get_timeout_signal(chn[chan].gchannel) )
	    {
              fprintf(stderr, "ERROR: (mgi_read) TIMEOUT reading Integer(s) \n" );
	      return READ_TIMEOUT;
	    }
	  else
	    {
              fprintf( stderr, "ERROR: (mgi_read) Bad data reading Integer(s)\n" );
	      return READ_ERROR;
	    }
	}
    }

  else if (*dtype == 'R')
    { /* float */

//      fprintf(stderr, "MGI_READ: 'Real', elts Nbr = %d, channel = '%s'\n", nelem, chn[chan].name);

      buffer = (float *)read_record(chn[chan].gchannel, (float *)buffer, &nelem, nelem, sizeof(int));

      if(buffer != NULL)
	{
	  get_ack_nack(chn[chan].gchannel);
	  return  nelem;

	}
      else
	{
	  if( get_timeout_signal( chn[chan].gchannel ) )
	    {
              fprintf(stderr, "ERROR: (mgi_read)  TIMEOUT for read 'Real' \n");
	      return READ_TIMEOUT;
	    }
	  else
	    {
              fprintf( stderr, "ERROR: (mgi_read) problem reading Float data\n" );
	      return READ_ERROR;
	    }
	}
    }
  else if (*dtype == 'D')
    { /* double */

//      fprintf(stderr, "MGI_READ: 'Double', Element's Nbr = %d, channel = '%s'\n", nelem, chn[chan].name);

      buffer = (double *)read_record(chn[chan].gchannel, (double *)buffer, &nelem, nelem, sizeof(double));

      if(buffer != NULL)
	{
	  get_ack_nack(chn[chan].gchannel);
	  return nelem;
	}
      else
	{
	  if( get_timeout_signal( chn[chan].gchannel ) )
	    {
              fprintf(stderr, "ERROR: (mgi_read) TIMEOUT reading Double(s)\n");
	      return READ_TIMEOUT;
	    }
	  else
	    {
              fprintf( stderr, "ERROR: (mgi_read) Bad data reading Double(s)\n" );
	      return READ_ERROR;
	    }
 	}
     }

  else if (*dtype == 'C')
    { /* character */
      int i;
      char *temp = (char *)buffer;

      for(i = 0; i < ltype ; i++ )
	{
	  temp[i] = ' ';
	}

      buffer = (char *)read_record(chn[chan].gchannel, (char *)buffer, &nelem, nelem, sizeof(char));

      for(i = nelem+1 ; i < ltype ; i++ )
	{
	  temp[i] = ' ';
	}

      if(buffer != NULL)
	{
	  get_ack_nack(chn[chan].gchannel);
	  return  nelem;
	}
      else
	{
	  if( get_timeout_signal( chn[chan].gchannel ) )
	    {
              fprintf(stderr, "ERROR: (mgi_read) TIMEOUT for read 'Character'\n");
	      return READ_TIMEOUT;
	    }
	  else
	    {
              fprintf( stderr, "ERROR: (mgi_read) Problem read Character data\n" );
	      return READ_ERROR;
	    }
	}
    }

  else
    {
      fprintf(stderr,"ERROR: (mgi_read) channel '%s': Unknown data type: %c\n", chn[chan].name, *dtype);
      return READ_TYPE_ERROR;
    }
#ifdef UNUSED_CODE
  if(ier == CLOSE)  /* cannot happen !! */
    {
      close_channel(chn[chan].gchannel, chn[chan].name);
    }
#endif
  return ier;
#endif
}
/*
 *added mgi_read_c Fortran callable routine because of extra character string
*/
ftnword f77_name(mgi_read_c) (ftnword *f_chan, void *buffer, ftnword *f_nelem, char *dtype, F2Cl lbuf, F2Cl ltype)
{
  return ( f77_name(mgi_read) (f_chan, buffer, f_nelem, dtype, lbuf) );
}
int C_mgi_read (int c_chan, void *buffer, int c_nelem, char c_dtype){
  ftnword f_chan=c_chan;
  ftnword f_nelem=c_nelem;
  F2Cl ltype=1;
  char dtype=c_dtype;
  return( f77_name (mgi_read) (&f_chan, buffer, &f_nelem, &dtype, ltype) );
}
#if defined(WITH_MAIN)
main(int argc, char**argv){
  printf("Hello World!!\n");
}
#endif
