/*
 * Copyright (C) 2014       ESCER center, UQAM
 *
 * This code is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this code; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*
Author: Michel Valin, UQAM, 2014/03/10

This program runs in the background and terminates when the parent process terminates.
It will also terminate if a specified number of "attaches" to the memory segment occurred.

Before it terminates, the program will mark the memory segment for destruction

It expects one or two input parameters:
Arg1: size of a shared memory segment to create and watch
Arg2: channel name for mgi

The program will wait 10 milliseconds between checks of parent process existence.

*/

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <mgi.h>

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

static char *pname="NoNe";

void usage(){
  fprintf(stderr,"ERROR: bad arguments\n");
  fprintf(stderr,"usage: %s [-v] shm_size[M] [mgi_channel_name]\n",pname);
  fprintf(stderr,"       %s --stop [mgi_channel_name]\n",pname);
  fprintf(stderr,"          shm_size = shared memory segment size in [K/M]Bytes\n");
  exit(1);
}

int mgi_stop(int argc, char **argv){

  char channel_filename[1024];
  char char_id[32];
  int shm_id;
  int items;
  int temp;
  FILE *FD;
  mgi_shm_buf *shm;
  struct shmid_ds shm_stat;
  size_t shm_size;
  int shm_size_k;

  if(*argv[1] == '-' ){
    usage();
  }
  snprintf(channel_filename,sizeof(channel_filename),"%s/.gossip/SHM/%s.id",getenv("HOME"),argv[1]);
  FD=fopen(channel_filename,"r");
  if(FD != NULL){
    items = fscanf(FD,"%d",&shm_id);  /* get shared memory segment ID */
    if(items != 1) {
      fprintf(stderr,"ERROR: (mgi_stop) invalid shared memory channel id \n");
      exit(1);
    }
    fprintf(stderr,"INFO: closing shared memory segment %d associated to channel %s\n",shm_id,channel_filename);
    fclose(FD);
  }

  shm = (mgi_shm_buf *) shmat(shm_id, NULL, 0);
  if(shm == (void *) -1){
    fprintf(stderr,"ERROR: (mgi_stop) cannot attach shared memory segment\n");
    exit(1);
  }

  temp = shmctl(shm_id,IPC_STAT,&shm_stat);
  shm_size = shm_stat.shm_segsz / 1024;
  shm_size_k = shm_size;
  fprintf(stderr,"INFO: Shared memory segment '%d' has size %d KBytes\n",shm_id,shm_size_k);
  temp = shm_stat.shm_nattch;
  fprintf(stderr,"INFO: %d process(es) currently attached\n",temp-1);

  if( shm->read_status == MGI_SHM_IDLE ) {
    shm->read_status = MGI_SHM_ACTIVE;
    sleep(1);
    shm->read_status = MGI_SHM_IDLE;
  }
  if( shm->write_status == MGI_SHM_IDLE ) {
    shm->write_status = MGI_SHM_ACTIVE;
    sleep(1);
    shm->write_status = MGI_SHM_IDLE;
  }
  temp = shmctl(shm_id,IPC_STAT,&shm_stat);
  shm_size = shm_stat.shm_segsz / 1024;
  shm_size_k = shm_size;
  temp = shm_stat.shm_nattch;
  if(temp==1) fprintf(stderr,"INFO: no process attaching %s (%d) any more\n",channel_filename,shm_id);
  exit(0);
}

/* unused as yet borrowed code */
int mgi_shm(int argc, char **argv){
int i, nc1,nc2;
pid_t pp=getppid();  /* get pid of my parent */
#if defined(DAEMONIZE)
pid_t sid;
#endif
useconds_t sleep_duration=10000;  /* 10 milliseconds */
int to_watch, to_create;
struct shmid_ds shm_stat;
size_t shm_size;
int shm_size_k;
//char *shm;
mgi_shm_buf *shm;
char channel_filename[1024];
char log_filename[1024];
int fd;
FILE *FD, *LOG;
int nints;
char multiplier='K';
int items=0;
int read_att=0;
int writ_att=0;
int attach=0;
int read_status=0;
int write_status=0;
volatile int in=0;
volatile int out=0;
long long time0;
long long time1;
double rtime;
struct timeval timetag;
size_t count, count1, count2;
int monitor=0;

pname = argv[0];

if (argc > 1){
  if(strcmp(argv[1],"-v") == 0) {    /* -v */
    monitor = 1;
    fprintf(stderr,"INFO: full logging will be enabled\n");
    argc--;
    argv++;
  }
}

if (argc > 1){
  if(strcmp(argv[1],"--stop") == 0) {    /* -v */
    if(argc != 3) {
      usage();
    }else{
      fprintf(stderr,"INFO: will send stop signal for channel %s\n",argv[2]);
      argc--;
      argv++;
      mgi_stop(argc,argv);
      exit(0);
    }
  }
}
if(argc < 2 || argc >3) {
  usage();
}
if(**argv == '-' ){
  usage();
}
//to_create=atoi(argv[1]);
items = sscanf(argv[1],"%d %1c",&to_create,&multiplier);
if(to_create <= 0){
  usage();
}

if(multiplier=='M' || multiplier=='m') to_create *= 1024;   /* megabytes, not kilobytes */
FD=NULL;
if(argc == 3){
  snprintf(channel_filename,sizeof(channel_filename),"%s/.gossip/SHM/%s.id",getenv("HOME"),argv[2]);
  unlink(channel_filename);   /* remove previously existing file, ignore errors */
  FD=fopen(channel_filename,"w");
  if(FD==NULL){
    fprintf(stderr,"ERROR: cannot open %s\n",channel_filename);
    exit(1);
  }
}
to_watch=shmget(IPC_PRIVATE,to_create*1024,0600);
if(to_watch == -1) {
  fprintf(stderr,"ERROR: shared memory segment creation failed\n");
  exit(1);
}
if(FD != NULL){
  fprintf(FD,"%d",to_watch);
  fclose(FD);
}
fprintf(stdout,"%d\n",to_watch);
fflush(stdout);
fclose(stdout);

shm = (mgi_shm_buf *) shmat(to_watch, NULL, 0);  /* attach segment just created */
if(shm == (void *) -1) {
  fprintf(stderr,"ERROR: memory segment %d cannot be attached\n",to_watch);
  exit(1);
}else{
  fprintf(stderr,"INFO: memory segment %d successfully attached\n",to_watch);
}

usleep(50000);  /* 50 milliseconds */
/* on non linux systems, it is not possible to attach a shared memory area that is marked for deletion */
#if defined(linux)
i = shmctl(to_watch,IPC_RMID,NULL);   /* immediately mark segment for removal if linux */
#endif
usleep(50000);  /* 50 milliseconds */

i = shmctl(to_watch,IPC_STAT,&shm_stat);
shm_size = shm_stat.shm_segsz / 1024;
shm_size_k = shm_size;
fprintf(stderr,"Monitoring shared memory id '%d' of size %d KBytes\n",to_watch,shm_size_k);

/* initialize memory segment for use by mgilib */
shm->read_lock = 0;
shm->write_lock = 0;
shm->read_status = MGI_SHM_IDLE;
shm->write_status = MGI_SHM_IDLE;
shm->first = 0;                     /* first position in buffer */
shm->in = 0;                        /* insertion position in buffer */
shm->out = 0;                       /* extraction position in buffer */
shm_size = shm_stat.shm_segsz;      /* size of memory area */
shm_size -= sizeof(mgi_shm_buf);    /* minus structure size */
shm_size /= sizeof(unsigned int);   /* convert to number of unsigned int */
shm->limit = shm_size;              /* last position in buffer */

snprintf(channel_filename,sizeof(channel_filename),"%s/.gossip/SHM/%s",getenv("HOME"),(argc>2) ? argv[2] : "default");
snprintf(log_filename,sizeof(log_filename),"%s/.gossip/SHM/%s.LOG",getenv("HOME"),(argc>2) ? argv[2] : "default");
fprintf(stderr,"INFO: trying to open %s\n",channel_filename);
fd=open(channel_filename,O_RDONLY);
if(fd > 0) {                        /* channel restart file exists, read it into buffer */
  nints = read(fd,shm->data,shm->limit*sizeof(unsigned int)) / sizeof(unsigned int) ;
  shm->in = nints;
  fprintf(stderr,"INFO: read %d elements from '%s'\n",nints,channel_filename);
  if(nints == shm->limit){
    fprintf(stderr,"WARNING: restart file possibly larger than shared memory buffer\n");
  }
  close(fd);
  unlink(channel_filename); /* remove what we have just read */
}

i = fork();
if(i > 0) exit(0);  /* parent exits */
if(i < 0) exit(1);  /* fork failed */

#if defined(DAEMONIZE)
sid = setsid();  /* daemonize */
if(sid < 0) exit(1) ;
i = fork();
if(i > 0) exit(0);  /* parent exits */
if(i < 0) exit(1);  /* fork failed */
#endif

LOG = stderr;
if(monitor){
  LOG = fopen(log_filename,"a");
  if(LOG == NULL){
    fprintf(stderr,"INFO: cannot open log file '%s'\n",log_filename);
    LOG = stderr;
  }
}
while(1){
  usleep(sleep_duration);            /* 10 milliseconds */
/* ----------------------------------------------------------------------------------------------------------------- */
  if(kill(pp,0)) {                   /* original parent no longer exists, time to save leftover data and quit */
    read_att=1 ;
    writ_att=1 ;
    shm->read_status=MGI_SHM_IDLE ;     /* simulate activate and close on read and write */
    shm->write_status=MGI_SHM_IDLE ;
  }
/* ----------------------------------------------------------------------------------------------------------------- */
  if(shm->read_status == MGI_SHM_ACTIVE) read_att = 1;     /* a process opend for reading */
  if(shm->write_status == MGI_SHM_ACTIVE) writ_att = 1;    /* a process opend for writing */
/* ----------------------------------------------------------------------------------------------------------------- */
//  if(read_att==1 && writ_att==1) break;  /* channel has been opened at both ends , we can quit */
  if(read_att==1 && writ_att==1 && shm->read_status==MGI_SHM_IDLE && shm->write_status==MGI_SHM_IDLE){ /* everybody opened, everything is closed */
    if(shm->in == shm->out) break; /* buffer is empty, quit */
    fd = open(channel_filename,O_CREAT+O_RDWR,0644);  /* save leftover data */
    if(fd >0) {
      if(shm->out < shm->in) {   /* write from out to in-1 */
        count = sizeof(shm->data[0])*(shm->in-shm->out);
        nc1 = write(fd, (void *) shm->data+shm->out, count);
        fprintf(stderr,"INFO: %d/%d bytes written to %s(%d)\n",nc1,count,channel_filename);
        if(monitor)fprintf(LOG,"INFO: %d/%d bytes written to %s\n",nc1,count,channel_filename);
      }else{                     /* write from out to limit, then first to in-1 */
        count1 = sizeof(shm->data[0])*(shm->limit-shm->out+1);
        nc1 = write(fd, shm->data+shm->out, count1);
        count2 = sizeof(shm->data[0])*(shm->in-shm->first);
        nc2 = write(fd, shm->data+shm->first, count2);
        count = count1 + count2;
        fprintf(stderr,"INFO: %d/%d bytes written to %s\n",nc1+nc2,count,channel_filename);
        if(monitor)fprintf(LOG,"INFO: %d/%d bytes written to %s\n",nc1+nc2,count,channel_filename);
      }
      fprintf(stderr,"INFO: leftover data saved to '%s'\n",channel_filename);
      if(monitor)fprintf(LOG,"INFO: leftover data saved to '%s'\n",channel_filename);
      close(fd);
    }else{
      fprintf(stderr,"ERROR: cannot open %s to save leftover data\n",channel_filename);
      if(monitor)fprintf(LOG,"ERROR: cannot open %s to save leftover data\n",channel_filename);
    }
    break;
  }
/* ----------------------------------------------------------------------------------------------------------------- */
  i = shmctl(to_watch,IPC_STAT,&shm_stat);
/* ----------------------------------------------------------------------------------------------------------------- */
  if(monitor){   /* if full monitoring is active, check if anything changed */
    if(in != shm->in || out != shm->out || read_status != shm->read_status || write_status != shm->write_status || attach != shm_stat.shm_nattch){
      gettimeofday(&timetag,NULL);
      time1 = timetag.tv_sec;
      time1 *= 1000000;
      time1 += timetag.tv_usec;
      rtime = time1 - time0;
      rtime /= 1000000.;
      if(monitor) {
        fprintf(LOG,"%12.3f: size = %d, attach count = %d, creator=%d, last attach=%d ",rtime,
                (int)shm_stat.shm_segsz,(int)shm_stat.shm_nattch,shm_stat.shm_cpid,shm_stat.shm_lpid);
        fprintf(LOG," first=%d, in=%d, out=%d, limit=%d, rs=%d, ws=%d\n",
                shm->first,shm->in,shm->out,shm->limit,shm->read_status,shm->write_status);
      }
    }
    in = shm->in;
    out = shm->out;
    read_status = shm->read_status;
    write_status = shm->write_status;
    attach = shm_stat.shm_nattch;
  }
/* ----------------------------------------------------------------------------------------------------------------- */
  if(i == -1) exit(0);               /* segment no longer accessible, quit */
  if(monitor) {
  }

#if defined(NOT_USED_ANYMORE)
  if(shm_stat.shm_nattch >= 4) {   /* segment attached by 3 other processes */
    break;                           /* job done, exit */
  }
#endif

}  /* while */
i = shmctl(to_watch,IPC_RMID,NULL);   /* mark for removal if not already done */
i = shmdt(shm);                       /* do not keep the memory area attached, we do not want to prevent release */
fclose(LOG);
exit (0);
}
