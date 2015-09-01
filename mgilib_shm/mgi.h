/* mgi.h */
#ifndef MGI_INCLUDE_ONCE

#define MGI_INCLUDE_ONCE

#define MAX_CHANNELS 24
#define MAX_NAME 125
#define MAX_STR 1024
#define BUFSIZE 40960

#ifndef FALSE
#define FALSE            0
#define TRUE        !FALSE
#endif

#define MGI_SHM_IDLE -2
#define MGI_SHM_ACTIVE 1

typedef struct
{
  int read_status;
  int write_status;
  int first;
  int in;
  int out;
  int limit;
  unsigned int data[1];
} mgi_shm_buf;

typedef struct
{
  mgi_shm_buf *shmbuf;
  int *buffer;
  int pos;
  int gchannel;
  int shmid;
  int fd_data;
  int fd_sig;
  int msgno_W;
  int msgno_R;
  int nblks;
  int timeout;
  char name[MAX_NAME];
  char mode;
} channel;
#endif
