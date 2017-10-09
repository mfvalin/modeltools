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

#define MGI_SHM_IDLE     -2
#define MGI_SHM_ACTIVE    1

#define WITHOUT_SHM_PORT -1
#define WITHOUT_MPI_PORT  1
#define WITH_MPI_PORT     0

typedef struct
{
  int read_lock;
  int write_lock;
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

#define MGI_CHAN_IDLE   -1
#define MGI_CHAN_ACTIVE  1
#define MGI_CHAN_CLOSED  2

typedef struct{
  int read_lock;
  int write_lock;
  int read_status;   // flags used by remote reader [ initialized as 0]
  int write_status;  // flags used by remote writer [ initialized as 0]
  int first;         // start of buffer [normally 0]
  int in;            // insertion index
  int out;           // extraction index
  int limit;         // end of buffer + 1  [index]
  unsigned int data[1];       // place holder, start of data buffer
} mgi_channel_buffer;         // this struct should be usable in the shm and MPI cases

#endif
