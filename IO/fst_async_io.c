#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <aio.h>
#include <signal.h>
#include <mpi.h>
#include <errno.h>

struct asyncio_info{
  struct aiocb *aio_item;
  struct aiocb **aio_list;
  off_t *aio_ret;
  int *aio_stat;
  int nreq;
  int fd;
};
typedef struct asyncio_info *asyncio_control;

asyncio_control wa_async_write(int fd, char **buffers, off_t *offsets, size_t *sizes, int nreq){
  struct aiocb *aio_item;
  struct aiocb **aio_list;
  int *aio_stat;
  off_t *aio_ret;
  asyncio_control answer;
  int i;

  aio_item = (struct aiocb *)malloc(sizeof(struct aiocb) * nreq);
  aio_list = (struct aiocb **) malloc(sizeof(struct aiocb*) * nreq);
  aio_ret = (off_t *)malloc(sizeof(off_t) * nreq);
  aio_stat = (int *)malloc(sizeof(int) * nreq);
  answer = (asyncio_control) malloc(sizeof(struct asyncio_info));
  answer->aio_item = aio_item;
  answer->aio_list = aio_list;
  answer->aio_ret = aio_ret;
  answer->aio_stat = aio_stat;
  answer->nreq = nreq;
  answer->fd = fd;
  
  for(i=0 ; i<nreq ; i++){
    aio_item[i].aio_fildes = fd ;
    aio_item[i].aio_offset = offsets[i] ;
    aio_item[i].aio_buf = (void *) buffers[i] ;
    aio_item[i].aio_nbytes = sizes[i] ;
    aio_item[i].aio_reqprio = 0 ;
    aio_item[i].aio_lio_opcode = LIO_WRITE ;
    aio_item[i].aio_sigevent.sigev_notify = SIGEV_NONE ;
    aio_list[i] = &aio_item[i];
    aio_stat[i] = EINPROGRESS ;
  }
  lio_listio(LIO_NOWAIT, aio_list, nreq, NULL);
  return answer;
}

int wa_async_wait(asyncio_control acb){
  int i;
  int nbusy = acb->nreq;

  while(nbusy > 0){
    for(i=0 ; i<32 ; i++){
      if(acb->aio_stat[i] == EINPROGRESS){
        acb->aio_stat[i] = aio_error(acb->aio_list[i]);
        if(acb->aio_stat[i] == 0) {
          nbusy--;
          acb->aio_ret[i] = aio_return(acb->aio_list[i]);
        }
      }
    }
    if(nbusy == 0) break;
    printf("nbusy = %d\n",nbusy);
    usleep(50000);
  }
  free(acb->aio_item);
  free(acb->aio_stat);
  free(acb->aio_list);
}

#if defined(SELF_TEST)
static int BLKSIZE = 1024*1024*8;

// static struct aiocb aio_item[32];
// static struct aiocb *aio_list[32];

main(int argc, char **argv){
  int i;
  int fd;
  char *buffers[32];
  off_t offsets[32];
  size_t sizes[32];
//   int aio_stat[32];
//   ssize_t aio_ret[32];
//   int nbusy;
  double t1, t2;
  void *acb;

  MPI_Init(&argc,&argv);
  fd = open(argv[1],O_WRONLY + O_CREAT + O_EXCL, 0777);
  if(fd < 0) exit(1);
  if(argc == 3) BLKSIZE = 1024*1024*atoi(argv[2]);

  for(i=0 ; i<32 ; i++){
    buffers[i] = (char *) malloc(BLKSIZE) ;
    offsets[i] = i * BLKSIZE ;
    sizes[i] = BLKSIZE ;
  }

  acb = wa_async_write(fd, buffers, offsets, sizes, 32);
  t1 = MPI_Wtime();
  wa_async_wait(acb);
  close(fd);
  t2 = MPI_Wtime();
//   for(i=0 ; i<32 ; i++) printf("%d ",aio_stat[i]);
  printf("\n time = %f\n",t2-t1);
  MPI_Finalize();
}
#endif
