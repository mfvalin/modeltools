#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <string.h>
#include <errno.h>

#define MEMBERS 10
#define OFFSET_FACTOR 4

main(int argc, char **argv)
{
  int member = 0;
  int fd;
  int buf[10];
  int i;
  ssize_t nc;
  size_t nw;
  off_t offset;
  useconds_t usec=100000;
  key_t key = 0;
  size_t shmsize = 1024*1024;
  int shmid;
  int *area;
  int *count;
  struct shmid_ds shctl;
  int totcnt = 0;
  int fetch_mem(void *);
  int timeout = 10;

  for (i=0 ; i<10 ; i++) buf[i] = i;
  fd = open(argv[1],O_RDWR|O_CREAT,0700);
  nw = 40;
  nc = write(fd,buf,nw);
  close(fd);

  shmid = shmget(0,shmsize,IPC_CREAT|S_IRUSR|S_IWUSR);
  area = shmat(shmid,NULL,0);
  if(area == (void *) -1) {
    printf("shmid = %d, %s\n",shmid,strerror(errno));
    exit(1);
  }
  shmctl(shmid,IPC_RMID,&shctl);

  while(member < MEMBERS - 1){
    if(fork() > 0) {
      member++;
      usleep(usec);
    }else{
      break;
    }
  }
  count = area+1000;
  fd = open(argv[1],O_RDONLY);
  offset = member * OFFSET_FACTOR;
  lseek(fd,offset,SEEK_SET);
  nw = 4;
  i = -1;
  nc = read(fd,&area[member],nw);
  close(fd);

  if(member != 0) {
    sleep(1);
  }
  count[member] = nc;

  usec = 1000;
  printf("I am member %d , processing %s, fd = %d, read %d, id = %d, area = %p\n",member,argv[1],fd,area[member],shmid,area);
  totcnt = 0;
  while(totcnt != MEMBERS) {
    totcnt = 0;
    for(i=0 ; i<MEMBERS ; i++) if(fetch_mem(count+i) > 0) totcnt++;
    if(totcnt != MEMBERS) usleep(usec);
  }
  printf("I am member %d , after barrier\n",member);
  if(member == 0) {
    for(i=0 ; i<MEMBERS ; i++) printf("%4.4d ",area[i]);
    printf("\n");
    for(i=0 ; i<MEMBERS ; i++) { printf("%4.4d ",fetch_mem(&count[i])); count[i] = 0; }
    printf("\n");
  }
  totcnt = 0;
  while(totcnt != MEMBERS) {
    totcnt = 0;
    for(i=0 ; i<MEMBERS ; i++) if(fetch_mem(count+i) == 0) totcnt++;
    if(totcnt != MEMBERS) usleep(usec);
//     timeout --;
  }
  if(member == 0) {
    for(i=0 ; i<MEMBERS ; i++) printf("%4.4d ",area[i]);
    printf("\n");
    for(i=0 ; i<MEMBERS ; i++) { printf("%4.4d ",fetch_mem(&count[i])); count[i] = 0; }
    printf("\n");
  }
}