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
#include <sys/time.h>

static struct timeval tv1, tv2;
static double t1, t2, sumt1, sum2t1, mint1, maxt1, sumt2, sum2t2, mint2, maxt2;
static double gt1[32768], gt2[32768];
static int status, nt1, nt2;

main(int argc, char **argv)
{
  int member = 0;
  int fd;
  int buf[10];
  char *block;
  int i;
  ssize_t nc;
  size_t nw, nr, ntot;
  off_t offset, eoffset, pos;
  useconds_t usec=100000;
  key_t key = 0;
  size_t shmsize, chunksize;
  int shmid;
  char *area;
  int *count;
  struct shmid_ds shctl;
  int totcnt, nblks;
  int fetch_mem(void *);
  int timeout = 10;
  int offset_factor;
  int tot_members;
  char *file_in = "DoesNotExist";
  int control_size = 1024 ;
//   int lblock[1024*1024];

  sumt1 = 0.0;
  sum2t1 = 0.0;
  mint1 = 999999.0;
  maxt1 = 0;
  nt1 = 0;
  sumt2 = 0.0;
  sum2t2 = 0.0;
  mint2 = 999999.0;
  maxt2 = 0;
  nt2 = 0;

  if(argc < 3) exit(1);
  file_in = argv[1];
  offset_factor = atoi(argv[2]);
  tot_members = atoi(argv[3]);
fprintf(stderr,"INFO: file '%s' read by %d processes in blocks of %d MBytes\n",file_in,tot_members,offset_factor);
  offset_factor = offset_factor * 1024 * 1024;
// exit(0);

//   for (i=0 ; i<10 ; i++) buf[i] = i;
//   fd = open(argv[1],O_RDWR|O_CREAT,0700);
//   nw = 40;
//   nc = write(fd,buf,nw);
//   close(fd);

  chunksize = 1024 * 1024 * 4;  // basic chunk is 4MBytes
  shmsize = chunksize * tot_members + control_size;
  shmid = shmget(0,shmsize,IPC_CREAT|S_IRUSR|S_IWUSR);
  area = shmat(shmid,NULL,0);
  if(area == (void *) -1) {
    fprintf(stderr,"shmid = %d, %s\n",shmid,strerror(errno));
    exit(1);
  }
  shmctl(shmid,IPC_RMID,&shctl);

  while(member < tot_members){      // we will have tot_members readers + 1 writer
    if(fork() > 0) {
      member++;
      usleep(usec);
    }else{
      break;
    }
  }
  if(member < tot_members)  {                          // reader
//     block = area + control_size + member * chunksize ; // 4MByte buffers for reading
    block = area + member * chunksize ; // 4MByte buffers for reading
    offset = member * offset_factor;
    fd = open(argv[1],O_RDONLY);
    eoffset = lseek(fd,offset,SEEK_SET);
    pos = block - area;
// fprintf(stderr,"Reader member %d , processing %s, fd = %d, read block %dMB, offset = %dMB, pos = %dMB, , id = %d, area = %p\n",
// member,argv[1],fd,offset_factor/1024/1024,eoffset/1024/1024,pos/1024/1024,shmid,block);
    nr = chunksize ;
    ntot = offset_factor;
    totcnt = 0;
    nblks = 0;
// fprintf(stderr,"Reader member %d , start reading\n",member);
    while(ntot > 0) {   // loop until the whole file is read
      if(nr > ntot) nr = ntot;
      nc = nr ;

      status = gettimeofday(&tv1, NULL);
      nc = read(fd,block,nr);  // read a block
      status = gettimeofday(&tv2, NULL);

      t1 = tv1.tv_sec + tv1.tv_usec * .000001 ;
      t2 = tv2.tv_sec + tv2.tv_usec * .000001 ;
      t1 = t2 - t1;
      gt1[nblks] = t1;
      sumt1 += t1;
      sum2t1 += t1*t1;
      nt1++;
      maxt1 = (maxt1 < t1) ? t1 : maxt1;
      mint1 = (mint1 > t1) ? t1 : mint1;

      totcnt = totcnt + nc;
      if(nc < 0) {
	offset = 0;
	pos = lseek(fd,offset,SEEK_CUR);
fprintf(stderr,"EOF member %d, pos = %dMB\n",member,pos/1024/1024);
	break;
      }
      ntot = ntot - nc;
//       block = block + nc;
      nblks++;
    }
// fprintf(stderr,"Reader member %d , done reading\n",member);
    close(fd);
fprintf(stderr,"Reader member %d , read %dMB in %d blocks ",member,totcnt/1024/1024,nblks);
fprintf(stderr,"  stats : min = %8.6f max = %8.6f avg = %8.6f\n",mint1,maxt1,sumt1/nt1);
  }else{                                               // writer
// sleep(1);
fprintf(stderr,"Writer member %d , processing %s\n",member,argv[1]);
  }
// fprintf(stderr,"Member %d , exiting\n",member);
exit(0);


  if(member != 0) {
//     sleep(1);
  }
  count[member] = nc;

  usec = 1000;
  printf("I am member %d , processing %s, fd = %d, read %d, id = %d, area = %p\n",member,argv[1],fd,area[member],shmid,area);
  totcnt = 0;
  while(totcnt != tot_members) {
    totcnt = 0;
    for(i=0 ; i<tot_members ; i++) if(fetch_mem(count+i) > 0) totcnt++;
    if(totcnt != tot_members) usleep(usec);
  }
  printf("I am member %d , after barrier\n",member);
  if(member == 0) {
    for(i=0 ; i<tot_members ; i++) printf("%4.4d ",area[i]);
    printf("\n");
    for(i=0 ; i<tot_members ; i++) { printf("%4.4d ",fetch_mem(&count[i])); count[i] = 0; }
    printf("\n");
  }
  totcnt = 0;
  while(totcnt != tot_members) {
    totcnt = 0;
    for(i=0 ; i<tot_members ; i++) if(fetch_mem(count+i) == 0) totcnt++;
    if(totcnt != tot_members) usleep(usec);
//     timeout --;
  }
  if(member == 0) {
    for(i=0 ; i<tot_members ; i++) printf("%4.4d ",area[i]);
    printf("\n");
    for(i=0 ; i<tot_members ; i++) { printf("%4.4d ",fetch_mem(&count[i])); count[i] = 0; }
    printf("\n");
  }
}
