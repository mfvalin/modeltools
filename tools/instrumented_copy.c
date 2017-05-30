#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>

#define BUFFER_SIZE 4096*1024

static char buf[BUFFER_SIZE];
static struct timeval tv1, tv2;
static double t1, t2, sumt1, sum2t1, mint1, maxt1, sumt2, sum2t2, mint2, maxt2;
int status, nt1, nt2;
double gt1[32768], gt2[32768];

main(int argc, char **argv)
{
  int fdi, fdo;
  ssize_t nread, nwritten;
  size_t toread = BUFFER_SIZE;
  size_t towrite;
  int i;
  
  if( (fdi = open(argv[1],O_RDONLY)) < 0) { printf("fdi=%d\n",fdi) ; exit(1); }
  if( (fdo = open(argv[2],O_RDWR+O_CREAT, 0777 )) < 0) { printf("fdo=%d\n",fdi) ; exit(1); }

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
  
  while(1){
    status = gettimeofday(&tv1, NULL);
    nread = read(fdi,buf,toread);
    if(nread <= 0) break;
    status = gettimeofday(&tv2, NULL);

    t1 = tv1.tv_sec + tv1.tv_usec * .000001 ;
    t2 = tv2.tv_sec + tv2.tv_usec * .000001 ;
    t1 = t2 - t1;
    gt1[nt1] = t1;
    sumt1 += t1;
    sum2t1 += t1*t1;
    nt1++;
    maxt1 = (maxt1 < t1) ? t1 : maxt1;
    mint1 = (mint1 > t1) ? t1 : mint1;

    towrite = nread;
    status = gettimeofday(&tv1, NULL);
    nwritten = write(fdo,buf,towrite);
    if(nwritten != towrite) exit(1);
    status = gettimeofday(&tv2, NULL);

    t1 = tv1.tv_sec + tv1.tv_usec * .000001 ;
    t2 = tv2.tv_sec + tv2.tv_usec * .000001 ;
    t2 = t2 - t1;
    gt2[nt2] = t2;
    sumt2 += t2;
    sum2t2 += t2*t2;
    nt2++;
    maxt2 = (maxt2 < t2) ? t2 : maxt2;
    mint2 = (mint2 > t2) ? t2 : mint2;
  }
  close(fdi);
  close(fdo);
  fprintf(stderr,"READ  stats : min = %8.6f max = %8.6f avg = %8.6f\n",mint1,maxt1,sumt1/nt1);
  fprintf(stderr,"WRITE stats : min = %8.6f max = %8.6f avg = %8.6f\n",mint2,maxt2,sumt2/nt2);
  for(i=0 ; i<nt1 ; i++){ fprintf(stdout,"%8.6f , %8.6f\n",gt1[i],gt2[i]); }
  exit(0);
}