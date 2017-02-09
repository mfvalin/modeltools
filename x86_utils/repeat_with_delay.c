#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>
#include <signal.h>

main(int argc, char **argv){
  int fd;
  int status;
  struct stat sbuf;
  time_t t;
  int delay;
  int tim;
  pid_t ppid = getppid();

  if(argc < 4) {
    printf("USAGE: %s flag_file age_limit command\n",argv[0]);
    exit(1);
  }
  delay = atoi(argv[2]);
  if(delay < 30) {
    printf("ERROR: %s age_linit must be >= 30 sec\n",argv[0]);
    exit(1);
  }

  if(fork()){
    printf("INFO: %s monitor started\n",argv[0]);
    exit(0);
  }
  umask(0111);
  while(1){
    if(kill(ppid, 0)) exit(0);
    fd = open(argv[1],O_RDWR|O_CREAT,0777);
    if(fd < 0 ) {
      printf("ERROR: %s cannot open flag file '%s'\n",argv[0],argv[1]);
      exit(1);
    }
    status = flock(fd,LOCK_EX);
    status = fstat(fd,&sbuf);
    t = time(NULL);
    tim = (t-sbuf.st_mtime);
    if( tim > delay ) {
//       printf("'%s'\n",argv[3]);
      write(fd,&tim,1);
      status = system(argv[3]);
      if(status){
        printf("ERROR: %s error in command'%s'\n",argv[0],argv[3]);
        exit(1);
      }
//     }else{
//       printf("waiting (%d) \n",tim);
    }
    status = flock(fd,LOCK_UN);
    close(fd);
    sleep(delay);
  }
}