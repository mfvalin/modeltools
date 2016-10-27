#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
// usage: makedirs ndirs n
// ndirs: number of directories to create
// n    : number of digits to use (1-9)
main(int argc, char **argv){
  char name[12];
  char format[6];
  int i, ndirs,nchars;
  if(argc != 3) exit(1);
  ndirs = atoi(argv[1]);
  nchars = atoi(argv[2]);
  if(ndirs<=0 || nchars <=0 || nchars>9) exit(1);
  sprintf(format,"%c%d.%dd",'%',nchars,nchars);
//printf("format='%s', ndirs=%d\n",format,ndirs);
  for (i=0 ; i<ndirs ; i++) { sprintf(name,format,i); mkdir(name,0777); }
}
