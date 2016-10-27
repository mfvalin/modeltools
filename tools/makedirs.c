#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
// usage: makedirs ndirs n
//        makedirs nx ny n
// ndirs: number of directories to create
// n    : number of digits to use (1-9)
// nx,ny  directories with names in nnx-nny (nnx = 0->nx-1 , nny = 0->ny-1) will be created
main(int argc, char **argv){
  char name[32];
  char format[32];
  int i, j, ndirs, nchars, nx, ny;
  if(argc == 3) {
    ndirs = atoi(argv[1]);
    nchars = atoi(argv[2]);
    if(ndirs<=0 || nchars <=0 || nchars>9) exit(1);
    sprintf(format,"%c%d.%dd",'%',nchars,nchars);
    printf("format='%s', ndirs=%d\n",format,ndirs);
    for (i=0 ; i<ndirs ; i++) { sprintf(name,format,i); mkdir(name,0777); }
    exit(0);
  }
  if(argc == 4) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nchars = atoi(argv[3]);
    if(nx<=0 || ny<=0 || nchars <=0 || nchars>9) exit(1);
    sprintf(format,"%c%d.%dd-%c%d.%dd",'%',nchars,nchars,'%',nchars,nchars);
    printf("format='%s', nx=%d, ny=%d\n",format,nx,ny);
    for (i=0 ; i<nx ; i++) {
      for (j=0 ; j<ny ; j++) { sprintf(name,format,i,j); mkdir(name,0777); }
    }
    exit(0);
  }
  exit (1);
}
