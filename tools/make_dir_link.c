#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
void usage(char *msg){
  printf("%s\n",msg);
  printf(" usage: makedirs [directory_name] ndirs n [link_to link_name]\n");
  printf("        makedirs [directory_name] nx ny n [link_to link_name]\n");
  printf(" ndirs: number of directories with name nnnn to create (0 -> ndirs-1)\n");
  printf(" n    : number of digits to use (1-9) for directory name elements\n");
  printf(" nx,ny  directories with names in nnx-nny (nnx = 0->nx-1 , nny = 0->ny-1) will be created\n");
  printf(" nx = number of process groups, ny = number of processes per group\n");
  printf(" if optional argument directory_name is an existing directory name, go to that directory first\n");
  exit(1);
}
int main(int argc, char **argv){
  char name[32];
  char format[32];
  char format2[64];
  char target[128];
  int i, j, ndirs, nchars, nx, ny;
  DIR *dirptr;

  if(argc < 3) usage("Invalid number of arguments"); // not enough arguments

  dirptr = opendir(argv[1]) ;
  if(NULL != dirptr){   // directory argv[1] exists
    closedir(dirptr);
    i = chdir(argv[1]); // cd argv[1]
    if(i != 0) usage("cannot cd to directory"); // cannot cd to argv[1]
    printf("cd %s\n",argv[1]);
    argc--;             // skip argument
    argv++;
  }

  if(argc == 3) {   // makedirs ndirs n
    ndirs  = atoi(argv[1]);
    nchars = atoi(argv[2]);
    if(ndirs<=0 || nchars <=0 || nchars>9) usage("ndirs, n must be >0 and n must be <10");
    sprintf(format,"%c%d.%dd",'%',nchars,nchars);
    printf("format='%s', ndirs=%d\n",format,ndirs);
    for (i=0 ; i<ndirs ; i++) { sprintf(name,format,i); mkdir(name,0777); }
    exit(0);
  }

  if(argc == 4) {   // makedirs nx ny n
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nchars = atoi(argv[3]);
    if(nx<=0 || ny<=0 || nchars <=0 || nchars>9) usage("nx, ny, n must be >0 and n must be <10");
    sprintf(format,"%c%d.%dd-%c%d.%dd",'%',nchars,nchars,'%',nchars,nchars);
    printf("format='%s', nx=%d, ny=%d\n",format,nx,ny);
    for (i=0 ; i<nx ; i++) {
      for (j=0 ; j<ny ; j++) { sprintf(name,format,i,j); mkdir(name,0777); }
    }
    exit(0);
  }

  if(argc == 5) {   // makedirs ndirs n link_to link_name
    ndirs  = atoi(argv[1]);
    nchars = atoi(argv[2]);
    if(ndirs<=0 || nchars <=0 || nchars>9) usage("ndirs, n must be >0 and n must be <10");
    sprintf(format,"%c%d.%dd",'%',nchars,nchars);
    printf("format='%s', ndirs=%d\n",format,ndirs);
    sprintf(format2,"%c%d.%dd/%s",'%',nchars,nchars,"%s");
    for (i=0 ; i<ndirs ; i++) { 
      sprintf(name,format,i); 
      mkdir(name,0777);
      sprintf(target,format2,i,argv[4]) ;
      symlink(argv[3],target);
    }
    exit(0) ;
  }

  if(argc == 6) {   // makedirs nx ny n link_to link_name
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nchars = atoi(argv[3]);
    if(nx<=0 || ny<=0 || nchars <=0 || nchars>9) usage("nx, ny, n must be >0 and n must be <10");
    sprintf(format,"%c%d.%dd-%c%d.%dd",'%',nchars,nchars,'%',nchars,nchars);
    printf("format='%s', nx=%d, ny=%d\n",format,nx,ny);
    sprintf(format2,"%c%d.%dd-%c%d.%dd/%s",'%',nchars,nchars,'%',nchars,nchars,"%s");
    printf("format2='%s', nx=%d, ny=%d, link name = '%s', target = '%s'\n",format2,nx,ny,argv[5], argv[4]);
    for (i=0 ; i<nx ; i++) {
      for (j=0 ; j<ny ; j++) { 
        sprintf(name,format,i,j); 
        mkdir(name,0777);
        sprintf(target,format2,i,j,argv[5]) ;
        symlink(argv[4],target);
      }
    }
    exit(0) ;
  }
  usage("Invalid number of arguments");
  return 1 ;
}
