//
// Copyright (C) 2021  Environnement et Changement climatique Canada
//
// This is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// Author:
//     M. Valin,   Recherche en Prevision Numerique, 2021
//
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

// return total number of directories to create
// string may contain ndirs | first,last | first,step,last
// i = first value
// j = increment
// k = last value
int get_ndirs(char *s, int *i, int *j, int*k){
  int ii, jj, kk, n1, n2, nx, ny, n;

  *i = 0 ;   // precondition for failure
  *j = 0 ;
  *k = 0 ;
  n1 = sscanf(s,"%dx%d",&nx,&ny);          // nx/ny ?
  if(n1 == 2){
    if(ny > 1){
      *i = nx;
      *j = ny;
      *k = -1;
      return(nx*ny);
    }
  }
  n2 = sscanf(s,"%d,%d,%d",&ii,&jj,&kk);    // ndirs | first,last | first,step,last
  if(n2 == 1) {               // ndirs
    *i = 0 ;
    *j = 1 ;
    *k = ii - 1 ;
    return (ii) ;            // return 0, 1, last into i, j, k
  }
  if(n2 == 2) {               // first,last
    *i = ii ;
    *j = 1 ;
    *k = jj ;
    return (jj - ii + 1) ;   // return first, 1, last into i, j, k
  }
  if(n2 == 3) {               // first,step,last
    *i = ii ;
    *j = jj ;
    *k = kk ;
    n = 0 ;
    while(ii <= kk) {        // compute number of directories
      n++ ;
      ii += jj ;
    }
    return (n) ;             // return first, step, last into i, j, k
  }
  return 0 ;                 // invalid string
}

void usage(char *msg, int status){
  fprintf(stderr,"%s\n",msg);
  fprintf(stderr," usage: makedirs [directory_name] ndirs n [link_target link_name]\n");
  fprintf(stderr,"        makedirs [directory_name] n1xn2 n [link_target link_name]\n");
  fprintf(stderr," ndirs: number of directories with name nnnn to create (0 -> ndirs-1)\n");
  fprintf(stderr,"        ndirs may also use a first,increment,last or  first,last syntax\n");
  fprintf(stderr," n    : number of digits to use (1-9) for directory name elements\n");
  fprintf(stderr," n1,n2  directories with names nny-nnx (nnx = 0->n1-1 , nny = 0->n2-1) will be created\n");
  fprintf(stderr," n1 = number of processes per group, n2 = number of process groups\n");
  fprintf(stderr," if optional argument directory_name is an existing directory name, cd into that directory first\n");
  exit(status);
}

int main(int argc, char **argv){
  char dirname[32];            // directory name
  char format1[32];         // format for directory name
  char linkname[128];       // link name
  int i, ndirs, ndigits, nx, ny;
  int first, delta, last ;
  DIR *dirptr;

  if(argc == 1) usage("", 0) ;

  if(argc < 3) usage("Invalid number of arguments", 1);

  dirptr = opendir(argv[1]) ;                   // open directory
  if(NULL != dirptr){                           // directory argv[1] exists
    closedir(dirptr);                           // close directory
    i = chdir(argv[1]);                         // try cd argv[1]
    if(i != 0) usage("cd directory failed", 1);
    fprintf(stderr, "cd %s\n", argv[1]);
    argc--;                                     // skip first argument
    argv++;
  }

  if(argc != 3 && argc != 5) usage("Invalid number of arguments", 1);  // 3 or 5 are the only 2 acceptable values

  nx = 0 ; 
  ny = 0 ;
  // makedirs [ndirs | first,step,last | nx/ny]  ndigits   [link_target   link_name]
  ndirs   = get_ndirs(argv[1], &first, &delta, &last) ;
  ndigits = atoi(argv[2]);
  if(ndirs<=0 || ndigits <=0 || ndigits>9) usage("ndirs must be > 0 ; ndigits must be >0 and <10", 1);

  if(last == -1){     // nx/ny
    nx = first ;
    ny = delta ;
    first = 0 ;
    delta = 1 ;
    last = nx * ny - 1 ;
    sprintf(format1,"%c%d.%dd-%c%d.%dd",'%',ndigits,ndigits,'%',ndigits,ndigits); // format for directory name %N.Nd-%N.Nd
  }else{
    nx = ndirs ;
    ny = 1 ;
    sprintf(format1,"%c%d.%dd",'%',ndigits,ndigits);               // format for directory name %N.Nd
  }
  fprintf(stderr,"format1='%s', ndirs=%d, nx=%d, ny=%d\n",format1,ndirs,nx,ny);

  for (i=first ; i<=last ; i+=delta) { 
    if(ny > 1){
      sprintf(dirname,format1,i/nx,i%nx);               // create directory name
    }else{
      sprintf(dirname,format1,i);                       // create directory name
    }
    mkdir(dirname,0755);                                // create directory
    if(argc == 5){
      sprintf(linkname,"%s/%s",dirname,argv[4]) ;       // build soft link name
      symlink(argv[3],linkname);                        // ln -s link_target link_name
    }
  }

  return 0 ;
}
