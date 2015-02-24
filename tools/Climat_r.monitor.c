/*
Author: Michel Valin, RPN, 16 Januar 2009

This script is always running in the background and gets killed
the moment the parent job finishes.

It expects pairs of two input parameters:
Arg1: Name of file (including path)  to be checked for
Arg2: Command to execute

When the program finds the file 'Arg1' it will execute the command:

  'Arg2' 'Arg1'

Like this it will execute all pairs of arguments it receives.
*/

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <glob.h>

main(int argc, char **argv){
int i, fd;
char buffer[32768];
pid_t pp=getppid();
char *actual_file;
glob_t globtab;
int matches;

if(argc-1 & 1) { printf("argument count must be even \n"); exit(1); }
if(fork()) exit(0);
while(1){
  if(kill(pp,0)) exit(0);
  matches = 0;
  for(i=1 ; i<argc ; i+=2){
    if( 0 == glob(argv[i],GLOB_NOSORT,NULL,&globtab) ) {     /* pattern match successful */
      matches ++;
      actual_file = globtab.gl_pathv[0];                     /* use name of first match  */
      if( (fd=open(actual_file,O_RDONLY )) >= 0 ) {          /* try to open file for reading */
        close(fd);
        printf("file=%s,cmd=%s\n",actual_file,argv[i+1]);
        snprintf(buffer,sizeof(buffer)-1,"%s %s",argv[i+1],actual_file);
        printf("Executing:%s\n",buffer);
        system(buffer);                                      /* execute command giving it the file name as an argument */
        printf("Deleting:%s\n",actual_file);
        unlink(actual_file);                                 /* after command execution, get rid of file */
        }
      }
      globfree(&globtab);
    }
  if(matches ==0) sleep(2);                                  /* nothing found this time around, let's sleep a bit */
  }
}
