/* execve.c */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char **argv, char **envp) {

  char *newargv[4];
  char buffer[32768];
  int bp=0;
  int i;
  pid_t pid=getpid();
  int pidst;
  int maxwall;

  if (argc < 3) {
     fprintf(stderr, "Usage: %s <max wall time> <command> [arguments] \n", argv[0]);
     exit(EXIT_FAILURE);
  }
  maxwall=atoi(argv[1]);

  if( fork() ) {   /* parent process */

     /* the parent process will execve to bash the requested command */

    newargv[0]="/bin/bash";
    newargv[1]="-c";
    newargv[2]=buffer;
    newargv[3]=NULL;

    for (i=2 ; i<argc ; i++ ){bp+=sprintf(buffer+bp, "%s ",argv[i]);} //fprintf(stderr,"%s\n",buffer);

    execve("/bin/bash", newargv, envp);
    perror("execve");   /* execve() only returns on error */
    exit(EXIT_FAILURE);

  }else{   /* child process */

    /* the child watches the parent */
    while( (maxwall-->0) && ((pidst=kill(pid,0))== 0)) sleep(1)  ;  /* while there is time and parent is alive */

    if(maxwall<=0) { pidst=kill(pid,9) ; }   /* time is up, kill parent */
    if (pidst==0)
      fprintf(stderr,"process %d killed, time left=%d\n",pid,++maxwall);
    else
      fprintf(stderr,"process %d finished, time left=%d\n",pid,++maxwall);
    exit (0);
  }

}

