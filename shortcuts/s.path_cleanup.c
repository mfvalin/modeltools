#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


static  char **eliminate=NULL;
static  int n_eliminate=0;

print_usage(char *name)
{
   fprintf(stderr,"USAGE: %s -pat1 -pat2 ... -patn [+supp1] item1 [+supp1] item2 ... [+suppn] itemn \n",name);
   fprintf(stderr,"       pat1..patn simple patterns to eliminate\n");
   fprintf(stderr,"       supp1..suppn if item/supp is a directory, add it in front of item\n");
   fprintf(stderr,"       item = EnvName \n");
   fprintf(stderr,"       EnvName is the name of an environment variable containing a path\n");
   fprintf(stderr,"ex:    %s PATH MANPATH EC_LD_LIBRARY_PATH\n",name);
   fprintf(stderr,"       %s +$BASE_ARCH PATH MANPATH +EC_ARCH EC_LD_LIBRARY_PATH\n",name);
   fprintf(stderr,"       %s LD_LIBRARY_PATH +EC_ARCH EC_INCLUDE_PATH\n",name);
   fprintf(stderr,"ex:    %s -/a/b/ -machin/patente PATH MANPATH EC_LD_LIBRARY_PATH\n",name);
}

int is_duplicate(char *what, char **table, int nentries)
{
  int i;
  struct stat stat_buf;

  for(i=0;i<n_eliminate;i++) { if ( strstr(what,1+eliminate[i]) ) return 1; }
  if( 0 == stat(what,&stat_buf)) {           /* what exists */
    if( S_ISDIR(stat_buf.st_mode) ) {        /* and is a directory */
      while(nentries>0){                     /* let's now find if it is a duplicate */
        if(table[nentries-1] != NULL) {
          if(strcmp(what,table[nentries-1]) == 0) return 1;
        }
        nentries--;
      }
      return 0;
    }
  }
  return 1;
}

main(int argc, char **argv)
{
  char *progname;
  char *env_var;
  char *supplement;
  char *sub_path[1204];
  char extra_path[1204];
  char env_value[32768];
  int nsubpath;
  int i;
  char separator=':';
  char *temp;
  char *TMPDIR=getenv("TMPDIR");
  char tmpdir_bin[1024];
  char *HOME=getenv("HOME");
  char ovbin[1024];
  int len_ovbin;
  char *varname;
  struct stat stat_buf;

  progname=argv[0];
  if(argc < 2 ) { print_usage(argv[0]) ; exit(1); }

  if(strcmp( "-h",argv[1])==0 || strcmp("--help",argv[1])==0 ) { print_usage(argv[0]) ; exit(1); }
  argc-- ; argv++;
  if(**argv == '-') {
    eliminate=argv;
    while(argc > 0 && **argv == '-') { n_eliminate++ ; argc-- ; argv++ ; }
  }
/* for ( i=0 ; i<n_eliminate ; i++ ) { fprintf(stderr,"eliminating '%s'\n",1+eliminate[i]); } */
  tmpdir_bin[0]='\0';  /* tempdir/bin, top override */
  if(TMPDIR != NULL) snprintf(tmpdir_bin,1023,"%s/bin%c",TMPDIR,'\0'); tmpdir_bin[1023]='\0';

  ovbin[0]='\0';  /* base of ovbin directories, ~/ovbin */
  if(HOME   != NULL) snprintf(ovbin     ,1023,"%s/ovbin%c",HOME,'\0'); ovbin[1023]='\0';
  len_ovbin=strlen(ovbin);

  supplement="";
  while(argc--){                /* loop over cEnvName arguments */
    varname=argv[0];
    if( *varname == '+' ) {
      supplement=varname+1;
      argv++;
      continue;
    }
    if( isalnum(*varname) || *varname=='_' ) {
      separator = '\0';
    }else{
      separator = *varname++;
    }
    env_var=getenv(varname);    /* get value of environment variable */
    if(env_var != NULL){        /* make sure it exists */
      strncpy(env_value,env_var,32767); env_value[32767]='\0';  /* take a copy */
      temp=env_value;
      if(separator == '\0')
        while(*temp != '\0' && *temp != ':' &&  *temp != ' ' && *temp != ';') temp++;  /* find separator, one of blank , colon, semicolon */
      separator=*temp;
      temp=env_value;
      while( *temp == separator && *temp != '\0') temp++;       /* get rid of leading separators */

      nsubpath=1;
      sub_path[0]=temp;

      while(1){                                                 /* split path into separate entries */
        while( *temp != separator && *temp != '\0') temp++;     /* skip to next separator */
        if(*temp == '\0') break;                                /* OOPS end */
        *temp='\0'; temp++;
        if(is_duplicate(sub_path[nsubpath-1],sub_path,nsubpath-1)) /* duplicate entry in path ? */
          nsubpath--;
        while( *temp == separator && *temp != '\0') temp++;     /* skip all consecutive separators */
        if(*temp == '\0') break;                                /* OOPS end */
        sub_path[nsubpath++]=temp;
      }
      if(nsubpath < 1) goto try_next;                            /* nothing was found */
      if(is_duplicate(sub_path[nsubpath-1],sub_path,nsubpath-1)) /* is last entry in path a duplicate ? */
        nsubpath--;
      if(nsubpath < 1) goto try_next;                            /* nothing was found */

      fprintf(stdout,"%s='",varname);                            /* variable_name= */

      for(i=0 ; i<nsubpath ; i++){
         if(0 == strcmp(sub_path[i],tmpdir_bin) ){               /* $TMPDIR/bin ? */
           fprintf(stdout,"%s%c",sub_path[i],separator);
           sub_path[i]=NULL;                                     /* yes, suppress entry */
         }
      }
      while(sub_path[nsubpath-1] == NULL) nsubpath--;            /* suppress NULLs at tail of list */
      for(i=0 ; i<nsubpath ; i++){
         if(sub_path[i] == NULL) continue;
         if(0 == strncmp(sub_path[i],ovbin,len_ovbin) ){         /* $HOME/ovbin... ? */
           fprintf(stdout,"%s%c",sub_path[i],separator);
           sub_path[i]=NULL;                                     /* yes, suppress entry */
         }
      }
      while(sub_path[nsubpath-1] == NULL) nsubpath--;            /* suppress NULLs at tail of list */
      for(i=0 ; i<nsubpath ; i++)                                /* output non override values */
        if(sub_path[i] != NULL) {
          if( *supplement != '\0' ) {
            snprintf(extra_path,1023,"%s/%s",sub_path[i],supplement);
            extra_path[1023]='\0';
            if ( ! is_duplicate(extra_path,sub_path,nsubpath ) ) {
              fprintf(stdout,"%s/%s%c",sub_path[i],supplement,separator);
            }
          }
          if( i == nsubpath-1 )
            fprintf(stdout,"%s'\n",sub_path[i]);
          else
            fprintf(stdout,"%s%c",sub_path[i],separator);
        }
    }else{   /* specified path environment variable does notexist, do nothing */
      fprintf(stderr,"%s: ERROR environment variable %s not found\n",progname,*argv);
    }
try_next:
    if(nsubpath == 0) {
      fprintf(stderr,"%s: ERROR no valid directory found in path variable %s\n",progname,*argv);
    }
    supplement="";
    argv++;
  }
}
