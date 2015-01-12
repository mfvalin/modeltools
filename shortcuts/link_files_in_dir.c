#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>

main(int argc, char **argv){

int file_count = 0;
DIR * dirp;
struct dirent * entry;
char old_path[4096], new_path[4096];
int status;
char *so_string=".so";
int look_for_so=0;

if(strcmp("-so",argv[1]) == 0) {
  argc--;
  argv++;
  look_for_so=1;
}
if(argc!=3)  {
  fprintf(stderr,"Usage: %s path_to_scan pat_to_link_into\n",argv[0]);
  exit(1);
}

// fprintf(stderr,"Scanning %s\n",argv[1]);
dirp = opendir(argv[1]); /* There should be error handling after this */
while ((entry = readdir(dirp)) != NULL) {

#ifdef __linux
    if (entry->d_type == DT_DIR) { /* If the entry is a directory */
      continue;
    }
#endif
    snprintf(old_path,sizeof(old_path),"%s/%s",argv[1],entry->d_name);
    snprintf(new_path,sizeof(new_path),"%s/%s",argv[2],entry->d_name);

//    if(NULL != strstr(entry->d_name,so_string)) fprintf(stderr,".so file located\n");
    if(look_for_so & (NULL == strstr(entry->d_name,so_string))) continue;

    status = symlink(old_path,new_path);
//    fprintf(stderr,"%d ln -s %s %s\n",status,old_path,new_path);
    if(status==0)file_count++;
}
fprintf(stderr,"Linked %d files from %s into %s\n",file_count,argv[1],argv[2]);
closedir(dirp);
}
