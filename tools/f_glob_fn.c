#if ! defined(FTEST)

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

/* following lines added for Portland Group C compiler */
#if defined(linux)
#ifndef __size_t
#define __size_t size_t
#endif
#endif

#include <glob.h>

/* modified strncpy that blank fills insteead of zero filling */
static char *strncpy_f(char *dest, const char *src, size_t n)
{
    size_t i;

   for (i = 0; i < n && src[i] != '\0'; i++)
        dest[i] = src[i];
    for ( ; i < n; i++)
        dest[i] = ' ';

   return dest;
}

static int in_use = -1;
static glob_t tglob = {0,NULL,0};

/* fill glob table using pattern name */
/* if called from FORTRAN, use trim(name)//achar(0) to make sure string is null terminated */
/* 0 is returned upon success, -1 if in use, glob status upon glob failure */
#pragma weak f_set_glob_=f_set_glob
#pragma weak f_set_glob__=f_set_glob
int f_set_glob(char *name)
{
  int status;

  if(in_use != -1) return -1;
  status = glob(name,0,NULL,&tglob);
  if(status==0) {
    in_use = 0;
  }else{
    globfree(&tglob);
    tglob.gl_pathc = 0;
    tglob.gl_offs = 0;
    in_use = -1;
  }
  return(status);
}

/* free glob table */
#pragma weak f_free_glob_=f_free_glob
#pragma weak f_free_glob__=f_free_glob
void f_free_glob()
{
  if(in_use == -1) return;
  globfree(&tglob);
  in_use = -1;
}

/* get next entry from glob table and store it into string name, *maxlen is size of name */
/* if called from C, make sure string will be null terminated (see test below) */
/* 0 is returned upon success, -1 upon failure */
#pragma weak f_next_glob_=f_next_glob
#pragma weak f_next_glob__=f_next_glob
int f_next_glob(char *name, int *maxlen)
{
  if(in_use == -1) return -1;
  if(in_use >= tglob.gl_pathc) return -1;
  strncpy_f(name,tglob.gl_pathv[in_use],*maxlen);
  in_use++;
  return 0;
}

/* is name an ordinary file ? */
#pragma weak name_is_a_file_=name_is_a_file
#pragma weak name_is_a_file__=name_is_a_file
int name_is_a_file(char *name)
{
  struct stat temp;
  int status;
  status = stat(name,&temp);
  if(status != 0) return(-1);
  if( S_IFREG & temp.st_mode ) return (0); /* it is a regular file */
  return (-1);                             /* not a regular file */
}

/* is name a text file ? */
/* if no character lower than \r is detected within the first 32 characters of a file, a text file is assumed */
#pragma weak name_is_txt_file_=name_is_txt_file
#pragma weak name_is_txt_file__=name_is_txt_file
int name_is_txt_file(char *name)
{
  int fd, status;
  unsigned char buffer[128];
  ssize_t nc;
  fd = open(name,0);  /* try to open file */
  if(fd < 0) return(-1);  /* failed */
  nc=read(fd,buffer,128);
  close(fd);
  if(nc<1) return(-1);
  nc--;
  while(nc>=0) 
  {
    if(buffer[nc] < 0x0a) return(-1);  /* control characters detected */
    nc--;
  }
  return (0);   /* is a text file */
}


/* is name a standard (RPN) file ?  (89 or 98 flavor) */
#define SIGN_STD89_RND  0x55555555
#pragma weak name_is_std_file_=name_is_std_file
#pragma weak name_is_std_file__=name_is_std_file
int name_is_std_file(char *nom) 
{
      FILE *pf;
      int buffer[4];
      char *cbuffer;
      int nc;

      pf = fopen(nom,"rb");
      if (pf == NULL) return(-1);
      cbuffer = (char *) &buffer[0];
      buffer[0]=0;
      buffer[1]=0;
      buffer[2]=0;
      buffer[3]=0;
      nc=fread(buffer,sizeof(int),4,pf);

     /* RANDOM89 */
      if (buffer[0] == SIGN_STD89_RND && buffer[1] == SIGN_STD89_RND)
	 return(0);
    /* STANDARD 98 RANDOM */
      if (cbuffer[8] == 'X' && cbuffer[9] == 'D' && cbuffer[12] == 'S' && cbuffer[10] == 'F' && 
          cbuffer[13] == 'T' && cbuffer[14] == 'D' && cbuffer[15] == 'R') 
             return(0);
      return(-1);  /* not recognized */
}

/* is name a directory ? */
#pragma weak name_is_a_dir_=name_is_a_dir
#pragma weak name_is_a_dir__=name_is_a_dir
int name_is_a_dir(char *name)
{
  struct stat temp;
  int status;
  status = stat(name,&temp);
  if(status != 0) return(-1);
  if( S_IFDIR & temp.st_mode ) return (0); /* it is a directory */
  return (-1);                             /* not a directory */
}

#if defined(TEST)
int main(int argc, char **argv)
{
  char buffer[33];
  int buflen=32;
  int status, size;
  char *name;
  
  buffer[32]='\0';
  status=f_set_glob(argv[1]);
  size = tglob.gl_pathc;
  printf("matches = %d, status=%d\n",size,status);
  while( (status=f_next_glob(buffer,&buflen)) == 0) {
    printf("%d:%d %s\n",status,in_use,buffer);
  }
  printf("%d:\n",status);
  f_free_glob();
  
  system("mkdir taratata");
  system("touch taratata2");

  name="taratata";
  status=name_is_a_dir(name);
  printf("%s is a dir = %d\n",name,status);
  status=name_is_a_file(name);
  printf("%s is a file = %d\n",name,status);
  
  name="taratata2";
  status=name_is_a_dir(name);
  printf("%s is a dir = %d\n",name,status);
  status=name_is_a_file(name);
  printf("%s is a file = %d\n",name,status);
  
  system("rmdir taratata");
  system("rm taratata2");
  
  return 0;
}
#endif
#endif
#if defined(FTEST) && ! defined(FULL)
program testglob
implicit none
character *32 buffer
character *32 target
integer f_set_glob, f_next_glob
external f_set_glob, f_free_glob, f_next_glob
integer status

target="raaffdfiudfh" ! look for non existent file
status=f_set_glob(trim(target)//achar(0))
print *,'status=',status,'for '//trim(target)

target="*.c"  ! look for C source files
status=f_set_glob(trim(target)//achar(0))
print *,'status=',status,'for '//trim(target)
do while(f_next_glob(buffer,32) == 0)
  print *,'"'//trim(buffer)//'"'
enddo

target="raaffdfiudfh"  ! in use, status should be -1
status=f_set_glob(trim(target)//achar(0))
print *,'status=',status,'for '//trim(target)

call f_free_glob()

stop
end
#endif
#if defined(FTEST) && defined(FULL)
program testglob
implicit none
character *32 buffer
character *32 target
integer f_set_glob, f_next_glob
external f_set_glob, f_free_glob, f_next_glob
integer name_is_a_dir, name_is_a_file
external name_is_a_dir, name_is_a_file
integer status
character *32 name

target="raaffdfiudfh" ! look for non existent file
status=f_set_glob(trim(target)//achar(0))
print *,'status=',status,'for '//trim(target)

target="*.[fc]"  ! look for .f or .c files
status=f_set_glob(trim(target)//achar(0))
print *,'status=',status,'for '//trim(target)
do while(f_next_glob(buffer,32) == 0)
  print *,'"'//trim(buffer)//'"'
enddo

target="raaffdfiudfh"  ! in use, status should be -1
status=f_set_glob(trim(target)//achar(0))
print *,'status=',status,'for '//trim(target)

call f_free_glob()

target="*.f??"  ! look for .ftn or .f90 files
status=f_set_glob(trim(target)//achar(0))
print *,'status=',status,'for '//trim(target)
do while(f_next_glob(buffer,32) == 0)
  print *,'"'//trim(buffer)//'"'
enddo
call f_free_glob()

call  system("mkdir taratata");
call  system("touch taratata2");

name="taratata";
status=name_is_a_dir(trim(name)//achar(0));
print *,name,' is a dir=',status
status=name_is_a_file(trim(name)//achar(0));
print *,name,' is a file=',status
name="taratata2";
status=name_is_a_dir(trim(name)//achar(0));
print *,name,' is a dir=',status
status=name_is_a_file(trim(name)//achar(0));
print *,name,' is a file=',status

call  system("rmdir taratata");
call  system("rm taratata2");


stop
end

#endif
