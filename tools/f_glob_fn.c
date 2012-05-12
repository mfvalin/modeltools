#if ! defined(FTEST)

#include <stdlib.h>
#include <stdio.h>

/* following 3 lines added for Portland Group C compiler */
#ifndef __size_t
#define __size_t size_t
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
#if defined(TEST)
int main(int argc, char **argv)
{
  char buffer[33];
  int buflen=32;
  int status, size;
  
  buffer[32]='\0';
  status=f_set_glob(argv[1]);
  size = tglob.gl_pathc
  printf("matches = %d, status=%d\n",size,status);
  while( (status=f_next_glob(buffer,&buflen)) == 0) {
    printf("%d:%d %s\n",status,in_use,buffer);
  }
  printf("%d:\n",status);
  f_free_glob();
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
integer status

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

stop
end

#endif
