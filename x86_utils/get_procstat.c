/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2015  Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

// struct timeval {
//    time_t      tv_sec;     /* seconds */
//    suseconds_t tv_usec;    /* microseconds */
// };

static long tickspersec;

#pragma weak get_my_procstat_ = get_my_procstat
#pragma weak get_my_procstat__ = get_my_procstat
int get_my_procstat_(int *verbose, long long int *pid, long long int *ppid, long long int *utime, long long int *stime, long long int *vsize, long long int *rss, long long int *rsslim);
int get_my_procstat_(int *verbose, long long int *pid, long long int *ppid, long long int *utime, long long int *stime, long long int *vsize, long long int *rss, long long int *rsslim);

int get_my_procstat(int *verbose, long long int *pid, long long int *ppid, long long int *utime, long long int *stime, long long int *vsize, long long int *rss, long long int *rsslim) {
  struct rusage usage;
  struct rlimit rlim;
  int status;

//   tickspersec = sysconf(_SC_CLK_TCK);
  tickspersec = 1000000;
  *rss = -1;  /* precondition outputs to miserable failure */
  *pid = getpid();
  *ppid = getppid();
  *utime = -1;
  *stime = -1;
  *vsize = 0;
  *rsslim = 0;

  status = getrlimit(RLIMIT_RSS,&rlim);
  *rsslim = rlim.rlim_cur;

  status = getrusage(RUSAGE_SELF,&usage);
  
  *rss = usage.ru_maxrss;
  *vsize = *rss * 1024;
  *rss = *rss / 4;

  *utime = usage.ru_utime.tv_sec;
  *utime = *utime * 1000000 + usage.ru_utime.tv_usec;

  *stime = usage.ru_stime.tv_sec;
  *stime = *stime * 1000000 + usage.ru_stime.tv_usec;
  

  return tickspersec;
}
#ifdef TEST
static void printone(char *name, long long int x) {  printf("%20s: %lld\n", name, x);}
static void printstr(char *name, char *x) {  printf("%20s: %s\n", name, x);}
static void printtime(char *name, long long int x) {  printf("%20s: %f\n", name, (((double)x) / tickspersec));}
main()
{
  int verbose = 0;
  int i;
  long scrap[10000000];
  long long int pid, ppid, utime, stime, vsize, rss, rsslim ;
  int tpr;

  for (i=0 ; i<10000000 ; i++) scrap[i]=i;
  tpr = get_my_procstat(&verbose, &pid, &ppid, &utime, &stime, &vsize, &rss, &rsslim);
  printstr("","------------------------------");
  printone("pid", pid);
  printone("ppid", ppid);
  printtime("utime", utime);
  printtime("stime", stime);
  printone("vsize", vsize);
  printone("rss", rss);
  printf("ticks per second: %d\n",tpr);
  printone("rsslim", rsslim);
  return 0;
}
#endif


