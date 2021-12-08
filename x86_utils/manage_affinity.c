/* RPN_COMM - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*
   this code is only applicable to Linux X86 for the time being

   Michel Valin , 2019
*/

#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static uint64_t time0 = 0;

uint64_t cpu_ticks_and_id(int *socket, int *core){  // get tsc/socket/core
#if defined(__x86_64__) &&  defined(__linux__)
   uint32_t lo, hi, c;
   // rdtscp instruction
   // EDX:EAX contain TimeStampCounter
   // ECX contains IA32_TSC_AUX[31:0] (MSR_TSC_AUX value set by OS, lower 32 bits contain socket+core)
   __asm__ volatile("rdtscp" : "=a" (lo), "=d" (hi), "=c" (c));
    *socket = (c>>12) & 0xFFF;
    *core   =  c      & 0xFFF;

   return time0 = (uint64_t)lo | (((uint64_t)hi) << 32) ;
#else
  time0   = cpu_real_time_clock();  // a tick will be a microsecond in this case
  *socket = -1;                     // no socket/core info available
  *core   = -1;
  return time0;
#endif
}


void rpn_comm_unbind_process(int nthreads)
{
#if defined(linux)

int atoi(const char *nptr);
char *getenv(const char *name);

cpu_set_t set;
int i;
int will_print=1;
int ncores=sysconf(_SC_NPROCESSORS_CONF);
int nbound=0;
// char *omp=getenv("OMP_NUM_THREADS");

CPU_ZERO(&set);
sched_getaffinity(0,sizeof(set),&set);
i=ncores;
while(--i >=0) { if (CPU_ISSET(i,&set)){ nbound++; } }  /* how many cores are we allowed to run on ? */
printf("nbound = %d\n",nbound);
i = CPU_COUNT(&set);
printf("CPUs in set = %d\n",i);

printf("_SC_PHYS_PAGES = %ld, _SC_AVPHYS_PAGES = %ld, _SC_PAGESIZE = %ld, POSIX2_VERSION = %ld\n",
       sysconf(_SC_PHYS_PAGES),
       sysconf(_SC_AVPHYS_PAGES),
       sysconf(_SC_PAGESIZE),
       sysconf(_SC_2_VERSION));
printf("MEM = %ld KBytes\n",sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGESIZE) / 1024 );
  
if(getenv("FULL_UNBIND") != NULL) nbound = 0;  /* FULL_UNBIND variable defined, unbind no matter what */
  
if(nthreads > nbound) {  /* need more threads than cores we can run on , unbind everything */
  if(will_print) printf("FULL unbinding will be done, cores=%d, threads needed=%d, usable cores=%d\n",ncores,nthreads,nbound);
  CPU_ZERO(&set);
  i=ncores;
  while(--i >=0) { CPU_SET(i,&set);}
  sched_setaffinity(0,sizeof(set),&set);  /* set affinity to all cores */
}else{
  if(will_print) printf("NO unbinding will be done\n");  /* enough resources available and no forced unbind */
}

#endif

return;
}

#include <mpi.h>
static cpu_set_t set;
static int lo_core=-1;
static int hi_core=-1;

int get_cpu_hyperthreads(){
  int fd=open("/sys/devices/system/cpu/smt/active",O_RDONLY);
  char ht = '0';
  read(fd,&ht,1l);
  close(fd);
  if(ht == '0') return 1;
  return 2;
}

void RebindBySocket(int init) {
  int rank, rank0, size, ncores, i, npersock;
  char *omp;
  MPI_Comm comm;
  int color;
  int get_cpu_cores();
  int get_cpu_hyperthreads();
  int *flag;

if(init){
  omp=getenv("OMP_NUM_THREADS");
  color=gethostid();
  color &= 0x7FFFFFFF;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank0);
  MPI_Comm_split(MPI_COMM_WORLD,color,rank0,&comm);
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  ncores=sysconf(_SC_NPROCESSORS_CONF);
  ncores /= get_cpu_hyperthreads();
  npersock=ncores/2;
  CPU_ZERO(&set);
  if(rank < size/2) {
    lo_core = 0 ;
    hi_core = npersock - 1 ;
  }else{
    lo_core = npersock ;
    hi_core = ncores - 1;
  }
  for(i=lo_core ; i < hi_core ;  i++) { CPU_SET(i,&set) ;}
}
  sched_setaffinity(0,sizeof(set),&set);
}

#include <mpi.h>

int main(int argc, char **argv){
  int nthreads = 1;
  int core, socket;
  uint64_t tmp;
  int ierr;
  int my_rank;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  rpn_comm_unbind_process(nthreads);
  tmp = cpu_ticks_and_id(&socket, &core);
  if(tmp) printf("rank = %3.3d, core = %3.3d, socket = %d\n",my_rank,core,socket);
  ierr = MPI_Finalize();
}
