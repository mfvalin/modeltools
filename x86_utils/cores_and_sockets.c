//     functions for C and FORTRAN programming
//     Copyright (C) 2020  Recherche en Prevision Numerique
// 
//     This software is free software; you can redistribute it and/or
//     modify it under the terms of the GNU Lesser General Public
//     License as published by the Free Software Foundation,
//     version 2.1 of the License.
// 
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//     Lesser General Public License for more details.
//
#define _GNU_SOURCE
#include<stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <sched.h>
#include <numa.h>

static cpu_set_t set;
static int cores = -1;
static int numanode = -1;
static int core = -1;
static int node_affinity = 0;
static int lo_core=-1;
static int hi_core=-1;
static uint8_t *flg;
static int8_t *nod;

int sched_getcpu(void) ;

static void init_core_numa(){
  int i;
  int last_numa = -1;
  cores = sysconf(_SC_NPROCESSORS_CONF);                 // get number of cores on node
  (void)sched_getaffinity(0, sizeof(set), &set);         // get processor affinity mask
  flg = (uint8_t *)malloc(sizeof(uint8_t) * cores);      // allocate flag table
  nod = (int8_t *)malloc(sizeof(uint8_t) * cores);      // allocate numa node table
  memset(flg, 0, cores);
  memset(nod, 0xFF, cores);
  for (i = 0; i < cores ; i++) { // translate affinity mask into flags
    if( CPU_ISSET(i, &set) ) {
      if(lo_core == -1) lo_core = i ;
      hi_core = i ;
      flg[i] = 1 ;
      nod[i] = numa_node_of_cpu(i);
      if(nod[i] != last_numa) { node_affinity++ ; last_numa = nod[i] ; }
    }
  }
}

void current_core_and_node(int *cpu, int *node){
  *cpu  = sched_getcpu()          ; core     = *cpu ;
  *node = numa_node_of_cpu(core)  ; numanode = *node ;
}

int32_t cores_on_node(){
  if(cores == -1) init_core_numa();
  return cores;
}

/* Derived from util-linux-2.13-pre7/schedutils/taskset.c */
char *cpuset_to_str(cpu_set_t *mask)    // cpu set to string conversion
{
  char *str = (char *) malloc(7 * CPU_SETSIZE);  // safe allocation, 7 chars per core
  char *ptr = str;
  int i, j, entry_made = 0;
  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, mask)) {
      int run = 0;
      entry_made = 1;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
      }
      if (!run)
        sprintf(ptr, "%d,", i);
      else if (run == 1) {
        sprintf(ptr, "%d,%d,", i, i + 1);
        i++;
      } else {
        sprintf(ptr, "%d-%d,", i, i + run);
        i += run;
      }
      while (*ptr != 0) ptr++;
    }
  }
  ptr -= entry_made;
  *ptr = 0;
  return str;
}
#if defined(SELF_TEST)
int main(int argc, char **argv){
  char *clbuf;
  int i;
  int mycore, mynuma;

  printf("cores on node = %d \n",cores_on_node());
  clbuf = cpuset_to_str(&set);
  printf("core affinity = %s (%d/%d), on %d numa node(s), ",clbuf, lo_core, hi_core, node_affinity);
  free(clbuf);
  printf("mask = "); for(i=0 ; i<cores ; i++) printf("%1d",flg[i]) ; printf("\n");
  printf("numa node mask = "); for(i=0 ; i<cores ; i++) printf("%1.1x",nod[i]&0xF) ; printf("\n");
  current_core_and_node(&mycore, &mynuma);
  printf("current core = %d in numa node %d\n",mycore, mynuma);
  return 0;
}
#endif

