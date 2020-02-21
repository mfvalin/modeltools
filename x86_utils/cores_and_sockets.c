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
static int core_affinity = 0;
static int lo_core=-1;
static int hi_core=-1;
static int8_t *flg;
static int8_t *nod;

int sched_getcpu(void) ;

// initialize internal tables
static void init_core_numa(){
  int i;
  int last_numa = -1;
  cores = sysconf(_SC_NPROCESSORS_CONF);                 // get number of cores on node
  (void)sched_getaffinity(0, sizeof(set), &set);         // get processor affinity mask
  flg = (int8_t *)malloc(sizeof(uint8_t) * cores);       // allocate flag table
  nod = (int8_t *)malloc(sizeof(uint8_t) * cores);       // allocate numa node table
  memset(flg,    0, cores);                              // core affinity
  memset(nod, 0xFF, cores);                              // numa node associated with allowed cores
  for (i = 0; i < cores ; i++) {       // translate affinity mask into flags
    if( CPU_ISSET(i, &set) ) {         // we can run on this core
      core_affinity++ ;                // number of cores this process can run on
      if(lo_core == -1) lo_core = i ;  // first allowed core encountered
      hi_core = i ;                    // last core allowed
      flg[i] = 1 ;
      nod[i] = numa_node_of_cpu(i);    // numa node this core belongs to
      if(nod[i] != last_numa) {        // new node number
	node_affinity++ ;              // number of numa nodes this process can run on (should be 1)
	last_numa = nod[i] ;
      }
    }
  }
}
//interface   !InTf!
// function GetCoreAffinity() result(p) bind(C,name='GetCoreAffinity')  !InTf!
//   import :: C_PTR                                     !InTf!
//   type(C_PTR) :: p                                    !InTf!
// end function GetCoreAffinity                          !InTf!

// get core (hyperthread) affinity table
int8_t *GetCoreAffinity(){
  if(cores == -1) init_core_numa();
  return flg;
}

// function GetNumaAffinity() result(p) bind(C,name='GetNumaAffinity')  !InTf!
//   import :: C_PTR                                     !InTf!
//   type(C_PTR) :: p                                    !InTf!
// end function GetNumaAffinity                          !InTf!

// get numa (socket) affinity table
int8_t *GetNumaAffinity(){
  if(cores == -1) init_core_numa();
  return nod;
}

// function CoresForProcess() result(n) bind(c,name='CoresForProcess')  !InTf!
//   import :: C_INT                                     !InTf!
//   integer(C_INT) :: n                                 !InTf!
// end function CoresForProcess                          !InTf!

// number of cores (and hyperthreads available to this process
int CoresForProcess(){
  if(cores == -1) init_core_numa();
  return core_affinity;   // number of cores this process can run on
}

// subroutine GetCurrentCoreAndNode(cpu, node) bind(C,name='GetCurrentCoreAndNode')  !InTf!
//   import :: C_INT                                     !InTf!
//   integer(C_INT), intent(OUT) :: cpu, node            !InTf!
// end subroutine GetCurrentCoreAndNode                  !InTf!

// return current core and numa node (socket in most cases)
void GetCurrentCoreAndNode(int *cpu, int *node){          // current core and numa node for this process
  *cpu  = sched_getcpu()          ; core     = *cpu ;     // code number (may be a hyperthread)
  *node = numa_node_of_cpu(core)  ; numanode = *node ;    // numa node
}

// function CoresOnNode() result(n) bind(c,name='CoresOnNode')  !InTf!
//   import :: C_INT                                     !InTf!
//   integer(C_INT) :: n                                 !InTf!
// end function CoresOnNode                              !InTf!

// return number of cores on node (including hyperthreads if any)
int32_t CoresOnNode(){
  if(cores == -1) init_core_numa();
  return cores;      // numer of cores
}
//end interface   !InTf!

/* Derived from util-linux-2.13-pre7/schedutils/taskset.c */
// convert affinity into printable string
// string is allocated internally and should be freed by user when no longer needed
// mask : affinity set returned by sched_getaffinity
char *cpuset_to_str(cpu_set_t *mask)    // cpu set to C string conversion
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
  int8_t *f, *n;

  printf("cores on node = %d \n",CoresOnNode());
  clbuf = cpuset_to_str(&set);
  printf("core affinity = %s (%d/%d)[%d], on %d numa node(s), ",clbuf, lo_core, hi_core, core_affinity, node_affinity);
  free(clbuf);
  f = GetCoreAffinity();
  n = GetNumaAffinity();
  printf("mask = "); for(i=0 ; i<cores ; i++) printf("%1d",f[i]) ; printf("\n");
  printf("numa node mask = "); for(i=0 ; i<cores ; i++) printf("%1.1x",n[i]&0xF) ; printf("\n");
  GetCurrentCoreAndNode(&mycore, &mynuma);
  printf("current core = %d in numa node %d\n",mycore, mynuma);
  return 0;
}
#endif

