#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define MAX_APIC 512
static short *APICids = NULL;
static int Maxcpus = -1;
static int Maxapicid = -1;
static int make_apicid_table(int *maxcpus){
  char line[4096];
  short apics[MAX_APIC];
  char *ptr;
  int nc;
  int apicid, napics;
  int i;
  int lines = 0;
  int maxapicid = -1;
  int cpunumber = -1;
  FILE *fd;

  if(APICids != NULL){ // already initialized
#if defined(DEBUG)
    fprintf(stderr,"make_apicid_table has already been called, returning previous values\n");
#endif
    *maxcpus = Maxcpus;
    return Maxapicid;
  }
  fd=fopen("/proc/cpuinfo","r");
  for(i=0 ; i<MAX_APIC ; i++) apics[i] = -1;
  if(fd == NULL) return 0;
  ptr = fgets(line,sizeof(line),fd);
  while(ptr != NULL) {
    lines++;
    ptr = fgets(line,sizeof(line),fd);
    if(strncmp(line,"apicid",6) == 0) {
      cpunumber++;
      while( *ptr != ':') ptr++; ptr++;
      apicid = atoi(ptr);
      maxapicid = (maxapicid < apicid) ? apicid : maxapicid;
      if(apicid < MAX_APIC) apics[apicid] = cpunumber;
    }
  }
  fclose(fd);
  napics = maxapicid + 1 ;
  APICids = malloc(napics*sizeof(short));
  for(i=0 ; i<=maxapicid ; i++) APICids[i] = apics[i];
#if defined(DEBUG)
  fprintf(stderr,"lines read = %d, number of cpus = %d\n",lines,cpunumber+1);
  fprintf(stderr,"allocated %d entry table for APIC ids\n",napics);
#endif
  *maxcpus = cpunumber+1;
  Maxcpus = cpunumber+1;
  Maxapicid = maxapicid;
  return maxapicid;
}
#if defined(SELF_TEST)
main(){
//  short apic_table[MAX_APIC];
  short *cpu_table;
  int maxcpus;
  int i;
  int maxapic;
//  int napics;

//  for(i=0 ; i<MAX_APIC ; i++) apic_table[i] = -1;
  maxapic = make_apicid_table(&maxcpus);
  maxapic = make_apicid_table(&maxcpus);
  maxapic = make_apicid_table(&maxcpus);
  fprintf(stderr,"maxcpus = %d\n",maxcpus);
  cpu_table = malloc(maxcpus*sizeof(short));
  for(i=0 ; i<maxcpus ; i++) cpu_table[i] = -1;
//  napics = maxapic + 1 ;
//  APICids = malloc(napics*sizeof(short));
//  fprintf(stderr,"allocated %d entry table for APIC ids\n",napics);
//  for(i=0 ; i<=maxapic ; i++) APICids[i] = apic_table[i];

  fprintf(stderr,"APICID -> CPU table\n");
  for(i=0 ; i<=maxapic ; i++) fprintf(stderr,"%4hd",APICids[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"CPU -> APICID table\n");
  for(i=0 ; i<=maxapic ; i++) if(APICids[i]>=0) cpu_table[APICids[i]] = i;
  for(i=0 ; i<maxcpus ; i++) if(cpu_table[i]>=0) fprintf(stderr,"%4hd",cpu_table[i]);
  fprintf(stderr,"\n");
  return(0);
}
#endif
