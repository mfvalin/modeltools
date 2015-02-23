#ifdef linux
#define _GNU_SOURCE
#define _LARGEFILE64_SOURCE
#define ioroutines64
#endif

/*
   jio.c

   Original Author: John Hague, ECMWF, 09-Jul-2010
   Thanks to Paul Burton and Will Weir for useful suggestions. 

   Routines below intercept C calls fopen, fclose, fread, fwrite, open, 
   close, read, write, and fgets. (No Fortran calls intercepted).

   At runtime:
   "export JIO_SUMMARY summarises overall I/O activity.
   "export JIO_DETAIL" plus I/O detail for each file.
   "export JIO_TRACE"  plus traces every I/O.

   Revision 001 : M.Valin Environment Canada 06-Dec-2010
          -added lseek/fseek/lseek64/open64/fopen64/kread intercepts
           (lseek64/open64/fopen64 for Linux only)
          -modified/tested with PGI Fortran compiler Linux 32/64 bits
          -gave most functions the static attribute
          -default made JIO_SUMMARY if no JIO_ENV environment variable is found
          -Linux version single thread at this time ( local omp_get_thread_num()
           always returns 0)
          -if the JIO_VERBOSE environment variable is defined, extra debug info is produced
          -can now append messages to a specific file (rather than stderr) via env variable
           export JIO_DIAG=file_name
   Revision 002 : M.Valin Environment Canada 10-Dec-2010
          -added usage of rdtsc instruction on x86 together with affinity setting 
           export JIO_RDTSC=cpuno
          -refactored code
*/

#define MAXF 1024  /* Default max filenames */
#define AVGL 100  /* Avg length of filenames */
#define NUME 16   /* Number of intercepted I/O routines */
#define MAXD 1024  /* Default Max File Desciptors */

#undef MIN
#define MIN(a,b) ( (a) < (b) ? (a) :  (b) )

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <sys/time.h>
#include <math.h>
#include <stdarg.h>

#ifdef linux
#include <sched.h>
#include <sys/syscall.h>
#endif

#ifdef _AIX
ssize_t kread(int , void *, size_t );
ssize_t kwrite___(int , const void *, size_t );
#include <omp.h>
extern int omp_get_thread_num();
/*
extern long long int irtc();
extern long long int get_thread_id_();
*/
#else
// threads ignored if not under AIX for the time being
static int omp_get_thread_num(){
return 0;
}
#endif

static int use_rdtsc=0;
static int iname = 0;
static double timj[NUME] ;
static double bytes[NUME];
static double mbps[NUME] ;
static long long int calls[NUME] ;
static char *namelist;
static int  *nameoff; 
static int  *calllist;
static double *timelist;
static double *bytelist;
static double *bytelsqr;
static int *indx;
static long long int jstart;  /* timestamp of first call */
static int jio_num=0;         /* diagnostic level 0,1,2,3 */
static int maxf = MAXF;
static int maxl = MAXF*AVGL;
static int maxd = MAXD;
static int jio_first = 0;   /* package is not initialized */
static int *filedescriptors;
static int verbose_mode=0;  /* extra debug flag */
static FILE *stdmsg;        /* normally standard error, can be redirected via JIO_DIAG env var */
static int stdmsg_fileno=2; /* file descriptor associated with stdmsg */
static float time_factor = 1.0;  /* time to true nanoseconds correction (1.0 / clock_in_GHz for X86, 1.0 otherwise) */

static  char XXX[NUME][8] = { "fopen: ","fclose:","fread: ","fwrite:",
                              "open:  ","close: ","read:  ","write: ",
                              "fgets: ","fseek: ","lseek: ","fopen64",
                              "open64:","lseek64","kread: ","kwrite:" };

typedef FILE * (*fopen_fn) (const char *, const char *);
static fopen_fn S_fopen_fn_ptr = fopen;
typedef int (*fclose_fn) (FILE *);
static fclose_fn S_fclose_fn_ptr = fclose;
typedef size_t (*fread_fn) (void *, size_t, size_t, FILE *);
static fread_fn S_fread_fn_ptr = fread;
typedef size_t (*fwrite_fn) (const void *, size_t, size_t, FILE *);
static fwrite_fn S_fwrite_fn_ptr = fwrite;
typedef int  (*open_fn) (const char *, int, ...);
static open_fn S_open_fn_ptr = open;
typedef int  (*close_fn) (int);
static close_fn S_close_fn_ptr = close;
typedef ssize_t  (*read_fn) (int, void *, size_t);
static read_fn S_read_fn_ptr = read;
typedef ssize_t  (*write_fn) (int, const void *, size_t);
static write_fn S_write_fn_ptr = write;
typedef char * (*fgets_fn) (char *, int, FILE *);
static fgets_fn S_fgets_fn_ptr = fgets;
typedef int  (*fseek_fn) (FILE *, long, int);
static fseek_fn S_fseek_fn_ptr = fseek;
typedef off_t  (*lseek_fn) (int, const off_t, int);
static lseek_fn S_lseek_fn_ptr = lseek;
#ifdef ioroutines64
typedef int  (*open64_fn) (const char *, int, ...);
static open64_fn S_open64_fn_ptr = open64;
typedef FILE * (*fopen64_fn) (const char *, const char *);
static fopen64_fn S_fopen64_fn_ptr = fopen64;
typedef off64_t  (*lseek64_fn) (int, const off64_t, int);
static lseek64_fn S_lseek64_fn_ptr = lseek64;
#endif
#ifdef _AIX
typedef ssize_t  (*kread_fn) (int, void *, size_t);
static kread_fn S_kread_fn_ptr = kread;
typedef ssize_t  (*kwrite_fn) (int, const void *, size_t);
static kwrite_fn S_kwrite_fn_ptr = kwrite___;
#endif

static void *my_ptr_init(char *text)
{
  void *addr = dlsym(RTLD_NEXT, text) ;
  if(verbose_mode)fprintf(stdmsg,"address of pointer to %s = %lx\n",text,addr);
  return addr;
}

void jio_init_extern()
{
  S_fopen_fn_ptr = (fopen_fn) my_ptr_init("fopen");
  S_fclose_fn_ptr = (fclose_fn) my_ptr_init("fclose");
  S_fread_fn_ptr = (fread_fn) my_ptr_init("fread");
  S_fseek_fn_ptr = (fseek_fn) my_ptr_init("fseek");
  S_fgets_fn_ptr = (fgets_fn) my_ptr_init("fgets");
  S_fwrite_fn_ptr = (fwrite_fn) my_ptr_init("fwrite");
  S_open_fn_ptr = (open_fn) my_ptr_init("open");
  S_close_fn_ptr = (close_fn) my_ptr_init("close");
  S_read_fn_ptr = (read_fn) my_ptr_init("read");
  S_write_fn_ptr = (write_fn) my_ptr_init("write");
  S_lseek_fn_ptr = (lseek_fn) my_ptr_init("lseek");
#ifdef ioroutines64
  S_fopen64_fn_ptr = (fopen64_fn) my_ptr_init("fopen64");
  S_open64_fn_ptr = (open64_fn) my_ptr_init("open64");
  S_lseek64_fn_ptr = (lseek64_fn) my_ptr_init("lseek64");
#endif
#ifdef _AIX
  S_kread_fn_ptr = (kread_fn) my_ptr_init("kread");
  S_kwrite_fn_ptr = (kwrite_fn) my_ptr_init("kwrite");
#endif
}

/* 
  jio_time function uses high res clock for IBM AIX and gettimeofday elsewhere
  rdtsc instruction under investigation for X86 platforms (coherence issues)
*/
#ifdef _AIX
static long long int jio_time()
{
  timebasestruct_t t_start;
  int              sec, n_sec;
  long long int    elapsed;
  read_real_time(&t_start, TIMEBASE_SZ);
  time_base_to_time(&t_start, TIMEBASE_SZ);
  sec   = t_start.tb_high;
  n_sec = t_start.tb_low; 
  elapsed = (long long int) sec*1000000000 + (long long int)n_sec;
  return elapsed;
}
#else
static long long int jio_time2()
{
  static long long last=0;  /* kept in nanoseconds, even if gettimeofday has at best microsec resolution */
  long long int elapsed;
  struct timeval tp;
  gettimeofday(&tp, NULL);
  elapsed=tp.tv_sec;
  elapsed=elapsed*1000000000;
  elapsed=elapsed+tp.tv_usec*1000;
  if (last == 0 ) { last = elapsed ;}  /* initialize last */
  if (elapsed<=last) elapsed=++last;   /* never return same value twice */
  return elapsed;  /* time in nanoseconds */
}
static unsigned long long jio_time1( void )
{
    unsigned long lo, hi;
    unsigned long long result;
    asm( "rdtsc" : "=a" (lo), "=d" (hi) ); 
    result = hi;
    result = lo | (result<<32);
    return( result );
}
static unsigned long long jio_time( void )
{
  if(use_rdtsc) return jio_time1();
  return jio_time2();
}
#endif

static void my_ptr_print( void *addr, char *text)
{
  if(verbose_mode)fprintf(stdmsg,"address of pointer to %s = %lx\n",text,addr);
}

/* this routine produces the JIO1 and JIO2 trace messages */
static void jio_trace(long long nsec, long long jstmp, const char * cc, int file, int len, const char *name) 
{
  int it;
  char buffer[1024];
  int buflen;

  it=omp_get_thread_num();
  buflen=sprintf(buffer,"JIO%s Trace ",name!=NULL?"2":"1");  /* JIO1 if no filename, JIO2 if there is one */
  buflen+=sprintf(buffer+buflen,"%s thrd=%2d, file =%3d, bytes =%10d",cc,it,file,len);
  buflen+=sprintf(buffer+buflen,", time =%12Ld nsec",nsec);
  buflen+=sprintf(buffer+buflen,", stamp=%20Ld",jstmp);
  if(name != NULL) buflen+=sprintf(buffer+buflen,", file= %s",name);
  buflen+=sprintf(buffer+buflen,"\n\0",name);
//  (*S_write_fn_ptr)(2,buffer,buflen);
  fprintf(stdmsg,"%s",buffer);
}

/* called at open time to open a new statistics session */
static void jio_detail_new(const char * name, int file)
{
  if(file < maxd) {
    if(indx[file] != 0) {
      if(verbose_mode)
        fprintf(stdmsg,"jio_detail_new: file=%d , oldname='%s' , newname='%s, indx = %d'\n",
                file,&namelist[nameoff[indx[file]]],name,indx[file]);
      if(strncmp(&namelist[nameoff[iname]],name,strlen(name))==0) {
        return ;  /* descriptor already in use with same file name */
      }
    }
  }
  iname=iname+1;
  if(verbose_mode)fprintf(stdmsg,"jio_detail_new: '%s' , %d, iname = %d \n",name,file,iname);
  /*--- Increasing Max files if necessary----*/
  if(iname >= maxf-1) {
    int ii;
    int    *nameoff1;
    int    *filedescriptors1;
    int    *calllist1;
    double *timelist1;
    double *bytelist1;
    double *bytelsqr1;
    fprintf(stdmsg,"JIO Detail increasing Max Number of Files to %d\n",2*maxf);
    nameoff1  = malloc(4*maxf);
    filedescriptors1  = malloc(4*maxf);
    calllist1 = malloc(4*maxf*NUME);
    timelist1 = malloc(8*maxf*NUME);
    bytelist1 = malloc(8*maxf*NUME);
    bytelsqr1 = malloc(8*maxf*NUME);
    for(ii=0; ii<maxf;      ii++) nameoff1[ii]=nameoff[ii];
    for(ii=0; ii<maxf;      ii++) filedescriptors1[ii]=filedescriptors[ii];
    for(ii=0; ii<maxf*NUME; ii++) {
      calllist1[ii] = calllist[ii];
      timelist1[ii] = timelist[ii];
      bytelist1[ii] = bytelist[ii];
      bytelsqr1[ii] = bytelsqr[ii];
    }
    free(nameoff);
    free(filedescriptors);
    free(calllist);
    free(timelist);
    free(bytelist);
    free(bytelsqr);
    nameoff  = malloc(4*2*maxf);
    filedescriptors  = malloc(4*2*maxf);
    calllist = malloc(4*2*maxf*NUME);
    timelist = malloc(8*2*maxf*NUME);
    bytelist = malloc(8*2*maxf*NUME);
    bytelsqr = malloc(8*2*maxf*NUME);
    for(ii=0; ii<maxf;      ii++) nameoff[ii]=nameoff1[ii];
    for(ii=0; ii<maxf;      ii++) filedescriptors[ii]=filedescriptors1[ii];
    for(ii=0; ii<maxf*NUME; ii++) {
      calllist[ii] = calllist1[ii];
      timelist[ii] = timelist1[ii];
      bytelist[ii] = bytelist1[ii];
      bytelsqr[ii] = bytelsqr1[ii];
    }
    for(ii=maxf;      ii<2*maxf;      ii++) nameoff[ii]=0;
    for(ii=maxf*NUME; ii<2*maxf*NUME; ii++) {
      calllist[ii] = 0;
      timelist[ii] = 0;
      bytelist[ii] = 0;
      bytelsqr[ii] = 0;
    }
    free(nameoff1);
    free(filedescriptors1);
    free(calllist1);
    free(timelist1);
    free(bytelist1);
    free(bytelsqr1);
    maxf=2*maxf;
  }
  /*--- Done Increasing Max files if necessary----*/

  size_t len = strlen(name);

  /*--- Increasing Max File Space if necessary----*/
  while ( nameoff[iname]+len+1 >= maxl) {
    int ii;
    char *namelist1;
    fprintf(stdmsg,"JIO Detail increasing Max File Space to %d bytes\n",2*maxl);
    namelist1 = malloc(maxl);
    for(ii=0; ii<maxl; ii++) namelist1[ii]=namelist[ii];
    free(namelist);
    namelist = malloc(2*maxl);
    for(ii=0; ii<maxl; ii++) namelist[ii]=namelist1[ii];
    for(ii=maxl; ii<2*maxl-1; ii++) namelist[ii]=(char)' ';
    free(namelist1);
    maxl=2*maxl;
  }
  /*--- Done Increasing Max File Space if necessary----*/

  strcpy(&namelist[nameoff[iname]],name);
  nameoff[iname+1]=nameoff[iname]+len+1;

  /*--- Incresase max File Descriptors if Necessary---*/ 
  while ( file >= maxd) {
    int i;
    int *indx1 = malloc(4*maxd);
    for (i=0; i<maxd; i++) indx1[i]=indx[i];
    free(indx);
    fprintf(stdmsg,"JIO Detail increasing size of File Desciptor Table to %d\n",maxd+100);
    indx=malloc(4*(maxd+100));
    for (i=0; i<maxd; i++) indx[i]=indx1[i];
    for (i=maxd; i<maxd+100; i++) indx[i]=indx1[i];
    free(indx1);
    maxd=maxd+100;
   }
  indx[file]=iname;
  /*--- Done Increase max File Descriptors if Necessary---*/ 
}

/*
  add information to statistics, SUMMARY, DETAILED. calls trace messages if appropriate
  argument 1 controls the level of detail
  1 summary
  2 detailed stats
  3 full trace messages
  name may be NULL if no filename present (any call other than an open)
  num if function number ( open, close, seek, etc ...) (see XXX table)
  file if file descriptor
  len is number of bytes (0 if not applicable as in open/close/seek)
*/
static void jio_add_info(long long t0, int jio_num, int num, int file, int len, const char *name)
{
  int ii;
  long long jstmp, nsec;

  if(jio_num == 0 ) return;   /* JIO_NONE */

  jstmp=t0;
  nsec =jio_time()-t0;
  calls[num]=calls[num]+1;
  timj[num]=timj[num]+(double)nsec;
  bytes[num]=bytes[num]+(double)len;
  if(jio_num < 2) return ;    /* JIO_SUMMARY , done */

  if(name != NULL) {
    jio_detail_new( name, file);                           /* JIO_DETAILED */
  }
  if(jio_num > 2) jio_trace(nsec,jstmp,&XXX[num][0],file,len,name);   /* JIO_TRACE active */
  ii=indx[file];                                           /* JIO_DETAILED follows */
  if(ii <= 0 ) return;
  if(ii < maxf) {
    filedescriptors[ii] = file; 
    calllist[ii*NUME+num]=calllist[ii*NUME+num]+1;
    timelist[ii*NUME+num]=timelist[ii*NUME+num]+(double)nsec;
    bytelist[ii*NUME+num]=bytelist[ii*NUME+num]+(double)len;
    bytelsqr[ii*NUME+num]=bytelsqr[ii*NUME+num]+(double)len * (double)len;
  }
}

/* print accumulated statistics */
/* "There are three kinds of lies: lies, damned lies, and statistics." */
static void jio_final()
{
  int i, ilist;
  double tbytes=0;
  double ttime=0;
  long long int tcalls=0;
  double tmbps;
  float Tbytes, Bytes, Ttime, Tmbps, Timj, Mbps;
  int ncalls;
  int filedescriptor;
  int j;
  int details_present;
  double avg,stdev;
  int jio_num_=jio_num;

  /* turn tracing off as we are printing final report */
  jio_num=0;
  verbose_mode=0;

  if(jio_num_ == 0) return;

  if(jio_num_ >= 1 )
  {
    fprintf(stdmsg,"JIO Summary Routine   Calls          MB        MSEC       MB/s\n");
    for(i=0; i<NUME; i++) {
      timj[i]=0.000001*timj[i]*time_factor;
      bytes[i]=bytes[i]/(1024.0*1024.0); 
      if( timj[i] > 0.0 ) mbps[i]=bytes[i]/(0.001*timj[i]);
      ttime=ttime+timj[i];
      tbytes=tbytes+bytes[i];
      tcalls=tcalls+calls[i];
      Bytes=bytes[i] ; Timj=timj[i] ; Mbps=mbps[i]; ncalls=calls[i];
//      fprintf(stdmsg,"JIO Summary %s  %6d  %10.3f  %10.3f %10.3f\n",&XXX[i][0],calls[i],bytes[i],timj[i],mbps[i]);
      if( calls[i] > 0 ) fprintf(stdmsg,"JIO Summary %s  %6Ld  %10.3f  %10.3f %10.3f\n",&XXX[i][0],calls[i],bytes[i],timj[i],mbps[i]);
    }
    tmbps = 0.0;
    if( ttime > 0.0 ) tmbps=tbytes/(0.001*ttime);
//    fprintf(stdmsg,"JIO Summary TOTAL:   %6d  %10.3f  %10.3f %10.3f\n",tcalls,tbytes,ttime,tmbps);
    Tbytes=tbytes ; Ttime=ttime ; Tmbps=tmbps;
    fprintf(stdmsg,"JIO Summary TOTAL:   %6Ld  %10.3f  %10.3f %10.3f\n",tcalls,Tbytes,Ttime,Tmbps);
    Ttime = (jio_time() - jstart) * time_factor;
    Ttime *= .001; Ttime *= .001; Ttime *= .001;
    fprintf(stdmsg,"JIO Summary Application Total Wall Clock time: %10.3f seconds\n",Ttime);
  }

  if (jio_num_ >= 2) {
    for (ilist=1; ilist<=MIN(iname,maxf-1); ilist++) {
      details_present=0;
      for(i=0; i<NUME; i++) {
         if( calllist[ilist*NUME+i] > 0) details_present++;
      }
      if(details_present == 0 && filedescriptors[ilist] < 0 ) continue;
      fprintf(stdmsg,"------------------------------%4.4d------------------------------\n",ilist);
      fprintf(stdmsg,"JIO Detail File=%s , fd= %d\n",&namelist[nameoff[ilist]],filedescriptors[ilist]);
      if(details_present == 0) continue;
      fprintf(stdmsg,"JIO Detail Routine   Calls          MB        MSEC       MB/s  (avg     stdev)\n");
      for(i=0; i<NUME; i++) {
        timelist[ilist*NUME+i]=0.000001*timelist[ilist*NUME+i]*time_factor;
        mbps[i]=0.0;
        bytelist[ilist*NUME+i]=bytelist[ilist*NUME+i]/(1024.0*1024.0); 
        bytelsqr[ilist*NUME+i]=bytelsqr[ilist*NUME+i]/(1024.0*1024.0*1024.0*1024.0);
        if( timelist[ilist*NUME+i] > 0.0 ) mbps[i]=bytelist[ilist*NUME+i]/(0.001*timelist[ilist*NUME+i]);
        if( calllist[ilist*NUME+i] > 0) { 
           fprintf(stdmsg,"JIO Detail %s   %6d  %10.3f  %10.3f %10.3f",
              &XXX[i][0],calllist[ilist*NUME+i],bytelist[ilist*NUME+i],timelist[ilist*NUME+i],mbps[i]); 
           avg=bytelist[ilist*NUME+i]/calllist[ilist*NUME+i] ;
           stdev=bytelsqr[ilist*NUME+i]/calllist[ilist*NUME+i] + (avg*avg)*(1.0/calllist[ilist*NUME+i] -2.0);
           if(stdev < 0.0 ) stdev = 0.0 ; stdev=sqrt(stdev);
           if(avg != 0.0) {
             if( avg > .9 ) fprintf(stdmsg," (%6.3f %6.3f)MB", avg, stdev);
             else           fprintf(stdmsg," (%6.3f %6.3f)KB", avg*1024., stdev*1024.);
           }
           fprintf(stdmsg,"\n");
        }
      }
    }
  }
  fprintf(stdmsg,"----------------------------------------------------------------\n");
}

/* initialize detailed stat tables */
static void jio_detail_init()
{
  int i;

  filedescriptors=malloc(4*maxf);
  namelist = malloc(maxf*AVGL);
  nameoff  = malloc(4*maxf);
  calllist = malloc(4*maxf*NUME);
  timelist = malloc(8*maxf*NUME);
  bytelist = malloc(8*maxf*NUME);
  bytelsqr = malloc(8*maxf*NUME);
  for (i=0; i<maxf; i++) nameoff[i]=0;
  for (i=0; i<maxf; i++) filedescriptors[i]=-1;
  for (i=0; i<maxf*NUME; i++) calllist[i]=0;
  for (i=0; i<maxf*NUME; i++) timelist[i]=0;
  for (i=0; i<maxf*NUME; i++) bytelist[i]=0;
  for (i=0; i<maxf*NUME; i++) bytelsqr[i]=0;
  for (i=0; i<maxd;      i++) indx[i]=0;
}

#ifdef linux
static pid_t gettid(void)
{
        return syscall(__NR_gettid);
}
#endif
/* initialize package at first call                                 */
/* return time stamp if jio_num is non zero, return 0 if jio_num is zero */
static long long jio_init()
{
  char * jio_env;
  int i;
  char *msgfile;
  FILE *temp;
  int affinity_cpu=0;
  int mpi_rank=-1;
  char *mpi_child;
  char msgfile_name[4096];

  if(jio_first) {
    return jio_num == 0 ? 0 : jio_time();
  }
  jio_init_extern();  /* initialize pointers to trapped functions */
  stdmsg=stderr;
  for (i=0 ; i<=NUME ; i++){
    calls[i] = 0;
    timj[i]  = 0.0;
    mbps[i]  = 0.0;
    bytes[i] = 0.0;
  }
  mpi_child = getenv("MP_CHILD");
  if(mpi_child != NULL){
    mpi_rank=atoi(mpi_child);
    fprintf(stderr,"MPI child %4.4d starting\n",mpi_rank);
  }else{
    mpi_child = getenv("PMI_RANK");
    if(mpi_child != NULL){
      mpi_rank=atoi(mpi_child);
      fprintf(stderr,"MPI child %4.4d starting\n",mpi_rank);
    }
  }
  jio_env = getenv("JIO_VERBOSE");
  if(jio_env != NULL) verbose_mode=1;
#ifdef linux
  jio_env = getenv("JIO_RDTSC");
  if(jio_env != NULL)
  {
    cpu_set_t mycpus;
    pid_t mytid;

    use_rdtsc=1;
    affinity_cpu=atoi(jio_env) ;
    fprintf(stderr,"Binding process to CPU %d\n",affinity_cpu);
    __CPU_ZERO((&mycpus)) ;
    __CPU_SET(affinity_cpu,&(mycpus)) ;
    sched_setaffinity(mytid,sizeof(cpu_set_t),&mycpus);
    time_factor=.4;
  }
#endif
  jio_env = getenv("JIO_ENV");
  if(jio_env == NULL) jio_env="JIO_SUMMARY";
  jio_num=1;
  if ( 0 == strcmp(jio_env,"JIO_NONE") )     jio_num=0; 
  if ( 0 == strcmp(jio_env,"JIO_SUMMARY") )  jio_num=1; 
  if ( 0 == strcmp(jio_env,"JIO_DETAILED") ) jio_num=2; 
  if ( 0 == strcmp(jio_env,"JIO_TRACE") )    jio_num=3; 
  jstart=jio_time();
  jio_first=1;
  indx = malloc(4*maxd);
  if(jio_num >=2 ) {
    jio_detail_init();
    jio_detail_new("STDIN",0);
    jio_detail_new("STOUT",1);
    jio_detail_new("STERR",2);
  }
  /* init done, time to switch trace file if requested */
  msgfile=getenv("JIO_DIAG");
  if(msgfile != NULL){
     if(mpi_rank < 0)
       snprintf(msgfile_name,sizeof(msgfile_name)-1,"%s",msgfile);
     else
       snprintf(msgfile_name,sizeof(msgfile_name)-1,"%s_%4.4d",msgfile,mpi_rank);
     stdmsg=fopen(msgfile_name,"a");  /* try to open file for append if one was specified */
  }
  if(stdmsg == NULL) {
    stdmsg = stderr ;
    fprintf(stderr,"ERROR: cannot open trace file %s, reverting to standard error\n",msgfile);
  }else{
    stdmsg_fileno=fileno(stdmsg);
    fprintf(stderr,"diagnostic file: %s, fd=%d\n",msgfile?msgfile:"Standard Error",stdmsg_fileno);
/*    fprintf(stdmsg,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"); */
  }
  fprintf(stdmsg,"JIO Diag JIO_ENV='%s'\n",jio_env);
  fprintf(stdmsg,"JIO Diag jio_num=%d\n",jio_num);
  atexit(jio_final);  /* call print routine automagically when program exits */
  return(jstart);
}

/*===================== interceptor functions start here =======================*/
/* fopen - fopen64 */
FILE *fopen(const char *name, const char *mode){
  FILE *result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_fopen_fn_ptr)(name, mode );
  if(result != NULL) {
    if(strcmp(name,"JIO_END") == 0) {
      fprintf(stdmsg,"JIO Calling jio_final in fopen, fd=%d\n",fileno(result));
      jio_final();
      return result;
    }
      jio_add_info(t0,jio_num,0,fileno(result),0,name);
  }
  return result;
}

#ifdef ioroutines64
FILE *fopen64(const char *name, const char *mode)
{
  FILE *result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_fopen64_fn_ptr)(name, mode );
  if(result != NULL) {
    if(strcmp(name,"JIO_END") == 0) {
      fprintf(stdmsg,"JIO Calling jio_final in fopen64, fd=%d\n",fileno(result));
      jio_final();
      return result;
    }
    jio_add_info(t0,jio_num,11,fileno(result),0,name);
  }
  return result;
}
#endif
/* fclose */
int fclose(FILE *fp)
{
  int iret;
  int fd;
  long long int t0;

  t0 = jio_init();
  fd=fileno(fp);
  if(verbose_mode)fprintf(stdmsg,"fclosing file = %d\n",fd);
  iret = (*S_fclose_fn_ptr)(fp );
  if(iret < 0) return iret;

  jio_add_info(t0,jio_num,1,fd,0,NULL);
  if(verbose_mode)fprintf(stdmsg,"fclosed file = %d\n",fd);
  indx[fd]=0;  /* file is now closed, no entry associated any more with this fd */
  return iret;
}
/* fread */
size_t fread(void *ptr, size_t size, size_t nobj, FILE *stream)
{
  size_t result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_fread_fn_ptr)(ptr, size, nobj, stream);
  jio_add_info(t0,jio_num,2,fileno(stream),size*nobj,NULL);
  return result;
}
/* fwrite */
size_t fwrite(const void *ptr, size_t size, size_t nobj, FILE *stream)
{
  size_t result;
  long long int t0;

  t0 = jio_init();
  result = (*S_fwrite_fn_ptr)(ptr, size, nobj, stream);
  jio_add_info(t0,jio_num,3,fileno(stream),size*nobj,NULL);
  return result;
}
/* open, open64 */
int open(const char * name, int flags, ...)
{
  int result;
  int perms=0777;  /* use umask value by default */
  long long int t0;

/*  get mode argument if O_CREAT */
  if (flags & O_CREAT)
    {
      va_list arg;
      va_start(arg, flags);
      perms = va_arg(arg, int);
      va_end(arg);
    }
  t0 = jio_init();

  result = (*S_open_fn_ptr)(name, flags, perms);
  if(result > 0) {
    if(strcmp(name,"JIO_END") == 0) {
      fprintf(stdmsg,"JIO Calling jio_final in open, fd=%d\n",result);
      jio_final();
      return result;
    }
    jio_add_info(t0,jio_num,4,result,0,name);
  }
  return result;
}
#ifdef ioroutines64
int open64(const char * name, int flags, ...)
{
  int result;
  int perms=0660;
  long long int t0;

  t0 =  jio_init();
  result = (*S_open64_fn_ptr)(name, flags, perms);
  if(result > 0) {
    if(strcmp(name,"JIO_END") == 0) {
      fprintf(stdmsg,"JIO Calling jio_final in open64, fd=%d\n",result);
      jio_final();
      return result;
    }
    jio_add_info(t0,jio_num,12,result,0,name);
  }
  return result;
}
#endif
/* close */
int close(int fd)
{
  int result;
  long long int t0;

  t0 = jio_init();
  if(verbose_mode)fprintf(stdmsg,"closing file = %d\n",fd);

  result = (*S_close_fn_ptr)(fd);
  if(result < 0) return result;

  jio_add_info(t0,jio_num,5,fd,0,NULL);
  indx[fd]=0;  /* file is now closed, no entry associated any more with this fd */
  if(verbose_mode)fprintf(stdmsg,"closed file = %d\n",fd);
  return result;
}
/* read */
ssize_t read(int fd, void * buf, size_t n)
{
  int result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_read_fn_ptr)(fd, buf, n);
  jio_add_info(t0,jio_num,6,fd,result,NULL);
  return result;
}

/* write */
ssize_t write(int fd, const void * buf, size_t n)
{
  int result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_write_fn_ptr)(fd, buf, n);

  jio_add_info(t0,jio_num,7,fd,result,NULL);
  return result;
}
/* fgets */
char *fgets(char *line, int maxline, FILE *fp)
{
  char *result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_fgets_fn_ptr)(line, maxline, fp);

  jio_add_info(t0,jio_num>2?2:jio_num, 8, fileno(fp), maxline, NULL);
  /* give at most 2 as first argument to fgets to keep verbosity under control */
  /* --------don't trace print because too much oputput
   if(jio_num >= 3) jio_trace2("&XXX[8][0],fp->_file,maxline); 
  -----------*/
  return result;
}

/* fseek */
int fseek(FILE *stream, long offset, int whence)
{
  int result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_fseek_fn_ptr)(stream, offset, whence);
  jio_add_info(t0,jio_num,9,fileno(stream),0,NULL);
  return result;
}

/* lseek, lseek64 */
off_t lseek(int fd, off_t offset, int whence)
{
  int result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_lseek_fn_ptr)(fd, offset, whence);
  jio_add_info(t0,jio_num,10,fd,0,NULL);
  return result;
}
#ifdef ioroutines64
off64_t lseek64(int fd, off64_t offset, int whence)
{
  int result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_lseek64_fn_ptr)(fd, offset, whence);

  jio_add_info(t0,jio_num,13,fd,0,NULL);
  return result;
}
#endif
#ifdef _AIX
/* kread */
ssize_t kread(int fd, void * buf, size_t n)
{
  int result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_kread_fn_ptr)(fd, buf, n);
  jio_add_info(t0,jio_num,14,fd,result,NULL);
  return result;
}

/* kwrite this is deactivated , it does not work */
ssize_t kwrite___(int fd, const void * buf, size_t n)
{
  int result;
  long long int t0;

  t0 =  jio_init();
  result = (*S_kwrite_fn_ptr)(fd, buf, n);

  jio_add_info(t0,jio_num,15,fd,result,NULL);
  return result;
}
#endif
