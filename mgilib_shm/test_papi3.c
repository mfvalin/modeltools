#include <papi.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_FLOPS 10000

void handle_error (int retval)
{
    /* print error to stderr and exit */
    PAPI_perror("ERROR: ");
    exit(retval);
}

void Handle_error(int code){
  printf("ERROR %d\n", code);
  exit(1);
}

void do_flops8(int n){
  double x=1.001;
  while(n-- > 0) x *= 1.000001;
  printf("x = %f\n",x);
}

void do_flops(int n){
  float x=1.001;
  float mult = 1.00001;
  while(n-- > 0) x *= mult;
  printf("x = %f\n",x);
}
#define NUM_EVENTS 3
main()
{
int retval, EventSet=PAPI_NULL;
int Events[NUM_EVENTS] = {PAPI_TOT_INS, PAPI_SP_OPS, PAPI_VEC_SP};
long_long values[11];
long_long values2[11];

/* Initialize the PAPI library */
retval = PAPI_library_init(PAPI_VER_CURRENT);
if (retval != PAPI_VER_CURRENT) {
  fprintf(stderr, "PAPI library init error!\n");
  exit(1);
}
fprintf(stderr, "PAPI has access to %d counters and %d components\n",PAPI_num_counters(),PAPI_num_components());
/* Create the Event Set */
if (PAPI_create_eventset(&EventSet) != PAPI_OK)
    handle_error(1);

/* Add Total Instructions Executed to our Event Set */
if (PAPI_add_event(EventSet, PAPI_TOT_INS) != PAPI_OK)
    handle_error(1);
if (PAPI_add_event(EventSet, PAPI_SP_OPS) != PAPI_OK)
    handle_error(1);
if (PAPI_add_event(EventSet, PAPI_VEC_SP) != PAPI_OK)
    handle_error(1);
//if (PAPI_add_event(EventSet, PAPI_DP_OPS) != PAPI_OK)
//    handle_error(1);

/* Start counting events in the Event Set */
if (PAPI_start(EventSet) != PAPI_OK)
    handle_error(1);
// if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK)
//     handle_error(1);

/* Defined in tests/do_loops.c in the PAPI source distribution */
do_flops(NUM_FLOPS);

/* Read the counting events in the Event Set */
if (PAPI_read(EventSet, values) != PAPI_OK)
    handle_error(1);
// if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
//     handle_error(1);

printf("After reading the counters: %lld %lld %lld\n",values[0],values[1],values[2]);

/* Reset the counting events in the Event Set */
if (PAPI_reset(EventSet) != PAPI_OK)
  handle_error(1);
if (PAPI_read(EventSet, values2) != PAPI_OK)
    handle_error(1);
printf("Read after resetting the counters: %lld %lld %lld\n",values2[0],values2[1],values2[2]);

do_flops(NUM_FLOPS);

/* Add the counters in the Event Set */
printf("Before adding the counters: %lld %lld %lld\n",values[0],values[1],values[2]);
// if (PAPI_read(EventSet, values) != PAPI_OK)
//    handle_error(1);
if (PAPI_accum(EventSet, values) != PAPI_OK)
   handle_error(1);
// if (PAPI_accum_counters(values, NUM_EVENTS) != PAPI_OK)
//    handle_error(1);
printf("After adding the counters: %lld %lld %lld\n",values[0],values[1],values[2]);

do_flops(NUM_FLOPS);

/* Stop the counting of events in the Event Set */
if (PAPI_stop(EventSet, values) != PAPI_OK)
    handle_error(1);
// if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
//     handle_error(1);

printf("After stopping the counters: %lld %lld %lld\n",values[0],values[1],values[2]);
if (PAPI_read(EventSet, values) != PAPI_OK)
    handle_error(1);
printf("After reading after stopping the counters: %lld %lld %lld\n",values[0],values[1],values[2]);
}

