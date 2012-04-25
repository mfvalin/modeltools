/* RMNLIB - Library of useful routines for C and FORTRAN programming
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
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "f77name.h"

/* dummy wall time clock, pseudo clock increases at one nanosecond per call */
static double dummy_time=0.0;
static double dummy_wtime()
{
   dummy_time += 1.0E-09;
   return dummy_time;
}

static double (*fn)() = dummy_wtime ;
static double time0 = 0.0;

/* returns wall time as a 64 bit real number */
/* example FORTRAN CALL:                     */
/* real *8 time                              */
/* time = RPN_COMM_Wtime()                   */
double f77_name(rpn_comm_wtime)()
{
  return (*fn)() - time0;
}

/* returns time of day as a 64 bit real number */
double f77name(rpn_comm_timeofday)()
{
  double temp;
  struct timeval tv;
  gettimeofday(&tv,NULL);
  temp = tv.tv_usec;
  temp *= .000001;
  temp += tv.tv_sec;
  return temp;
}

/* set wall clock function used by rpn_comm_wtime to arbitrary clock function */
/* that function MUST return a 64 bit real number (e.g MPI clock MPI_Wtime) */
/* example FORTRAN call : call RPN_COM_Wtime_set(MPI_Wtime) */
/* example FORTRAN call : call RPN_COM_Wtime_set(RPN_COMM_Timeofday)   */
/* special case :                                                      */
/*   call RPN_COM_Wtime_set(RPN_COM_Wtime_set)                         */
/*   resets the clock to the very fast dummy internal clock            */
void f77_name(rpn_comm_wtime_set)(double (*function)())
{
  fn = function;
  if(fn == (void *)f77_name(rpn_comm_wtime_set)) fn = (void *)dummy_wtime;
  time0 = (*fn)();
}

#ifdef TEST
#include <stdio.h>
#include <mpi.h>
int main(int argc,char **argv)
{
  int i;
  int ierr = MPI_Init(&argc, &argv);
  
  fprintf(stdout,"Phase 1, dummy timing function\n");
  for (i=0 ; i<5 ; i++){
    double x = f77_name(rpn_comm_wtime)();
    fprintf(stdout,"TIME1= %G\n",x);
  }
  fprintf(stdout,"Phase 2a, using MPI function MPI_Wtime\n");
  f77_name(rpn_comm_wtime_set)(MPI_Wtime);
  for (i=0 ; i<5 ; i++){
    double x = f77_name(rpn_comm_wtime)();
    fprintf(stdout,"TIME2= %G\n",x);
  }
  fprintf(stdout,"Phase 2b, using MPI function MPI_Wtime\n");
  f77_name(rpn_comm_wtime_set)(MPI_Wtime);
  for (i=0 ; i<5 ; i++){
    double x = f77_name(rpn_comm_wtime)();
    fprintf(stdout,"TIME2= %G\n",x);
  }
  fprintf(stdout,"Phase 3, using get_time_of_day\n");
  f77_name(rpn_comm_wtime_set)(f77name(rpn_comm_timeofday));
  for (i=0 ; i<5 ; i++){
    double x = f77_name(rpn_comm_wtime)();
    fprintf(stdout,"TIME2= %G\n",x);
  }
  fprintf(stdout,"Phase 4, dummy timing function\n");
  f77_name(rpn_comm_wtime_set)((void *)f77_name(rpn_comm_wtime_set));
  for (i=0 ; i<5 ; i++){
    double x = f77_name(rpn_comm_wtime)();
    fprintf(stdout,"TIME1= %G\n",x);
  }
  MPI_Finalize();
}
#endif