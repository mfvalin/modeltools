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
#if defined(linux)
#define _GNU_SOURCE
#endif
#include <sched.h>
#include <stdlib.h>
#include <stdio.h>

/*
   unset any processor affinity for the current process
   if the FULL_UNBIND environment variable is defined
   ( FULL_UNBIND would normally be set by r.mpirun / ord_soumet )
   in the case of an MPI launch, only process 0 will print a message
   this code is only applicable to Linux for the time being
   and has been written to counter openmpi behaviour in some cases.
   (written specifically for use in the RPN_COMM library)

   this code is Linux only

   Michel Valin , 2011 / 11 / 02 UQAM
*/

/* use pragma weak to create the alternate FORTRAN entry points */
#pragma weak rpn_comm_unbind_process_=rpn_comm_unbind_process
#pragma weak rpn_comm_unbind_process__=rpn_comm_unbind_process


/* WARNING: some old versions of gcc may not generate the weak entry points correctly */

void rpn_comm_unbind_process(void)
{
#if defined(linux)

cpu_set_t set;
int i;
int will_print=1;
char *ompi=getenv("OMPI_COMM_WORLD_RANK");  /* openmpi */

if(ompi == NULL) ompi=getenv("PMI_RANK");   /* mpich family */

if(ompi != NULL) if(0 != atoi(ompi)) will_print=0;  /* not process 0, no message */

if(getenv("FULL_UNBIND") == NULL) { if(will_print) printf("NO unbinding will be done\n"); return ; }

if(will_print) printf("FULL unbinding will be done\n");

CPU_ZERO(&set);
for (i=0;i<=11;i++) CPU_SET(i,&set);
sched_setaffinity(0,sizeof(set),&set);

#endif

return;
}
