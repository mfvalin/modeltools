/*
 * RMNLIB 
 * Copyright (C) 1995-2018 Recherche en Prevision Numerique
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
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

// generate a quick and dirty NULL terminated "unique" character string
// the string including the NULL character must fit within sz characters
// the function returns the number of usable characters in name

// the string is "locally" unique (at least within a cluster)
// it uses hostid, process id, and the time of day in microseconds

// generate the most "popular" Fortran name manglings
#pragma weak generate_qad_unique_name__=generate_qad_unique_name
int32_t generate_qad_unique_name__(char *str, int32_t *sz);
#pragma weak generate_qad_unique_name_=generate_qad_unique_name
int32_t generate_qad_unique_name_(char *str, int32_t *sz);

int32_t generate_qad_unique_name(char *str, int32_t *sz){
  int32_t hid = gethostid();
  int32_t pid = getpid();
  int32_t t1, nc, i;
  int32_t ts, tu;
  struct timeval t;
  size_t szt = *sz;
  char lsz[33];

  t1 = gettimeofday(&t, NULL);
  ts = t.tv_sec ;
  tu = t.tv_usec ;
  t1 = snprintf(lsz, 33, "%8.8x%8.8x%8.8x%8.8x" , hid, pid, ts, tu); // host pid sec usec
  nc = t1 > szt ? szt - 1 : t1;
  for(i=0 ; i<nc ; i++) { str[i] = lsz[t1-1-i] ; }  // reverse string to get max changes at beginning
  str[nc] = '\0';                                   // make sure the string is NULL terminated
  return nc;
}
/*
 Fortran usage :
 character(len=33) :: buffer
 integer, external :: generate_qad_unique_name
 integer :: number_of_chars
 number_of_chars = generate_qad_unique_name(buffer,len(buffer))
*/
#if defined(WITH_MAIN)
int main(int argc,char **argv){
  char str[256];
  int nc, sz;

  sz = sizeof(str);
  if(argc > 1) sz = atoi(argv[1]);
  nc = generate_qad_unique_name(str,&sz);
  str[255] = '\0';
  fprintf(stdout,"%s",str);
#if defined(SELF_TEST)
  fprintf(stdout,"\nnc = %d\n",nc);
#endif
  return 0;
}
#endif
