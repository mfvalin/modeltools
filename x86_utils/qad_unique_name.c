/*
 * Copyright (C) 1995-2018 Recherche en Prevision Numerique
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
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

// base 64 encoding table
static char *b64="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789@=";

// generate a quick and dirty NULL terminated "unique" character string
// the string including the NULL character must fit within sz characters
// the function returns the number of usable characters in name

// the string is "locally" unique (at least within a cluster)
// it uses hostid, process id, and the time of day in microseconds

// generate the most "popular" Fortran name manglings
#pragma weak qad_unique_name__=qad_unique_name
int32_t qad_unique_name__(char *str, int32_t *sz);
#pragma weak qad_unique_name_=qad_unique_name
int32_t qad_unique_name_(char *str, int32_t *sz);

int32_t qad_unique_name(unsigned char *str, int32_t *sz){
  uint32_t hid = gethostid();
  uint32_t pid = getpid();
  int32_t t1, nc, i;
  uint64_t ts;
  struct timeval t;
  size_t szt = *sz;
  unsigned char lsz[33];

  t1 = gettimeofday(&t, NULL);
  ts = t.tv_sec ;
  ts = ts * 10000000 ;       // convert seconds to microseconds
  ts = ts + t.tv_usec ;      // 32 bits for seconds, 20 bits for microseconds, we will encode 9*6 =  54 bits
  for(i=0 ; i<9 ; i++) {lsz[i] = b64[ts & 0x3F] ; ts = ts >> 6; } ;  // encode 6 bits at a time (9 chars)

  ts = pid ;
  ts = ts << 32;
  ts = ts | hid;
  for(i=0 ; i<11 ; i++) {lsz[i+9] = b64[ts & 0x3F] ; ts = ts >> 6; } ;  // encode 6 bits at a time (11 chars)

  lsz[20] = '\0' ;
  t1 = 20 ; // need 20 chars to encode the whole thing

  nc = t1 > szt - 1 ? szt - 1 : t1;                 // there is room for at most szt - 1 characters
  for(i=0 ; i<nc ; i++) { str[i] = lsz[i] ; } ;     // copy nc characters
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
  nc = qad_unique_name(str,&sz);
  str[255] = '\0';
  fprintf(stdout,"%s",str);
#if defined(SELF_TEST)
  fprintf(stdout,"\nnc = %d\n",nc);
#endif
  return 0;
}
#endif
