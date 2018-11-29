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

#pragma weak generate_unique_name__=generate_unique_name
int32_t generate_unique_name__(char *str, int32_t *sz);
#pragma weak generate_unique_name_=generate_unique_name
int32_t generate_unique_name_(char *str, int32_t *sz);
int32_t generate_unique_name(char *str, int32_t *sz){
  int32_t hid = gethostid();
  int32_t pid = getpid();
  int32_t t1;
  int64_t tf;
  struct timeval t;
  size_t szt = *sz;

  t1 = gettimeofday(&t, NULL);
  tf = t.tv_sec;
  tf *= 1000000;
  tf += t.tv_usec;
  t1 = snprintf(str, szt, "%8.8x%8.8x%16.16x" , hid, pid, tf);
  return t1 > *sz ? *sz : t1;
}
#if defined(SELF_TEST)
main(){
  char str[33];
  int szt = sizeof(str);
  int nc = generate_unique_name(str, &szt);
  fprintf(stderr,"nc = %d\n",nc);
  fprintf(stderr,"'%s'\n",str);
}
#endif
/*
 Fortran calling sequence
 character(len=33) :: buffer
 integer, external :: generate_unique_name
 integer :: number_of_chars
 number_of_chars = generate_unique_name(buffer,len(buffer))
*/