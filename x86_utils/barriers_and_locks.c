/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2018  Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

static int volatile __attribute__ ((aligned (64))) spins[64];
static int volatile __attribute__ ((aligned (64))) barrs[64];
static int volatile __attribute__ ((aligned (64))) locks[64];

void wait_barrier(int barrier, int maxcount){
 int count = __sync_fetch_and_add (barrs+barrier, 1);
 if(count == maxcount-1) {
  spins[barrier] = maxcount;
 }else{
  while(spins[barrier] != maxcount);
 }
}

void reset_barrier(int barrier){
 spins[barrier] = 0;
 barrs[barrier] = 0;
}

int test_lock(int lock){
 return locks[lock];
}

int release_lock(int lock, int me){
 return (__sync_val_compare_and_swap(locks+lock, me, 0) == me);
}

void acquire_lock(int lock, int me){
 while(__sync_val_compare_and_swap(locks+lock, 0, me) != 0) ;
}
