/* 
 * Copyright (C) 2018  Environnement Canada
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#define FIOL_VERSION 0x1BAD
// in == out means buffer is empty
// in == out-1 (or in=limit-1 && out==0) means buffer is full
typedef struct{        // circular buffer management variables
  int32_t version;     // version marker
  int32_t first;       // should be 0 (assumed to be 0 in circular_buffer.c)
  int32_t in;          // start inserting data at data[in]
  int32_t out;         // start extracting data at data[out]
  int32_t limit;       // size of data buffer (last available index + 1)
} fiol_management;
typedef fiol_management *fiol_management_p;

typedef struct{        // skeleton for circular buffer
  fiol_management m;   // management variables
  int32_t data[1];     // data buffer (contains at most limit -1 useful data elements)
} circular_buffer;
typedef circular_buffer *circular_buffer_p;

int circular_buffer_space_available(circular_buffer_p p);
int circular_buffer_wait_space_available(circular_buffer_p p, int n);
int circular_buffer_data_available(circular_buffer_p p);
int circular_buffer_wait_data_available(circular_buffer_p p, int n);
int circular_buffer_atomic_get(circular_buffer_p p, int *dst, int n);
int circular_buffer_atomic_put(circular_buffer_p p, int *src, int n);
