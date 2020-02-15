//     functions for C and FORTRAN programming
//     Copyright (C) 2020  Recherche en Prevision Numerique
// 
//     This software is free software; you can redistribute it and/or
//     modify it under the terms of the GNU Lesser General Public
//     License as published by the Free Software Foundation,
//     version 2.1 of the License.
// 
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//     Lesser General Public License for more details.
//
static void inline lissage_121_1D(float *dst, float *src, int n){  // 1 2 1 smoothing pass along i
  int i;
  dst[0] = .5f * (src[0] + src[1]);
  for(i=1 ; i<n-1 ; i++) dst[i] = .25f * src[i-1] + .5f * src[i] + .25f * src[i+1];
  dst[n-1] = .5f * (src[n-2] + src[n-1]);
}

void lissage_121_2D(float *z, int ni, int nj, int li){
  float t1[li], t2[li], *bot, *row, *top, *new, *old, *t;
  int i, j;
  old = t1 ; new = t2;
  row = z ; top = z + ni;                                  // rows 0 and 1
  for(i=0 ; i<li ; i++) old[i] = .5f * (row[i] + top[i]);  // row 0, smooth along j

  for(j=1; j<nj-1 ; j++){              // middle rows, smoothing along j, then i
    bot = row; row = top ; top += ni;  // current set of rows (bot, row, top : j-1 , j , j+1)
    for(i=0 ; i<li ; i++) new[i] = .25f * bot[i] + .5f * row[i] + .25f * top[i];  // row j, smooth along j
    lissage_121_1D(bot, old, li);      // row j-1, smooth along i and store (row j-1 no longer needed as a source)
    t = old ; old = new ; new = t;     // swap old and new
  }
  for(i=0 ; i<li ; i++) new[i] = .5f * row[i] + .5f * top[i];   // row nj-1, smooth along j
  lissage_121_1D(row, old, li);        // row nj-2 , smooth along i and store
  lissage_121_1D(top, new, li);        // row nj-1 , smooth along i and store
}

