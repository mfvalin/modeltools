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
void avg_2x2(float *s, int ni, int nj, int li, float *d, int nj2, int li2){
  int i, j, k;
  float t[4];
  float *s1, *s2, *d0;
  s1 = s ; s2 = s1 + li;
  for(j=0 ; j<nj-1 ; j+=2){   // pairs of rows
    d0 = d;
//     for(i=0 ; i<ni-7 ; i+=8){
//       for(k=0 ; k<8 ; k++) t[k] = s1[k] + s2[k];
//       d0[0] = t[0] + t[1];
//       d0[1] = t[2] + t[3];
//       d0[2] = t[4] + t[5];
//       d0[3] = t[6] + t[7];
//       d0 += 4;
//     }
    for(i=0 ; i<ni-3 ; i+=4){
      for(k=0 ; k<4 ; k++) t[k] = s1[k] + s2[k];
      d0[0] = t[0] + t[1];
      d0[1] = t[2] + t[3];
      d0 += 2;
    }
    s1 = s2 + li; s2 = s1 + li;
    d = d + li2;
  }
}

// in place expansion, source MUST start at z[ni - ni2]
void expand_2x2_1D(float *z, int ni, int ni2){
  float *s = z + ni - ni2;
  int lim = ni2 - (ni & 1);
  int i, j;
  float t;

  z[0] = 1.25f * s[0] - .25f * s[1];
  t = 1.25f * s[ni2-1] - .25f * s[ni2-2];;
  for(i=1, j=2 ; i < lim ; i++, j+=2){
    z[j-1] = .75f * s[i-1] + .25f * s[i];
    z[j  ] = .25f * s[i-1] + .75f * s[i];
printf("+++ %d %d %5.2f %5.2f\n",j-1,j,z[j-1],z[j]);
  }
  if(ni & 1) { // ni is odd
    z[ni-2] = (2.0f/3.0f) * s[ni2-2] + (1.0f/3.0f) * s[ni2-1];
    z[ni-1] = s[ni2-1];   // effectively a no-op
  }else{       // ni is even
    z[ni-1] = t;
printf("=== %5.2f \n",t);
  }
}
#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define NI 10
#define NI2 5
main(int argc, char **argv){
  float s1[NI];
  float s2[NI2];
  float s3[NI];
  float s4[NI2];
  int i;

  for(i=0 ; i<NI ; i++) s3[i] = 0;
  for(i=0 ; i<NI ; i++) s1[i] = i + .01f * i * i;
  for(i=0 ; i<NI ; i++) printf(" %5.2f",s1[i]); printf("\n");
  for(i=0 ; i<NI2 ; i++) s2[i] =  .5f * (s1[2*i] + s1[2*i+1]);
  for(i=0 ; i<NI2 ; i++) printf(" %5.2f",s2[i]); printf("\n");
  for(i=0 ; i<NI2 ; i++) s3[i+NI-NI2] = s2[i];
  for(i=0 ; i<NI ; i++) printf(" %5.2f",s3[i]); printf("\n");
  expand_2x2_1D(s3, NI, NI2);
  for(i=0 ; i<NI ; i++) printf(" %5.2f",s3[i]); printf("\n");
  for(i=0 ; i<NI2 ; i++) s4[i] = .5f * (s3[2*i] + s3[2*i+1]);
  for(i=0 ; i<NI2 ; i++) printf(" %5.2f",s4[i]); printf("\n");
}
#endif
