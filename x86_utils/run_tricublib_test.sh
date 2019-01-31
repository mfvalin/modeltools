#!/bin/bash
set -ex
rm -f tricublin_beta_test.Abs tricublin_beta.o
s.cc -c -O3 tricublin_beta.c -march=core-avx2   ; s.f90 tricublin_beta_test.F90 -o tricublin_beta_test.Abs tricublin_beta.o
./tricublin_beta_test.Abs
