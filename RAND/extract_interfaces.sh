#!/bin/bash
cat randomgeneric.c random_r250.c random_shr3.c random_mt19937.c random_gaussian.c | grep InTf | sed 's/[ ]*!InTf!.*//' >randomfunctions.inc
cat randomgeneric.c random_r250.c random_shr3.c random_mt19937.c random_gaussian.c | grep InTc | sed -e 's://[ ]*!InTc!.*:;:' >randomfunctions.h
