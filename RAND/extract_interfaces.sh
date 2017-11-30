#!/bin/bash
cat random_r250.c random_gaussian.c | grep InTf | sed 's/[ ]*!InTf!.*//' >randomfunctions.inc
cat random_r250.c random_gaussian.c | grep InTc | sed -e 's://[ ]*!InTc!.*:;:' >randomfunctions.h
