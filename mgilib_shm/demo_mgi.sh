#!/bin/bash
export MGI_SHM_CFG=" 2 : oce2atm 32 : atm2oce 32 "
mpirun -n 1 ./oce.Abs >./oce.Abs.lst 2>&1 &
mpirun -n 1 ./atm.Abs | tee ./atm.Abs.lst
