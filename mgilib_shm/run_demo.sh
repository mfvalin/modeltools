#!/bin/bash
set -x
mpicc -o server_1 server_1.c
mpicc -o client_1 client_1.c
rm -f DemoClientServerPort
mpirun -n 1 ./server_1 &
sleep 2
mpirun -n 2 ./client_1 "$(cat DemoClientServerPort)"
wait
rm -f DemoClientServerPort
