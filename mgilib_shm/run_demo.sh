#!/bin/bash
set -x
rm -f DemoClientServerPort
mpirun -n 1 server_1 &
sleep 2
mpirun -n 2 ./client_1 "$(cat DemoClientServerPort)"
wait
