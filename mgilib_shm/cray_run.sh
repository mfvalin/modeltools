#!/bin/bash
set -x
rm -f DemoClientServerPort
apmgr pdomain -c mfv_domain
aprun -n 1 -p mfv_domain ./server_1 &
sleep 5
aprun -n 2 -p mfv_domain ./client_1 "$(cat DemoClientServerPort)"
wait
sleep 2
apmgr pdomain -r mfv_domain
rm -f DemoClientServerPort
