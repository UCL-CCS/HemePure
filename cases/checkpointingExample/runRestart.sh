#!/bin/bash

rm -rf results-restart
mpirun -np 3 ./hemepure -in input-restart.xml -out results-restart
