#!/bin/bash

rm -rf results
mpirun -np 3 ./hemepure -in input.xml -out results
