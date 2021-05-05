#!/bin/bash
## Compilation/build script for HEMELB
## Run from found location

## MODULE loads
##GCC compilers
MODULES(){

#Module environment on ARCHER2
module restore PrgEnv-gnu 

module list
#Export compiler shortcuts as named on given machine
export FC=ftn && export CC=cc && export CXX=CC
}

## HEMELB build
# 1) Dependencies
DEPbuild(){
cd dep
rm -rf build
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} ..
make -j  && echo "Done HemeLB Dependencies"

cd ../..
}


SRCbuild(){
cd src
rm -rf build
mkdir build
cd build

cmake -DCMAKE_Fortran_COMPILER=${FC} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCTEMPLATE_LIBRARY=<PathToRepoOnARCHER2>/dep/install/lib/libctemplate.a -DHEMELB_USE_MPI_WIN=OFF ..

make -j && echo "Done HemeLB Source"

cd ../..
}

MODULES
DEPbuild
SRCbuild
