#!/bin/bash
## Compilation/build script for HEMELB
## Run from found location

## MODULE loads
##GCC compilers
MODULES(){

#Module environment on ARCHER2
#module restore PrgEnv-gnu 

#Module environment on SuperMUC-NG, default compilers are fine
#module load cmake

#module list
#Export compiler shortcuts as named on given machine
export CC=mpicc
export CXX=mpicxx

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
FOLDER=build_WKFtest
rm -rf $FOLDER
mkdir $FOLDER
cd $FOLDER

cmake -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DHEMELB_USE_GMYPLUS=OFF -DHEMELB_USE_MPI_WIN=OFF -DHEMELB_USE_SSE3=ON -DHEMELB_USE_AVX2=OFF -DHEMELB_OUTLET_BOUNDARY=GRINBERGKARNIADAKISWKIOLET -DHEMELB_WALL_OUTLET_BOUNDARY=GRINBERGKARNIADAKISWKBFL -DHEMELB_USE_VELOCITY_WEIGHTS_FILE=ON ..

make -j && echo "Done HemeLB Source"

cd ../..
}

SRCbuild_ARCHER2(){
cd src
FOLDER=build
rm -rf $FOLDER
mkdir $FOLDER
cd $FOLDER

cmake -DCMAKE_Fortran_COMPILER=${FC} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCTEMPLATE_LIBRARY=<PathToRepoOnARCHER2>/dep/install/lib/libctemplate.a -DHEMELB_USE_MPI_WIN=OFF ..


make -j && echo "Done HemeLB Source"

cd ../..
}


MODULES
#DEPbuild
SRCbuild
