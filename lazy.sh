#!/bin/bash

cd bin
rm -r i*
cd ..
make clean
#--with-eos=isothermal
#--enable-smr
./configure --with-problem=taylor-green --with-eos=isothermal --with-order=3 --enable-mpi
make all MACHINE=hpcc

cd bin
mpirun -np 64 athena -i ../tst/3D-mhd/athinput.new time/tlim=4.0 --parallel 
