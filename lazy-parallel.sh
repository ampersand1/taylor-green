#!/bin/bash

rm -r i*
cd ~/taylor-green4/ 
make clean
#--with-eos=isothermal
#--enable-smr
./configure --with-problem=taylor-green --with-eos=isothermal --with-order=3 --enable-mpi
make all MACHINE=hpcc

cd /mnt/scratch/shinjasm/newrun3/
mpirun -np 64 ~/taylor-green4/bin/athena -i athinput.new time/tlim=25.0 --parallel  
