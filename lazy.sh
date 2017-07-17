#!/bin/bash

cd bin
rm -r *.vtk
cd ..
make clean
#--with-eos=isothermal
#--enable-smr
#--with-eos=isothermal
./configure --with-problem=taylor-green --with-eos=isothermal --with-order=3
make all

cd bin
./athena -i ../tst/3D-mhd/athinput.taylor-green
