#!/bin/bash

make clean
./configure --with-problem=taylor-green --with-order=3
make all

cd bin
./athena -i ../tst/3D-mhd/athinput.taylor-green
