#!/bin/bash -x

#-Wall -W -o join_vtk join_vtk.c -lm
#./join_vtk -o <outfile.vtk> infile1.vtk infile2.vtk ...

declare -a FILENUMBERS
declare -a THREADS


echo "enter number of files"
read filecount

echo "number of threads"
read threadcount

tcount="$(($threadcount - 1))"
#vals={1..$tcount}

#for i in $(seq -f "%04g" 0 $filecount)
#do
#  echo $i
#done
FILENUMBERS=($(seq -f "%04g" 460 $filecount))
THREADS=($(seq 1 1 $tcount))

#printf '%s\n' "${THREADS[@]}"

# current_x

#mkdir current_x
#change the 63 depending on the number of threads.  it is threads -1.  so 64-1.i
#Couldn't do $tcount because brace expansion happens first and then $, but I guess that got confused somewhere.  for the computer.   
for i in "${FILENUMBERS[@]}"
  do
       ./join_vtk -o combined_newrun3-4_"$i"_current_z.vtk /mnt/scratch/shinjasm/newrun3/newrun3-4/id0/TaylorGreen3."$i".current_z.vtk /mnt/scratch/shinjasm/newrun3/newrun3-4/id{1..63}/*."$i".current_z.vtk
  done

#mv *_current_x.vtk current_x
