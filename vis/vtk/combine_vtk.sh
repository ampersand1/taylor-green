#!/bin/bash -x

#gcc -Wall -W -o join_vtk join_vtk.c -lm
#./join_vtk -o <outfile.vtk> infile1.vtk infile2.vtk ...

declare -a FILENUMBERS
declare -a THREADS


echo "enter number of files"
read filecount

echo "number of threads"
read threadcount

threadcount = threadcount-1

for i in $(seq -f "%04g" 0 $filecount)
do
  echo $i
done
FILENUMBERS= $(seq -f "%04g" 0 $filecount)
THREADS=($(seq 0 1 $threadcount))


# current_x

mkdir current_x
for i in "${FILENUMBERS}"
  do
    for j in "${THREADS}"
      do
       ./join_vtk -o combined_newrun3-4_${i}_current_x.vtk /mnt/scratch/shinjasm/newrun3/newrun3-4/id0/TaylorGreen3.${i}.current_x.vtk /mnt/scratch/shinjasm/newrun3/newrun3-4/id[1-$threadcount]/TaylorGreen3-id$j.${i}.current_x.vtk
      done
   done

mv *_current_x.vtk current_x


mkdir current_y
# current_y
for i in "${FILENUMBERS}"
   do
     for j in "${THREADS}"
      do
        ./join_vtk -o combined_newrun3-4_$i_current_y.vtk /mnt/scratch/shinjasm/newrun3/newrun3-4/id$j/TaylorGreen3-id$j.00$i.current_y.vtk
      done
    done

mv *_current_y.vtk current_y


mkdir current_z
# current_z
for i in "${FILENUMBERS}"
    do
     for j in "${THREADS}"
      do
        ./join_vtk -o combined_newrun3-4_$i_current_z.vtk /mnt/scratch/shinjasm/newrun3/newrun3-4/id$j/TaylorGreen3-id$j.00$i.current_z.vtk
        done
    done

mv *_current_z.vtk current_z

#./join_vtk -o combined_newrun3-4.vtk /mnt/scratch/shinjasm/newrun3/newrun3-4/id1/TaylorGreen3-id1.0000.current_x.vtk /mnt/scratch/shinjasm/newrun3/newrun3-4/id2/TaylorGreen3-id2.0000.current_x.vtk
