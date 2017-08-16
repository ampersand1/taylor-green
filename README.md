# taylor-green
Taylor Green Simulation with Athena4.2

Jasmin Shin
jshin156@g.ucla.edu
-----------------------------------------------------------
Taylor Green Simulation was done following the paper called "Paradigmatic flow for
small-scale magnetohydrodynamics: Properties of the ideal case the collision of current sheets" by
E. Lee et. al.  

In src/prob you will see the code called taylor-green.c that implements equations (6),(7), and (8) (the velocity field)
and equations (9),(10),(11) (the insulating box case of the magnetic field). The code can be run serially by using the commands in the lazy.sh file.   

In tst/3D-mhd there is the input file for this called athinput.taylor-green.  I added output[2-4].  These are the user defined outputs for the current.  (This code was a heavily modified version of the Orzag-Tang code).  They will come out as separate files like TaylorGreen.XXXX.current_x.vtk and there will also be the main fields outputed in TaylorGreen.XXXX.vtk, and they will appear in bin/ unless you specify otherwise.

Parallelizing on the HPCC
-----------------------------------------------------------
The code is parallelizable using mpirun and all the files have already been set in the Makefile.in file.  If you would like to run it, check out the lazy-parallel.sh file.      

I have also added the qsub file, so you can just qsub taylor-green.qsub right away after downloading the repository.  

The files will output into folders with names "id*" id0/ is supposed to hold the combined picture.  

--> Beware of the fact that when you run this with mpirun, the current_(x,y,z).vtk files will only give you a fraction of the picture while the TaylorGreen.XXXX.vtk will give you the full fields.  

--->It will probably be up to you out of convenience to figure out how to get the file to combine all the processes into the id0 file for the user defined output, but for now, in vis/vtk, you will have to combine them using join_vtk.c provided in the original Athena package.  This is unfortunately, where I left off.  Either combine_vtk.sh or combine.sh will combine all the current files, but you must set them by manually setting them inside the code.  I left off in making this work better.

Lastly, analysis files are in analysis/   

-----------------------------------------------------------
See https://princetonuniversity.github.io/Athena-Cversion/ for more information on the Athena code and how it works.
