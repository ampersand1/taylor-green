import yt
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ts = yt.load("/mnt/scratch/shinjasm/newrun3/id0/TaylorGreen3.0???.vtk")
times = []

K_en= []
B_en= []

for ds in ts:
    dd = ds.all_data()
    K_en.append(dd.quantities.total_quantity("kinetic_energy"))
    B_en.append(dd.quantities.total_quantity("magnetic_energy"))
    times.append(ds.current_time)
B_en = np.array(B_en)
K_en = np.array(K_en)

times=np.array(times)

n=590
plt.plot(times[0:n], K_en[0:n], '-g',times[0:n], B_en[0:n],'-k')     
plt.savefig('energy.png')

