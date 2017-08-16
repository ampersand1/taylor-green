
#get_ipython().magic(u'matplotlib inline')
import yt
import numpy as np
from matplotlib import pylab
#combined_newrun3-4_0000_current_x.vtk
a = yt.load("/mnt/home/shinjasm/taylor-green4/vis/vtk/current_x/combined_newrun3-4_0???_current_x.vtk")
b = yt.load("/mnt/home/shinjasm/taylor-green4/vis/vtk/current_y/combined_newrun3-4_0???_current_y.vtk")
c = yt.load("/mnt/home/shinjasm/taylor-green4/vis/vtk/current_z/combined_newrun3-4_0???_current_z.vtk")


#a = yt.load("/mnt/scratch/shinjasm/newrun2/newrun2-1/id0/TaylorGreen3.00??.current_x.vtk")
#b = yt.load("/mnt/scratch/shinjasm/newrun2/newrun2-1/id0/TaylorGreen3.00??.current_y.vtk")
#c = yt.load("/mnt/scratch/shinjasm/newrun2/newrun2-1/id0/TaylorGreen3.00??.current_z.vtk")

j1=[]
j2=[]
j3=[]
times = []

for ds in a:
    #print (ds)
    dd = ds.all_data()
    j1.append(dd.quantities.extrema("current_x"))
    #j1.append(dd.quantities.weighted_average_quantity("current_x",'ones'))
    times.append(ds.current_time)
j1 = np.array(j1)

for ds in b:
    #print (ds)
    dd =ds.all_data()
    #j2.append(dd.quantities.weighted_average_quantity("current_y",'ones'))
    j2.append(dd.quantities.extrema("current_y"))
j2 = np.array(j2)

for ds in c:
    ##print (ds)
    dd = ds.all_data()
    #j3.append(dd.quantities.weighted_average_quantity("current_z",'ones'))
    j3.append(dd.quantities.extrema("current_z"))
j3 = np.array(j3)



#for ds in ts:
#    print (ds)
#    dd = ds.all_data()
#    j1.append(dd.quantities.extrema("current_x"))
#    times.append(ds.current_time)
#j1 = np.array(j)

#pylab.plot(times, j1a, '-xk')
#pylab.title("J3")

#ts = yt.load("TaylorGreen.0000.J3.vtk")
#for ds in ts:
#   yt.ProjectionPlot(ds,"x","cell_centered_B_y").save()


# In[21]:

times=np.array(times)

n=300
#pylab.plot(times[0:n], j1[0:n,1], '-g',times[0:n], j3[0:n,1],'-k',times[0:n], j3[0:n,1],'-r')
#pylab.savefig('jmax-nolog.png')

# In[14]:
#pylab.plot(times[0:n],np.sqrt(j1[0:n,1]*j1[0:n,1]+ j2[0:n,1]*j2[0:n,1]+j3[0:n,1]*j3[0:n,1]),'-r') # for averaged
#pylab.yscale('log')

#this is just taking the minima.  



#maxima
pylab.plot(times[0:n],np.sqrt(j1[0:n,1]*j1[0:n,1]+ j2[0:n,1]*j2[0:n,1]+j3[0:n,1]*j3[0:n,1]),'-r')  #averaged
pylab.yscale('log')


pylab.savefig('jmax.png')









