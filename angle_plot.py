import numpy as np
import matplotlib.pyplot as plt

a1 = np.loadtxt("Alpha.txt",usecols = 0)
a2 = np.loadtxt("Alpha.txt",usecols = 1)
#a3 = np.loadtxt("Alpha.txt",usecols = 2)

t = np.arange(len(a1))

#plt.subplot(311)
line1,= plt.plot(t,a1, label = "Alpha1" ,linestyle = '-.') 
line2, = plt.plot(t, a2, label="Alpha2",  linestyle = '-')
plt.ylabel('Alpha')
plt.xlabel('Time')

plt.legend(handles=[line1, line2], loc =1)

#plt.subplot(312)
#plt.plot(t,a2, '-')

#plt.subplot(313)
#plt.plot(t,a3)

plt.show()


