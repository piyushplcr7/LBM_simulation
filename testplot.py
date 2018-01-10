import numpy as np
import matplotlib.pyplot as plt

a1 = np.loadtxt("test.txt",usecols = 0)
a2 = np.loadtxt("test.txt",usecols = 1)
#a3 = np.loadtxt("test.txt",usecols = 2)

t = np.arange(len(a1))

plt.subplot(311)
plt.plot(t,a1)

plt.subplot(312)
plt.plot(t,a2)

#plt.subplot(313)
#plt.plot(t,a3)

plt.show()

