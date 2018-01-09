import numpy as np
import matplotlib.pyplot as plt

a1 = np.loadtxt("test.txt",usecols = 0)
a2 = np.loadtxt("test.txt",usecols = 1)
t = np.arange(len(a1))

plt.subplot(211)
plt.plot(t,a1)

plt.subplot(212)
plt.plot(t,a2)

plt.show()

