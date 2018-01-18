import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import math

Fx = np.loadtxt("Result/Force.txt",usecols = 0)
Fy = np.loadtxt("Result/Force.txt",usecols = 1)
F = Fx
for i in range(len(Fx)):
    F[i] = math.sqrt(Fx[i]**2 + Fy[i]**2)
t = np.arange(len(Fx))
plt.subplot(221)
plt.plot(t,F)
plt.xlabel("Time")
plt.ylabel("Force")
plt.title("|F| vs time")
#plt.show()

plt.subplot(222)
plt.plot(t,Fx)
plt.xlabel("Time")
plt.ylabel("Force")
plt.title("Fx vs time")
#plt.show()

plt.subplot(223)
plt.plot(t,Fy)
plt.xlabel("Time")
plt.ylabel("Force")
plt.title("Fy vs time")
plt.show()
