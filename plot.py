import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import math

Fx = np.loadtxt("Force.txt",usecols = 0)
Fy = np.loadtxt("Force.txt",usecols = 1)
print "check: " , Fx[-10:]
"""F = Fx
for i in range(len(Fx)):
    F[i] = math.sqrt(Fx[i]**2 + Fy[i]**2)"""
t = np.arange(len(Fx))
"""plt.subplot(221)
plt.plot(t,F)
plt.xlabel("Time")
plt.ylabel("Force")
plt.title("|F| vs time")"""
#plt.show()
rho = 1
uinf = 0.05
D = 2 * 10
#plt.subplot(222)
print .5*rho*uinf**2*D
plt.plot(t,Fx/(.5*rho*uinf**2*D))
plt.xlabel("Time")
plt.ylabel("Drag coefficient")
plt.title("Drag coefficient vs time")
#plt.show()

"""plt.subplot(223)
plt.plot(t,Fy)
plt.xlabel("Time")
plt.ylabel("Force")
plt.title("Fy vs time")"""
plt.show()
