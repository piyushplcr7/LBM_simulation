import numpy as np
import matplotlib.pyplot as plt
import math


y = np.loadtxt("Nodes_scatter.txt",usecols = 1)
X = np.loadtxt("Nodes_scatter.txt",usecols = 0)

ys = np.loadtxt("solid_nodes.txt",usecols = 1)
xs = np.loadtxt("solid_nodes.txt",usecols = 0)
"""X=[]
y=[]
with open("Nodes_scatter.txt") as f:
	lines = f.readlines()
lines = lines[:-1]
for line in lines:
	print line.split()[0]
	X += [int(line.split()[0])]
	y += [int(line.split()[1])]"""

R = 50
nx = 800
ny = 400
circle = np.array([[R*math.cos(theta)+nx/4+.5, R*math.sin(theta)+ny/2+.5] for theta in np.linspace(0,2*np.pi,100)])
#X = np.array([X])
#y = np.array([y])
plt.figure()
flagella = np.array([[x,ny/2+.5] for x in np.linspace(R+nx/4+.5,3*R+nx/4+.5)])
plt.plot(flagella[:,0],flagella[:,1],c="k")
plt.plot(circle[:,0],circle[:,1],c="k")
#print "x is",x
#print y
plt.scatter(X,y,c="b")
plt.scatter(xs,ys,c="k")
plt.show()
