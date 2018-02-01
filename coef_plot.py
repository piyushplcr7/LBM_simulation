import numpy as np
import matplotlib.pyplot as plt
import math

cyl_r = 10;
u_inf = 0.05;
rho_inf = 1;
A = 2*cyl_r;
#norm=1;
norm = 0.5*u_inf*u_inf * rho_inf*A;

F_x = np.loadtxt("Force.txt",usecols = 0)
F_y = np.loadtxt("Force.txt",usecols = 1)
F_x_2 = np.loadtxt("Force.txt",usecols = 2)
F_y_2 = np.loadtxt("Force.txt",usecols = 3)
F_x_3 = np.loadtxt("Force.txt",usecols = 4)
F_y_3 = np.loadtxt("Force.txt",usecols = 5)

F_x = F_x_2 + F_x_3;
F_y = F_y_2 + F_y_3;

t = np.arange(len(F_x))

#plt.subplot(211)
line1,= plt.plot(t,F_x/norm, label = "Cd" ,linestyle = '-.') 
line2, = plt.plot(t, F_y/norm, label= "Cl",  linestyle = '-')
#line1_a,= plt.plot(t,F_x_added/norm, label = "Fx_added" ,linestyle = '-.') 
#line2_a, = plt.plot(t, F_y_added/norm, label= "Fy_added",  linestyle = '-')
plt.ylabel('Coefficients [-]')
plt.xlabel('Time')

plt.legend(handles=[line1, line2 ], loc =1)


plt.show()
