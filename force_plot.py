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
cyl_x = np.loadtxt("Force.txt",usecols = 2)
cyl_y = np.loadtxt("Force.txt",usecols = 3)
flag_x = np.loadtxt("Force.txt",usecols = 4)
flag_y = np.loadtxt("Force.txt",usecols = 5)

nmin = 500;
Fx_mean= np.mean(F_x[nmin:])
Fy_mean = np.mean(F_y[nmin:])
np.mean(cyl_x[nmin:])
np.mean(cyl_y[nmin:])
np.mean(flag_x[nmin:])
np.mean(flag_y[nmin:])

F_x_added = cyl_x + flag_x;
F_y_added = cyl_y + flag_y;

t = np.arange(len(F_x))



plt.subplot(311)
line1,= plt.plot(t,F_x/norm, label = "Fx" ,linestyle = '-.') 
line2, = plt.plot(t, F_y/norm, label= "Fy",  linestyle = '-')
line1_a,= plt.plot(t,F_x_added, label = "Fx_added" ,linestyle = '-.') 
line2_a, = plt.plot(t, F_y_added, label= "Fy_added",  linestyle = '-')
plt.ylabel('Forces')
plt.xlabel('Time')

plt.legend(handles=[line1, line2, line1_a, line2_a ], loc =1)

plt.subplot(312)
line3,= plt.plot(t, cyl_x, label = "F_cyl_x" ,linestyle = '-.') 
line4, = plt.plot(t, cyl_y, label= "F_cyl_y",  linestyle = '-')
line5,= plt.plot(t, flag_x, label = "F_flag_x" ,linestyle = '-') 
line6, = plt.plot(t, flag_y, label= "F_flag_y",  linestyle = '--')

plt.ylabel('Forces')
plt.xlabel('Time')
plt.legend(handles=[line3, line4, line5, line6], loc =1)

plt.subplot(313)
line3b,= plt.plot([1, 2], [Fx_mean/norm, Fx_mean/norm], label = "F_x" ,linestyle = '-.') 
line4b, = plt.plot([1, 2], [Fy_mean/norm, Fy_mean/norm], label= "F_y",  linestyle = '-')
#line5b,= plt.plot([1, 2], [np.mean(flag_x)/norm, np.mean(flag_x)/norm], label = "F_flag_x" ,linestyle = '-') 
#line6b, = plt.plot([1, 2], [np.mean(flag_y)/norm, np.mean(flag_y)/norm], label= "F_flag_y",  linestyle = '--')
plt.ylabel('Mean Coefficients')
plt.xlabel('')
plt.legend(handles=[line3b, line4b], loc =1)

plt.show()


