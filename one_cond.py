import numpy as np

#sample locations
x = np.arange(-25, 26)
y = 3*np.ones((len(x),))

#conductor properties
f_cond = 60     #frequency (degrees)
x_cond = 0.    #x coordinate (feet)
y_cond = 4.    #y coordinate (feet)
subconds = 1.   #number of subconductors
d_cond = .89   #diameter (inches)
d_bund = .89
V_cond = 664.   #phase-phase voltage (kV)
I_cond = 220.   #phase current (amp)
p_cond = 0.    #phase angle (degrees)

#time interval
T = 1./f_cond
t = np.linspace(0, T, 1001)

#convenient variables/constants
epsilon = 8.854e-12
C = 1./(2.*np.pi*epsilon)
L = len(t)          #number of time steps
N = 1     #number of conductors
Z = len(x)          #number of sample points or x,y pairs

#calculate the effective conductor diameters
d_cond  = d_bund*((subconds*d_cond/d_bund)**(1./subconds))

#conversions
w_cond = 2*np.pi*f_cond    #convert to radians
x_cond *= 0.3048                #convert to meters
y_cond *= 0.3048                #convert to meters
d_cond *= 0.0254                #convert to meters
V_cond *= 1./np.sqrt(3)         #convert to ground reference from line-line
                                #leave in kV
                                #reference
p_cond *= 2*np.pi/360.          #convert to radians
x       = x*0.3048              #convert to meters
y       = y*0.3048              #convert to meters

#do the calcs
V = np.cos(w_cond*t + p_cond)

P = C*np.log(4*y_cond)

Q = V/P

E_x = np.zeros((L, Z))
E_y = np.zeros((L, Z))
for i in range(Z):
    nx = C*Q*(x[i] - x_cond)
    d1 = (x[i] - x_cond)**2 + (y[i] - y_cond)**2
    d2 = (x[i] - x_cond)**2 + (y[i] + y_cond)**2
    E_x[:,i] = nx/d1 - nx/d2

    ny1 = C*Q*(y[i] - y_cond)
    ny2 = C*Q*(y[i] + y_cond)
    E_y[:,i] = ny1/d1 - ny2/d2

E = np.sqrt(E_x**2 + E_y**2)

E_max = np.zeros((Z,))
for i in range(Z):
    E_max[i] = max(E[:,i])

print E_max

with open('one_cond_out.csv','w') as ofile:
    for z in zip(x,y,E_max):
        ofile.write('%s,%s,%s\n' % (str(z[0]), str(z[1]), str(z[2])))
