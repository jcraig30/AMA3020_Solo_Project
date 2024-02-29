#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:04:01 2024

@author: jamescraig
"""

import numpy as np
import math as mth
import matplotlib.pyplot as plt

# Constants
m = 1.3
R = 0.06
d = 0.005
mu = 0.032
g = 9.81

# Moment of Inertia
I_cm = 0.4*m*R**2

I_0 = I_cm +m*R**2


def traj_plotter(v_0,v_f,d):
    """
    traj_plotter - determine trajectory of a lawn bowl

    Parameters
    ----------
    v_0 : initial velocity.
    v_f : final velocity
    d : eccentricity of centre of mass.

    Returns
    -------
    x_array : x-values of trajectory.
    y_array : y-values of trajectory.

    """
    x_array = []
    y_array = []
    for t in np.arange(0,200,0.1):
        v = v_0-mu*g*t
        if abs(v)>v_f:
            a = (2*I_0*mu)/(m*d*R)
            r_0 = (a*v_0**2)/(2*mu*g)
            
            phi = 2/a *np.log(v_0/v)
        
            b = np.exp(-a*phi)
            
            x = r_0/(1+a**2) * (a-a*b*np.cos(phi)+b*np.sin(phi))

            y = r_0/(1+a**2) * (1-b*np.cos(phi)-a*b*np.sin(phi)) 
            
            if np.isnan(x):
                break
            else:
                x_array.append(x)
                y_array.append(y)
        else:
            break
    return x_array, y_array


# Varying v_0 and plotting trajectory
fig = plt.figure()
end_x = []
end_x.append(0)
end_y = []
end_y.append(0)
for v_0 in np.arange(1,5,1):
    traj = traj_plotter(v_0,0, d)
    plt.plot(traj[0], traj[1], label = str(v_0) + '$ms^{-1}$')
    end_x.append(traj[0][-1])
    end_y.append(traj[1][-1])
    
plt.plot(end_x, end_y, 'k:')
plt.legend()
plt.xlabel('x/m')
plt.ylabel('y/m')
plt.savefig('trajectory_v_0.pdf', bbox_inches = 'tight')
plt.show()
theta_1 = np.arctan(np.divide(end_y,end_x))

# Distance from origin as function of v_0
end_x = []
end_x.append(0)
end_y = []
end_y.append(0)
for v_0 in np.arange(0.25,5,0.25):
    traj = traj_plotter(v_0,0, d)
    # plt.plot(traj[0], traj[1], label = str(v_0) + '$ms^{-1}$')
    end_x.append(traj[0][-1])
    end_y.append(traj[1][-1])
D = np.sqrt(np.square(end_x)+np.square(end_y))
v_0 = np.arange(0,5,0.25)

fig = plt.figure()
plt.plot(v_0, D)
plt.xlabel('$v_0/ms^{-1}$')
plt.ylabel('$D/m$')
plt.savefig('D_v_0.pdf', bbox_inches = 'tight')
plt.show()

## Approximate square relation so took log plot to verify
# Log plot
fig = plt.figure()
# Remove initial zeros
v_0 = v_0[1:]
D = D[1:]
ln_v_0 = np.log(v_0)
ln_D = np.log(D)
plt.plot(ln_v_0, ln_D)
plt.xlabel('$\ln{v_0}$')
plt.ylabel('$\ln{D}$')
plt.savefig('ln_D_plot.pdf', bbox_inches = 'tight')
plt.show()


# Equation of line
def Least_Squares(x,y):
    """
    Least_Squares: Function to perform a least squares fit on the given data 
        
    Arguments
    ----------
    x : array of x-values (numpy array).
    y : array of y-values (numpy array).

    Returns
    -------
    m : Best-fit gradient value.
    c : Best-Fit y-intercept value.

    """
    N = np.shape(x)[0]
    S_x = np.sum(x)
    S_y = np.sum(y)
    S_xy = np.sum(x*y)
    S_xx = np.sum(x*x)
    m = (N*S_xy-S_x*S_y)/(N*S_xx-S_x**2)
    c = (S_xx*S_y-S_x*S_xy)/(N*S_xx-S_x**2)
    return m, c

grad,c = Least_Squares(ln_v_0, ln_D)
print('grad = ' + str(grad))
print('intercept = ' + str(c))
a = (2*I_0*mu)/(m*d*R)
c_theory = np.log(a/(2*mu*g*np.sqrt(a**2+1)))
print('theoretical grad = ' +str(2))
print('theoretical intercept= ' +str(c_theory))

# Calculation of required angle of release
a = (2*I_0*mu)/(m*d*R)
theta = np.arctan(1/a)
print('angle of release = '+str(mth.degrees(theta)))

# Varying eccentricity of COM and plotting trajectory
fig = plt.figure()
for d in np.arange(0.005,0.025,0.005):
    traj = traj_plotter(2,0, d)
    plt.plot(traj[0], traj[1], label = str(d) + '$m$')
    
plt.legend()
plt.xlabel('x/m')
plt.ylabel('y/m')
plt.savefig('trajectory_d.pdf', bbox_inches = 'tight')
plt.show()


# Consider needing to hit another bowl at 0.25ms^-1
d = 0.005
fig = plt.figure()
end_x = []
end_x.append(0)
end_y = []
end_y.append(0)
for v_0 in np.arange(2,7,1):
    traj = traj_plotter(v_0,1.75, d)
    plt.plot(traj[0], traj[1], label = str(v_0) + '$ms^{-1}$')
    end_x.append(traj[0][-1])
    end_y.append(traj[1][-1])
    
plt.plot(end_x, end_y, 'k:')
plt.legend()
plt.xlabel('x/m')
plt.ylabel('y/m')
plt.show()
theta_2 = np.arctan(np.divide(end_y,end_x))
# angle converges to ~ 0.759rad on increasing v_0



#Velocity required to move target by 0.5m
vel_target = np.sqrt(2*mu*g*0.5)
print(vel_target)
d = 0.005
fig = plt.figure()
end_x = []
end_x.append(0)
end_y = []
end_y.append(0)
for v_0 in np.arange(2,7,1):
    traj = traj_plotter(v_0,vel_target, d)
    plt.plot(traj[0], traj[1], label = str(v_0) + '$ms^{-1}$')
    end_x.append(traj[0][-1])
    end_y.append(traj[1][-1])
    
plt.plot(end_x, end_y, 'k:')
plt.legend()
plt.xlabel('x/m')
plt.ylabel('y/m')
plt.show()


# Plot for same v_0 the different target velocities to show variation in impact direction

d = 0.005
v_0 = 4
fig = plt.figure()
# impact angle to x-axis
angle_imp = []
end_x = []
end_x.append(0)
end_y = []
end_y.append(0)

# plot last 5 entries in trajectory values, assume straight line
for v_target in np.arange(0.5,3,0.5):
    traj = traj_plotter(v_0, v_target,d)
    plt.plot(traj[0][-5:], traj[1][-5:], label = str(v_target) + '$ms^{-1}$')
    end_x.append(traj[0][-1])
    end_y.append(traj[1][-1])
    # approximate each line segment as a straight line
    # use trig to find angle to x axis
    angle_imp.append(np.arctan((traj[1][-1]-traj[1][-5])/((traj[0][-1]-traj[0][-5]))))
    

plt.legend()
plt.xlabel('x/m')
plt.ylabel('y/m')
plt.ylim([0, 14])
plt.show()

# plot of impact angle to x-axis as a function of impact velocity
fig = plt.figure()
plt.plot(np.arange(0.5,3,0.5), np.degrees(angle_imp), 'ko')
plt.xlabel('$v_{impact}/ms^{-1}$')
plt.ylabel('$\phi_{impact}/ \degree$')
plt.show()


   

# Say jack is 21m away - take max distance as 21m and show angle of impact for varying v_0
    
    
def traj_plotter_dist(v_0,v_f,d, D):
    """
    traj_plotter - determine trajectory of a lawn bowl given the distance to a target bowl

    Parameters
    ----------
    v_0 : initial velocity.
    v_f : final velocity
    d : eccentricity of centre of mass.
    D : distance to target

    Returns
    -------
    x_array : x-values of trajectory.
    y_array : y-values of trajectory.

    """
    x_array = []
    y_array = []
    phi_array = []
    for t in np.arange(0,200,0.1):
        v = v_0-mu*g*t

        a = (2*I_0*mu)/(m*d*R)
        r_0 = (a*v_0**2)/(2*mu*g)
            
        phi = 2/a *np.log(v_0/v)
        
        b = np.exp(-a*phi)
            
        x = r_0/(1+a**2) * (a-a*b*np.cos(phi)+b*np.sin(phi))

        y = r_0/(1+a**2) * (1-b*np.cos(phi)-a*b*np.sin(phi)) 
            
        if np.isnan(x):
            break
        elif np.sqrt(x**2+y**2)<=D:
            x_array.append(x)
            y_array.append(y)
            phi_array.append(phi)               
        else:
            break
    return x_array, y_array, phi_array

# fig = plt.figure()
end_x = []
end_y = []
end_phi = []
for v_0 in np.arange(2,10,0.01):
    traj = traj_plotter_dist(v_0,0, d, 21)
    end_x.append(traj[0][-1])
    end_y.append(traj[1][-1])
    if np.degrees(traj[2][-1])<180:
        end_phi.append(traj[2][-1])
    
theta_3 = np.degrees(np.arctan(np.divide(end_y,end_x)))    

# v_0 as a function of theta
fig = plt.figure()
v_0 = np.arange(2,10,0.01)
plt.plot(v_0,theta_3, 'k.', markersize = 2)
plt.xlabel('$v_0/ms^{-1}$')
plt.ylabel('$\Theta/\degree$')
plt.show()

# theta remains constant for bowls which come to rest before the 21m
# theta decreases for further increasing v_0

# log plot
fig = plt.figure()
plt.plot(np.log(v_0), np.log(theta_3), 'k.', markersize = 2)
plt.xlabel('$ln{v_0}$')
plt.ylabel('$ln{\Theta}$')
plt.show()

# log plot not a straight line after theta starts to decrease
# therefore no exact exponential relation

# Angle of projection of target bowl
id_cut_off = len(theta_3)-len(end_phi)
sigma = np.degrees(end_phi)-theta_3[id_cut_off:]
# initial velocity for each collision angle
v_0_col = v_0[id_cut_off:]

# plot of angle of projection of target bowl from line of sight as a function of v_0
fig = plt.figure()
plt.plot(v_0_col, sigma, 'k.', markersize = 2)
plt.xlabel('$v_0/ms^{-1}$')
plt.ylabel('$\sigma/\degree$')
plt.savefig('impact angle.pdf', bbox_inches = 'tight')
plt.show()

