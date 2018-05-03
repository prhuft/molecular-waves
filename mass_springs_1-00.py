""""
Mechanical Vibration Transfer Simulation, v1.05
	
Preston Huft, Spring 2018. 

Numerical simulation of _name here_, solved iteratively with
the Runge-Kutta (4th order) method. 

Version notes: This first version is in the spirit of (and frankly, the flesh
of) double_pendulum_1-05.py. 


- Generalize get_initial_states for any type of system; i.e. pass in a function
	to handle system-specific things, such as getting the initial mth derivs. 
"""

## LIBRARIES

import helpfunclib as hfl
from rk4_multicoord import rk4_update as rk4
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from random import random as rn
import numpy as np
from copy import copy as cp
import math as m
from math import cos
from math import sin
from math import sqrt

## GLOBAL CONSTANTS

g = 9.8 # [m/s]
DEBUG = 1 # set to 1 if we're debugging
iter = 0;

## METHODS

def derivs(state,tau,params):
	""" Returns a list of the first and second derivatives of the generalized 
		positions of each mass in the system.
		
		'state' is a list of the lists of the derivatives of generalized coordinates
		for each object in the system, grouped by derivative order. There are lists
		of these derivative list for each coordinate need to describe the system:
		
		state = [coordinate1,coordinate2,...,coordinateq], where
			coordinate1 = [[f1^(0),f2^(0),...,fn^(0)],[f1^(1),f2^(1),...fn^(1)],
					...,[f1^(m),f2^(m),...fn^(m)]]
			coordinate2 = likewise, and so on for the other coordinates
		tau = the timestep. It does not get used here.
	"""
	
	
	iter = iter + 1
	print("iter: ",iter)
	
	# The mass positions. state[i][j] is the jth deriv group for coord. i
	x_list = state[0][0] # <-- on second call, this is a scalar, 
						 # and state is length 3. wth?
	y_list = state[1][0]
	vx_list = state[0][1]
	vy_list = state[1][1]
	
	if DEBUG:
		print("x_list ", x_list)
		print("len(state) ",len(state))
		print("state[0] ",state[0])
		print("state[0][0] ",state[0][0])
	
	# The spring constant,mass, and network dimensions
	kx,ky,m,r0,xlen,ylen = params

	o2_x = kx/m
	o2_y = ky/m
	
	ax_list = xlen*ylen*[0]
	for i in range(0,xlen):
		for j in range(0,ylen):
			n_ij = xlen*i + j # map the grid location to a 1D index
			# ax[n_ij] = 0 # ax for the (ij)th mass (thinking in 2D)
			
			# Add acceleration term if we're not on the...
			if (i > 0): # ... top edge
				n_hj = xlen*(i-1) + j
				dx = x_list[n_hj]-x_list[n_ij]
				dy = y_list[n_hj]-y_list[n_ij]
				ax_list[n_ij] = o2_x*(1-r0/sqrt(dx**2+dy**2))*dx
			if (i < xlen-1): # ... bottom edge
				n_jj = xlen*(i+1) + j
				dx = x_list[n_jj]-x_list[n_ij]
				dy = y_list[n_jj]-y_list[n_ij]
				ax_list[n_ij] = o2_x*(1-r0/sqrt(dx**2+dy**2))*dx
			if (j > 0): # ... left edge
				n_ii = xlen*i + (j-1)
				dx = x_list[n_ii]-x_list[n_ij]
				ax_list[n_ij] = -o2_x*(dx+r0)
			if (j < ylen-1): # ... right edge
				n_ik = xlen*i + (j+1)
				dx = x_list[n_ik]-x_list[n_ij]
				ax_list[n_ij] = o2_x*(dx-r0)
	
	ay_list = xlen*ylen*[0]
	for i in range(0,xlen):
		for j in range(0,ylen):
			n_ij = xlen*i + j # map the grid location to a 1D index
			# ay[n_ij] = 0 # ax for the (ij)th mass (thinking in 2D)
		
			# Add acceleration term if we're not on the...
			if (j > 0): # ... left edge
				n_ii = xlen*i + (j-1)
				dx = x_list[n_ii]-x_list[n_ij]
				dy = y_list[n_ii]-y_list[n_ij]
				ay_list[n_ij] = o2_y*(1-r0/sqrt(dx**2+dy**2))*dy
			if (j < ylen-1): # ... right edge 
				n_ik = xlen*i + (j+1)
				dx = x_list[n_ik]-x_list[n_ij]
				dy = y_list[n_ik]-y_list[n_ij]
				ay_list[n_ij] = o2_y*(1-r0/sqrt(dx**2+dy**2))*dy
			if (i > 0): # ... top edge
				n_hj = xlen*(i-1) + j
				dy = y_list[n_hj]-y_list[n_ij]
				ay_list[n_ij] = -o2_y*(dy+r0)
			if (i < xlen-1): # ... bottom edge
				n_jj = xlen*(i+1) + j
				dy = y_list[n_ik]-y_list[n_ij]
				ay_list[n_ij] = o2_y*(dy-r0)
	
	return [vx_list,ax_list],[vy_list,ay_list]
		
def get_data(state,tau,steps,params,num_update):
	""" Returns a list of the states of the system at each timestep, for steps
		number of iterations, tau [s] apart. num_update is thenumerical method 
		function to be used. 
		
		state = [q1,q2,...,qz] where each q is a necessary parameter for 
			describing the generalized position of a specific subsystem:
			q = [[f1^(0),f2^(0),...,fn^(0)],[f1^(1),f2^(1),...fn^(1)],
				  ...,[f1^(m),f2^(m),...fn^(m)]]
		params = system-specific constants. E.g., mass, length, density, etc.
	"""
	
	# Extract system parameters
	kx,ky,m,r0,xlen,ylen = params
		
	# Timestep
	dt = tau
	
	xdata = [] # the x positions
	ydata = [] # ditto, but in y
	
	# Forward feed the solver method for i = 0 to i = steps
	for i in range(0,steps): 
		try:
			# Update the state of the network
			new_state = num_update(state,dt,params,derivs)
			
			xdata.append(new_state[0][0]) 
			ydata.append(new_state[1][0])
			
		except ValueError:		
			print('value error at iteration ',i)
			break
			
	return xdata,ydata
	
def get_initial_state(params,tau):
	""" Returns a list describing the state, in the form:
		state = [q1,q2,...,qz] where each q is a necessary parameter for 
			describing the generalized position of a specific subsystem:
			q = [[f1^(0),f2^(0),...,fn^(0)],[f1^(1),f2^(1),...fn^(1)],
				  ...,[f1^(m),f2^(m),...fn^(m)]]
		params = system-specific constants. E.g., mass, length, density, etc.
		
		Assume that the particles, initially stretched random distances out of 
		equilibrium, are starting from rest. 
	"""
	
	# The spring constant,mass, and network dimensions
	kx,ky,m,r0,xlen,ylen = params
	total =xlen*ylen # number of subsystems
	
	rx_list,ry_list = total*[0],total*[0]
	vx_list,vy_list = total*[0],total*[0] # replace later with derivs()
	ax_list,ay_list = total*[0],total*[0] # replace later with derivs()
	
	# Build the part of the state which is projected on x
	state_x = []
	for i in range(0,xlen): # iterate over the columns
		for j in range(0,ylen): # iterate over the rows
			rx_list[xlen*i + (j-1)] = i*xlen + r0*(rn()-.5) # the x coord of mass ij
	state_x.append(rx_list)
	state_x.append(vx_list)
	state_x.append(ax_list)
	
	state_y = []
	for i in range(0,xlen): # iterate over the columns
		for j in range(0,ylen): # iterate over the rows
			ry_list[xlen*i + (j-1)] = j*ylen + r0*(rn()-.5) # the y coord of mass ij
	state_y.append(ry_list)
	state_y.append(vy_list)
	state_y.append(ay_list)
		
	state = [state_x,state_y]
	print("state: ", state)
	# Replace ax_list and ay_list using derivs()
	deriv_list = derivs(state,tau,params)
	state[0][2] = deriv_list[0][1]
	state[1][2] = deriv_list[1][1]
	print("state: ", state)
	
	return state
	
## SYSTEM INITIALIZATION

# Simulation parameters
m = .1 # [kg] these are massive particles lol
kx = 1 # [N/m] Spring constant in x
ky = 1 # [N/m] Spring constant in y
r0 = .1 # [m] the spring equilibrium length
x_num =  2 # number of columns of masses
y_num = 2 # number of rows of masses
params = [kx,ky,m,r0,x_num,y_num]

dt = 0.01 # [s]
iters = 10000 # times to update the systems

# Generate the initial state
state_0 = get_initial_state(params,dt)

# # Generate the data
xdata,ydata = get_data(state_0,dt,iters,params,rk4)

## SIMULATION SETUP

# Initialize the figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_facecolor('black')
ax.set_aspect(aspect='equal')
fig.patch.set_facecolor('black')

# # Initialize the lines; actually, it will be more of a scatter plot
scatter = []
scatter, = ax.plot([],[],color='purple',marker='o')

def init():
	""" Set the axes limits with global values. """
	ax.set_ylim(-r0,xlen*(r0+1))
	ax.set_xlim(-r0,ylen*(r0+1))
	return scatter,
	
def update(i):
	""" Set the ith data points with global values."""
	
	# The points describing the network at the ith step
	scatter.set_data(xdata[i],ydata[i])
	return scatter,

# Run the animation
anim = animation.FuncAnimation(fig, update, frames=range(0,iters+1), 
	init_func=init, blit=True, interval=1000*dt, repeat=False)
plt.show()

