"""
Material Heat Transfer (Thermal Gradient) Simulation, v1.00 
(based on mass_springs_1-00.py)
	
Preston Huft, Spring 2018. 

Numerical simulation of masses connected by invisible, massless springs, solved 
iteratively with the Runge-Kutta (4th order) method. 

Version notes: Particles on leftmost and rightmost are fixed in position, but 
those on the left have a non-decaying sinuisoidal noise (thermal energy), while
those on the left have 

To-do:
- Generalize get_initial_states for any type of system; i.e. pass in a function
	to handle system-specific things, such as getting the initial mth derivs. 
- Color particles according to mean frequency from doing a FFT.
"""

## LIBRARIES

import helpfunclib as hfl
from rk4_multicoord import rk4_update as rk4
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import msvcrt as ms
from random import random as rn
import numpy as np
from copy import copy as cp
from math import cos,sin,exp,sqrt,pi

## GLOBAL CONSTANTS

g = 9.8 # [m/s]
DEBUG = 0 # set to 1 if we're debugging

## METHODS - NECESSARY

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
	
	# The mass positions. state[i][j] is the jth deriv group for coord. i
	x_list = state[0][0] 
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
	
	# Sum the acceleration in x
	ax_list = (xlen*ylen)*[0]
	for j in range(0,xlen):
		for i in range(0,ylen):
			n_ij = xlen*i + j # map the grid location to a 1D index
			ax_list[n_ij] = 0 # ax for the (ij)th mass (thinking in 2D)
			# Add acceleration term if we're not on the...
			if (i > 0): # ... top edge; F_hj,ij
				n_hj = xlen*(i-1) + j
				dx = x_list[n_hj]-x_list[n_ij]
				dy = y_list[n_hj]-y_list[n_ij]
				ax_list[n_ij] += o2_x*(1-r0/sqrt(dx**2+dy**2))*dx
			if (i < ylen-1): # ... bottom edge; F_jj,ij
				n_jj = xlen*(i+1) + j
				dx = x_list[n_jj]-x_list[n_ij]
				dy = y_list[n_jj]-y_list[n_ij]
				ax_list[n_ij] += o2_x*(1-r0/sqrt(dx**2+dy**2))*dx
			if (j > 0): # ... left edge; F_ii,ij
				n_ii = xlen*i + (j-1)
				dx = x_list[n_ii]-x_list[n_ij]
				ax_list[n_ij] += o2_x*(dx+r0)
			if (j < xlen-1): # ... right edge; F_ik,ij
				n_ik = xlen*i + (j+1)
				dx = x_list[n_ik]-x_list[n_ij]
				ax_list[n_ij] += o2_x*(dx-r0)
				
	# Sum the acceleration in y
	ay_list = xlen*ylen*[0]
	for j in range(0,xlen):
		for i in range(0,ylen):
			n_ij = xlen*i + j # map the grid location to a 1D index
			ay_list[n_ij] = 0 # ax for the (ij)th mass (thinking in 2D)
			# Add acceleration term if we're not on the...
			if (j > 0): # ... left edge; F_ii,ij
				n_ii = xlen*i + (j-1)
				dx = x_list[n_ii]-x_list[n_ij]
				dy = y_list[n_ii]-y_list[n_ij]
				ay_list[n_ij] += o2_y*(1-r0/sqrt(dx**2+dy**2))*dy
			if (j < xlen-1): # ... right edge; F_ik,ij 
				n_ik = xlen*i + (j+1)
				dx = x_list[n_ik]-x_list[n_ij]
				dy = y_list[n_ik]-y_list[n_ij]
				ay_list[n_ij] += o2_y*(1-r0/sqrt(dx**2+dy**2))*dy
			if (i > 0): # ... bottom edge; F_hj,ij
				n_hj = xlen*(i-1) + j
				dy = y_list[n_hj]-y_list[n_ij]
				ay_list[n_ij] += o2_y*(dy+r0)
			if (i < ylen-1): # ... top edge; F_jj,ij
				n_jj = xlen*(i+1) + j
				dy = y_list[n_jj]-y_list[n_ij]
				ay_list[n_ij] += o2_y*(dy-r0)
	
	return [vx_list,ax_list],[vy_list,ay_list]
		
def get_data(state,tau,steps,params,noise_func,num_update,other):
	""" Returns a list of the states of the system at each timestep, for steps
		number of iterations, tau [s] apart. num_update is the numerical method 
		function to be used. 
		
		state = [q1,q2,...,qz] where each q is a necessary parameter for 
			describing the generalized position of a specific subsystem:
			q = [[f1^(0),f2^(0),...,fn^(0)],[f1^(1),f2^(1),...fn^(1)],
				  ...,[f1^(m),f2^(m),...fn^(m)]]
		params = system-specific constants. E.g., mass, length, density, etc.
		num_update = the numerical method function to be used, with signature 
			num_update(state,dt,params,derivs).
		noise_func = the optional noise function with signature noise(t)
		other = a list of any other parameters, e.g. arguments needs by noise
	"""
	# Fix all outermost masses?
	# fix_outer = True
	
	# Extract system parameters
	kx,ky,m,r0,xlen,ylen = params
		
	# Timestep
	dt = tau
	
	# Other parameters
	if other is not None:
		n_params = other
		res_params = list(n_params)
		res_params[0] *= 10
	
	xdata = [] # the x positions
	ydata = [] # ditto, but in y
	
	# Define a noise function
	if noise_func is not None:
		def noise(t,phi):
			return noise_func(t,phi,n_params)
		def res_noise(t,phi):
			res_params[1] = -t # so that -t/a = 1 always
			return noise_func(t,phi,res_params)
		xphases = hfl.rand_list(len(state[0][0]),0,1)
		yphases = hfl.rand_list(len(state[1][0]),0,1)
	else: 
		def noise(t):
			return 0
		xphases = len(state[0][0])*[0]
		yphases = len(state[1][0])*[0]
	
	# Forward feed the solver method
	for i in range(0,steps): 
		try:
			# Create a "thermal" gradient by making the columns at opposite ends 
			# effective hot and cold reservoirs
			# xlen_m1 = xlen-1
			# state[0][0] = [(x+res_noise(i*tau,phi_x))*(i%xlen==0)+
				# (x+noise(i*tau,phi_x))*(i%xlen!=0) for (x,phi_x) in 
				# zip(state[0][0],xphases)]
			# state[1][0] = [(y+res_noise(i*tau,phi_y))*(i%xlen==0)+
				# (y+noise(i*tau,phi_y))*(i%xlen!=0) for (y,phi_y) in 
				# zip(state[1][0],yphases)]
			
			# print(len(state[0][0]))
			
			# Add noise to particles' positions
			state[0][0] = [x + noise(i*tau,phi_x) 
						   for (x,phi_x) in zip(state[0][0],xphases)]
			state[1][0] = [y + noise(i*tau,phi_y) 
						   for (y,phi_y) in zip(state[1][0],yphases)]
		
			# Update the state of the network
			state = num_update(state,dt,params,derivs)
			
			xdata.append(state[0][0]) 
			ydata.append(state[1][0])

		except ZeroDivisionError:		
			print('error at iteration ',i)
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
	state_x,state_y = [],[]
	for j in range(0,xlen): # iterate over the columns
		for i in range(0,ylen): # iterate over the rows
			rx_list[xlen*i + j] = j*r0 #+ r0*(rn()-.5)/10. # the x coord of mass ij
			ry_list[xlen*i + j] = i*r0 #+ r0*(rn()-.5)/10. # the y coord of mass ij
	state_x.append(rx_list)
	state_x.append(vx_list)
	state_x.append(ax_list)
	state_y.append(ry_list)
	state_y.append(vy_list)
	state_y.append(ay_list)
	
	state = [state_x,state_y]
	# Replace ax_list and ay_list using derivs()
	deriv_list = derivs(state,tau,params)
	state[0][2] = deriv_list[0][1]
	state[1][2] = deriv_list[1][1]
	# print("state: ", state)
	
	return state

## METHODS - MISCELLANEOUS

def sin_noise(t,phi,consts):
	""" A time-dependent system noise that can be added on to generalized 
		coordinates each iteration. 

		t = current time
		consts = a list containing:
			o = noise frequency
			a = decay constant
			c = scaling constant
		
		For a non-decaying noise, set a = -t.
	"""
	o,a,c = consts
	
	return c*sin(o*t+phi)*exp(-t/a)

## SYSTEM INITIALIZATION

# Simulation parameters
m = .1 # [kg] these are massive particles lol
kx = .01 # [N/m] Spring constant in x
ky = .01 # [N/m] Spring constant in y
r0 = .1 # [m] the spring equilibrium length
x_num =  6 # number of columns of masses
y_num = 6 # number of rows of masses
params = [kx,ky,m,r0,x_num,y_num]

# Noise (i.e. a frequency added atop the primary motion) parameters
# frequency [rad/s], decay const. [s], scaling []
n_params = [2*pi*10,(1/3.),r0/15.] 

dt = 0.01 # [s]
iters = 100#1500 # times to update the systems

# Generate the initial state
state_0 = get_initial_state(params,dt)
				 
# Generate the data
xdata,ydata = get_data(state_0,dt,iters,params,sin_noise,rk4,n_params)

## SIMULATION SETUP

# Initialize the figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_axis_off() # gif background shows up white when we do this. We're only
				  # drawing the point now... But the figure is as expected.
ax.set_aspect(aspect='equal')
ax.set_facecolor('black')
fig.patch.set_facecolor('black')
fig.add_axes(ax)

# Initialize the lines; actually, it will be more of a scatter plot
scatter = [] # this line is probably unnecessary
scatter, = ax.plot(xdata[0],ydata[0],color='blue',marker='o',linestyle='None')
#scatter, = ax.plot([],[],color='blue',marker='o',linestyle='None')

def init():
	""" Set the axes limits with global values. """
	ax.set_ylim(min(ydata[0])-r0,max(ydata[0])+r0)
	ax.set_xlim(min(xdata[0])-r0,max(xdata[0])+r0)
	return scatter,
	
def update(i):
	""" Set the ith data points with global values."""
	# The points describing the network at the ith step
	scatter.set_data(xdata[i],ydata[i])
	return scatter,

# Run the animation
anim = animation.FuncAnimation(fig, update, frames=range(0,iters), 
	init_func=init, blit=True, interval=1000*dt, repeat=True)

# # Save the animation as a gif. Note this only draws the points we plot, not the
# # background, so the gif background is white. Weird.
# plt.rcParams["animation.convert_path"] = "C:\Program Files\ImageMagick-7.0.7-Q16\magick.exe"
# anim.save('py_fast_cooling_6x6_2.gif', writer="imagemagick",extra_args="convert",fps=30)
plt.show()

