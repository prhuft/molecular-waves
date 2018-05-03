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
from rk4 import rk4_update as rk4
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np
from copy import copy as cp
import math as m
from math import cos
from math import sin
from math import sqrt

## GLOBAL CONSTANTS

g = 9.8 # [m/s]
DEBUG = 0 # set to 1 if we're debugging

## METHODS

def derivs(state,params):
	"""Returns a list of the first and second derivatives of the generalized 
	positions of each mass in the system.
	
	'state' is a list of the lists of the derivatives of generalized coordinates
	for each object in the system, grouped by derivative order. There are lists
	of these derivative list for each coordinate need to describe the system:
	
	state = [coordinate1,coordinate2,...,coordinateq], where
		coordinate1 = [[f1^(0),f2^(0),...,fn^(0)],[f1^(1),f2^(1),...fn^(1)],
				...,[f1^(m),f2^(m),...fn^(m)]]
		coordinate2 = likewise, and so on for the other coordinates
	"""
	# The mass positions. state[i][j] is the jth deriv group for coord. i
	x_list = state[0][0]
	y_list = state[1][0]
	vx_list = state[0][1]
	vy_list = state[1][1]
	
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
				ax[n_ij] += o2_x*(1-r0/sqrt(dx**2+dy**2))*dx
			if (i < xlen-1): # ... bottom edge
				n_jj = xlen*(i+1) + j
				dx = x_list[n_jj]-x_list[n_ij]
				dy = y_list[n_jj]-y_list[n_ij]
				ax[n_ij] += o2_x*(1-r0/sqrt(dx**2+dy**2))*dx
			if (j > 0): # ... left edge
				n_ii = xlen*i + (j-1)
				dx = x_list[n_ii]-x_list[n_ij]
				ax[n_ij] += -o2_x*(dx+r0)
			if (j < ylen-1); # ... right edge
				n_ik = xlen*i + (j+1))
				dx = x_list[n_ik]-x_list[n_ij]
				ax[n_ij] += o2_x*(dx-r0)
	
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
				ay[n_ij] += o2_y*(1-r0/sqrt(dx**2+dy**2))*dy
			if (j < ylen-1); # ... right edge 
				n_ik = xlen*i + (j+1))
				dx = x_list[n_ik]-x_list[n_ij]
				dy = y_list[n_ik]-y_list[n_ij]
				ay[n_ij] += o2_y*(1-r0/sqrt(dx**2+dy**2))*dy
			if (i > 0): # ... top edge
				n_hj = xlen*(i-1) + j
				dy = y_list[n_hj]-y_list[n_ij]
				ay[n_ij] += -o2_y*(dy+r0)
			if (i < xlen-1): # ... bottom edge
				n_jj = xlen*(i+1) + j
				dy = y_list[n_ik]-y_list[n_ij]
				ay[n_ij] += o2_y*(dy-r0)
	
	return [vx_list,ax_list],[vy_list,ay_list]
		
# def get_data(state,tau,steps,params,num_update):
def get_data(states,tau,steps,params,num_update):
	""" Returns a list of the states of the double pendulum systems at each
	timestep, for steps number of iterations, dt [s] apart. num_update is the
	numerical method function to be used. """
	
	# Get pendulum parameters; assume same params for each system
	# m1,m2,l1,l2 = params 
	
	# # The arrays of endpoints of arm1 and arm2 for each system. I.e., 
	# # arm1_data = [[[x1_1,y1_1],[x1_2,y1_2],...,[x1_n,y1_n]],
	# #				 [[x2_1,y2_1],[x2_2,y2_2],...,[x2_n,y2_n]],...,
	# #				 [[xm_1,y2_1],[xm_2,y2_2],...,[xm_n,ym_n]]]
	# # where [xi_j,yi_j] is the endpoint coordinate of arm1 in the ith 
	# # double pendulum system at step j. Likewise for arm2.
	# arm1_data,arm2_data = [],[] 
	
	# Generate the states at each timestep over the specified number of
	# steps for each double pendulum system
	for state in states:
	
		# The initial state of the nth double pendulum
		# t1,t2 = state[0]
		# o1,o2 = state[1]
		# a1,a2 = state[2]
		
		# if DEBUG:
			# print('Iter 0',': t1,t2= ',t1,t2)
			# print('o1,o2= ',o1,o2,' a1,a2= ',a1,a2)
			# print()
		
		# Timestep
		dt = tau
		
		# # Initialize the endpoints of each arm, in Cartesian coordinates
		# xData1,yData1 = [toXY(t1,l1)[0]],[toXY(t1,l1)[1]]
		# xData2,yData2 = [toXY(t2,l2)[0]+xData1[0]],[toXY(t2,l2)[1]+yData1[0]]
		
		# Forward feed the solver method for i = 0 to i = steps
		for i in range(0,steps): 
			try:
				# Update each variable
				# new_state = num_update([[t1,t2],[o1,o2],[a1,a2]],dt,params,derivs)
				# t1,t2 = new_state[0]
				# o1,o2 = new_state[1]
				# a1,a2 = new_state[2]
				
				# xData1 += [toXY(t1,l1)[0]]
				# yData1 += [toXY(t1,l1)[1]]
				# xData2 += [toXY(t2,l2)[0]+xData1[i+1]]
				# yData2 += [toXY(t2,l2)[1]+yData1[i+1]]
				
				# if DEBUG:
					# print('Iter ',i,': t1,t2= ',t1,t2)
					# print('o1,o2= ',o1,o2,' a1,a2= ',a1,a2)
					# print()
				
			except ValueError:		
				print('value error at iteration ',i)
				break
			
		# Append the coordinate lists for the nth system  
		# arm1_data.append([xData1,yData1]) # endpoints of arm 1
		# arm2_data.append([xData2,yData2]) # endpoints of arm 2
			
	return arm1_data,arm2_data
	
def get_initial_state(params):
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
			rx_list[xlen*i + (j-1)] = j*ylen + r0*(rn()-.5) # the y coord of mass ij
	state_y.append(ry_list)
	state_y.append(vy_list)
	state_y.append(ay_list)
		
	state = [state_x,state_y]
	# Replace ax_list and ay_list using derivs()
	deriv_list = derivs(state,params)
	state[0][2] = deriv_list[0][1]
	state[1][2] = deriv_list[1][1]
	
	return states
	
## SYSTEM INITIALIZATION

# Simulation parameters
m = .1 # [kg] these are massive particles lol
kx = 1 # [N/m] Spring constant in x
ky = 1 # [N/m] Spring constant in y
r0 = .1 # [m] the spring equilibrium length
x_num =  10 # number of columns of masses
y_num = 10 # number of rows of masses
params = [kx,ky,m,r0,x_num,y_num]

# Generate the initial state
state_0 = get_initial_state(params)

dt = 0.01 # [s]
iters = 10000 # times to update the systems

# # Generate the data
# data = get_data(states_0,dt,iters,params,rk4)

## SIMULATION SETUP

# Initialize the figure
# xpts,ypts = [],[]
# for armdata in data: # iterates twice
	# for data_n in armdata: # iterate through the data for all n systems
		# xpts.append(data_n[0])
		# ypts.append(data_n[1])
# # xpts and ypts should now each be of length 2n
# n = int(len(xpts)/2)
# x1pts,x2pts = xpts[:n],xpts[n:]
# y1pts,y2pts = ypts[:n],ypts[n:]

# if DEBUG: print("pt. lengths: ",len(x1pts),len(x2pts),len(y1pts),len(y2pts))

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_facecolor('black')
# ax.set_aspect(aspect='equal')
# fig.patch.set_facecolor('black')

# # Initialize the lines to be plotted 
# pen_lines = []
# trail1_lines = []
# trail2_lines = []
# for k in range(0,total):
	# trail1_line, = ax.plot([],[],lw=1)
	# trail2_line, = ax.plot([],[],lw=1)
	# trail1_lines.append(cp(trail1_line))
	# trail2_lines.append(cp(trail2_line))
	# pen_line, = ax.plot([],[],color='white',lw=3)
	# pen_lines.append(cp(pen_line))

def init():
	""" Set the axes limits. """
	# l = 0
	# if (params[2] > params[3]):
		# l = params[2]
	# else:
		# l = params[3]

	# ax.set_ylim(-2*l*1.1,2*l*1.1)
	# ax.set_xlim(-2*l*1.1,2*l*1.1)
	# # Concatenate the lists, then return as a tuple
	return tuple(trail1_lines + trail2_lines + pen_lines)
	
def update(i):
	""" Uses values established previously as globals."""
	# j = i + 1;
	# # Set the lines to plot
	# for k in range(0,total):
		# # The lines from the arm endpoints
		# trail1_lines[k].set_data(x1pts[k][:j],y1pts[k][:j])
		# trail2_lines[k].set_data(x2pts[k][:j],y2pts[k][:j])
		
		# # The line describing the kth double pendulum
		# pen_lines[k].set_data([0,x1pts[k][i],x2pts[k][i]],
							  # [0,y1pts[k][i],y2pts[k][i]])
	# # Concatenate the lists, then return as a tuple
	return tuple(trail1_lines + trail2_lines + pen_lines)

# Run the animation
anim = animation.FuncAnimation(fig, update, frames=range(0,iters+1), 
	init_func=init, blit=True, interval=1000*dt, repeat=False)
plt.show()

