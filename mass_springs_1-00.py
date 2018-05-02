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

def derivs(state,tau,params):
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
	
	# # The spring constant,mass, and network dimensions
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
	
def get_initial_states(params,state_template,dstate,sys_size):
	""" Returns a list of length 'sys_size', the elements of which are the 
		initial states for each double pendulum system. That is, 
		states = [[state1_0],[state2_0],...,[statem_0
		where statei_0 is the initial state of the ith double pendulum, and 
		m = sys_size. 
		dstate = the initial difference between adjacent systems
		state_template = the state of one system; each other state is built
		by adding adding dstate scaled by an integer in (1,sys_size).
		sys_size = the number of systems to generate."""
	# l1,l2 = params[2:]
	
	# state = hfl.copy_nested_list(state_template)
	# states_0 = [state_template]
	# for i in range(0,sys_size-1): # append n-1 one times for system of size n
		# state = hfl.add_nested_lists(state,dstate) # add dstate each iteration
		# # overwrite [a1,a2] which depend on the initial angle
		# state[2] = [alpha_init(state[0][0],l1),alpha_init(state[0][1],l2)]
		# states_0.append(state)
	# if DEBUG: [print(x) for x in states_0]
	return states_0
	
def alpha_init(theta,length):
	""" Returns the initial angular acceleration due only to gravitational 
	torque exerted on the center of mass of a uniform arm of length 'length'
	about one end, held at angle theta from the vertical."""
	return -6*g*sin(theta)/length
	
def toXY(theta,length):
	""" Returns the (x,y) coordinate of the end of a single pendulum arm."""
	return length*sin(theta), -length*cos(theta)

## SYSTEM INITIALIZATION

# # Simulation parameters - assume each double pendulum identical
# m1 = 1 #1 # [kg]
# m2 = .5 # [kg]
# l1 = 1 # [m]
# l2 = 1 # [m]
# params = [m1,m2,l1,l2]

# # The state variables for one double pendulum
# t1_0 = m.pi/2 # [rad] from vertical
# t2_0 = m.pi/2 # ditto
# o1_0 = 0 # [rad/s] 
# o2_0 = 0 # ditto
# a1_0 = alpha_init(t1_0,l1) # [rad/s^2]
# a2_0 = alpha_init(t2_0,l2) # ditto

# # The difference in initial variables between "adjacent" systems
# dt1 = m.pi/360
# dt2 = m.pi/360
# do1 = 0
# do2 = 0
# da1 = 0 # a1_0 has to be calculated for each system in get_initial_states()
# da2 = 0 # same. just leave as 0 for now

# # Initial variables of one double_pendulum, grouped by derivative order. 
# state1 = [[t1_0,t2_0],[o1_0,o2_0],[a1_0,a2_0]]

# # The initial state difference between "adjacent" systems
# delta_state = [[dt1,dt2],[do1,do2],[da1,da2]]

# # Generate the initial state for each double pendulum system
# states_0 = get_initial_states(params,state1,delta_state,total)

# total = 10 # number of double pendulum systems
# dt = 0.01 # [s]
# iters = 10000 # times to update the systems

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

