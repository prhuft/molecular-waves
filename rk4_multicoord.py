""" 
Runge-Kutta Numerical ODE Solving Method - Multiple Coordinates Version

Preston Huft, Spring 2018.

Notes: This solver is for systems whose motion is singly-parameterizeable, i.e. 
one variable is sufficient for describing the system's generalized position. 
To use this solver for systems which are not singly-parameterizeable, this
method can be called for updating each generalized coordinate separately. 

To-do: use better nomenclature. 'state' should be replaced with 'coord' or 
'state1D' throughout to specify that it is a state in only one parameter.
"""

## LIBRARIES

## GLOBAL VARIABLES

DEBUG = 1

## THE METHOD

def rk4_update(state,h,params,derivatives):
	""" state = [q1,q2,...,qz] where each q is a necessary parameter for 
			describing the generalized position of a specific subsystem:
			q = [[f1^(0),f2^(0),...,fn^(0)],[f1^(1),f2^(1),...fn^(1)],
				  ...,[f1^(m),f2^(m),...fn^(m)]]
				
				where fi^(j) denotes the jth derivative of the generalized
				position of the ith object in the system.
		h = step size to be used in derivative calculations. 
		params = a list of system parameters which required by the 
			derivatives method; see the specific derivatives function for more
			information.
		derivatives = derivatives(state,h,params)
			a method which returns a list of the m derivatives of 
			the generalized position of the n bodies in the system, i.e.:
			[[[fx1^(1),fx2^(1),...,fxn^(1)], ...,[fx1^(m),fx2^(m),...,fxn^(m)]],
			 [[fy1^(1),fy2^(1),...,fyn^(1)], ...,[fy1^(m),fy2^(m),...,fyn^(m)],
			 ...] where lists of lists of ordered derivatives are grouped by 
			 coordinate.
			
		This method returns the updated state of the system, in the exact form 
		of the 'state' list which is passed in. 
	"""
	
	def k_mat(dh_mat,u): 
		""" The output is of the form: 
			k_mat = [[f1^(0)(q+dh),f2^(0)(q+sh),...,fn^(0)(q+dh)],[f1^(1)(q+dh),
					f2^(1)(q+dh),...fn^(1)(q+dh)],...,[f1^(m-1)(q+dh),f2^(m-1)(q+dh),
					...fn^(m-1)(q+h)]]
					
			where fi^(j)(q+h) is the jth derivative of the generalized position
			of the ith object, evaluated at q+h.
			u = the current coordinate we're dealing with; the uth coordinate.
			"""
		def copy_state(state):
			temp = []
			for substate in state:
				temp.append(copy_nested_list(substate))
			return temp
			
		# if DEBUG: 
			# print("state[u]: ",u)
		
		k_mat = []
		for i in range(0,len(state[u][:-1])): # iter over orders of derivs up to m-1, inclusive
			k_list = []
			for j in range(0,len(state[u][i])): # iter through the object indices
				temp_s = copy_state(state)
				temp_s[u][i][j] = state[u][i][j]+dh_mat[i][j] # add dh to the ith f^(j)

				k_list.append(h*derivatives(temp_s,h,params)[u][i][j])
			k_mat.append(list(k_list)) # append a copy of l; this is the mth k list
		return k_mat
				
	def dh_mat(k_mat,c):
		""" Returns the matrix of steps dh, given the previous matrix of k vals. 
			The output is in the format of input expected by k_mat(), i.e.:
			[[h1^1,h2^1,...,hn^1],[h1^2,h2^2,...,hn^2],...,[h1^m,h2^m,...,hn^m]]
			
			where hi^j = dh is the step in fi^(j)(q+dh)."""
		dh_mat = []
		for k_list in k_mat:
			dh_list = []
			for k in k_list:
				dh_list.append(k*c)
			dh_mat.append(list(dh_list))
		return dh_mat
		
	def copy_nested_list(mylist):
		""" This assumes a depth of 1, i.e. a list of list of non-lists."""
		copy_list = []
		for l in mylist:
			copy_list.append(list(l))
		return copy_list
		
	# Run the single-parameter solver of rk4.py for each parameter
	new_state = []
	for u in range(0,len(state)):
		# Create the initial list of zeros to pass into k()
		dh_mat_0 = []
		for i in range(0,len(state[u][:-1])): # up to to m-1, inclusive
			dh_mat_0.append([0]*len(state[u][i]))
			
		# Generate the k matrices, where k(dh) def: f(q_i+dh) for systems which do 
		# not depend explicitly on the parameter (e.g. time) with respect to which 
		# we are differentiating.
		k1 = k_mat(dh_mat_0,u)
		k2 = k_mat(dh_mat(k1,1/2.),u)
		k3 = k_mat(dh_mat(k2,1/2.),u)
		k4 = k_mat(k3,u) # passing in k3 is equivalent to passing in dh_mat(k3,1)

		new_state_u = [] # the part of the state projected on the uth coordinate
		for i in range(0,len(state[u][:-1])): # iter over deriv orders up to m-1, inclusive
			new_derivs = []
			for j in range(0,len(state[u][0])): # iter over coords up to q_n, inclusive
				new_derivs.append(state[u][i][j]+(k1[i][j]+2*(k2[i][j]+k3[i][j])+
									k4[i][j])/6.)
			new_state_u.append(list(new_derivs))
			
		# Get the mth derivatives with the analytical method passed in
		new_state_u.append(derivatives(state,h,params)[u][-1])
		new_state.append(new_state_u)
			
	return new_state