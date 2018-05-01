""" 
Runge-Kutta Numerical ODE Solving Method, v1.00

Preston Huft, Spring 2018.

Known bugs:
- I think there may be an error in the multiplying of k values somewhere, or in
adding the appropriate step values to the state variables in calculating k_mat.
The double_pendulum_1-04 simulation gives a substantially different (and severely
unrealistic, no less) result from double_pendulum_1-03 which uses rk4_two_bodies. 

"""

## LIBRARIES

## GLOBAL VARIABLES

DEBUG = 0

## THE METHOD

def rk4_update(state,h,params,derivatives):
	""" state = [[f1^(0),f2^(0),...,fn^(0)],[f1^(1),f2^(1),...fn^(1)],
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
			[[f1^(1),f2^(1),...fn^(1)], ...,[f1^(m),f2^(m),...fn^(m)]]
			
		This method returns the updated state of the system, in the exact form 
		of the 'state' list which is passed in. 
	"""
		
	def k_mat(dh_mat):
		""" The output is of the form: 
			k_mat = [[f1^(0)(q+dh),f2^(0)(q+sh),...,fn^(0)(q+dh)],[f1^(1)(q+dh),
					f2^(1)(q+dh),...fn^(1)(q+dh)],...,[f1^(m-1)(q+dh),f2^(m-1)(q+dh),
					...fn^(m-1)(q+h)]]
					
			where fi^(j)(q+h) is the jth derivative of the generalized position
			of the ith object, evaluated at q+h."""
		k_mat = []
		for i in range(0,len(state[:-1])): # iter over orders of derivs up to m-1, inclusive
			k_list = []
			for j in range(0,len(state[i])): # iter through the object indices
				temp_s = copy_nested_list(state)
				temp_s[i][j] = state[i][j]+dh_mat[i][j] # add dh to the ith f^(j)
				k_list.append(h*derivatives(temp_s,h,params)[i][j])
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
	
	# Create the initial list of zeros to pass into k()
	dh_mat_0 = []
	for i in range(0,len(state[:-1])): # up to to m-1, inclusive
		dh_mat_0.append([0]*len(state[i]))
		
	# Generate the k matrices, where k(dh) def: f(q_i+dh) for systems which do 
	# not depend explicitly on the parameter (e.g. time) with respect to which 
	# we are differentiating.
	k1 = k_mat(dh_mat_0)
	k2 = k_mat(dh_mat(k1,1/2.))
	k3 = k_mat(dh_mat(k2,1/2.))
	k4 = k_mat(k3) # passing in k3 is equivalent to passing in dh_mat(k3,1)
	
	if DEBUG: 
		print('dh1: ',dh_mat_0)
		print('k1: ',k1)
		print('dh2: ',dh_mat(k1,1/2.))
		print('k2: ',k2)
		print('dh3: ',dh_mat(k2,1/2.))
		print('k3: ',k3)
		print('dh4: ',k3)
		print('k4: ',k4)
	
	new_state = []
	for i in range(0,len(state[:-1])): # iter over deriv orders up to m-1, inclusive
		new_derivs = []
		for j in range(0,len(state[0])): # iter over coords up to q_n, inclusive
			new_derivs.append(state[i][j]+(k1[i][j]+2*(k2[i][j]+k3[i][j])+
								k4[i][j])/6.)
		new_state.append(list(new_derivs))
		
	# Get the mth derivatives with the analytical method passed in
	new_state.append(derivatives(state,h,params)[-1])
		
	return new_state