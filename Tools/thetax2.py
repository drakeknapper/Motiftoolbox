#! /usr/bin/python

import ctypes as ct
import numpy as np
import tools as tl
import os

try:
	lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/thetax2_cuda.so')
	CUDA_ENABLED = True

except:
	lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_thetax2.so')
	CUDA_ENABLED = False

print '# CUDA_ENABLED', CUDA_ENABLED

PI2 = tl.PI2


params = dict(omega_0=1.15, beta_0=0.07)

params['omega_b'] = params['omega_0'] #
params['beta_b'] = params['beta_0'] #

params['omega_g'] = params['omega_0'] #
params['beta_g'] = params['beta_0'] #

params['omega_r'] = params['omega_0'] #
params['beta_r'] = params['beta_0'] #

params['omega_y'] = params['omega_0'] #
params['beta_y'] = params['beta_0'] #


def parameters_one(color='0'):
	return np.array([params['omega_'+color], params['beta_'+color]])


def parameters_three():
	return np.concatenate((parameters_one('b'), parameters_one('g'), parameters_one('r')))


def parameters_four():
	return np.concatenate((parameters_one('b'), parameters_one('g'), parameters_one('r'), parameters_one('y')))


#===

N_EQ1 = 1
N_EQ3 = 3*N_EQ1
N_EQ4 = 4*N_EQ1


if CUDA_ENABLED:
	lib.cuda_integrate_three.argtype = [ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.c_double, ct.c_uint, ct.c_uint]
	def cuda_integrate_three_rk4(initial_states, coupling, dt, N_integrate, stride=1):
		initial_states = np.asarray(initial_states)
		N_initials = initial_states.size/N_EQ3 # initial states / state variables = initial conditions
	
		coupling = np.asarray(coupling)
	
		X_out = np.zeros((3*N_initials*N_integrate), float) # save only one per oscillator per time
		p = parameters_three()
	
		lib.cuda_integrate_three(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(initial_states.size),
					X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
					coupling.ctypes.data_as(ct.POINTER(ct.c_double)),
					p.ctypes.data_as(ct.POINTER(ct.c_double)),
					ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))
	
		return np.reshape(X_out, (N_initials, N_integrate, 3), 'C')



#===


lib.integrate_one_rk4.argtypes = [ct.POINTER(ct.c_double), 
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.c_double, ct.c_uint, ct.c_uint]
def integrate_one_rk4(initial_state, dt, N_integrate, stride=42):
	initial_state = np.asarray(initial_state)
	assert initial_state.size == N_EQ1

	X_out = np.zeros((N_EQ1*N_integrate), float)
	p = parameters_one()

	lib.integrate_one_rk4(initial_state.ctypes.data_as(ct.POINTER(ct.c_double)),
				p.ctypes.data_as(ct.POINTER(ct.c_double)),
				X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
				ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))
	return np.reshape(X_out, (N_EQ1, N_integrate), 'F')


INITIAL_ORBIT = [0.0]

def single_orbit(DT_ORBIT=0.05, N_ORBIT=5*10**4, STRIDE_ORBIT=10, V_threshold=0., verbose=0):

	X = integrate_one_rk4(INITIAL_ORBIT, DT_ORBIT/float(STRIDE_ORBIT), N_ORBIT, STRIDE_ORBIT)
	x_raw, y = np.cos(X[0]), np.sin(X[0])
	x_m, y_m = tl.splineLS1D(),  tl.splineLS1D()

	try:
		ni = tl.crossings(x_raw, V_threshold) # convert to millivolts
		x, y = x_raw[ni[-2]:ni[-1]], y[ni[-2]:ni[-1]]
		t = tl.PI2*np.arange(x.size)/float(x.size-1)
		x_m.makeModel(x, t)
		y_m.makeModel(y, t)
	
	except:
		print '# single_orbit:  No closed orbit found!'
		T = DT_ORBIT*x_raw.size	 # in msec.
		x = x_raw
		raise ValueError

	T = DT_ORBIT*x.size	 # in msec.

	return x_m, y_m, T


#===


lib.integrate_three_rk4.argtypes = [ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.c_double, ct.c_uint, ct.c_uint]
def integrate_three_rk4(initial_states, coupling, dt, N_integrate, stride=1):
	initial_states = np.asarray(initial_states) #
	assert initial_states.size == N_EQ3

	coupling = np.asarray(coupling)
	assert coupling.size == 6

	X_out = np.zeros((3*N_integrate), float)
	p = parameters_three()

	lib.integrate_three_rk4(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)),
		p.ctypes.data_as(ct.POINTER(ct.c_double)),
		coupling.ctypes.data_as(ct.POINTER(ct.c_double)),
		X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
		ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))

	return np.reshape(X_out, (N_integrate, 3), 'C')


lib.integrate_four_rk4.argtypes = [ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.c_double, ct.c_uint, ct.c_uint]
def integrate_four_rk4(initial_states, coupling, dt, N_integrate, stride=1):
	initial_states = np.asarray(initial_states) #
	assert initial_states.size == N_EQ4

	coupling = np.asarray(coupling)
	assert coupling.size == 12

	X_out = np.zeros((4*N_integrate), float)
	p = parameters_four()

	lib.integrate_four_rk4(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)),
		p.ctypes.data_as(ct.POINTER(ct.c_double)),
		coupling.ctypes.data_as(ct.POINTER(ct.c_double)),
		X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
		ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))

	return np.reshape(X_out, (N_integrate, 4), 'C')



def right_hand_side(phase, color='0'):
	p = parameters_one(color=color)
	return p[0]-p[1]*np.cos(phase)-np.cos(2.*phase)


#===


if __name__ == '__main__':

	from pylab import *
	import tools
	import time

	dt = 0.02
	stride = 10
	N = 10**5
	N_initials = 100

	coupling = 0.01*ones((6), float)

	"""
	X = integrate_one_rk4(INITIAL_ORBIT, dt/float(stride), N, stride)
	print X
	plot(X[0])
	show()

	phase = arange(0, 2.*pi, 0.01)
	x_m, y_m, T = single_orbit(verbose=0)
	print "period", T
	ax = subplot(211)
	plot(x_m(phase), 'k-')
	subplot(212, sharex=ax)
	plot(y_m(phase), 'k-')
	show()
	#"""

	X = integrate_three_rk4(2.*pi*rand(N_EQ3), coupling, dt/float(stride), N, stride)


	plot(X[:, 0], 'b-')
	plot(X[:, 1]-pi, 'g-')
	plot(X[:, 2]-2.*pi, 'r-')
	show()
