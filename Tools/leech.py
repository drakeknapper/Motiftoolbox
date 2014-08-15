#! /usr/bin/python

import ctypes as ct
import time
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
import tools as tl
import numpy as np
import fork_master as fm
import scipy.optimize as opt
import scipy.interpolate as trp

lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_leech.so')


CUDA_ENABLED = False
try:
	lib_cuda = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_leech_cuda.so')
	CUDA_ENABLED = True
except:
	pass

#===

PI2 = 2.*np.pi

params = {}
params['C'] = 0.5	  # nF
params['I_ext'] = 0.006	  # nA

params['g_K2'] = 30.	  # nS
params['g_Na'] = 160.	  # nS 	(oder 200 nS ?)
params['g_L'] = 8.	  # nS

params['tau_K2'] = 0.9 	  # ms
params['tau_Na'] = 0.0405  # ms

params['E_Na'] = 0.045	  # V
params['E_K2'] = -0.07	  # V
params['E_L'] = -0.046	  # V
params['E_syn'] = -0.0625

params['V_Na'] = -0.0305
params['V_h'] = -0.0325
params['V_K2'] = -0.018
params['V_shift'] = -0.021 # V

params['s_Na'] = 150.
params['s_h'] = -500. 
params['s_K2'] = 83.

params['Theta'] = -0.04	  # V

N_EQ1 = 3
N_EQ3 = 3*N_EQ1
N_EQ4 = 4*N_EQ1


def params_one():
	return np.array([params['E_L'], params['E_K2'], params['E_Na'], params['g_L'] , params['g_K2'], params['g_Na'], params['C'], params['tau_Na'], params['tau_K2'], params['I_ext'], params['V_shift']])


def params_three():
	return np.array([params['E_L'], params['E_K2'], params['E_Na'], params['g_L'] , params['g_K2'], params['g_Na'], params['C'], params['tau_Na'], params['tau_K2'], params['I_ext'], params['V_shift'], params['Theta']])


#===

if CUDA_ENABLED:
	
	lib_cuda.cuda_integrate_three.argtype = [ct.POINTER(ct.c_double), ct.c_uint,
				    	ct.POINTER(ct.c_double), ct.c_uint,
				    	ct.POINTER(ct.c_double), ct.c_uint,
				    	ct.POINTER(ct.c_double), ct.c_uint,
				    	ct.c_double, ct.c_uint, ct.c_uint]
	def cuda_integrate_three_rk4(initial_states, coupling, dt, N_integrate, stride=1):
		coup = np.zeros((9), float)
		coup[:coupling.size] = coupling
	
		initial_states = np.asarray(initial_states)
		N_initials = initial_states.size/N_EQ3 # initial states / state variables
		X_out = np.zeros((3*N_initials*N_integrate), float)
		p = params_three()
	
		lib_cuda.cuda_integrate_three(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(initial_states.size),
					X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
					coup.ctypes.data_as(ct.POINTER(ct.c_double)),
					p.ctypes.data_as(ct.POINTER(ct.c_double)),
					ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))
	
		return np.reshape(X_out, (N_initials, N_integrate, 3), 'C')
	

	
	lib_cuda.cuda_relax_one.argtype = [ct.POINTER(ct.c_double), ct.c_uint,
				    	    ct.POINTER(ct.c_double),
                                            ct.c_double, ct.c_uint]
	def cuda_relax_one(initial_states, dt, N_integrate):
	
		initial_states = np.asarray(initial_states)
		N_initials = initial_states.size/N_EQ1 # initial states / state variables
		p = params_one()
	
		lib_cuda.cuda_relax_one(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(initial_states.size),
					p.ctypes.data_as(ct.POINTER(ct.c_double)),
					ct.c_double(dt), ct.c_uint(N_integrate))
	
		return np.reshape(initial_states, (N_initials, 3), 'C')
	
#===


lib.derivs_one.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]
def derivs_one(state):
	state = np.array(state, dtype=float)
	derivative = np.zeros((3), float)
	parameters = params_one()
	lib.derivs_one(state.ctypes.data_as(ct.POINTER(ct.c_double)),
			derivative.ctypes.data_as(ct.POINTER(ct.c_double)),
			parameters.ctypes.data_as(ct.POINTER(ct.c_double)))

	return derivative



lib.integrate_one_rk4.argtypes = [ct.POINTER(ct.c_double), ct.c_double, ct.c_uint, ct.c_uint, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]
def integrate_one_rk4(V0, h0, m0, dt, N_integrate, stride=1):
	X_out = np.zeros((N_EQ1*N_integrate), float)
	X_in = np.array([V0, h0, m0])
	parameters = params_one()
	lib.integrate_one_rk4(X_in.ctypes.data_as(ct.POINTER(ct.c_double)),
		ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride),
		parameters.ctypes.data_as(ct.POINTER(ct.c_double)),
		X_out.ctypes.data_as(ct.POINTER(ct.c_double)))

	return np.reshape(X_out, (N_EQ1, N_integrate), 'F')

#===

INITIAL_ORBIT = [-0.04346306, 0.99599451, 0.02006609]
DT_ORBIT = 0.001
N_ORBIT = 10**4
STRIDE_ORBIT = 1
THRESHOLD = -0.04

def single_orbit(DT_ORBIT=DT_ORBIT, N_ORBIT=N_ORBIT, STRIDE_ORBIT=STRIDE_ORBIT, V_threshold=THRESHOLD, verbose=0):
	X = integrate_one_rk4(-0.04346306, 0.99599451, 0.02006609, dt=DT_ORBIT/float(STRIDE_ORBIT), N_integrate=N_ORBIT, stride=STRIDE_ORBIT)
	V, h, m = X[0], X[1], X[2]
	Vm, Hm, Mm = tl.splineLS1D(), tl.splineLS1D(), tl.splineLS1D()

	try:
		ni = tl.crossings(V, V_threshold) # convert to millivolts
		V, h, m = V[ni[-2]:ni[-1]], h[ni[-2]:ni[-1]], m[ni[-2]:ni[-1]]
		t = PI2*np.arange(V.size)/float(V.size-1)
		Vm.makeModel(V, t); Hm.makeModel(h, t); Mm.makeModel(m, t)

	except:
		print '# single_orbit:  No closed orbit found!'
		raise ValueError

	
	T = DT_ORBIT*V.size

	return Vm, Hm, Mm, T


#===


lib.integrate_three_rk4.argtypes = [ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.c_double, ct.c_uint, ct.c_uint]
def integrate_three_rk4(initial_state, coupling, dt, N_integrate, stride=1):
	coup = np.zeros((9), float)
	coup[:coupling.size] = coupling

	initial_state = np.asarray(initial_state)
	X_out = np.zeros((3*N_integrate), float)
	parameters = params_three()

	lib.integrate_three_rk4(initial_state.ctypes.data_as(ct.POINTER(ct.c_double)),
				parameters.ctypes.data_as(ct.POINTER(ct.c_double)),
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

	coup = np.zeros((18), float)
	coup[:coupling.size] = np.asarray(coupling);

	X_out = np.zeros((4*N_integrate), float)
	p = params_three()

	lib.integrate_four_rk4(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)),
		p.ctypes.data_as(ct.POINTER(ct.c_double)),
		coup.ctypes.data_as(ct.POINTER(ct.c_double)),
		X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
		ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))

	return np.reshape(X_out, (N_integrate, 4), 'C')


#===

def activation(V, s, V_0):
	return 1./(1.+np.exp(s*(V_0-V)))


def dactivation(V, s, V_0):
	expf = np.exp(s*(V_0-V))
	return s*expf/(1.+expf)**2


def m_Na(V):
	return activation(V, params['s_Na'], params['V_Na'])

def dm_Na(V):
	return dactivation(V, params['s_Na'], params['V_Na'])

def h_Na(V):
	return activation(V, params['s_h'], params['V_h'])

def dh_Na(V):
	return dactivation(V, params['s_h'], params['V_h'])

def m_K2(V):
	return activation(V, params['s_K2'], params['V_K2']-params['V_shift'])

def dm_K2(V):
	return dactivation(V, params['s_K2'], params['V_K2']-params['V_shift'])


def sum_I_dI(X):
	(V, g_syn) = X

	I_L = params['g_L']*(V-params['E_L'])
	I_syn = g_syn*(V-params['E_syn'])

	I_K2 = params['g_K2']*m_K2(V)**2*(V-params['E_K2'])
	dI_K2 = params['g_K2']*(2.*m_K2(V)*dm_K2(V)*(V-params['E_K2']) + m_K2(V)**2)

	I_Na = params['g_Na']*m_Na(V)**3*h_Na(V)*(V-params['E_Na'])
	dI_Na = params['g_Na']*(3.*m_Na(V)**2*dm_Na(V)*h_Na(V)*(V-params['E_Na']) + m_Na(V)**3*(dh_Na(V)*(V-params['E_Na']) + h_Na(V)))

	I_sum = I_Na   + I_K2  + I_L + params['I_ext'] + I_syn
	dI_sum = dI_Na + dI_K2 + params['g_L']         + g_syn

	return np.array([I_sum, dI_sum])


def nullcline_m(V): # nullcline_m, voltages in mV
	return m_K2(V)


def nullcline_V(V, g_syn=0.): # m_K2(V), nullcline V'=0, V in volts
	I_Na = params['g_Na']*m_Na(V)**3*h_Na(V)*(V-params['E_Na'])
	I_L = params['g_L']*(V-params['E_L'])
	I_syn = g_syn*(V-params['E_syn'])
	return np.sqrt((I_Na+I_L+params['I_ext']+I_syn)/(params['g_K2']*(params['E_K2']-V)))

#===

lib.step_three_rk4.argtypes = [ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.c_double, ct.c_uint]
def step_three_rk4(state, parameters, coupling, dt, stride=1):
	assert state.size == N_EQ3
	assert coupling.size == 9

	lib.step_three_rk4(state.ctypes.data_as(ct.POINTER(ct.c_double)),
				parameters.ctypes.data_as(ct.POINTER(ct.c_double)),
				coupling.ctypes.data_as(ct.POINTER(ct.c_double)),
				ct.c_double(dt), ct.c_uint(stride))

	return state




if __name__ == '__main__':

	from pylab import *
	import time

	dt = 0.03
	stride = 1
	N = 10**3
        N_osci = 100
        q = 0.001

	coupling = q*ones((9), float)
	parameters = params_three()

	print step_three_rk4(array([0., 1., 0., 0., 1., 0., 0., 1., 0.]), parameters, coupling=coupling, dt=dt/float(stride), stride=stride)
	exit(0)


	X = integrate_one_rk4(0., 1., 0., dt=dt/float(stride), N_integrate=N, stride=stride)
	ti = tl.crossings(X[0], -0.04)
	#plot(X[0])
	#for t in ti:
	#	axvline(x=t)
	#show()
	#exit(0)
        initial_states = X[:, :ti[0]]
	#plot(initial_states[0])
        #show()

        initial_states = transpose(initial_states).flatten()
        plot(initial_states[::3])
        show()


        X = cuda_relax_one(initial_states, dt/float(stride), N)
        plot(X[:, 0])
        show()
        subplot(211)
        plot(X[:, 1], X[:, 0], 'k.')
        subplot(212)
        plot(X[:, 2], X[:, 0], 'k.')
        show()
        exit(0)

        #X = cuda_integrate_three_rk4(0.1*randn(N_osci*N_EQ3), q*ones((6), float), dt/float(stride), N, stride)
        #plot(X[0, :, 0], 'b-')
        #plot(X[0, :, 1]-0.06, 'g-')
        #plot(X[0, :, 2]-0.12, 'r-')
        #show()
        #exit(0)


	X = 1000.*integrate_four_rk4(0.1*randn(12), coupling=zeros((12), float), dt=dt/float(stride), N_integrate=N, stride=stride)
	t = dt*arange(X.shape[0])

	plot(t, X[:, 0], 'b-')
	plot(t, X[:, 1]-60., 'g-')
	plot(t, X[:, 2]-120., 'r-')
	plot(t, X[:, 3]-180., 'y-')
	ylim(-250., 50)
	show()

