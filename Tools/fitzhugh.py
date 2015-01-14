#! /usr/bin/python

import ctypes as ct
import numpy as np
import tools as tl
import scipy.optimize as opt
import autoAdapter
import os

try:
	lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_fitzhugh_cuda.so')
	CUDA_ENABLED = True

except:
	lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_fitzhugh.so')
	CUDA_ENABLED = False

print '# CUDA_ENABLED', CUDA_ENABLED

PI2 = tl.PI2


# set hard:
V_0 = 0.                # coupling threshold
threshold_slope = 100.          # coupling threshold slope

params = dict(I_0=0.4, k_0=10., x_0=0., epsilon_0=0.3, E_0=-1.5, m_0=1., sigma_0=0.)

description = dict(I="External current", k="K+ activation slope",
			x="K+ activation shift", epsilon="time scale", E="Inh. reversal potential",
                        m="K+ current scale", sigma="noise amplitude")


params['I_b'] = params['I_0'] #
params['epsilon_b'] = params['epsilon_0'] #
params['x_b'] = params['x_0'] #
params['k_b'] = params['k_0'] #
params['E_b'] = params['E_0'] #
params['m_b'] = params['m_0'] #
params['sigma_b'] = params['sigma_0'] #

params['I_g'] = params['I_0'] #
params['epsilon_g'] = params['epsilon_0'] #
params['x_g'] = params['x_0'] #
params['k_g'] = params['k_0'] #
params['E_g'] = params['E_0'] #
params['m_g'] = params['m_0'] #
params['sigma_g'] = params['sigma_0'] #

params['I_r'] = params['I_0'] #
params['epsilon_r'] = params['epsilon_0'] #
params['x_r'] = params['x_0'] #
params['k_r'] = params['k_0'] #
params['E_r'] = params['E_0'] #
params['m_r'] = params['m_0'] #
params['sigma_r'] = params['sigma_0'] #

params['I_y'] = params['I_0'] #
params['epsilon_y'] = params['epsilon_0'] #
params['x_y'] = params['x_0'] #
params['k_y'] = params['k_0'] #
params['E_y'] = params['E_0'] #
params['m_y'] = params['m_0'] #
params['sigma_y'] = params['sigma_0'] #


def parameters_one(color='0'):
	return np.array([params['I_'+color], params['epsilon_'+color], params['x_'+color], params['k_'+color], params['E_'+color], params['m_'+color]])


def parameters_one_noise(color='0'):
	return np.array([params['I_'+color], params['epsilon_'+color], params['x_'+color], params['k_'+color], params['E_'+color], params['m_'+color], params['sigma_'+color]])



def parameters_three():
	return np.concatenate((parameters_one('b'), parameters_one('g'), parameters_one('r')))


def parameters_four():
	return np.concatenate((parameters_one('b'), parameters_one('g'), parameters_one('r'), parameters_one('y')))


def parameters_n():
        return parameters_one_noise('0')


def setParams(**kwargs):
	for k in kwargs.keys():
		for c in ['0', 'b', 'g', 'r']:
			params[k+'_'+c] = kwargs[k]

	print '# update params', kwargs

#===

N_EQ1 = 2
N_EQ3 = 3*N_EQ1
N_EQ4 = 4*N_EQ1


#=== CUDA ===#


if CUDA_ENABLED:

	lib.cuda_integrate_one_nosave.argtype = [ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.c_double, ct.c_uint, ct.c_uint]
	def cuda_integrate_one_rk4_nosave(initial_states, dt, N_integrate):
		initial_states = np.array(initial_states).flatten()
		N_initials = initial_states.size/N_EQ1 # initial states / state variables = initial conditions

		p = parameters_one()
	
		lib.cuda_integrate_one_nosave(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(initial_states.size),
					p.ctypes.data_as(ct.POINTER(ct.c_double)),
					ct.c_double(dt), ct.c_uint(N_integrate))

		return np.reshape(initial_states, (N_initials, N_EQ1), 'C')



	lib.cuda_integrate_three.argtype = [ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.c_double, ct.c_uint, ct.c_uint]
	def cuda_integrate_three_rk4(initial_states, coupling, dt, N_integrate, stride=1):
		initial_states = np.asarray(initial_states)
		N_initials = initial_states.size/N_EQ3 # initial states / state variables = initial conditions
	
		coupling = np.asarray(coupling)
		# coupling[0] : 2 -> 1
		# coupling[1] : 3 -> 1
		# coupling[2] : 1 -> 2
		# coupling[3] : 3 -> 2
		# coupling[4] : 1 -> 3
		# coupling[5] : 2 -> 3
	
		X_out = np.zeros((3*N_initials*N_integrate), float) # save only one per oscillator per time
		p = parameters_three()
	
		lib.cuda_integrate_three(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(initial_states.size),
					X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
					coupling.ctypes.data_as(ct.POINTER(ct.c_double)),
					p.ctypes.data_as(ct.POINTER(ct.c_double)),
					ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))
	
		return np.reshape(X_out, (N_initials, N_integrate, 3), 'C')



	lib.cuda_crossing_three.argtype = [ct.POINTER(ct.c_double), ct.c_uint,
                                            ct.POINTER(ct.c_double), 
                                            ct.POINTER(ct.c_double), 
                                            ct.POINTER(ct.c_double),
                                            ct.c_double, ct.c_uint,
                                            ct.c_double, ct.c_uint]
	def cuda_crossing_three(initial_states, coupling, threshold, num_crossings, dt, stride=1):
		initial_states = np.asarray(initial_states)
		N_initials = initial_states.size/N_EQ3 # initial states / state variables = initial conditions
	
		coupling = np.asarray(coupling)
	
		X_out = np.zeros((3*N_initials*num_crossings), float) # save only one per oscillator per time
		p = parameters_three()
	
		lib.cuda_crossing_three(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(initial_states.size),
					X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
					coupling.ctypes.data_as(ct.POINTER(ct.c_double)),
					p.ctypes.data_as(ct.POINTER(ct.c_double)),
                                        ct.c_double(threshold), ct.c_uint(num_crossings),
					ct.c_double(dt), ct.c_uint(stride))

                print X_out

		return np.reshape(X_out, (N_initials, num_crossings, 3), 'C')



	lib.cuda_integrate_four.argtype = [ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.POINTER(ct.c_double), ct.c_uint, ct.c_double, ct.c_uint, ct.c_uint]
	def cuda_integrate_four_rk4(initial_states, coupling, dt, N_integrate, stride=1):
		initial_states = np.asarray(initial_states)
		N_initials = initial_states.size/N_EQ4 # initial states / state variables = initial conditions
	
		coupling = np.asarray(coupling)
	
		X_out = np.zeros((4*N_initials*N_integrate), float) # save only one per oscillator per time
		p = parameters_four()
	
		lib.cuda_integrate_four(initial_states.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(initial_states.size),
					X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
					coupling.ctypes.data_as(ct.POINTER(ct.c_double)),
					p.ctypes.data_as(ct.POINTER(ct.c_double)),
					ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))
	
		return np.reshape(X_out, (N_initials, N_integrate, 4), 'C')




#=== CUDA ===#



lib.integrate_one_rk4_nosave.argtypes = [ct.POINTER(ct.c_double), 
					ct.POINTER(ct.c_double),
					ct.c_double, ct.c_uint]
def integrate_one_rk4_nosave(initial_state, dt, N_integrate, stride=42):
	initial_state = np.array(initial_state)
	assert initial_state.size == N_EQ1

	p = parameters_one()

	lib.integrate_one_rk4_nosave(initial_state.ctypes.data_as(ct.POINTER(ct.c_double)),
				p.ctypes.data_as(ct.POINTER(ct.c_double)),
				ct.c_double(dt), ct.c_uint(N_integrate))
	return initial_state



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



INITIAL_ORBIT = np.array([-0.62376542, 0.00650901])
dt = 0.05
stride = 50
N_integrate = 2*10**5
IDX_THRESHOLD = 0
THRESHOLD = 0.

def single_orbit(DT_ORBIT=dt, N_ORBIT=N_integrate, STRIDE_ORBIT=stride, V_threshold=THRESHOLD, verbose=0):

	X = integrate_one_rk4(INITIAL_ORBIT, DT_ORBIT/float(STRIDE_ORBIT), N_ORBIT, STRIDE_ORBIT)
	x_raw, y = X[0], X[1]
	x_m, y_m = tl.splineLS1D(), tl.splineLS1D()

	try:
		ni = np.asarray(tl.crossings(x_raw, V_threshold), dtype=int) # convert to millivolts
		x, y = x_raw[ni[-2]:ni[-1]], y[ni[-2]:ni[-1]]
		t = tl.PI2*np.arange(x.size)/float(x.size-1)
		x_m.makeModel(x, t); y_m.makeModel(y, t)
	
	except:
		print '# single_orbit:  No closed orbit found!'
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



lib.step_n_em.argtypes = [ct.c_uint,
                            ct.POINTER(ct.c_double),
                            ct.POINTER(ct.c_double),
			    ct.c_double,
			    ct.c_double, ct.c_uint,
                            ct.POINTER(ct.c_double)]
def step_n_em(initial_state, p, coupling_strength, dt, stride):

        num_osci = initial_state.size/N_EQ1
        
        np.random.seed()
        noise = np.sqrt(dt)*p[-1]*np.random.randn(stride*num_osci)

	lib.step_n_em(ct.c_uint(num_osci),
                initial_state.ctypes.data_as(ct.POINTER(ct.c_double)),
		p.ctypes.data_as(ct.POINTER(ct.c_double)),
                ct.c_double(coupling_strength),
		ct.c_double(dt), ct.c_uint(stride),
                noise.ctypes.data_as(ct.POINTER(ct.c_double)))

        return initial_state



#===

IDX_COUPLING = 0
def coupling_function(state, other, THRESHOLD_SLOPE=100., COUPLING_THRESHOLD=0.):
	return (params['E_0']-state)/(1.+np.exp(-THRESHOLD_SLOPE*(other-COUPLING_THRESHOLD)))



def nullcline_x(x, I, m=1., E=0., g=0.):
	return m*(x - x**3) + I + g*(E - x)


def nullcline_y(x, V_0, k):
	return 1./(1.+np.exp(-k*(x-V_0)))



def g_critical(I):

        k, E, m = params['k_0'], params['E_0'], params['m_0']

        def g_of_V(V):
                expfct = np.exp(-k*(V-V_0))
                return m*(1.-3.*V**2)-k*expfct/(1.+expfct)**2
        
        def to_minimize(V):
                return m*(V-V**3)+I+g_of_V(V)*(E-V)-nullcline_y(V, V_0, k)

        V_opt = opt.fsolve(to_minimize, x0=-1.)
        print V_opt

        return g_of_V(V_opt)



def createAutoCode(filename, period=1.):

	diffEqs_txt = """
	f[0] = par[5]*(u[0]-u[0]*u[0]*u[0])-u[1]+par[0];
	f[1] = par[1]*(1./(1.+exp(-par[3]*(u[0]-par[2])))-u[1]);
	"""

	params_txt = """
	par[0] = %lf;	// I		:	shift of fast nullcline
	par[1] = %lf;	// epsilon	:	time scale separation
	par[2] = %lf;	// x		:	slow nullcline elevation
	par[3] = %lf;	// k		:	slow nullcline inclination
	par[4] = %lf;	// E		:	unused
	par[5] = %lf;	// m		:	width of the 'Z'
	par[10] = %lf;	// period
	""" % (params['I_0'], params['epsilon_0'], params['x_0'], params['k_0'], params['E_0'], params['m_0'], period)
	
	autoAdapter.createAutoCode(filename, diffEqs_txt, params_txt)


def createAutoCode2(filename, period=1.):

	diffEqs_txt = """
	doublereal x, y, xh, yh, m, I, eps, k, v0;
	doublereal fxdx, fxdy, fydy, fydx, x2, E, Einv;

	x = u[0];
	y = u[1];
	xh = u[2];
	yh = u[3];

	I = par[0];
	eps = par[1];
	v0 = par[2];
	k = par[3];
	m = par[5];
	
	x2 = x*x;
	E = exp(-k*(x-v0));
	Einv = 1./(1.+E);

	fxdx = m*(1.-3.*x2);
	fxdy = -1.;
	fydy = -eps;
	fydx = eps* k*E * Einv*Einv;

	f[0] = m*(x-x2*x)-y+I;
	f[1] = eps*(Einv-y);
	f[2] = - fxdx*xh - fydx*yh;
	f[3] = - fxdy*xh - fydy*yh;
	"""

	params_txt = """
	par[0] = %lf;	// I		:	shift of fast nullcline
	par[1] = %lf;	// epsilon	:	time scale separation
	par[2] = %lf;	// x		:	slow nullcline elevation
	par[3] = %lf;	// k		:	slow nullcline inclination
	par[4] = %lf;	// E		:	unused
	par[5] = %lf;	// m		:	width of the 'Z'
	par[10] = %lf;	// period
	""" % (params['I_0'], params['epsilon_0'], params['x_0'], params['k_0'], params['E_0'], params['m_0'], period)
	
	autoAdapter.createAutoCode(filename, diffEqs_txt, params_txt, 'adjoint')




if __name__ == '__main__':

	from pylab import *
	import tools
	import time

	#createAutoCode('test.c')
	#exit(0)

        #print g_critical(I=0.5)
        #exit(0)

	dt = 0.02
	stride = 10
	N = 10**5
	N_initials = 100

	#coupling = 0.001*ones((12), float)
	coupling = 0.001*ones((6), float)

	X = array([INITIAL_ORBIT+0.3*randn(2) for i in xrange(100)])
	plot(X[:, 1], X[:, 0], 'ko')
	dtx, Nx = dt/float(stride), 3010*stride
	#X = array([integrate_one_rk4_nosave(X[i], dtx, Nx) for i in xrange(100)])
	X = cuda_integrate_one_rk4(X, dtx, Nx)
	plot(X[:, 1], X[:, 0], 'ro')
	show()
	exit(0)
	#print X[:, -1]
	#plot(X[0], X[1])
	#show()

	#single_orbit(verbose=1)
	#show()

	#params.update(shift=2.)
	#X = integrate_one_rk4(INITIAL_ORBIT, dt/float(stride), N, stride)
	#print X[:, -1]
	#plot(X[0], X[1])
	#show()

	X = cuda_crossing_three(randn(N_initials*N_EQ4), coupling, 0.0, 3, dt/float(stride), stride)
	#X = cuda_integrate_four_rk4(randn(N_initials*N_EQ4), coupling, dt/float(stride), N, stride)
        for i in xrange(X.shape[0]):
            ti, d = tl.compute_phase_difference(transpose(X[i]))
            tl.plot_phase_2D(d[:, 0], d[:, 1])
        show()
	exit(0)
	#X = integrate_four_rk4(randn(N_EQ4), coupling, dt/float(stride), N, stride)


	plot(X[:, 0], 'b-')
	plot(X[:, 1]-1., 'g-')
	plot(X[:, 2]-2., 'r-')
	plot(X[:, 3]-3., 'y-')
	show()
