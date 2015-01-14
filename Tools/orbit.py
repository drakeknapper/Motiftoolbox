#! /usr/bin/python

import numpy as np
import tools as tl
import scipy.optimize as opt
import autoAdapter

try:
	import auto
	import parseS
	AUTO_ENABLED = True
except:
	AUTO_ENABLED = False

print "# AUTO_ENABLED", AUTO_ENABLED


class orbit(object):

	def __init__(self, model):
		self.set_model(model)


	def evaluate_orbit(self, phase):
		return np.array([self.trajectory[j](phase) for j in xrange(self.dimensions)])


	def set_model(self, model):
		self.model = model
		self.integrate = model.integrate_one_rk4	# generates a single trace of all variables
		self.dimensions = model.N_EQ1
		self.trajectory = [tl.splineLS1D(isphase=False) for i in xrange(self.dimensions)]
		self.ORBIT_COMPUTED = False

		self.dt = model.dt
		self.stride = model.stride
		self.N_integrate = model.N_integrate


	def setParams(self, **kwargs):
		self.model.setParams(**kwargs)
		self.ORBIT_COMPUTED = False


	def find_orbit(self):

		if self.ORBIT_COMPUTED:
			return

		X = self.integrate(self.model.INITIAL_ORBIT, self.dt/float(self.stride), self.N_integrate, self.stride)
		assert X.shape[0] == self.dimensions
	
		try:
			ni = tl.crossings(X[self.model.IDX_THRESHOLD], self.model.THRESHOLD) # convert to millivolts
			X = X[:, ni[-2]:ni[-1]]
	
			self.period = self.dt*(ni[-1]-ni[-2])

			phase = tl.PI2*np.arange(X.shape[1])/float(X.shape[1]-1)
	
			if AUTO_ENABLED: # if enabled, find the exact solution
				t = self.period*np.arange(X.shape[1])/float(X.shape[1]-1)
				self.write_solution(t, X, mode=0)
				phase, X = self.auto_orbit()
	
			# splinefit the trajectory
			for i in xrange(self.dimensions):
				self.trajectory[i].makeModel(X[i], phase)
		
			self.ORBIT_COMPUTED = True
			
		except:
			print '# single_orbit:  No closed orbit found!'
			raise ValueError
	


	def write_solution(self, t, X, mode=0): # mode=0: orbit, mode=1: adjoint
		f = open('orbit.dat', 'w+')
		(dim, N) = X.shape

		for i in xrange(N):
			string = repr(t[i])
			for j in xrange(dim):
				string += '\t'+repr(X[j, i])

			if mode:
				for j in xrange(dim):
					string += '\t1.'

			f.write(string+'\n')

		f.close()

		if mode:	mode = 'adjoint'
		else:		mode = 'orbit'

		autoAdapter.writeConstantsFile('c.orbit', NDIM=dim, NTST=512, mode=mode)



	def auto_orbit(self):
		self.model.createAutoCode('orbit.c', period=self.period)
		solution = auto.run('orbit')(2)	# the last label
		self.period = solution["p"](11)	# save new period
		tX = np.array(solution.toArray())
		phase = tl.PI2*tX[:, 0]
		X = tX[:, 1:] # save new solution
		return phase, np.transpose(X)



	def shifted_orbit(self, kick=0.01):

		self.find_orbit()
		phase = tl.PI2*np.arange(self.N_PRC)/float(self.N_PRC)
		O = self.evaluate_orbit(phase)

		shifted_states = []
		for dim in xrange(self.dimensions):	# shift trajectories in each direction (dimension)
			x_dim = O.copy()
			x_dim[dim, :] = x_dim[dim, :]+kick
			shifted_states.append(x_dim.transpose())

		return phase, shifted_states



	def closestPhase(self, state, phaseGuess=0.):

		def distance(phase):
			return np.max( np.abs(self.evaluate_orbit(phase[0])-state) )

		return opt.fmin(distance, [phaseGuess], disp=0)[0]




	def show(self, ax=None):
		from pylab import subplot, plot, show

		if ax == None:
			ax = subplot(111)

		phase = np.arange(0., 2.*np.pi+0.01, 0.01)

		X = self.evaluate_orbit(phase)
		ax.plot(phase, X[0])
		ax.plot(phase, X[1])
		




if __name__ == '__main__':

	import fitzhugh as model
	from pylab import *

	orb = orbit(model)
	orb.find_orbit()

	phase = 1.0
	state = orb.evaluate_orbit(phase)
	kickedstate = state.copy()
	kickedstate[1] = kickedstate[1]+0.1

	bestphase = orb.closestPhase(kickedstate, phase-0.3)

	phi = np.arange(0., 2.*np.pi+0.01, 0.01)
	X = orb.evaluate_orbit(phi)
	#subplot(211)
	#plot(phi, X[0], 'k-')
	#subplot(212)
	#plot(phi, X[1], 'k-')
	#show()
	#exit(0)
	plot(X[1], X[0])
	plot([state[1]], [state[0]], 'ko')
	plot([kickedstate[1]], [kickedstate[0]], 'ko')
	x = orb.evaluate_orbit(bestphase)
	plot([x[1]], [x[0]], 'rs')

	show()
