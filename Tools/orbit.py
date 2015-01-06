#! /usr/bin/python

import numpy as np
import tools as tl
import scipy.optimize as opt


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


	def find_orbit(self):

		if self.ORBIT_COMPUTED:
			return

		X = self.integrate(self.model.INITIAL_ORBIT, self.dt/float(self.stride), self.N_integrate, self.stride)
		assert X.shape[0] == self.dimensions
	
		try:
			ni = tl.crossings(X[self.model.IDX_THRESHOLD], self.model.THRESHOLD) # convert to millivolts
			X = X[:, ni[-2]:ni[-1]]
			phase = tl.PI2*np.arange(X.shape[1])/float(X.shape[1]-1)
	
			for i in xrange(self.dimensions):
				self.trajectory[i].makeModel(X[i], phase)

			self.period = self.dt*(ni[-1]-ni[-2])

			self.ORBIT_COMPUTED = True
		
		except:
			print '# single_orbit:  No closed orbit found!'
			raise ValueError



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
		self.find_orbit()

		def distance(phase):
			return np.sqrt( ((self.evaluate_orbit(phase[0])-state)**2.).sum() )

		return opt.fmin(distance, [phaseGuess])[0]




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
	state = orb(phase)
	kickedstate = state.copy()
	kickedstate[1] = kickedstate[1]+0.1

	bestphase = orb.closestPhase(kickedstate, phase-0.3)

	phi = np.arange(0., 2.*np.pi+0.01, 0.01)
	X = orb(phi)
	plot(X[1], X[0], 'k-')
	plot([state[1]], [state[0]], 'ko')
	plot([kickedstate[1]], [kickedstate[0]], 'ko')
	x = orb(bestphase)
	plot([x[1]], [x[0]], 'rs')

	show()
