#! /usr/bin/python

import numpy as np
import tools as tl


class orbit(object):

	def __init__(self, model):
		self.set_model(model)


	def __call__(self, phase):
		return np.array([self.trajectory[j](phase) for j in xrange(self.dimensions)])


	def set_model(self, model):
		self.model = model
		self.integrate = model.integrate_one_rk4
		self.dimensions = model.N_EQ1
		self.trajectory = [tl.splineLS1D(isphase=False) for i in xrange(self.dimensions)]

		self.initial_state = model.initial_state
		self.dt = model.dt
		self.stride = model.stride
		self.N_integrate = model.N_integrate


	def find_orbit(self):
		X = self.integrate(self.initial_state, self.dt/float(self.stride), self.N_integrate, self.stride)
		assert X.shape[0] == self.dimensions
	
		try:
			ni = tl.crossings(X[self.model.IDX_THRESHOLD], self.model.THRESHOLD) # convert to millivolts
			X = X[:, ni[-2]:ni[-1]]
			phase = tl.PI2*np.arange(X.shape[1])/float(X.shape[1]-1)
	
			for i in xrange(self.dimensions):
				self.trajectory[i].makeModel(X[i], phase)
		
		except:
			print '# single_orbit:  No closed orbit found!'
			raise ValueError
	
		self.period = self.dt*X.shape[1]	 # in msec.


	def show(self):
		from pylab import plot, show
		phase = np.arange(0., 2.*np.pi, 0.01)
		X = self(phase)
		plot(phase, X[0])
		plot(phase, X[1])
		show()




if __name__ == '__main__':

	import plant as model
	from pylab import *

	orb = orbit(model)
	orb.find_orbit()
	phase = 2.*pi*arange(0., 1., 0.001)
	X = orb(phase)
	plot(X[0])
	plot(X[1])
	show()

