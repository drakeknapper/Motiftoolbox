#! /usr/bin/python

import numpy as np
import tools as tl
import orbit as orb


class phaseResettingCurve(orb.orbit):

	N_PRC = 256

	def __init__(self, model):
		orb.orbit.__init__(self, model)


	def set_model(self, model):
		orb.orbit.set_model(self, model)
		self.PRC_COMPUTED = False
		self.prcurve = [tl.splineLS1D(isphase=False) for i in xrange(self.dimensions)]



	def compute_prc(self, kick=0.005, N_integrate=10**4):

		if self.PRC_COMPUTED:
			return

		phase, shiftedStates = self.shifted_orbit(kick=kick)	# copies of shifted orbits separate for each dimension
		shiftedStates = np.concatenate(shiftedStates)	# make one long set of initial values

		dt = self.period/float(N_integrate)
	
		# integrate the shifted trajectories!
		if self.model.CUDA_ENABLED:	iterStates = self.model.cuda_integrate_one_rk4_nosave(shiftedStates, dt, N_integrate)
		else:				iterStates = np.array([self.model.integrate_one_rk4_nosave(shiftedStates[i], dt, N_integrate) for i in xrange(O.shape[1])])

		iterStates = [iterStates[i*self.N_PRC:(i+1)*self.N_PRC] for i in xrange(self.dimensions)]

		kickedPhases = []
		for dim in xrange(self.dimensions): # compute phases from iterated kicked states
			iSdim = iterStates[dim]
			kickedPhase = np.array([self.closestPhase(iSdim[i], phase[i]) for i in xrange(iSdim.shape[0])])
			self.prcurve[dim].makeModel((kickedPhase-phase)/kick, phase)

		self.PRC_COMPUTED = True



	def evaluate_prc(self, phase):
		return np.array([self.prcurve[j](phase) for j in xrange(self.dimensions)])


	def show(self):
		from pylab import figure, subplot, plot, show, tight_layout

		phase = np.arange(0., 2.*np.pi+0.01, 0.01)

		fig = figure()
		X = self.evaluate_orbit(phase)
		PRC = self.evaluate_prc(phase)
		ax = fig.add_subplot(self.dimensions, 1, 1)

		for i in xrange(self.dimensions):
			if i > 0: ax = fig.add_subplot(self.dimensions, 1, i+1, sharex=ax)
			ax.plot(phase, X[i], 'k-', lw=2.)
			ax2 = ax.twinx()
			ax2.plot(phase, PRC[i], 'r-', lw=2.)
			ax2.axhline(y=0., ls='--')

			ax.set_xlim(0., 2.*np.pi)

		tight_layout()
		show()




if __name__ == '__main__':

	import fitzhugh as model
	from pylab import *

	prc = phaseResettingCurve(model)
	prc.compute_prc()
	prc.show()
