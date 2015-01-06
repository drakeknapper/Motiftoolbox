#! /usr/bin/python

import numpy as np
import tools as tl
import scipy.integrate as itg
import phaseResettingCurve as prc


class prcNetwork(prc.phaseResettingCurve):

	def __init__(self, model):
		prc.phaseResettingCurve.__init__(self, model)


	def set_model(self, model):
		prc.phaseResettingCurve.set_model(self, model)
		self.coupling = model.coupling_function	# bivariate coupling function
		self.coupling_idx = model.IDX_COUPLING	# variable idx to which coupling is added



	def twoCellCoupling(self):

		if not self.PRC_COMPUTED: self.compute_prc()

		phase =  np.arange(0., 2.*np.pi, 0.01)
		Q = self.prcurve[self.coupling_idx](phase)			# Q[i] = Q(phi_i)
		K = [self.coupling(phi, phase) for phi in phase]	# K[i, j] = K(phi_i, phi_j)
		K = np.array(K)

		integral = []
		idx = np.arange(phase.size)
		for i in idx:
			idx_q = np.mod(idx+i, phase.size)
			idx_1 = [[idx[k], idx_q[k]] for k in xrange(idx.size)]
			idx_2 = [[idx_q[k], idx_q[k]] for k in xrange(idx.size)]
			K_1 = [K[idx_i] for idx_i in idx_1]
			K_2 = [K[idx_i] for idx_i in idx_2]
			integrand = Q*K_1-Q[idx_q]*K_2
			
			integral.append(itg.trapz(integrand, phase))

		return phase, np.array(integral)



	def show(self):
		from pylab import figure, subplot, plot, show, tight_layout

		phase = np.arange(0., 2.*np.pi+0.01, 0.01)

		fig = figure()
		ax = fig.add_subplot(1, 1, 1)
		tight_layout()
		show()




if __name__ == '__main__':

	from pylab import *
	import fitzhugh as model

	net = prcNetwork(model)
	phase, Delta = net.twoCellCoupling()

	plot(phase, Delta, 'k-', lw=2)
	show()
