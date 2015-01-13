#! /usr/bin/python

import numpy as np
import tools as tl
import phaseResettingCurve as prc
import distribute

import scipy.interpolate
import scipy.optimize as opt
import scipy.linalg

import ctypes as ct
import os
lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_prcNetwork.so')


lib.trapz_twoCell.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_int, ct.POINTER(ct.c_double), ct.c_double, ct.POINTER(ct.c_double)]
def trapz_twoCell(Q, K, dphi, strength):
	# strenght[0] : 2 -> 1
	# strenght[1] : 1 -> 2
	Q, K, strength = np.array(Q),  np.array(K), np.array(strength)
	result = np.zeros((Q.size), float)

	lib.trapz_twoCell(Q.ctypes.data_as(ct.POINTER(ct.c_double)), K.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_int(Q.size),
			strength.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_double(dphi), result.ctypes.data_as(ct.POINTER(ct.c_double)))
	
	return result


lib.trapz_threeCell.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_int, ct.POINTER(ct.c_double), ct.c_double, ct.POINTER(ct.c_double)]
def trapz_threeCell(Q, K, dphi, strength):
	Q, K, strength = np.array(Q), np.array(K), np.array(strength)
	result = np.zeros((2, Q.size, Q.size), float)

	lib.trapz_threeCell(Q.ctypes.data_as(ct.POINTER(ct.c_double)), K.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_int(Q.size),
			strength.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_double(dphi), result.ctypes.data_as(ct.POINTER(ct.c_double)))

	return result



def phase_difference(x, y, period):
	period2 = 0.5*period
	return np.mod(x-y+period2, period)-period2


def phase_distance(x, y, period):
	return np.sqrt(np.sum(phase_difference(x, y, period)**2))


class fixedPoint_torus(object):

	stability = [['sink', 'o'], ['source', 'x'], ['saddle', 'd']]

	def __init__(self, x, eigenvalues):
		self.fp = np.mod(np.asarray(x), tl.PI2)
		self.eigenValues = np.asarray(eigenvalues)
		self.color = tl.clmap(self.fp[1], self.fp[0])

		if np.prod(self.eigenValues) < 0.:
			self.stability_idx = 2 # saddle

		elif all(self.eigenValues > 0.):
			self.stability_idx = 1 # source

		elif all(self.eigenValues < 0.):
			self.stability_idx = 0 # sink

		else:
			self.stability_idx = 'undefined'


	def distance2(self, other):
		return phase_distance(self.fp, other.fp, period=tl.PI2)


	def plot(self, axis, *args, **kwargs):

		if self.stability_idx == 'undefined': return

		if 'period' in kwargs:
			period = kwargs.pop("period")
			self.fp = period/tl.PI2*self.fp
		else:
			period = tl.PI2

		if not 'ms' in kwargs:	kwargs['ms'] = 15.
		if not 'mfc' in kwargs:	kwargs['mfc'] = self.color

		args = (self.stability[self.stability_idx][1])

		axis.plot([self.fp[0]-period, self.fp[0], self.fp[0]+period, self.fp[0]-period, self.fp[0], self.fp[0]+period, self.fp[0]-period, self.fp[0], self.fp[0]+period],
		  	[self.fp[1]+period, self.fp[1]+period, self.fp[1]+period, self.fp[1], self.fp[1], self.fp[1], self.fp[1]-period, self.fp[1]-period, self.fp[1]-period],
				*args, **kwargs)




class interp_torus_vec(object):


	def __init__(self, X, Y, F, **kwargs):
		self.dimensions = len(F)
		self.fixedPoints = []
		self.f = [scipy.interpolate.RectBivariateSpline(X, Y, F[i], **kwargs)
				for i in xrange(self.dimensions)]

	def __call__(self, XY):
		X, Y = np.mod(XY[0], tl.PI2), np.mod(XY[1], tl.PI2)	# X : dphi12, Y : dphi13
		return [self.f[i](X, Y)[0, 0] for i in xrange(self.dimensions)]
	


	def findRoots(self, GRID=10):

		phase = tl.PI2*np.arange(GRID)/float(GRID)
		dphi = phase[1]
		phase += dphi/2.				# shift into the open field (0, 2pi)x(0, 2pi)

		### find roots from a grid ###

		def findRoot(n):

			i, j = n/GRID, np.mod(n, GRID)

			try:
				fp = [[0., 0.], [0., 0.]]
				sol = opt.root(self.__call__, x0=[phase[i], phase[j]], tol=10.**-11, method='broyden2')

				if sol.success:
					jac = np.array([scipy.optimize.approx_fprime(sol.x, lambda x: self.f[i](x[0], x[1]), epsilon=0.25*dphi)
								for i in xrange(self.dimensions)])
					eigenvalues =  scipy.linalg.eigvals(jac)

					if not any([eigen == 0. for eigen in eigenvalues]):
						fp = [sol.x, np.real(eigenvalues)]


			except:
				pass

			return fp

		fps = distribute.distribute(findRoot, 'n', np.arange(GRID**2))

		fixedPoints = [fixedPoint_torus(x=fp[0], eigenvalues=fp[1]) for fp in fps]


		### discard doubles ###

		self.fixedPoints = [fixedPoints[0]]
		for newfp in fixedPoints:
			conditions = [oldfp.distance2(newfp) > 0.1*dphi for oldfp in self.fixedPoints]

			if all(conditions):
				self.fixedPoints.append(newfp)



	def plot(self, GRID=10, ax=None, **kwargs):

		NEWAXIS = False
		if ax == None:
			from pylab import figure, show
			fig = figure()
			ax = fig.add_subplot(111)
			NEWAXIS = True

		if "period" in kwargs:
			period = kwargs.pop("period")
		else:
			period = tl.PI2


		phase = tl.PI2*np.arange(GRID)/float(GRID-1)
		phase += phase[1]/2.
		UV = np.asarray([ [self([phase[i], phase[j]]) for i in xrange(GRID)]	# Spalten
								for j in xrange(GRID)])	# Zeilen

		X, Y = np.meshgrid(phase, phase)
		U, V = UV[:, :, 0], UV[:, :, 1]
	
		Q = ax.quiver(period*X/tl.PI2, period*Y/tl.PI2, U, V, units='width')

		for fp in self.fixedPoints:
			fp.plot(axis=ax, period=period)
			
		if NEWAXIS:
			ax.set_xlim(0., period)
			ax.set_ylim(0., period)
			fig.tight_layout()
			show()

		return Q




class prcNetwork(prc.phaseResettingCurve):

	def __init__(self, model):
		prc.phaseResettingCurve.__init__(self, model)


	def set_model(self, model):
		prc.phaseResettingCurve.set_model(self, model)
		self.coupling = model.coupling_function	# bivariate coupling function
		self.coupling_idx = model.IDX_COUPLING	# variable idx to which coupling is added


	def getCouplingFunctions(self, dphi):

		if not self.PRC_COMPUTED: self.compute_prc()
		phase =  np.arange(0., 2.*np.pi+dphi, dphi)	# [ 0, 2pi ]
		traj_n = self.trajectory[self.coupling_idx](phase)
		Q = self.prcurve[self.coupling_idx](phase)			# Q[i] = Q(phi_i)
		K = np.array([self.coupling(state, traj_n) for state in traj_n])	# K[i, j] = K(phi_i, phi_j)

		return phase, Q, K



	def twoCellCoupling(self, dphi, strength = np.ones((2), float)):
		phase, Q, K = self.getCouplingFunctions(dphi)

		integral = trapz_twoCell(Q, K, dphi, strength)

		return phase, integral



	def threeCellCoupling(self, dphi, strength=np.ones((6), float)):
		phase, Q, K = self.getCouplingFunctions(dphi)	# phase range, Q[j]=Q(phi_j), K[i, j]=K(phi_i, phi_j)

		integral = trapz_threeCell(Q, K, dphi, strength)

		return phase, integral



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
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	import scipy.interpolate as interp

	model.setParams(I=0.42, epsilon=0.3);
	net = prcNetwork(model)

	#figure()
	#phase, coupling = net.twoCellCoupling(0.05, strength=[1., ratio])
	#plot(phase, coupling, 'ko')
	#tight_layout()
	
	phase, coupling = net.threeCellCoupling(0.01)
	coupling_function = interp_torus_vec(phase, phase, coupling) # dphi12, dphi13, coupling=[q12, q13]
	coupling_function.findRoots(GRID=15)
	coupling_function.plot()

	




