#! /usr/bin/python

import numpy as np
import tools as tl
import phaseResettingCurve as prc
import attractor as att

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




class fixed_point_2d(object):

	def __init__(self, x, eigenvalues):
		self.fixedPoint = asarray(x)
		self.eigenValues = asarray(eigenvalues)
		self.color = tl.clmap(self.fixedPoint[1], self.fixedPoint[0])
		if np.prod(self.eigenValues) < 0.:
			self.symbol = 'd'
		elif all(self.eigenValues > 0.):
			self.symbol = 's'
		else:
			self.symbol = 'o'


	def distance2fixedPoint(self, other):
		return att.phase_distance(self.fixedPoint, other)


	def plot(self, axis, *args, **kwargs):

		if "period" in kwargs:
			period = kwargs.pop("period")
			self.fixedPoint = period/tl.PI2*self.fixedPoint
		else:
			period = tl.PI2

		kwargs['mfc'] = self.color
		kwargs['ms'] = 15.
		args = (self.symbol)
		#axis.plot([self.fixedPoint[0]],
		  	#[self.fixedPoint[1]],
				#*args, **kwargs)
		axis.plot([self.fixedPoint[0]-period, self.fixedPoint[0], self.fixedPoint[0]+period, self.fixedPoint[0]-period, self.fixedPoint[0], self.fixedPoint[0]+period, self.fixedPoint[0]-period, self.fixedPoint[0], self.fixedPoint[0]+period],
		  	[self.fixedPoint[1]+period, self.fixedPoint[1]+period, self.fixedPoint[1]+period, self.fixedPoint[1], self.fixedPoint[1], self.fixedPoint[1], self.fixedPoint[1]-period, self.fixedPoint[1]-period, self.fixedPoint[1]-period],
				*args, **kwargs)




class interp_torus_vec(object):

	pi2 = np.pi*2

	def __init__(self, X, Y, F, **kwargs):
		self.dimensions = len(F)
		self.fixedPoints = []
		self.f = [scipy.interpolate.RectBivariateSpline(X, Y, F[i], **kwargs)
				for i in xrange(self.dimensions)]

	def __call__(self, XY):
		X, Y = np.mod(XY[0], self.pi2), np.mod(XY[1], self.pi2)	# X : dphi12, Y : dphi13
		return [self.f[i](X, Y)[0, 0] for i in xrange(self.dimensions)]
	


	def findRoots(self, GRID=10):

		phase = self.pi2*np.arange(GRID)/float(GRID)
		phase += phase[1]/2.				# shift into the open field (0, 2pi)x(0, 2pi)

		for i in xrange(GRID):

			for j in xrange(GRID):

				try:
					sol = opt.root(self.__call__, x0=[phase[i], phase[j]], tol=10.**-7, method='broyden2')
			
					if sol.success:
						jac = np.array([scipy.optimize.approx_fprime(sol.x, lambda x: self.f[i](x[0], x[1]), epsilon=0.05)
									for i in xrange(self.dimensions)])
						eigenvalues =  scipy.linalg.eigvals(jac)
						self.fixedPoints.append( fixed_point_2d(x=sol.x, eigenvalues=np.real(eigenvalues)) )
				except:
					print "findRoots :\tphi =", [phase[i], phase[j]], "failed."



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
	
	phase, coupling = net.threeCellCoupling(0.03)
	coupling_function = interp_torus_vec(phase, phase, coupling) # dphi12, dphi13, coupling=[q12, q13]
	coupling_function.findRoots(GRID=10)
	coupling_function.plot()

	




