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


lib.trapz_twoCell.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_int, ct.c_double, ct.POINTER(ct.c_double)]
def trapz_twoCell(Q, K, dphi):
	Q, K = np.array(Q),  np.array(K)
	result = np.zeros((Q.size), float)

	lib.trapz_twoCell(Q.ctypes.data_as(ct.POINTER(ct.c_double)), K.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_int(Q.size),
			ct.c_double(dphi), result.ctypes.data_as(ct.POINTER(ct.c_double)))
	
	return result


lib.trapz_threeCell.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_int, ct.c_double, ct.POINTER(ct.c_double)]
def trapz_threeCell(Q, K, dphi):
	Q, K = np.array(Q), np.array(K)
	result = np.zeros((2, Q.size, Q.size), float)

	lib.trapz_threeCell(Q.ctypes.data_as(ct.POINTER(ct.c_double)), K.ctypes.data_as(ct.POINTER(ct.c_double)),
				ct.c_int(Q.size), ct.c_double(dphi), result.ctypes.data_as(ct.POINTER(ct.c_double)))

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
		kwargs['mfc'] = self.color
		kwargs['ms'] = 15.
		args = (self.symbol)
		#axis.plot([self.fixedPoint[0]],
		  	#[self.fixedPoint[1]],
				#*args, **kwargs)
		axis.plot([self.fixedPoint[0]-tl.PI2, self.fixedPoint[0], self.fixedPoint[0]+tl.PI2, self.fixedPoint[0]-tl.PI2, self.fixedPoint[0], self.fixedPoint[0]+tl.PI2, self.fixedPoint[0]-tl.PI2, self.fixedPoint[0], self.fixedPoint[0]+tl.PI2],
		  	[self.fixedPoint[1]+tl.PI2, self.fixedPoint[1]+tl.PI2, self.fixedPoint[1]+tl.PI2, self.fixedPoint[1], self.fixedPoint[1], self.fixedPoint[1], self.fixedPoint[1]-tl.PI2, self.fixedPoint[1]-tl.PI2, self.fixedPoint[1]-tl.PI2],
				*args, **kwargs)




class interp_torus_vec(object):

	pi2 = np.pi*2

	def __init__(self, X, Y, F, **kwargs):
		self.dimensions = len(F)
		self.fixedPoints = []
		self.f = [scipy.interpolate.RectBivariateSpline(X, Y, F[i], **kwargs)
				for i in xrange(self.dimensions)]

	def __call__(self, XY):
		X, Y = np.mod(XY[0], self.pi2), np.mod(XY[1], self.pi2)
		return [self.f[i](X, Y)[0, 0] for i in xrange(self.dimensions)]
	


	def findRoots(self, GRID=10):

		phase = self.pi2*np.arange(GRID)/float(GRID)

		for i in xrange(GRID):

			for j in xrange(GRID):

				try:
					sol = opt.root(self.__call__, x0=[phase[i], phase[j]], tol=10**-13, method='broyden2')
			
					if sol.success:
						jac = np.array([scipy.optimize.approx_fprime(sol.x, lambda x: self.f[i](x[0], x[1]), epsilon=0.05)
									for i in xrange(self.dimensions)])
						eigenvalues =  scipy.linalg.eigvals(jac)
						self.fixedPoints.append( fixed_point_2d(x=sol.x, eigenvalues=np.real(eigenvalues)) )
				except:
					print "phi =", [phase[i], phase[j]], "failed"



	def plot(self, GRID=10):
		phase = self.pi2*np.arange(GRID)/float(GRID-1)
		UV = np.asarray([[self([phase[j], phase[i]])
					for i in xrange(GRID)]
					for j in xrange(GRID)])

		X, Y = meshgrid(phase, phase)
		U = UV[:, :, 1]
		V = UV[:, :, 0]
	
		figure()
		ax = subplot(111)
		Q = quiver(X, Y, U, V, units='width')

		for fp in self.fixedPoints:
			fp.plot(axis=ax)
			
		xlim(0, self.pi2)
		ylim(0, self.pi2)
		tight_layout()
		show()





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



	def twoCellCoupling(self, dphi):
		phase, Q, K = self.getCouplingFunctions(dphi)

		integral = trapz_twoCell(Q, K, dphi)

		return phase, integral



	def threeCellCoupling(self, dphi):
		phase, Q, K = self.getCouplingFunctions(dphi)

		integral = trapz_threeCell(Q, K, dphi)

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

	model.setParams(I=0.6, epsilon=0.1);
	net = prcNetwork(model)

	figure()
	phase, coupling = net.twoCellCoupling(0.05)
	plot(phase, coupling, 'ko')
	tight_layout()
	
	phase, coupling = net.threeCellCoupling(0.05)
	#phase, coupling = phase[::50], coupling[:, ::50, ::50]
	coupling_function = interp_torus_vec(phase, phase, coupling)
	coupling_function.findRoots(GRID=20)
	coupling_function.plot()

	




