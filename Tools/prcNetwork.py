#! /usr/bin/python

import numpy as np
import tools as tl
import scipy.integrate as itg
import phaseResettingCurve as prc

import ctypes as ct
import os
lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_prcNetwork.so')


lib.trapz_twoCell.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_uint, ct.c_double, ct.POINTER(ct.c_double)]
def trapz_twoCell(Q, K, dphi):
	Q, K = np.array(Q),  np.array(K)
	result = np.zeros((Q.size), float)

	lib.trapz_twoCell(Q.ctypes.data_as(ct.POINTER(ct.c_double)), K.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_uint(Q.size),
			ct.c_double(dphi), result.ctypes.data_as(ct.POINTER(ct.c_double)))
	
	return result


lib.trapz_threeCell.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_uint, ct.c_double, ct.POINTER(ct.c_double)]
def trapz_threeCell(Q, K, dphi):
	Q, K = np.array(Q), np.array(K)
	result = np.zeros((2, Q.size, Q.size), float)

	lib.trapz_threeCell(Q.ctypes.data_as(ct.POINTER(ct.c_double)), K.ctypes.data_as(ct.POINTER(ct.c_double)),
				ct.c_uint(Q.size), ct.c_double(dphi), result.ctypes.data_as(ct.POINTER(ct.c_double)))

	return result




class prcNetwork(prc.phaseResettingCurve):

	def __init__(self, model):
		prc.phaseResettingCurve.__init__(self, model)


	def set_model(self, model):
		prc.phaseResettingCurve.set_model(self, model)
		self.coupling = model.coupling_function	# bivariate coupling function
		self.coupling_idx = model.IDX_COUPLING	# variable idx to which coupling is added


	def getPhaseCoupling(self, dphi):

		if not self.PRC_COMPUTED: self.compute_prc()
		phase =  np.arange(0., 2.*np.pi, dphi)	# [0, 2pi], including both ends
		traj_n = self.trajectory[self.coupling_idx](phase)
		Q = self.prcurve[self.coupling_idx](phase)			# Q[i] = Q(phi_i)
		K = np.array([self.coupling(state, traj_n) for state in traj_n])	# K[i, j] = K(phi_i, phi_j)

		return phase, Q, K



	def twoCellCoupling(self, dphi=0.01):
		phase, Q, K = self.getPhaseCoupling(dphi)

		integral = trapz_twoCell(Q, K, dphi)

		return phase, integral



	def threeCellCoupling(self, dphi=0.01):
		phase, Q, K = self.getPhaseCoupling(dphi)

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

	model.setParams(I=0.6);
	net = prcNetwork(model)

	figure()
	phase, coupling = net.twoCellCoupling(0.01)
	plot(phase/(2.*pi), coupling)
	tight_layout()
	
	phase, coupling = net.threeCellCoupling(0.05)
	phase, coupling = phase[::10]/(2.*pi), coupling[:, ::10, ::10]
	X, Y = meshgrid(phase, phase)
	U, V = coupling[1], coupling[0]
	figure()
	Q = quiver(X, Y, U, V, units='width')
	tight_layout()
	show()
	
