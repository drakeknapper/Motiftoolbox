#! /usr/bin/python

import window as win


class system(win.window, object):

	title = "System"
	figsize = (5, 5)
	params = {}


	def __init__(self, position=None):
		win.window.__init__(self, position)



	def inputParams(self):

		new_params = dict()

		print "\n=======SET PARAMS======="
		print "Enter new parameter value. (If you don't want to change a parameter, simply hit enter.)"
		for p in type(self).params.keys():
			[p_name, n_osci] = p.split('_')

			if n_osci == '0':

				try:
					p_new = float(raw_input("%s = " % p_name))
					new_params[p_name] = p_new

				except:
					pass
		print "=======SET PARAMS=======\n"
		
		return new_params



	def setParams(self, **kwargs):

		if len(kwargs) == 0:
			kwargs = self.inputParams()

		# setParams(**kwargs)
		# self.refresh_nullclines()
		# self.refresh_orbit()
	

	def load_initial_condition(self, x, y): # only one phase:  everything's square!
		"""
		Outdated
		"""
		X = np.zeros((model.N_EQ3), float)
		phi_x, phi_y = tl.PI2*(1.-x), tl.PI2*(1.-y)
		X[::model.N_EQ1] = self.x_orbit([0., phi_x, phi_y])
		X[1::model.N_EQ1] = self.y_orbit([0., phi_x, phi_y])
		return X



	def loadInitialCondition(self, dphi): # only one phase:  everything's square!

		dphi = np.asarray(dphi)
		NODES = dphi.size

		phi = np.concatenate(([0.], tl.PI2*(1.-dphi)))

		X = np.zeros((model.N_EQ1*NODES), float)

		for i in xrange(NODES):
			X[i::model.N_EQ1] = self.orbit[i](phi)

		return X



	def load_initial_conditions(self, initial_phase): # only one phase:  everything's square!
		initial_phase = np.asarray(initial_phase)

		n = initial_phase.size
		X = np.zeros((n**2, model.N_EQ3), float)
		phi = tl.PI2*(1.-initial_phase)
		X[:, 0], X[:, 1] = self.x_orbit(0.), self.y_orbit(0.)
	
		for i in xrange(n):
			V_i, H_i = self.x_orbit(phi[i]), self.y_orbit(phi[i])
	
			for j in xrange(n):
				X[i*n+j, model.N_EQ1:] = np.array([V_i, H_i, self.x_orbit(phi[j]), self.y_orbit(phi[j])])
	
		return X






if __name__ == "__main__":

	import pylab as pl

	mysys = system()
	pl.show()


