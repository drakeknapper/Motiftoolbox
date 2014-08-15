#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import window as win
import leech as model
import tools as tl

import numpy as np
import pylab as pl


def set_params(**kwargs):

	if len(kwargs) == 0:
		kwargs = default_params

	model.params.update(kwargs)
	#for k in kwargs.keys():
		#for c in ['0', 'b', 'g', 'r']:
		#	model.params[k+'_'+c] = kwargs[k]

	print '# update params', kwargs



class system(win.window):

	title = "Hodgkin-Huxley Dynamics (orbit+nullclines)"
	figsize = (6, 5)

	def __init__(self, info=None, position=None):
		win.window.__init__(self, position)

		self.x_orbit, self.y_orbit, self.z_orbit, self.orbit_period = None, None, None, None
		self.dt = 0.03
		self.stride = 10
		self.THRESHOLD = model.THRESHOLD
		self.info = info

		self.ax = self.fig.add_subplot(111, yticks=[-50, -25, 0, 25])
		self.ax.set_xlim(0., 0.6)
		self.V_min, self.V_max = -60.0, 40.0
		self.ax.set_ylim(self.V_min, self.V_max)

		self.ax.set_xlabel(r'K$^{+}$ gating $m_{K2}$ (a.u.)', fontsize=17)
		self.ax.set_ylabel(r'membrane voltage $V$ (mV)', fontsize=17)

		self.li_traj, = self.ax.plot([], [], 'k-', lw=1., label='Bursting Orbit')
		self.li_ncl_x, = self.ax.plot([], [], 'r--', lw=2., label='$\dot{V}=0$')
		self.li_ncl_z, = self.ax.plot([], [], 'g--', lw=2., label='$\dot{m}_{K2}=0$')
		pl.legend(loc=0)
		self.li_threshold = self.ax.axhline(y=1000.*model.THRESHOLD)
		self.tx_state_space = self.ax.text(0.2, -0.01, '', color='r')

		self.refresh_nullclines()
		self.refresh_orbit()

		self.fig.canvas.mpl_connect('button_press_event', self.on_button)
		self.fig.canvas.mpl_connect('button_release_event', self.off_button)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focus_in)


	def focus_in(self, event=None):
		descriptor = "System Parameters :\n V_shift = %lf \n tau_Na = %lf \n I_ext = %lf \n g_L = %lf" % (model.params['V_shift'], model.params['tau_Na'], model.params['I_ext'], model.params['g_L'])

		if self.info == None:
			print descriptor

		else:
			self.info.set(descriptor)


	def load_initial_condition(self, x, y): # only one phase:  everything's square!
		X = np.zeros((model.N_EQ3), float)
		phi_x, phi_y = tl.PI2*(1.-x), tl.PI2*(1.-y)
		X[::model.N_EQ1] = self.x_orbit([0., phi_x, phi_y])
		X[1::model.N_EQ1] = self.y_orbit([0., phi_x, phi_y])
		X[2::model.N_EQ1] = self.z_orbit([0., phi_x, phi_y])
		return X


	def load_initial_conditions(self, initial_phase): # only one phase:  everything's square!
		initial_phase = np.asarray(initial_phase)

		n = initial_phase.size
		X = np.zeros((n**2, model.N_EQ3), float)
		phi = tl.PI2*(1.-initial_phase)
		X[:, 0], X[:, 1], X[:, 2] = self.x_orbit(0.), self.y_orbit(0.), self.z_orbit(0.)
	
		for i in xrange(n):
			x_i, y_i, z_i = self.x_orbit(phi[i]), self.y_orbit(phi[i]), self.z_orbit(phi[i])
	
			for j in xrange(n):
				X[i*n+j, model.N_EQ1:] = np.array([x_i, y_i, z_i, self.x_orbit(phi[j]), self.y_orbit(phi[j]), self.z_orbit(phi[j])])
	
		return X


	def refresh_nullclines(self):
		V_i = np.arange(self.V_min, self.V_max, 0.1)
		self.li_ncl_x.set_data(model.nullcline_V(V_i/1000.), V_i)
		self.li_ncl_z.set_data(model.nullcline_m(V_i/1000.), V_i)
		self.li_threshold.set_ydata(1000.*self.THRESHOLD)
		self.fig.canvas.draw()


	def refresh_orbit(self):
	
		new_x, new_y, new_z, new_period = model.single_orbit(N_ORBIT=10**5, V_threshold=self.THRESHOLD)
		try:
			new_x, new_y, new_z, new_period = model.single_orbit(N_ORBIT=10**5, V_threshold=self.THRESHOLD)
			self.x_orbit, self.y_orbit, self.z_orbit, self.orbit_period = new_x, new_y, new_z, new_period
			self.tx_state_space.set_text("")
	
		except:
			self.tx_state_space.set_text("No closed orbit found!")
			pass
		
		phi = np.arange(2000.)/float(1999.)
		self.li_traj.set_data(self.z_orbit(tl.PI2*phi), 1000.*self.x_orbit(tl.PI2*phi))
		self.fig.canvas.draw()


	def on_button(self, event):
		self.event_start = np.array([event.xdata, event.ydata/1000.])


	def N_output(self, CYCLES):
		return int(CYCLES*self.orbit_period/self.dt)


	def off_button(self, event):
		delta_params = np.array([event.xdata, event.ydata/1000.])-self.event_start

		if event.button == 1:
			set_params(I_ext=model.params['I_ext']-delta_params[0]/10., V_shift=model.params['V_shift']-delta_params[1])

		elif event.button == 2:
			self.THRESHOLD = event.ydata/1000.

		elif event.button == 3:
			new_eps = model.params['tau_Na']+delta_params[0]
			set_params(tau_Na=model.params['tau_Na']*(new_eps<0.)+(new_eps>0.)*new_eps, g_L=model.params['g_L']+delta_params[1]*3.)
			
		self.refresh_nullclines()
		self.refresh_orbit()
		self.focus_in()



if __name__ == "__main__":

	sys = system()
	pl.show()






