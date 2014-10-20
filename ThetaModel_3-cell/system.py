#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import thetax2 as th2
import tools as tl
import window as win

import numpy as np
import pylab as pl


def set_params(**kwargs):

	if len(kwargs) == 0:
		kwargs = default_params

	for k in kwargs.keys():

		for c in ['0', 'b', 'g', 'r']:
			th2.params[k+'_'+c] = kwargs[k]

	print '# update params', kwargs



class system(win.window):

	title = "System"
	figsize = (6, 5)

	def __init__(self, info=None, position=None, network=None, traces=None):
		win.window.__init__(self, position)

		self.x_orbit, self.y_orbit, self.orbit_period = None, None, None
		self.dt = 0.1
		self.stride = 30
		self.info = info
		self.network = network
		self.traces = traces

		self.ax = self.fig.add_subplot(111, yticks=[-1.5, -0.5, 0.5, 1.5])
		self.ax.set_xlim(0., 1.)
		self.x_min, self.x_max = -1.5, 1.5
		self.ax.set_ylim(self.x_min, self.x_max)
		self.ax.set_xlabel(r'phase $\theta$', fontsize=15)
		self.ax.set_ylabel(r'inst. frequency', fontsize=15)
		self.ax.set_xticks([0., np.pi/2., np.pi, 3*np.pi/2., 2.*np.pi])
		self.ax.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2.$', r'$2\pi$'], fontsize=15)
		self.ax.axvline(x=2.*np.pi, color='k')

		self.ax.axhline(y=0, color='k')
		self.li_rhs, = self.ax.plot([], [], 'r-', lw=2.)
		self.tx_state_space = self.ax.text(0.2, -0.5, '', color='r')
		self.ax.set_xlim(0., 3.*np.pi)
		self.ax.set_ylim(-0.2, 3.)

		self.refresh_rhs()
		self.refresh_orbit()

		if not position == None:
			try:
				self.fig.canvas.manager.window.wm_geometry(position)
			except:
				pass

		self.fig.canvas.mpl_connect('button_press_event', self.on_button)
		self.fig.canvas.mpl_connect('button_release_event', self.off_button)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focus_in)


	def focus_in(self, event=None):
		descriptor = "System Parameters :\n omega_0 = %lf \n beta_0 = %lf" % (th2.params['omega_0'], th2.params['beta_0'])

		if self.info == None:
			print descriptor

		else:
			self.info.set(descriptor)


	def load_initial_condition(self, x, y): # only one phase:  everything's square!
		phi_x, phi_y = tl.PI2*(1.-x), tl.PI2*(1.-y)
		X = np.arctan2(self.y_orbit([0., phi_x, phi_y]), self.x_orbit([0., phi_x, phi_y]))
		return X


	def load_initial_conditions(self, initial_phase): # only one phase:  everything's square!
		initial_phase = np.asarray(initial_phase)

		n = initial_phase.size
		X = np.zeros((n**2, th2.N_EQ3), float)
		phi = tl.PI2*(1.-initial_phase)
		X[:, 0] = np.arctan2(self.y_orbit(0.), self.x_orbit(0.))
	
		for i in xrange(n):
			theta_i = np.arctan2(self.y_orbit(phi[i]), self.x_orbit(phi[i]))
	
			for j in xrange(n):
				theta_j = np.arctan2(self.y_orbit(phi[j]), self.x_orbit(phi[j]))
				X[i*n+j, th2.N_EQ1:] = np.array([theta_i, theta_j])
	
		return X


	def refresh_rhs(self):
		phase = np.arange(0., 3.*np.pi+0.01, 0.01)
		self.li_rhs.set_data(phase, th2.right_hand_side(phase))
		self.fig.canvas.draw()


	def refresh_orbit(self):
	
		try:
			new_x, new_y, new_period = th2.single_orbit(N_ORBIT=10**4)
			self.x_orbit, self.y_orbit, self.orbit_period = new_x, new_y, new_period
			self.tx_state_space.set_text("")
	
		except:
			self.tx_state_space.set_text("No closed orbit found!")
			pass
		
		self.fig.canvas.draw()

		try:
			self.traces.computeTraces()
		except:
			pass


	def on_button(self, event):
		self.event_start = np.array([event.xdata, event.ydata])


	def N_output(self, CYCLES):
		return int(CYCLES*self.orbit_period/self.dt)


	def off_button(self, event):
		delta_params = np.array([event.xdata, event.ydata])-self.event_start

		if event.button == 1:
			set_params(omega=th2.params['omega_0']+delta_params[1], beta=th2.params['beta_0']+delta_params[0])

		self.refresh_rhs()
		self.refresh_orbit()
		self.focus_in()



if __name__ == "__main__":

	sys = system()
	pl.show()






