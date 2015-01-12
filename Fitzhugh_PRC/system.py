#!/usr/bin/env python

import sys
sys.path.insert(0, '../Tools')
import window as win
import prcNetwork
import fitzhugh as model
import tools as tl
import numpy as np
import pylab as pl



class system(win.window, prcNetwork.prcNetwork):

	title = "System"
	figsize = (6, 12)
	debug = False

	def __init__(self, info=None, position=None, network=None, traces=None, torus=None):
		prcNetwork.prcNetwork.__init__(self, model)
		win.window.__init__(self, position)

		self.dt = 0.1
		self.stride = 30
		self.info = info
		self.network = network
		self.traces = traces
		self.torus = torus

		### state space

		self.ax = self.fig.add_subplot(311, yticks=[-1.5, -0.5, 0.5, 1.5], xticks=[0, 0.5, 1.0])
		self.ax.set_xlim(-0.05, 1.05)
		self.x_min, self.x_max = -1.5, 1.5
		self.ax.set_ylim(self.x_min, self.x_max)
		self.ax.set_ylabel(r'Membrane Voltage $V$', fontsize=18)
		self.ax.set_xlabel(r'Inactivation $x$', fontsize=18)
                self.ax.set_title('Single-Cell Dynamics (+coupling)', fontsize=18)

		self.li_traj, = self.ax.plot([], [], 'k-', lw=1., label='Oscillation Cycle')
		self.li_ncl_y, = self.ax.plot([], [], 'b--', lw=2., label='nullcline $\dot{V}=0$')
		self.li_ncl_x, = self.ax.plot([], [], 'g--', lw=2., label='nullcline $\dot{x}=0$')
		self.li_sft_ncl_x, = self.ax.plot([], [], 'r-.', lw=2., label='with coupling')
		self.tx_state_space = self.ax.text(0.2, -0.5, '', color='r')

		### phase resetting curves

		self.ax_prc1 = self.fig.add_subplot(312, xticks=[0, 0.5, 1.0])
		self.ax_prc1.set_xlim(0., 1.)
		self.ax_prc1.set_ylim(-5., 5.)
		self.ax_prc1.set_ylabel(r'PRC$_V$', fontsize=18)
		self.ax_prc1.axhline(y=0, ls='--')

		self.li_prc1_traj, = self.ax_prc1.plot([], [], 'k-', lw=1.)
		self.li_prc1, = self.ax_prc1.plot([], [], 'r-', lw=2.)

		self.ax_prc2 = self.fig.add_subplot(313, xticks=[0, 0.5, 1.0], sharex=self.ax_prc1)
		self.ax_prc2.set_xlim(0., 1.)
		self.ax_prc2.set_ylim(-5., 5.)
		self.ax_prc2.set_ylabel(r'PRC$_x$', fontsize=18)
		self.ax_prc2.set_xlabel(r'phase', fontsize=18)
		self.ax_prc2.axhline(y=0, ls='--')

		self.li_prc2_traj, = self.ax_prc2.plot([], [], 'k-', lw=1.)
		self.li_prc2, = self.ax_prc2.plot([], [], 'r-', lw=2.)

                #pl.legend(loc=0, prop=dict(size=18))

		self.key_func_dict.update(dict(C=self.setParams))
		
                self.fig.canvas.mpl_connect('button_press_event', self.on_button)
		self.fig.canvas.mpl_connect('button_release_event', self.off_button)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focus_in)

		self.setParams(epsilon=0.3)
		self.refresh_nullclines()
		self.refresh_orbit()


	def focus_in(self, event=None):

		descriptor = "System Parameters :\n I_0 = %lf \n x_0 = %lf \n epsilon_0 = %lf \n k_0 = %lf \n m_0 = %lf" % (model.params['I_0'], model.params['x_0'], model.params['epsilon_0'], model.params['k_0'], model.params['m_0'])
                descriptor += "\n\n'C': change parameters"
		if self.info == None:
			print descriptor

		else:
			self.info.set(descriptor)



	def inputParams(self):

		new_params = dict()

		print "\n=======SET PARAMS======="
		print "Enter new parameter value. (If you don't want to change a parameter, simply hit enter.)"
		for p in model.params.keys():
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

		self.ORBIT_COMPUTED = False
		self.PRC_COMPUTED = False

		if len(kwargs) == 0:
			kwargs = self.inputParams()


		for k in kwargs.keys():
	
			for c in ['0', 'b', 'g', 'r']:
				model.params[k+'_'+c] = kwargs[k]

		print '# update params', kwargs

		self.refresh_nullclines()
		self.refresh_orbit()

	
	
	def load_initial_condition(self, x, y): # only one phase:  everything's square!
		X = np.zeros((model.N_EQ3), float)
		phi_x, phi_y = tl.PI2*(1.-x), tl.PI2*(1.-y)
		XY = self.evaluate_orbit([0., phi_x, phi_y])
		X[::model.N_EQ1] = XY[0]
		X[1::model.N_EQ1] = XY[1]
		return X


	def load_initial_conditions(self, initial_phase): # only one phase:  everything's square!
		initial_phase = np.asarray(initial_phase)

		n = initial_phase.size
		X = np.zeros((n**2, model.N_EQ3), float)
		phi = tl.PI2*(1.-initial_phase)
		XY = self.evaluate_orbit(0.)
		X[:, 0], X[:, 1] = XY[0], XY[1]
	
		for i in xrange(n):
			XY_i = self.evaluate_orbit(phi[i])
	
			for j in xrange(n):
				XY_j = self.evaluate_orbit(phi[j])
				X[i*n+j, model.N_EQ1:] = np.array([XY_i[0], XY_i[1], XY_j[0], XY_j[1]])
	
		return X


	def refresh_nullclines(self):
		x_i = np.arange(self.x_min, self.x_max, 0.01)

		nullcline_x = model.nullcline_x(x_i, model.params['I_0'], model.params['m_0'])
		nullcline_y = model.nullcline_y(x_i, model.params['x_0'], model.params['k_0'])

		xscale, yscale = 1., 3.
		nullcline_x, X_i = tl.adjustForPlotting(nullcline_x, x_i, ratio=xscale/yscale, threshold=0.03*xscale)
		nullcline_y, x_i = tl.adjustForPlotting(nullcline_y, x_i, ratio=xscale/yscale, threshold=0.03*xscale)

		self.li_ncl_x.set_data(nullcline_x, X_i)
		self.li_ncl_y.set_data(nullcline_y, x_i)

		if not self.network == None:
			g = self.network.coupling_strength

			if all( (g[0] == g_i for g_i in g) ):
				nullcline_xg = model.nullcline_x(X_i, model.params['I_0'], model.params['m_0'], model.params['E_0'], g[0])
				self.li_sft_ncl_x.set_data(nullcline_xg, X_i)

			else:
				self.li_sft_ncl_x.set_data([], [])
			

		self.fig.canvas.draw()


	def refresh_orbit(self):
	
		try:
			self.compute_prc(kick=0.01, N_integrate=10**4) # also computes orbit

			self.tx_state_space.set_text("")
	
		except:
			self.tx_state_space.set_text("No closed orbit found!")
			pass
		
		phi = np.arange(500)/float(499.)

		X = self.evaluate_orbit(tl.PI2*phi)
		xscale, yscale = 1., 3.
		y, x = tl.adjustForPlotting(X[1], X[0], ratio=xscale/yscale, threshold=0.03*xscale)
		y[-1], x[-1] = y[0], x[0]
		self.li_traj.set_data(y, x)

		PRC = self.evaluate_prc(tl.PI2*phi)	
		self.li_prc1.set_data(tl.adjustForPlotting(phi, PRC[0], ratio=xscale/yscale, threshold=0.03*xscale))
		self.li_prc1_traj.set_data(tl.adjustForPlotting(phi, X[0], ratio=xscale/yscale, threshold=0.03*xscale))
		self.li_prc2.set_data(tl.adjustForPlotting(phi, PRC[1], ratio=xscale/yscale, threshold=0.03*xscale))
		self.li_prc2_traj.set_data(tl.adjustForPlotting(phi, X[1], ratio=xscale/yscale, threshold=0.03*xscale))
		

		self.fig.canvas.draw()

		try:	self.torus.vectorField_prc()
		except: pass

		try:	self.traces.computeTraces()
		except: pass


	def on_button(self, event):

		if not event.inaxes == self.ax:
			return

		self.event_start = np.array([event.xdata, event.ydata])


	def N_output(self, CYCLES):
		return int(CYCLES*self.orbit_period/self.dt)


	def off_button(self, event):

		if not event.inaxes == self.ax:
			return

		delta_params = np.array([event.xdata, event.ydata])-self.event_start

		if event.button == 1:
			self.setParams(I=model.params['I_0']+delta_params[0], x=model.params['x_0']+delta_params[1])

		elif event.button == 3:
			new_m = model.params['m_0']+delta_params[0]
			self.setParams(m=model.params['m_0']*(new_m<0.)+(new_m>0.)*new_m, k=model.params['k_0']+delta_params[1]*3.)
			
		elif event.button == 2:
			new_eps = model.params['epsilon_0']+delta_params[0]
			self.setParams(epsilon=model.params['epsilon_0']*(new_eps<0.)+(new_eps>0.)*new_eps)
			
		self.refresh_nullclines()
		self.refresh_orbit()
		self.focus_in()



if __name__ == "__main__":
	
	import pylab as pl

	sys = system()
	pl.show()






