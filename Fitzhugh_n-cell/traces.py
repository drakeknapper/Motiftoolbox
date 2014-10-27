#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import fitzhugh as model
import tools as tl
import window as win
from matplotlib import animation

import numpy as np
import pylab as pl

class traces(win.window):

	title = 'Dynamics traces (Oscilloscope)'
	figsize = (16, 4)

	def __init__(self, phase_potrait, info=None, position=None):
		win.window.__init__(self, position)
		self.system = phase_potrait
		self.info = info
		self.CYCLES = 8
		self.running = False
		self.pulsed = 0
                self.num_osci = 8
		self.state = np.random.randn(self.num_osci*model.N_EQ1)
		self.initial_condition = self.system.load_initial_condition(pl.rand(), pl.rand())

		self.ax = self.fig.add_subplot(111, frameon=False, yticks=[])

                self.li = [self.ax.plot([], [], 'k-', lw=2.)[0] for i in xrange(self.num_osci)]

		self.ax.set_xlabel(r'time (sec.)', fontsize=20)

		self.ax.set_xticklabels(np.arange(0., 1., 0.1), fontsize=15)
		self.ax.set_yticklabels(np.arange(0., 1., 0.1), fontsize=15)
		
		self.ax.set_xlim(0., 100.)
		self.ax.set_ylim(-1.5-self.num_osci*2, 1.5)

		self.key_func_dict.update(dict(u=traces.increase_cycles, i=traces.decrease_cycles))
		self.fig.canvas.mpl_connect('button_press_event', self.on_click)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focus_in)



	def adjust_cycles(self, adjustment):
		self.CYCLES = adjustment


	def increase_cycles(self):
                self.adjust_cycles(self.CYCLES+1)
		self.focus_in()


	def decrease_cycles(self):
                self.adjust_cycles(self.CYCLES-1*(self.CYCLES>0))
		self.focus_in()



	def focus_in(self, event=None):
		descriptor = "CYCLES : "+str(self.CYCLES)+" ('u' > 'i')"
		if self.info == None: print descriptor
		else: self.info.set(descriptor)

		


	def on_click(self, event):

		if event.inaxes == self.ax:

		        self.params = model.parameters_n()
		        self.g_inh = model.params['g_inh_0']

		        length = self.system.N_output(self.CYCLES)
		        self.t = self.system.dt*np.arange(length)
		        self.trajectory = np.zeros((length, self.num_osci), float)

                        for i in xrange(self.num_osci):
		            self.li[i].set_data(self.t, self.trajectory[:, i]-i*2.)

			ticks = np.asarray(self.t[::self.t.size/10], dtype=int)
			self.ax.set_xticks(ticks)
			self.ax.set_xticklabels(ticks)
                        self.ax.set_xlim(self.t[0], self.t[-1])

                        self.fig.canvas.draw()

                        self.anim = animation.FuncAnimation(self.fig, self.compute_step, init_func=self.init, frames=self.trajectory.shape[0], interval=0, blit=True, repeat=False)



        def init(self):
                for i in xrange(self.num_osci):
	            self.li[i].set_data([], [])
                
                return self.li



	def compute_step(self, idx):

		self.params = model.parameters_n()
		self.g_inh = model.params['g_inh_0']

		model.step_n_em(
			self.state, self.params, self.g_inh,
			self.system.dt/float(self.system.stride),
			self.system.stride)

		self.trajectory[idx, :] = self.state[::model.N_EQ1]

		if self.pulsed:
			self.trajectory[idx-1, self.pulsed-1] = -0.1
			self.trajectory[idx, self.pulsed-1] = 0.1
			self.pulsed = 0

                for i in xrange(self.num_osci):
		    self.li[i].set_data(self.t, self.trajectory[:, i]-2.*i)

                return self.li




	def compute_traces(self, initial_condition=None, plotit=True):

		if initial_condition == None:
			initial_condition = self.initial_condition

		V_i = model.integrate_three_rk4(
				initial_condition,
				model.params['g_inh_0']*np.ones((6), float),
				self.system.dt/float(self.system.stride),
				self.system.N_output(self.CYCLES),
				self.system.stride)

		t = self.system.dt*np.arange(V_i.shape[0])

		if plotit:
			ticks = np.asarray(t[::t.size/10], dtype=int)

                        for i in xrange(self.num_osci):
			    self.li[i].set_data(t, V_i[:, i]-i*2.)

			self.ax.set_xticks(ticks)
			self.ax.set_xticklabels(ticks)
			self.ax.set_xlim(t[0], t[-1])
			self.fig.canvas.draw()

		return t, V_i






	


if __name__ == "__main__":

	import system as sys
		
	s = sys.system()
	

	tra = traces(s)


	pl.show()



