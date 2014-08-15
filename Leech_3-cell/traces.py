#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import leech as model
import tools as tl
import window as win
import fork_master as fm
import gtk
import multiprocessing as mp
import numpy as np
import pylab as pl
from matplotlib import animation




class traces(win.window):

	title = 'Membrane Voltage Traces (Oscilloscope)'
	figsize = (16, 4)

	def __init__(self, phase_potrait, network, info=None, position=None):
		win.window.__init__(self, position)
		self.system = phase_potrait
		self.network = network
		self.info = info
		self.CYCLES = 8
		self.state = np.array([0., 1., 0., 0., 1., 0., 0., 1., 0.])
		self.running = False
		self.pause = False
		self.pulsed = 0
		self.initial_condition = self.system.load_initial_condition(pl.rand(), pl.rand())

		self.ax = self.fig.add_subplot(111, frameon=False, yticks=[])

		self.li_b, = self.ax.plot([], [], 'b-', lw=1.)
		self.li_g, = self.ax.plot([], [], 'g-', lw=1.)
		self.li_r, = self.ax.plot([], [], 'r-', lw=1.)

		self.ax.set_xlabel(r'time (sec.)', fontsize=20)

		self.ax.set_xticklabels(np.arange(0., 1., 0.1), fontsize=15)
		self.ax.set_yticklabels(np.arange(0., 1., 0.1), fontsize=15)
		
		self.ax.set_xlim(0., 100.)
		self.ax.set_ylim(-0.06-0.12, 0.04)

		self.key_func_dict = dict(u=traces.increase_cycles, i=traces.decrease_cycles)
		self.fig.canvas.mpl_connect('button_press_event', self.on_click)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focus_in)


	def adjust_cycles(self, adjustment):
		self.CYCLES = adjustment
		self.compute_traces()
		self.focus_in()

	
	def increase_cycles(self):
		self.adjust_cycles(self.CYCLES+1)


	def decrease_cycles(self):
		self.adjust_cycles(self.CYCLES-1*(self.CYCLES>0))


	def focus_in(self, event=None):
		descriptor = "CYCLES : "+str(self.CYCLES)+" ('u' > 'i')"

		if self.info == None: print descriptor
		else: self.info.set(descriptor)
		


	def on_click(self, event):

		if event.inaxes == self.ax:

		        self.params = model.params_three()
		        self.coupling = np.zeros((9), float)
		        self.coupling[:6] = self.network.coupling_strength
		        length = self.system.N_output(self.CYCLES)
		        self.t = self.system.dt*np.arange(length)
		        self.trajectory = np.zeros((length, 3), float)

		        self.li_b.set_data(self.t, self.trajectory[:, 0])
		        self.li_g.set_data(self.t, self.trajectory[:, 1]-0.06)
		        self.li_r.set_data(self.t, self.trajectory[:, 2]-0.12)

			ticks = np.asarray(self.t[::self.t.size/10], dtype=int)
			self.ax.set_xticks(ticks)
			self.ax.set_xticklabels(ticks)
                        self.ax.set_xlim(self.t[0], self.t[-1])

                        self.fig.canvas.draw()

                        self.anim = animation.FuncAnimation(self.fig, self.compute_step, init_func=self.init, frames=self.trajectory.shape[0], interval=0, blit=True, repeat=False)



        def init(self):
	        self.li_b.set_data([], [])
	        self.li_g.set_data([], [])
	        self.li_r.set_data([], [])
                
                return self.li_b, self.li_g, self.li_r



	def compute_step(self, idx):



		model.step_three_rk4(
			self.state, self.params, self.coupling,
			self.system.dt/float(self.system.stride),
			self.system.stride)

		self.trajectory[idx, :] = self.state[::model.N_EQ1]

		if self.pulsed:
			self.trajectory[idx-1, self.pulsed-1] = -0.1
			self.trajectory[idx, self.pulsed-1] = 0.1
			self.pulsed = 0

		self.li_b.set_data(self.t, self.trajectory[:, 0])
		self.li_g.set_data(self.t, self.trajectory[:, 1]-0.06)
		self.li_r.set_data(self.t, self.trajectory[:, 2]-0.12)

                return self.li_b, self.li_g, self.li_r






	def compute_traces(self, initial_condition=None, plotit=True):

		if initial_condition == None:
			initial_condition = self.initial_condition

		V_i = model.integrate_three_rk4(
				initial_condition,
				self.network.coupling_strength,
				self.system.dt/float(self.system.stride),
				self.system.N_output(self.CYCLES),
				self.system.stride)

		t = self.system.dt*np.arange(V_i.shape[0])

		if plotit:
			ticks = np.asarray(t[::t.size/10], dtype=int)
			self.li_b.set_data(t, V_i[:, 0])
			self.li_g.set_data(t, V_i[:, 1]-0.06)
			self.li_r.set_data(t, V_i[:, 2]-0.12)
			self.ax.set_xticks(ticks)
			self.ax.set_xticklabels(ticks)
			self.ax.set_xlim(t[0], t[-1])
			self.fig.canvas.draw()

		return t, V_i


	



if __name__ == "__main__":

	import system as sys
	import network as netw
		
	s = sys.system()
	n = netw.network()
	

	tra = traces(s, n)


	pl.show()



