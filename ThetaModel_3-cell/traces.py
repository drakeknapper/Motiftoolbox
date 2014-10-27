#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import thetax2 as th2
import tools as tl
import window as win
import fork_master as fm

import numpy as np
import pylab as pl

class traces(win.window):

	title = 'Amplitude Traces'
	figsize = (9, 2)

	def __init__(self, phase_potrait, network, info=None, position=None):
		win.window.__init__(self, position)

		self.system = phase_potrait
		self.network = network
		self.CYCLES = 10
		self.info = info
		self.initial_condition = self.system.load_initial_condition(np.random.rand(), np.random.rand())

		self.ax = self.fig.add_subplot(111, frameon=False, yticks=[])

		self.li_b, = self.ax.plot([], [], 'b-', lw=2.)
		self.li_g, = self.ax.plot([], [], 'g-', lw=2.)
		self.li_r, = self.ax.plot([], [], 'r-', lw=2.)

		self.ax.set_xlabel(r'time (sec.)', fontsize=20)

		self.ax.set_xticklabels(np.arange(0., 1., 0.1), fontsize=15)
		self.ax.set_yticklabels(np.arange(0., 1., 0.1), fontsize=15)
		
		self.ax.set_xlim(0., 100.)
		self.ax.set_ylim(-5.5, 1.5)

		self.key_func_dict = dict(u=traces.increaseCycles, i=traces.decreaseCycles)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focusIn)

		self.computeTraces()



	def adjustCycles(self, adjustment):
		self.CYCLES = adjustment
		self.computeTraces()
		self.focusIn()


	def increaseCycles(self):
		self.adjustCycles(self.CYCLES+1)

	def decreaseCycles(self):
		self.adjustCycles(self.CYCLES-1*(self.CYCLES>0))


	def focusIn(self, event=None):
		descriptor = "CYCLES : "+str(self.CYCLES)+" ('u' > 'i')"

		if self.info == None:
			print descriptor

		else:
			self.info.set(descriptor)


	def on_key(self, event):

		try:
			self.key_func_dict[event.key](self)

		except:
			self.key_func_dict[event.key] = lambda x: None


	def computeTraces(self, initial_condition=None, plotit=True):

		if initial_condition == None:
			initial_condition = self.initial_condition

		V_i = th2.integrate_three_rk4(
				initial_condition,
				self.network.coupling_strength,
				self.system.dt/float(self.system.stride),
				self.system.N_output(self.CYCLES),
				self.system.stride)

		t = self.system.dt*np.arange(V_i.shape[0])

		if plotit:
			ticks = np.asarray(t[::t.size/10], dtype=int)
			self.li_b.set_data(t, np.cos(V_i[:, 0]))
			self.li_g.set_data(t, np.cos(V_i[:, 1])-2.)
			self.li_r.set_data(t, np.cos(V_i[:, 2])-4.)
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



