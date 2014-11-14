#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import fitzhugh as fh
import tools as tl
import window as win

import numpy as np
import pylab as pl

class traces(win.window):

	title = 'Voltage Traces'
	figsize = (9, 2)

	def __init__(self, system, network, info=None, position=None):
		win.window.__init__(self, position)
		self.system = system
		self.network = network
		self.info = info
		self.CYCLES = 10

		self.ax = self.fig.add_subplot(111, frameon=False, yticks=[])

		self.li_b, = self.ax.plot([], [], 'b-', lw=2.)
		self.li_g, = self.ax.plot([], [], 'g-', lw=2.)
		self.li_r, = self.ax.plot([], [], 'r-', lw=2.)

		self.ax.set_xlabel(r'time (sec.)', fontsize=20)

		self.ax.set_xticklabels(np.arange(0., 1., 0.1), fontsize=15)
		self.ax.set_yticklabels(np.arange(0., 1., 0.1), fontsize=15)
		
		self.ax.set_xlim(0., 100.)
		self.ax.set_ylim(-5.5, 1.5)

		#self.fig.tight_layout()

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
			initial_condition = self.system.load_initial_condition(pl.rand(), pl.rand())

		V_i = fh.integrate_three_rk4(
				initial_condition,
				self.network.coupling_strength,
				self.system.dt/float(self.system.stride),
				self.system.N_output(self.CYCLES),
				self.system.stride)

		t = self.system.dt*np.arange(V_i.shape[0])

		if plotit:
			ticks = np.asarray(t[::t.size/10], dtype=int)

			xscale, yscale = t[-1], 2.
			for (i, li) in enumerate([self.li_b, self.li_g, self.li_r]):
				tj, Vj = tl.adjustForPlotting(t, V_i[:, i], ratio=xscale/yscale, threshold=0.05*xscale)
				li.set_data(tj, Vj-i*2)

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



