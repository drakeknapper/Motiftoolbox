#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import network3N as netw
import tools as tl

import numpy as np






class network(netw.network):

	title = "Neural Network Motif"

	def __init__(self, g_inh=0.01, info=None, position=None):
		netw.network.__init__(self, g_inh, info, position)

                self.ax.text(0.2, 0.1, 'inhibitory coupling strength in nS', fontsize=14)

		

	def on_button(self, event):
		self.event_start = np.array([event.xdata, event.ydata])
		
		if event.button == 1:
			dist = np.zeros((netw.N_COUPLING+3), float)

			for n in xrange(netw.N_COUPLING+3):
				dist[n] = max(abs(np.array(self.ax.texts[n].get_position())-self.event_start))

			self.i_text = np.argmin(dist)

			self.ax.texts[self.i_text].set_color('r')
			self.fig.canvas.draw()

			if not self.i_text < netw.N_COUPLING:
				self.traces.state[3*(self.i_text-netw.N_COUPLING)] += 0.001 # voltage pulse
				self.traces.pulsed = (self.i_text-netw.N_COUPLING)+1


	def off_button(self, event):

		delta = (event.ydata-self.event_start[1])/50.

		if event.button == 1:

			if self.i_text < netw.N_COUPLING:
				i_coupling = text2coupling[self.i_text]
				new_coupling = self.coupling_strength[i_coupling]+delta
				self.coupling_strength[i_coupling] = (new_coupling>0.)*new_coupling
				self.refresh_coupling(self.i_text)

			self.ax.texts[self.i_text].set_color('k')

		else:
			new_coupling = self.coupling_strength+delta
			self.coupling_strength = (new_coupling>0.)*new_coupling
			self.show_coupling()



		self.fig.canvas.draw()




if __name__ == "__main__":
		
	import pylab as pl

	net = network()

	pl.show()




