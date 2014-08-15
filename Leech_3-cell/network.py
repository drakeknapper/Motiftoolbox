#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import window as win
import fitzhugh as fh
import tools as tl

import numpy as np
import pylab as pl


win_width, win_height, margin = 700, 600, 10


text2coupling = {}
text2coupling[0] = 2
text2coupling[1] = 0
text2coupling[2] = 5
text2coupling[3] = 3
text2coupling[4] = 1
text2coupling[5] = 4


coupling2text = {}
coupling2text[2] = 0
coupling2text[0] = 1
coupling2text[5] = 2
coupling2text[3] = 3
coupling2text[1] = 4
coupling2text[4] = 5
N_COUPLING = len(text2coupling)


def set_params(**kwargs):

	if len(kwargs) == 0:
		kwargs = default_params

	for k in kwargs.keys():

		for c in ['0', 'b', 'g', 'r']:
			fh.params[k+'_'+c] = kwargs[k]

	print '# update params', kwargs



class network(win.window):

	title = "Neural Network Motif"
	figize = (5, 3)

	def __init__(self, g_inh=0.01, info=None, position=None):
		win.window.__init__(self, position)
		self.coupling_strength = g_inh*np.ones((6), float)
		self.COUPLING_LOCKED = True
		self.info = info
		self.traces = None

		self.ax = self.fig.add_axes([-0.12, -0.1, 1.1, 1.33])
		tl.three_cells_alt(1000.*self.coupling_strength, ax=self.ax)

                self.ax.text(0.2, 0.1, 'inhibitory coupling strength in pS', fontsize=14)

		
		self.show_coupling()
		self.fig.canvas.mpl_connect('button_press_event', self.on_button)
		self.fig.canvas.mpl_connect('button_release_event', self.off_button)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focus_in)


	def echo(self, string):
		if self.info == None: 	print string
		else: 			self.info.set(string)


	def focus_in(self, event=None):
		descriptor = "('left : change one strength)\n"
		descriptor += "('right : change all strengths)"
		self.echo(descriptor)



	def refresh_coupling(self, i):
		c = 1000.*self.coupling_strength[text2coupling[i]]
		self.ax.texts[i].set_text('%.2f' % c)


	def show_coupling(self):
		c = 1000.*self.coupling_strength

		for i in xrange(N_COUPLING):
			self.ax.texts[i].set_text('%.2f' % c[text2coupling[i]])

		self.fig.canvas.draw()


	def on_button(self, event):
		self.event_start = np.array([event.xdata, event.ydata])
		
		if event.button == 1:
			dist = np.zeros((N_COUPLING+3), float)

			for n in xrange(N_COUPLING+3):
				dist[n] = max(abs(np.array(self.ax.texts[n].get_position())-self.event_start))

			self.i_text = np.argmin(dist)

			self.ax.texts[self.i_text].set_color('r')
			self.fig.canvas.draw()

			if not self.i_text < N_COUPLING:
				self.traces.state[3*(self.i_text-N_COUPLING)] += 0.001 # voltage pulse
				self.traces.pulsed = (self.i_text-N_COUPLING)+1


	def off_button(self, event):

		delta = (event.ydata-self.event_start[1])/50.

		if event.button == 1:

			if self.i_text < N_COUPLING:
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


	def unlock_coupling(self, switch, data=None):
	
		if switch.get_active():
			self.COUPLING_LOCKED = False
	
        	else:
			self.COUPLING_LOCKED = True
	








if __name__ == "__main__":
		

	net = network()

	pl.show()




