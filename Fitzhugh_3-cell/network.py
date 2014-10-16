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
N_TEXTS = len(text2coupling)

def set_params(**kwargs):

	if len(kwargs) == 0:
		kwargs = default_params

	for k in kwargs.keys():

		for c in ['0', 'b', 'g', 'r']:
			fh.params[k+'_'+c] = kwargs[k]

	print '# update params', kwargs


class network(win.window):

	title = "Network"
	figize = (5, 3)

	def __init__(self, g_inh=0.005, info=None, position=None, system=None):
		win.window.__init__(self, position)
		self.coupling_strength = g_inh*np.ones((6), float)
		self.COUPLING_LOCKED = True
		self.info = info
		self.system = system

		self.ax = self.fig.add_axes([-0.12, -0.1, 1.1, 1.33])
		tl.three_cells_alt(self.coupling_strength, ax=self.ax)

		
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
		c = self.coupling_strength[text2coupling[i]]+0.00001
		self.ax.texts[i].set_text('%.4f' % c)


	def show_coupling(self):
		c = self.coupling_strength+0.00001

		for i in xrange(N_TEXTS):
			self.ax.texts[i].set_text('%.4f' % c[text2coupling[i]]) # 1 -> 2

		if not self.system == None:
			self.system.refresh_nullclines()
		self.fig.canvas.draw()


	def on_button(self, event):
		self.event_start = np.array([event.xdata, event.ydata])
		
		if event.button == 1:
			dist = np.zeros((N_TEXTS), float)

			for n in xrange(N_TEXTS):
				dist[n] = max(abs(np.array(self.ax.texts[n].get_position())-self.event_start))

			self.i_text = np.argmin(dist)
			self.ax.texts[self.i_text].set_color('r')
			self.fig.canvas.draw()


	def off_button(self, event):
		delta = (event.ydata-self.event_start[1])/50.

		if event.button == 1:
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




