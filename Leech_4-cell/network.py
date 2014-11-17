#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import tools as tl
import window as win

import numpy as np
import pylab as pl
import neuro_net as nn


coupling_dict = {}
coupling_dict['2-o1'] = 0
#coupling_dict['3-o1'] = 1
coupling_dict['4-o1'] = 2
coupling_dict['1-o2'] = 3
coupling_dict['3-o2'] = 4
#coupling_dict['4-o2'] = 5
coupling_dict['1-o3'] = 6
#coupling_dict['2-o3'] = 7
coupling_dict['4-o3'] = 8
#coupling_dict['1-o4'] = 9
coupling_dict['2-o4'] = 10
coupling_dict['3-o4'] = 11
coupling_dict['1w2'] = 12

COUPLINGS = len(coupling_dict)
N_COUPLINGS = np.max(coupling_dict.values())+1


class network(win.window):

	title='Network'
	figsize=(6, 5)
	testMode = True

	def __init__(self, g_inh=0.003, info=None, position=None):

		win.window.__init__(self, position, info)

		self.coupling_strength = {}
		for k in coupling_dict.keys():
			self.coupling_strength[k] = g_inh


		self.COUPLING_LOCKED = True
		self.info = info
		
		self.ax = self.fig.add_axes([0., 0., 1., 1.], xticks=[], yticks=[], frameon=False)
		self.texts = nn.neuro_net_four(self.coupling_strength, self.ax)

		self.key_func_dict.update({'0' : type(self).resetCoupling})

		self.fig.canvas.mpl_connect('button_press_event', self.on_button)
		self.fig.canvas.mpl_connect('button_release_event', self.off_button)




	def get_coupling(self):
		coupling = np.zeros((N_COUPLINGS), float)

		for k in coupling_dict.keys():
			coupling[coupling_dict[k]] = self.coupling_strength[k]

		return coupling
	


	def focusIn(self, event=None):
		descriptor = "('left : change one strength)\n"
		descriptor += "('right : change all strengths)"
		self.echo(descriptor)



	def setCoupling(self, key, value):

		self.coupling_strength[key] = value
		self.texts[key].set_text('%.4f' % (value))



	def resetCoupling(self):

		for key in coupling_dict.keys():
			self.setCoupling(key=key, value=0.)

		self.fig.canvas.draw()



	def on_button(self, event):
		self.event_start = np.array([event.xdata, event.ydata])
		
		if event.button == 1:
			dist = np.zeros((COUPLINGS), float)
			keys = coupling_dict.keys()
			for (i, k) in enumerate(keys):
				dist[i] = max(abs(np.array(self.texts[k].get_position())-self.event_start))

			self.selected_key = keys[np.argmin(dist)]
			self.texts[self.selected_key].set_color('r')
			self.fig.canvas.draw()


	def off_button(self, event):
		delta = (event.ydata-self.event_start[1])/20.

		if event.button == 1:
			self.setCoupling(self.selected_key, self.coupling_strength[self.selected_key]+delta)
			self.texts[self.selected_key].set_color('k')

		else:
			for k in coupling_dict.keys():
				self.setCoupling(k, self.coupling_strength[k]+delta)


		self.fig.canvas.draw()


	








if __name__ == "__main__":
		

	net = network()

	pl.show()




