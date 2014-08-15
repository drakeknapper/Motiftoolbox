#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import fitzhugh as fh
import tools as tl

import numpy as np
import pylab as pl
import neuro_net as nn


coupling_dict = {}
coupling_dict['2-o1'] = 0
coupling_dict['3-o1'] = 1
coupling_dict['4-o1'] = 2
coupling_dict['1-o2'] = 3
coupling_dict['3-o2'] = 4
coupling_dict['4-o2'] = 5
coupling_dict['1-o3'] = 6
coupling_dict['2-o3'] = 7
coupling_dict['4-o3'] = 8
coupling_dict['1-o4'] = 9
coupling_dict['2-o4'] = 10
coupling_dict['3-o4'] = 11

N_COUPLINGS = len(coupling_dict)


class network:

	def __init__(self, g_inh=0.003, info=None, position=None):

		self.coupling_strength = {}
		for k in coupling_dict.keys():
			self.coupling_strength[k] = g_inh


		self.COUPLING_LOCKED = True
		self.info = info
		
		self.fig = pl.figure('Network', figsize=(6, 5), facecolor='#EEEEEE')

		self.ax = self.fig.add_axes([0., 0., 1., 1.], xticks=[], yticks=[], frameon=False)
		self.texts = nn.neuro_net_four(self.coupling_strength, self.ax)

		self.show_coupling()
		self.fig.canvas.mpl_connect('button_press_event', self.on_button)
		self.fig.canvas.mpl_connect('button_release_event', self.off_button)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focus_in)

		if not position == None:
			try:
				self.fig.canvas.manager.window.wm_geometry(position)
			except:
				pass


	def get_coupling(self):
		coupling = np.zeros((N_COUPLINGS), float)

		for k in coupling_dict.keys():
			coupling[coupling_dict[k]] = self.coupling_strength[k]

		return coupling
	

	def echo(self, string):
		if self.info == None: 	print string
		else: 			self.info.set(string)


	def focus_in(self, event=None):
		descriptor = "('left : change one strength)\n"
		descriptor += "('right : change all strengths)"
		self.echo(descriptor)



	def refresh_coupling(self, k):
		c = self.coupling_strength[k]
		self.texts[k].set_text('%.4f' % (c))


	def show_coupling(self):
		c = self.coupling_strength

		for k in coupling_dict.keys():
			self.texts[k].set_text('%.4f' % (c[k]))

		self.fig.canvas.draw()


	def on_button(self, event):
		self.event_start = np.array([event.xdata, event.ydata])
		
		if event.button == 1:
			dist = np.zeros((N_COUPLINGS), float)
			keys = coupling_dict.keys()
			for (i, k) in enumerate(keys):
				dist[i] = max(abs(np.array(self.texts[k].get_position())-self.event_start))

			self.selected_key = keys[np.argmin(dist)]
			self.texts[self.selected_key].set_color('r')
			self.fig.canvas.draw()


	def off_button(self, event):
		delta = (event.ydata-self.event_start[1])/20.

		if event.button == 1:
			new_coupling = self.coupling_strength[self.selected_key]+delta
			self.coupling_strength[self.selected_key] = (new_coupling>0.)*new_coupling
			self.refresh_coupling(self.selected_key)
			self.texts[self.selected_key].set_color('k')

		else:

			for k in coupling_dict.keys():
				new_coupling = self.coupling_strength[k]+delta
				self.coupling_strength[k] = (new_coupling>0.)*new_coupling

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




