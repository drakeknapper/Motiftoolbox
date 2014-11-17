#!/usr/bin/env python

import sys
sys.path.insert(0, '../Tools')
import window as win
import numpy as np
import pylab as pl
import matplotlib.patches as mpatches


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


xb, yb, xg, yg, xr, yr, width, radius = 0.23, 0.65, 0.77, 0.65, 0.47, 0.25, 0.03, 0.1


class network(win.window):

	title = "3-Node Network"
	figize = (5, 3)

	def __init__(self, g_inh=0.005, info=None, position=None, system=None):
		win.window.__init__(self, position)
		self.coupling_strength = g_inh*np.ones((6), float)
		self.COUPLING_LOCKED = True
		self.info = info
		self.system = system

		self.ax = self.fig.add_axes([-0.12, -0.1, 1.1, 1.33])
		self.three_cells_alt(self.coupling_strength, ax=self.ax)

		
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

		for i in xrange(N_COUPLING):
			self.ax.texts[i].set_text('%.4f' % c[text2coupling[i]]) # 1 -> 2

		if not self.system == None:
			self.system.refresh_nullclines()
		self.fig.canvas.draw()



	def on_button(self, event):
		self.event_start = np.array([event.xdata, event.ydata])
		
		if event.button == 1:
			dist = np.zeros((N_COUPLING), float)

			for n in xrange(N_COUPLING):
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




	def three_cells_alt(self, coupling_strength, ax):

		c = coupling_strength+0.00001
		patches = []
		art = mpatches.Circle(np.array([xb, yb]), radius, fc='b')
		patches.append(art)
		art = mpatches.Circle(np.array([xg, yg]), radius, fc='g')
		patches.append(art)
		art = mpatches.Circle(np.array([xr, yr]), radius, fc='r')
		patches.append(art)

		sh = 0.02
		shr = 0.25
		alpha = 1.3

		# 1 blue -> green
		xbg = [xb+alpha*radius, xg-alpha*radius]
		ybgreen = np.array([yb, yg])+sh
		ax.add_line(pl.Line2D(xbg, ybgreen, lw=5., c='k'))
		ax.add_patch(mpatches.Circle(np.array([xbg[1], ybgreen[1]]), radius/4., fc='k'))
		ax.text((xb+xg)/2.-0.05, (yb+yg)/2.+0.05, repr(c[2])[:6], fontsize=15)

		# 2 green -> blue
		xgb = xbg
		ygb = ybgreen-2.*sh
		ax.add_line(pl.Line2D(xgb, ygb, lw=5., c='k'))
	
		ax.add_patch(mpatches.Circle(np.array([xgb[0], ygb[0]]), radius/4., fc='k'))

		ax.text((xb+xg)/2.-0.05, (yb+yg)/2.-0.1, repr(c[0])[:6], fontsize=15)

		# 3 green -> red
		xrgreen = np.array([xr+0., xg+0.])*0.45+0.375
		yrg = np.array([yr+0., yg+0.])*0.45+0.2185
		ax.add_line(pl.Line2D(xrgreen, yrg, lw=5., c='k'))
		ax.add_patch(mpatches.Circle(np.array([xrgreen[0], yrg[0]]), radius/4., fc='k'))
		ax.text((xg+xr)/2.+0.015, (yg+yr)/2.-0.015, repr(c[5])[:6], fontsize=15, rotation=45)
	
		# 4 red -> green
		xrgreen = np.array([xr+0., xg+0.])*0.45+0.33
		yrg = np.array([yr+0., yg+0.])*0.45+0.255
		ax.add_line(pl.Line2D(xrgreen, yrg, lw=5., c='k'))
		ax.add_patch(mpatches.Circle(np.array([xrgreen[1], yrg[1]]), radius/4., fc='k'))
		ax.text((xg+xr)/2.-0.1, (yg+yr)/2.+0.07, repr(c[3])[:6], fontsize=15, rotation=45)
	
		# 5 red -> blue
		xbr = np.array([xb+0., xr+0.])*0.45+0.17
		ybr = np.array([yb+0., yr+0.])*0.45+0.23
		ax.add_line(pl.Line2D(xbr, ybr, lw=5., c='k'))
		ax.add_patch(mpatches.Circle(np.array([xbr[0], ybr[0]]), radius/4., fc='k'))
		ax.text((xb+xr)/2.-0.1, (yb+yr)/2.-0.0, repr(c[1])[:6], fontsize=15, rotation=-55)
	
		# 6 blue -> red
		xbr = np.array([xb+0., xr+0.])*0.45+0.2
		ybr = np.array([yb+0., yr+0.])*0.45+0.265
		ax.add_line(pl.Line2D(xbr, ybr, lw=5., c='k'))
		ax.add_patch(mpatches.Circle(np.array([xbr[1], ybr[1]]), radius/4., fc='k'))
		ax.text((xb+xr)/2.-0.02, (yb+yr)/2.+0.1, repr(c[4])[:6], fontsize=15, rotation=-55)
	
		ax.add_patch(patches[0])
		ax.add_patch(patches[1])
		ax.add_patch(patches[2])
	
		ax.text(xb-0.04, yb-0.05, r'$1$', fontsize=40)
		ax.text(xg-0.04, yg-0.05, r'$2$', fontsize=40)
		ax.text(xr-0.04, yr-0.05, r'$3$', fontsize=40)
	
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_axis_off()


	def moveText(self, nText, vector):
		self.ax.texts[nText].set_position( np.asarray(self.ax.texts[nText].get_position())
									+ np.asarray(vector) )

	





if __name__ == "__main__":
		
	import pylab as pl

	net = network()

	pl.show()




