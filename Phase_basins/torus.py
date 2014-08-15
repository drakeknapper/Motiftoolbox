#!/usr/bin/env python

import sys
sys.path.insert(0, '../Tools')
import fitzhugh as model
import torus_2D
import numpy as np



class torus(torus_2D.torus_2D):

        model = model
        V_trigger = 0.


	def __init__(self, system, network, traces, info=None, position=None):
                torus_2D.torus_2D.__init__(self, system, network, traces, info, position)


	def erase_traces(self):

		for j in xrange(len(self.ax_traces.lines)):
			self.ax_traces.lines.pop(0)

		for j in xrange(len(self.ax_basins.lines)):
			self.ax_basins.lines.pop(0)

		self.ax_basins.imshow(np.ones((10, 10, 4)), interpolation='nearest', aspect='equal', extent=(0., 1., 0., 1.), origin='lower')

		self.fig.canvas.draw()
		


if __name__ == "__main__":

	import system as sys
	import network as netw
	import traces as tra
	import info as nf
		
	i = nf.info()
	s = sys.system(info=i)
	n = netw.network(info=i)
	t = tra.traces(s, n, info=i)
	tor = torus(s, n, t, info=i)

	pl.show()





