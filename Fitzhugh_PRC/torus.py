#!/usr/bin/env python

import sys
sys.path.insert(0, '../Tools')
import fitzhugh as model
import prcNetwork
import torus_2D
import numpy as np



class torus(torus_2D.torus_2D):

        model = model
        V_trigger = 0.
	figsize = (13, 6.5)

	def __init__(self, system, network, traces, info=None, position=None):

                torus_2D.torus_2D.__init__(self, system, network, traces, info, position)
		self.quiver = None


	def erase_traces(self):
		torus_2D.torus_2D.erase_traces(self)

		if self.quiver:
			self.quiver.remove()
			self.quiver = None
                

	def vectorField_prc(self):

		self.erase_traces()

		phase, coupling = self.system.threeCellCoupling(0.05, self.network.coupling_strength)
		coupling_function = prcNetwork.interp_torus_vec(phase, phase, coupling)
		self.quiver = coupling_function.plot(self.GRID, self.ax_traces, period=1.)

		self.fig.canvas.draw()




	



if __name__ == "__main__":

	import pylab as pl
	import system as sys
	import network3N as netw
	import traces as tra
	import info as nf
		
	info = nf.info()
	system = sys.system(info=info)
	network = netw.network(info=info)
	traces = tra.traces(system, network, info=info)
	t = torus(system, network, traces, info=info)
	system.torus = t
	t.vectorField_prc()

	pl.show()





