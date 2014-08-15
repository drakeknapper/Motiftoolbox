#!/usr/bin/env python

import sys
sys.path.insert(0, '../Tools')
import leech as model
import torus_2D
import numpy as np



class torus(torus_2D.torus_2D):

        model = model
        V_trigger = -0.04

	def __init__(self, system, network, traces, info=None, position=None):
                torus_2D.torus_2D.__init__(self, system, network, traces, info, position)
		self.ax_traces.set_xlabel(r'phase difference $\Delta\theta_{12}$', fontsize=20)
		self.ax_traces.set_ylabel(r'phase difference $\Delta\theta_{13}$', fontsize=20)
		self.ax_basins.set_xlabel(r'phase difference $\Delta\theta_{12}$', fontsize=20)
                


	



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





