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

                

	def vectorField_prc(self):

		phase, coupling = self.system.threeCellCoupling(0.05)
		coupling_function = prcNetwork.interp_torus_vec(phase, phase, coupling)
		coupling_function.plot(self.GRID)




	



if __name__ == "__main__":

	import system as sys
	import network3N as netw
	import traces as tra
	import info as nf
		
	info = nf.info()
	system = sys.system()
	network = netw.network(info=info)
	traces = tra.traces(system, network, info=info)
	t = torus(system, network, traces, info=info)
	system.torus=t

	pl.show()





