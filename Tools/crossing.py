#!/usr/bin/python

import scipy.weave as weave















if __name__ == "__main__":

	import numpy as np
	import pylab as pl


	dt = 0.001
	f = 1.
	T = 10.
	t = np.arange(0., T, dt)
	V = np.sin(2.*np.pi*f*t)


	pl.plot(t, V, 'k-')

	Cr = dt*crossing(V, threshold=0.)
	for t_i in Cr:
		pl.axvline(x=t_i, color='r')

	pl.show()
