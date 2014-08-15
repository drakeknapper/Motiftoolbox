#! /usr/bin/env python

import multiprocessing as mp


__func = None

def distribute(f, parameter, param_values, kwargs={}, max_threads=None, nice=19, verbose=0):
	global __func
	
	def __func(value):
		kwargs[parameter] = value
		return f(**kwargs)

	p = mp.Pool(max_threads)
	res = p.map(__func, param_values)

	p.close()
	p.join()

	return res




if __name__== "__main__":

	from pylab import *
	import time

	def test_func(a, q=5):

		for i in xrange(10):
			a += q*a

		return a


	for i in xrange(100000):
		distribute(test_func, 'a', ones((100), float))

	time.sleep(10)
