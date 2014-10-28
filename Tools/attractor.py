#! /usr/bin/python

import tools as tl
import numpy as np
import pylab as pl



def phase_difference(x, y):
	return np.mod(x-y+0.5, 1.)-0.5


def phase_distance(x, y):
	return np.sqrt(np.sum(phase_difference(x, y)**2))




class attractor(object):

	def __init__(self, trajectory, initial_state=None):	 # trajectory[time, space]
		self.trajectory = np.asarray(trajectory)
		self.initial_state = initial_state

		if len(self.trajectory.shape) == 2:

			if initial_state == None:
				initial_state = self.trajectory[0, :]

			
		else:
			raise ValueError


	def distance2state(self, state):
		pass


	def distance2attractor(self, att):

		if isinstance(att, point_attractor):
			return self.distance2state(att.attractor)

		else:
			raise ValueError


	def classify(self):
		return point_attractor(self.trajectory, self.initial_state)

			
	def plot(self, axis, *args, **kwargs):
		pass




class point_attractor(attractor, object):

	def __init__(self, trajectory, initial_state=None):	 # trajectory[time, space]
		attractor.__init__(self, trajectory, initial_state)
		
		self.attractor = self.trajectory[-1, :]
		self.color = tl.clmap(tl.PI2*self.attractor[1], tl.PI2*self.attractor[0])


	def distance2state(self, state):
		return phase_distance(self.attractor, state)


	def plot(self, axis, *args, **kwargs):
		kwargs['mfc'] = self.color
		axis.plot([self.attractor[0]-1., self.attractor[0], self.attractor[0]+1., self.attractor[0]-1., self.attractor[0], self.attractor[0]+1., self.attractor[0]-1., self.attractor[0], self.attractor[0]+1.],
		  	[self.attractor[1]+1., self.attractor[1]+1., self.attractor[1]+1., self.attractor[1], self.attractor[1], self.attractor[1], self.attractor[1]-1., self.attractor[1]-1., self.attractor[1]-1.],
				*args, **kwargs)


	
if __name__  == "__main__":
	from pylab import *

	trajectory = [[1., 0.]]
	trajectory2 = [[1., 1./5.]]
	att = attractor(trajectory)
	att = att.classify()
	print type(att)


	ax = subplot(111)
	att.plot(ax, 'o')
	show()

