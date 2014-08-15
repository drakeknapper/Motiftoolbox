#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import fitzhugh as fh
import tools as tl

import numpy as np
import pylab as pl


class info:

	def __init__(self, position=None):
		self.fig = pl.figure('Info', figsize=(3, 2), facecolor='#EEEEEE')
		self.description = self.fig.text(0.1, 0.1, 'here comes the description')

		if not position == None:

			try:
				self.fig.canvas.manager.window.wm_geometry(position)

			except:
				pass


	
	def set(self, text):
		self.description.set_text(text)
		self.fig.canvas.draw()







if __name__ == "__main__":
	

	net = info(position="+200+0")
	pl.show()




