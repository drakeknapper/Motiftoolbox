#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import tools as tl
import window as win

import numpy as np


class info(win.window):

	title = 'Info'
	figsize = (3, 2)

	def __init__(self, position=None):
		win.window.__init__(self, position)
		self.description = self.fig.text(0.1, 0.1, 'here comes the description')


	def set(self, text):
		self.description.set_text(text)
		self.fig.canvas.draw()







if __name__ == "__main__":
	
	import pylab as pl

	net = info()
	pl.show("+200+0")




