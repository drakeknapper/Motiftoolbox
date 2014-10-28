#! /usr/bin/python

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

try:
    from flask import Flask
    print '# web server available'

except:
    print '# web server not available'


for k in plt.rcParams.keys():

	if k.split('.')[0] == 'keymap':
		plt.rcParams[k] = ''


class window(object):

	title = "Window"
	figsize = (5, 3)

	def __init__(self, position=None):
		self.fig = plt.figure(type(self).title, figsize=type(self).figsize)

                self.key_func_dict = dict(Q=window.quit)
		self.fig.canvas.mpl_connect('key_press_event', self.on_key)


		if not position == None:

			try:
				self.fig.canvas.manager.window.wm_geometry(position)

			except:
				pass


	def on_key(self, event):
                #self.key_func_dict[event.key](self)
		#"""
                if event.key == 'Q':
		        self.key_func_dict[event.key](self)
                        exit(0)

                else:
		        try:
                                self.key_func_dict[event.key](self)
		        except:
                                self.key_func_dict[event.key] = lambda x: None
		#"""



        def quit(self):
            print "exiting ..."
            exit(0)


if __name__ == "__main__":
	win = window()
	plt.show()
