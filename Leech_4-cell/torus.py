#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import leech as model
import tools as tl
import distribute
import window as win

import numpy as np
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as mp3d


def guess_coord(self, xd, yd):

	if self.M is None:
	    return ''

	p = (xd, yd)
	edges = self.tunit_edges()
	ldists = [(mp3d.proj3d.line2d_seg_dist(p0, p1, p), i) for \
	        	i, (p0, p1) in enumerate(edges)]
	ldists.sort()
	edgei = ldists[0][1]
	p0, p1 = edges[edgei]
	x0, y0, z0 = p0
	x1, y1, z1 = p1
	d0 = np.hypot(x0-xd, y0-yd)
	d1 = np.hypot(x1-xd, y1-yd)
	dt = d0+d1
	z = d1/dt * z0 + d0/dt * z1
	return mp3d.proj3d.inv_transform(xd, yd, z, self.M)

Axes3D.guess_coord = guess_coord


PROCESSOR = ['CPU', 'GPU']
ticks = np.arange(0., 1.05, 0.1)


class torus(win.window):

	title = 'Phase Torus'
	figsize = (9, 8)
        V_trigger = -0.04
	testMode = False

	def __init__(self, phase_potrait, network, traces, info=None, position=None):
		win.window.__init__(self, position, infoWindow=info)
		self.system = phase_potrait
		self.network = network
		self.traces = traces
		self.GRID = 6
		self.USE_GPU = False

		self.ax = self.fig.add_subplot(111, projection='3d', xticks=ticks, yticks=ticks, zticks=ticks)
		self.ax.mouse_init(zoom_btn=None)
		self.ax.set_xlabel(r'$\Delta\theta_{12}$', fontsize=15)
		self.ax.set_ylabel(r'$\Delta\theta_{13}$', fontsize=15)
		self.ax.set_zlabel(r'$\Delta\theta_{14}$', fontsize=15)
		self.ax.set_xlim(0., 1.)
		self.ax.set_ylim(0., 1.)
		self.ax.set_zlim(0., 1.)
		self.fig.tight_layout()
		self.ax.autoscale(False)
		
		self.key_func_dict.update( dict( u=torus.increase_grid,		i=torus.decrease_grid,
						R=torus.show_torus,		E=torus.resetTorus,	v=torus.reset_view ) )
						#g=torus.switchProcessor))


		self.fig.canvas.mpl_connect('button_press_event', self.clickTorus)



	def clickTorus(self, event):

		if event.button == 1:
			return

		x, y, z = self.ax.guess_coord(event.xdata, event.ydata)

		for z in np.arange(0., 1., 1./float(self.GRID)) + 1./float(2*self.GRID):
			t, V_j = self.traces.compute_traces(self.system.load_initial_condition(x, y, z))
			V_j = np.transpose(V_j)
			ti, d = tl.phase_difference(V_j, V_trigger=torus.V_trigger)
			last = d[-1, :]
			tl.plot_phase_3D(d[:, 0], d[:, 1], d[:, 2], self.ax, c=tl.clmap_patterns(last[1], last[0]), PI=0.5)

		self.fig.canvas.draw()



	def increase_grid(self):
		self.GRID += 1
		self.focusIn()



	def decrease_grid(self):
		self.GRID -= 1*(self.GRID > 0)
		self.focusIn()


	

	def focusIn(self, event=None):
		descriptor = "GRID : "+str(self.GRID)+" ('u' > 'i')\n"
		descriptor += "'R' : Run torus\n"
		descriptor += "'v' : Reset view\n"
		descriptor += "'E' : Reset torus\n"
		if model.CUDA_ENABLED:	descriptor += "'g' : switch processor (%s)" % (PROCESSOR[self.USE_GPU])
		self.echo(descriptor)
	


	def phaseTrace_gpu(i, V):
		ti, d = tl.phase_difference(V[i, :, :], V_trigger=torus.V_trigger)
		return d
	
	

	def phaseTrace_cpu(self, i, X):
		V_i = np.transpose(
			model.integrate_four_rk4(
				X[i, :],
				self.network.get_coupling(),
				self.dt/float(self.stride),
				self.N_output, self.stride))
		ti, d = tl.phase_difference(V_i, V_trigger=torus.V_trigger)
		return d



	def computePhaseTrace_cpu(self, X):
		return distribute.distribute(self.phaseTrace_cpu, 'i', range(X.shape[0]), kwargs={'X': X})



	def computePhaseTrace_gpu(self, X):
		X = X.flatten()
		V_all = model.cuda_integrate_four_rk4(X,
				self.network.get_coupling(),
				self.dt/float(self.stride),
				self.N_output, self.stride)

		V_all = np.swapaxes(V_all, 1, 2)	 # V_all[initial condition, voltage trace, time index]
		return distribute.distribute(self(type).phaseTrace_gpu, 'i', range(V_all.shape[0]), kwargs={'V': V_all})



	def reset_view(self, widget=None, event=None):
		self.ax.view_init(elev=90, azim=0)
		self.fig.canvas.draw()


	
	def resetTorus(self, widget=None, event=None):

		for j in xrange(len(self.ax.lines)):
			self.ax.lines.pop(0)

		self.reset_view()
		self.fig.canvas.draw()
	


	def show_torus(self, widget=None, event=None):
		self.dt, self.stride, self.N_output = self.system.dt, self.system.stride, self.system.N_output(self.traces.CYCLES)
		self.echo('computing using '+PROCESSOR[self.USE_GPU]+' ...')
		self.fig.canvas.draw()
		initial_phases = np.arange(0., 1., 1./float(self.GRID)) + 1./float(2*self.GRID)
		X = self.system.load_initial_conditions(initial_phases) # X[initial_condition, variable]
		t0 = time.time()
		D = computePhaseTrace[self.USE_GPU](self, X)
		self.echo("used %.2g sec" % (time.time()-t0))
	
		for j in xrange(len(self.ax.lines)):
			self.ax.lines.pop(0)
	
		for d in D:
	
			try:
				last = d[-1, :]
				tl.plot_phase_3D(d[:, 0], d[:, 1], d[:, 2], self.ax, c=tl.clmap_patterns(last[1], last[0]), PI=0.5, lw=0.25)
	
			except:
				pass
		
		self.fig.canvas.draw()


	def switchProcessor(self):

		if model.CUDA_ENABLED:
			self.USE_GPU = not self.USE_GPU
			self.focusIn()
	
        	else:
			print 'GPU disabled'
			pass



computePhaseTrace = [torus.computePhaseTrace_cpu, torus.computePhaseTrace_gpu]


	



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





