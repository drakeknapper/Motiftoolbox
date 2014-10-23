#!/usr/bin/env python

import time
import sys
sys.path.insert(0, '../Tools')
import tools as tl
import distribute as fm
import window as win
import attractor as att
import numpy as np



class torus_2D(win.window):

        model = None
	PROCESSOR = ['CPU', 'GPU']
	V_trigger = 0.
	title = 'Phase Torus'
	figsize = (10, 5)

	def __init__(self, system, network, traces, info=None, position=None):
		win.window.__init__(self, position)
		self.system = system
		self.network = network
		self.traces = traces
		self.info = info
		self.GRID = 10
		self.USE_GPU = False
                self.basin_image = None

		ticks = 0.1*np.arange(11)
		self.ax_traces = self.fig.add_subplot(121, xticks=ticks, yticks=ticks[1:])
		self.ax_traces.set_xlabel(r'$\Delta\theta_{12}$', fontsize=20)
		self.ax_traces.set_ylabel(r'$\Delta\theta_{13}$', fontsize=20)
		self.ax_traces.set_xlim(0., 1.)
		self.ax_traces.set_ylim(0., 1.)

		self.ax_basins = self.fig.add_subplot(122, xticks=ticks, yticks=ticks[1:])
		self.ax_basins.set_xlabel(r'$\Delta\theta_{12}$', fontsize=20)

		self.key_func_dict.update(dict(u=torus_2D.increase_grid, i=torus_2D.decrease_grid, E=type(self).erase_traces, g=torus_2D.switch_processor))

		self.fig.canvas.mpl_connect('button_press_event', self.on_click)
		self.fig.canvas.mpl_connect('axes_enter_event', self.focus_in)



	def on_click(self, event):
		if event.inaxes == self.ax_traces:
			self.click_traces(event.xdata, event.ydata)

		if event.inaxes == self.ax_basins:
			self.sweep_phase_space()



	def click_traces(self, x, y):
		t, V_j = self.traces.computeTraces(self.system.load_initial_condition(x, y))
		V_j = np.transpose(V_j)
		ti, d = tl.phase_difference(V_j, V_trigger=type(self).V_trigger)
		last = d[-1, :]
		tl.plot_phase_2D(d[:, 0], d[:, 1], axes=self.ax_traces, c=tl.clmap(tl.PI2*last[1], tl.PI2*last[0]), PI=0.5)
		#tl.plot_phase_2D(d[:, 0], d[:, 1], axes=self.ax_traces, c=tl.clmap_patterns(last[1], last[0]), PI=0.5)
		self.fig.canvas.draw()



	def erase_traces(self):

		for j in xrange(len(self.ax_traces.lines)):
			self.ax_traces.lines.pop(0)

		self.fig.canvas.draw()



	def erase_basins(self):

		for j in xrange(len(self.ax_basins.lines)):
			self.ax_basins.lines.pop(0)

                if not self.basin_image == None:
		    self.basin_image.set_data(np.ones((10, 10, 3), dtype=float))

		self.fig.canvas.draw()



	def switch_processor(self):

		if type(self).model.CUDA_ENABLED:
			self.USE_GPU = not self.USE_GPU

		self.focus_in()



	def increase_grid(self):
		self.GRID += 1
		self.focus_in()



	def decrease_grid(self):
		self.GRID -= 1*(self.GRID > 0)
		self.focus_in()



	def echo(self, string):
		if self.info == None: 	print string
		else: 			self.info.set(string)



	def focus_in(self, event=None):
		descriptor = "GRID : "+str(self.GRID)+" ('u' > 'i')\n"
		descriptor += "'E' : erase traces\n"
		descriptor += "'g' : switch PU (now: %s)\n\n" % (torus_2D.PROCESSOR[self.USE_GPU])
		descriptor += "'Q' : quit program"
		self.echo(descriptor)
	


	def phase_trace_cpu(self, i, X):
		V_i = np.transpose(
			type(self).model.integrate_three_rk4(
				X[i, :],
				self.network.coupling_strength,
				self.dt/float(self.stride),
				self.N_output,
				self.stride))
		ti, d = tl.phase_difference(V_i, V_trigger=type(self).V_trigger)
		return d



	def compute_phase_trace_cpu(self, X):
		return fm.distribute(self.phase_trace_cpu, 'i', range(X.shape[0]), kwargs={'X': X})



	def phase_trace_gpu(self, i, V):
		ti, d = tl.phase_difference(V[i, :, :], V_trigger=type(self).V_trigger)
		return d


	
	def compute_phase_trace_gpu(self, X):
		X = X.flatten()
		V_all = type(self).model.cuda_integrate_three_rk4(X,
				self.network.coupling_strength,
				self.dt/float(self.stride),
				self.N_output, self.stride)
		V_all = np.swapaxes(V_all, 1, 2)	 # V_all[initial condition, voltage trace, time index]
		return fm.distribute(self.phase_trace_gpu, 'i', range(V_all.shape[0]), kwargs={'V': V_all})
	


	def sweep_phase_space(self, widget=None, event=None):
		self.dt, self.stride, self.N_output = self.system.dt, self.system.stride, self.system.N_output(self.traces.CYCLES)

		self.echo('computing using '+torus_2D.PROCESSOR[self.USE_GPU]+' ...')

		torus_2D.erase_traces(self)
		self.erase_basins()

		initial_phases = np.arange(0., 1., 1./float(self.GRID))
	
		X = self.system.load_initial_conditions(initial_phases) # X[initial_condition, variable]

		t0 = time.time()
		D = compute_phase_trace[self.USE_GPU](self, X)
                self.echo('traces computed.')
	
		# find attractors
		attractors = [att.attractor( [[0., 0.]], initial_state=[0., 0.] ).classify()]
		basin_hash = -np.ones((self.GRID, self.GRID), int)
		for i in xrange(self.GRID):

			for j in xrange(self.GRID):
				trajectory = D[i*self.GRID+j]
				initial_state = initial_phases[[j, i]]

				try:
					new_attractor = att.attractor( trajectory, initial_state=initial_state ).classify()

				except:
					basin_hash[i, j] = 0
					continue
				
				conditions = [attr_i.distance2attractor(new_attractor) > 0.1 for attr_i in attractors]
				for (c, condition) in enumerate(conditions):

					if not condition:
						basin_hash[i, j] = c

				if basin_hash[i, j] < 0:
                			attractors.append(new_attractor)
					basin_hash[i, j] = len(attractors)-1



		# refine attractor
		#"""
		def refine_attractor(initial_condition):
			t, V_j = self.traces.computeTraces(initial_condition, plotit=False)
			V_j = np.transpose(V_j)
			ti, trajectory = tl.phase_difference(V_j, V_trigger=type(self).V_trigger)

			try:
				new_attr = att.attractor(trajectory).classify()
				return new_attr

			except:
				return att.attractor( [[0., 0.]], initial_state=[0., 0.] ).classify()

			

		initial_conditions = []
		for i in xrange(1, len(attractors), 1):
			attr = attractors[i]
			initial_conditions.append( self.system.load_initial_condition(attr.initial_state[1], attr.initial_state[0]))

		if len(initial_conditions):
			self.traces.CYCLES *= 4
			new_attractors = [att.attractor( [[0., 0.]], initial_state=[0., 0.] ).classify()]
			new_attractors.extend( fm.distribute(refine_attractor, 'initial_condition', initial_conditions) )
			self.traces.CYCLES /= 4
		#"""

		attractors = [att.attractor( [[0., 0.]], initial_state=[0., 0.] ).classify()]
                for (j, new_att_j) in enumerate(new_attractors):
                    conditions = [new_att_j.distance2attractor(attr_i) < 0.1 for attr_i in attractors]
		    
                    if any(conditions):
                        which = conditions.index(True)

		    else:
			which = len(attractors)
			attractors.append(new_att_j)

                
		    for l in xrange(basin_hash.shape[0]):

			for m in xrange(basin_hash.shape[1]):

			    if basin_hash[l, m] == j:
				basin_hash[l, m] = which








		basins = np.ones((self.GRID, self.GRID, 4), float)
		for i in xrange(self.GRID):

			for j in xrange(self.GRID):
				basins[j, i, :3] = attractors[basin_hash[i, j]].color[:3]


                # plot attractors
		for attractor in attractors:
                    attractor.plot(self.ax_basins, 'o', mec='w', mew=1.0, ms=7.)

		
                # plot basins
                if self.basin_image == None:
		    self.basin_image = self.ax_basins.imshow(basins, interpolation='nearest', aspect='equal', extent=(0., 1., 0., 1.), origin='lower')

                else:
		    self.basin_image.set_data(basins)
                    self.ax_basins.set_xlim(0., 1.)
                    self.ax_basins.set_ylim(0., 1.)


                # plot traces
		for i in xrange(self.GRID):
			for j in xrange(self.GRID):
				d = D[i*self.GRID+j]
				color = basins[j, i, :3]
				#"""
				try:
					tl.plot_phase_2D(d[:, 0], d[:, 1], axes=self.ax_traces, c=color, PI=0.5)
					#last = d[-1, :]
					#tl.plot_phase_2D(d[:, 0], d[:, 1], axes=self.ax_traces, c=tl.clmap(tl.PI2*last[1], tl.PI2*last[0]), PI=0.5)
					#tl.plot_phase_2D(d[:, 0], d[:, 1], axes=self.ax_traces, c=tl.clmap_patterns(last[1], last[0]), PI=0.5)
	
				except:
					pass
				#"""

		for attractor in attractors:
                    attractor.plot(self.ax_traces, 'o', mec='w', mew=1.0, ms=7.)

		self.fig.canvas.draw()


		self.echo("used %.2g sec" % (time.time()-t0))

                return basin_hash, attractors


compute_phase_trace = [torus_2D.compute_phase_trace_cpu, torus_2D.compute_phase_trace_gpu]


	



if __name__ == "__main__":

	import system as sys
	import network as netw
	import traces as tra
	import info as nf
	import pylab as pl
        import fitzhugh as model
		
	i = nf.info()
	s = sys.system(info=i)
	n = netw.network(info=i)
	t = tra.traces(s, n, info=i)
	tor = torus_2D(s, n, t, info=i)
        tor.model = model

	pl.show()





