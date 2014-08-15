#! /usr/bin/python

import numpy as np
import pylab as pl
import matplotlib.patches as mpatches

np.random.seed()
jitter = 0.2
rnd = lambda x: jitter*(2*np.random.rand(x)-1)

colors = ['b', 'g', 'r', 'y']
width, radius = 0.03, 0.1
three_cells = [[radius+0.01, 0.97-radius], [0.97-radius, 0.97-radius], [0.47, radius+0.01]]
four_cells = [[radius+0.1, radius+0.1], [0.9-radius, radius+0.1], [radius+0.1, 0.9-radius], [0.9-radius, 0.9-radius]]

four_cells = np.array(four_cells)+0.02*np.reshape(rnd(8), (4, 2))


def neuro_net_three(coupling_strength, ax):
	[blue, green, red] = three_cells

	patches = []
	for (i, c) in enumerate(three_cells):
		art = mpatches.Circle(c, radius, fc=colors[i])
		patches.append(art)

	sh = 0.02
	shr = 0.25
	alpha = 1.3
	bar_color = '#888888'

	#idx = [0, 1]

	# 1 blue -> green
	xbluegreen = [blue[0]+alpha*radius, green[0]-alpha*radius]
	ybluegreen = np.array([blue[1], green[1]])+sh
	ax.add_line(pl.Line2D(xbluegreen, ybluegreen, lw=5., c=bar_color))
	ax.add_patch(mpatches.Circle(np.array([xbluegreen[1], ybluegreen[1]]), radius/4., fc=bar_color))
	ax.text((blue[0]+4.*green[0])/5., (blue[1]+green[1])/2.+0.05, repr(coupling_strength[2])[:6], fontsize=15)

	# 2 green -> blue
	xgreenblue = xbluegreen
	ygreenblue = ybluegreen-2.*sh
	ax.add_line(pl.Line2D(xgreenblue, ygreenblue, lw=5., c=bar_color))
	ax.add_patch(mpatches.Circle(np.array([xgreenblue[0], ygreenblue[0]]), radius/4., fc=bar_color))
	ax.text((4.*blue[0]+green[0])/5., (blue[1]+green[1])/2.+0.05, repr(coupling_strength[0])[:6], fontsize=15)

	# 3 green -> red
	xredgreen = np.array([red[0], green[0]])*0.45+0.375
	yredgreen = np.array([red[1], green[1]])*0.45+0.2185
	ax.add_line(pl.Line2D(xredgreen, yredgreen, lw=5., c=bar_color))
	ax.add_patch(mpatches.Circle(np.array([xredgreen[0], yredgreen[0]]), radius/4., fc=bar_color))
	ax.text((green[0]+4.*red[0])/5., (green[1]+4.*red[1])/5., repr(coupling_strength[5])[:6], fontsize=15)

	# 4 red -> green
	xredgreen = np.array([red[0], green[0]])*0.45+0.33
	yredgreen = np.array([red[1], green[1]])*0.45+0.255
	ax.add_line(pl.Line2D(xredgreen, yredgreen, lw=5., c=bar_color))
	ax.add_patch(mpatches.Circle(np.array([xredgreen[1], yredgreen[1]]), radius/4., fc=bar_color))
	ax.text((4.*green[0]+red[0])/5., (4.*green[1]+red[1])/5., repr(coupling_strength[3])[:6], fontsize=15)

	# 5 red -> blue
	xbluered = np.array([blue[0], red[0]])*0.45+0.17
	ybluered = np.array([blue[1], red[1]])*0.45+0.23
	ax.add_line(pl.Line2D(xbluered, ybluered, lw=5., c=bar_color))
	ax.add_patch(mpatches.Circle(np.array([xbluered[0], ybluered[0]]), radius/4., fc=bar_color))
	ax.text((4.*blue[0]+red[0])/5., (4.*blue[1]+red[1])/5., repr(coupling_strength[1])[:6], fontsize=15)

	# 6 blue -> red
	xbluered = np.array([blue[0], red[0]])*0.45+0.2
	ybluered = np.array([blue[1], red[1]])*0.45+0.265
	ax.add_line(pl.Line2D(xbluered, ybluered, lw=5., c=bar_color))
	ax.add_patch(mpatches.Circle(np.array([xbluered[1], ybluered[1]]), radius/4., fc=bar_color))
	ax.text((blue[0]+4.*red[0])/5., (blue[1]+4.*red[1])/5., repr(coupling_strength[4])[:6], fontsize=15)

	for (i, c) in enumerate(three_cells):
		ax.add_patch(patches[i])
		ax.text(c[0]-0.02, c[1]-0.03, r'$'+str(i+1)+'$', fontsize=40)

	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_axis_off()


#=== four cells


net_keys = ['2-o1','3-o1', '4-o1', '1-o2', '3-o2', '4-o2', '1-o3', '2-o3', '4-o3', '1-o4', '2-o4', '3-o4',
		'1w2', '1w3', '2w4', '3w4']


def neuro_net_four(coupling_strength, ax):
	[blue, green, red, yellow] = four_cells

	patches = []
	for (i, c) in enumerate(four_cells):
		art = mpatches.Circle(c, radius, fc=colors[i])
		patches.append(art)

	sh = 0.02
	shr = 0.25
	alpha = 1.3
	bar_color = '#888888'

	texts = {}
	for k in coupling_strength.keys():

		try:
			[i, j] = np.array(k.split('-o'), dtype=int)-1 # i -o j
			c = coupling_strength[k]
			cell_i, cell_j = four_cells[i]+jitter*rnd(2), four_cells[j]+jitter*rnd(2)
			ax.annotate('', xytext=cell_i, xy=cell_j, arrowprops=dict(frac=0.05, width=3., shrink=0.2, color=bar_color))
			ax.text((cell_i[0]+3.*cell_j[0])/4.-0.1*(i==0)-0.02, (cell_i[1]+3.*cell_j[1])/4., 'i', color='r', fontsize=19)
			ax.text((cell_i[0]+3.*cell_j[0])/4.-0.1*(i==0), (cell_i[1]+3.*cell_j[1])/4., '%.3g' % (c), fontsize=19)
			texts[k] = ax.texts[-1]

		except:
			pass

		try:
			[i, j] = np.array(k.split('w'), dtype=int)-1 # i -o j
			c = coupling_strength[k]
			cell_i, cell_j = four_cells[i]+jitter*rnd(2), four_cells[j]+jitter*rnd(2)
			ax.annotate('', xytext=cell_i, xy=cell_j, arrowprops=dict(frac=0.0, width=3., shrink=0.2, color='y'))
			ax.annotate('', xytext=cell_j, xy=cell_i, arrowprops=dict(frac=0.0, width=3., shrink=0.2, color='y'))
			ax.text((cell_i[0]+cell_j[0])/2.-(1.-(i*j))*0.02-0.02, (cell_i[1]+cell_j[1])/2-(1.-(i*j))*0.02, 'E', color='r', fontsize=19)
			ax.text((cell_i[0]+cell_j[0])/2.-(1.-(i*j))*0.02, (cell_i[1]+cell_j[1])/2-(1.-(i*j))*0.02, '%.3g' % (c), fontsize=19)
			texts[k] = ax.texts[-1]

		except:
			pass


	for (i, c) in enumerate(four_cells):
		ax.add_patch(patches[i])
		ax.text(c[0]-0.02, c[1]-0.03, r'$'+str(i+1)+'$', fontsize=40)

	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_axis_off()

	return texts




if __name__ == '__main__':
	
	K = {}
	for k in net_keys:
		K[k] = 0.001


	fig = pl.figure()
	ax = fig.add_axes([0., 0., 1., 1.])
	neuro_net_four(K, ax=ax)
	fig.show()

	pl.show()

