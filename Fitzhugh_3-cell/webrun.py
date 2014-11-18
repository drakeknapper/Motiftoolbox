#!/usr/bin/env python


import mpld3 as m
from mpld3 import plugins
from flask import Flask,request
import pylab as pl
import numpy as np
import system as sys
import info as nf
import network3N as netw
import torus as tor
import traces as tra
import fitzhugh as model
from WebSupport.Plugins.clickPlugin import ClickPlugin
from WebSupport.Plugins.dragPlugin import DragPlugin
from matplotlib.backend_bases import Event

from functools import wraps
from flask import request, Response



app = Flask(__name__)


def check_auth(username, password):
	"""
	This function is called to check if a username /
	password combination is valid.
	"""
	try:	keypair = open('/home/server/web_password', 'r').readline().split(',')
	except:	keypair = ['user', 'password']

	return username == keypair[0] and password == keypair[1]


def authenticate():
	"""Sends a 401 response that enables basic auth"""
	return Response(
	'Could not verify your access level for that URL.\n'
	'You have to login with proper credentials', 401,
	{'WWW-Authenticate': 'Basic realm="Login Required"'})



def requires_auth(f):
	@wraps(f)
	def decorated(*args, **kwargs):
		auth = request.authorization
		if not auth or not check_auth(auth.username, auth.password):
			return authenticate()
		return f(*args, **kwargs)
	return decorated







info, network, system, traces, torus = None, None, None, None, None
sweepingPhasespace = False;
def initialize():
	pl.close('all')
	global info, network, system, traces, torus, sweepingPhasespace
	reload(model)

	info = nf.info()
	network = netw.network(info=info)
	system = sys.system(info=info, network=network)
	traces = tra.traces(system, network, info=info)
	torus = tor.torus(system, network, traces, info=info)
	network.system = system
	system.traces = traces


	# customize system for web
	system.setParams(epsilon=0.3)
	system.ax.set_xlabel(r'Inactivation Variable')
	system.ax.set_ylabel(r'Voltage Variable')
	system.ax.set_title('')
	system.fig.tight_layout()
	plugins.connect(system.fig, DragPlugin(eventHandlerURL="updatesystem",radioButtonID="systemRadio"))

	# customize network
	network.ax.patch.set_facecolor('#777777')
	network.moveText(2, [0.02, -0.1])
	network.moveText(3, [0.02, -0.1])
	network.ax.texts[6].set_text('1')
	network.ax.texts[7].set_text('2')
	network.ax.texts[8].set_text('3')
	plugins.connect(network.fig, DragPlugin(eventHandlerURL="updatenetwork",radioButtonID="networkRadio"))

	# customize traces
	traces.ax.patch.set_facecolor('#777777')
	traces.fig.tight_layout()

	# customize torus
	torus.ax_traces.set_xlabel(r'phase lag: 2-1')
	torus.ax_basins.set_xlabel(r'phase lag: 2-1')
	torus.ax_traces.set_ylabel(r'phase lag: 3-1')
	torus.fig.tight_layout()
	torus.switch_processor()	# switches on the gpu if available
	if torus.USE_GPU: torus.setGridsize(24)
	plugins.connect(torus.fig, ClickPlugin(eventHandlerURL="updatetorus", radioButtonID="torusRadio"))

	# reload timing variable
	sweepingPhasespace = False;





@app.route("/reset")
@requires_auth
def reset():
	initialize()
	return ""




@app.route("/")
@requires_auth
def loadLayout():
	f = open('fitzhugh_3cell.html','r')
	HTML = f.read();
	return HTML



@app.route("/system")
@requires_auth
def loadSystem():
	return m.fig_to_html(system.fig)




@app.route("/updatesystem")
@requires_auth
def systemOnclick():
	event = Event("mockEvent",system.fig.canvas)
	event.xdata = float(request.args['startX'])
	event.ydata = float(request.args['startY'])
	system.on_button(event)
	event.xdata = float(request.args['endX'])
	event.ydata = float(request.args['endY'])
	event.button = int(request.args['type'])
	system.off_button(event)
	return ""




@app.route("/torus")
@requires_auth
def loadTorus():
	return m.fig_to_html(torus.fig)




@app.route("/updatetorus")
@requires_auth
def torusOnclick():
	global sweepingPhasespace

	if request.args['type'] == 'sweep':

		if sweepingPhasespace == False :
			sweepingPhasespace = True
			torus.sweep_phase_space()
			sweepingPhasespace = False

		return ""

	elif request.args['type'] == 'trace':
		torus.click_traces(float(request.args['xval']), float(request.args['yval']))
		return ""

	return ""


@app.route("/torusClear")
@requires_auth
def torusClear():
	torus.erase_traces()
	torus.erase_basins()
	return ""




@app.route("/network")
@requires_auth
def loadNetwork():
	return m.fig_to_html(network.fig)




@app.route("/updatenetwork")
@requires_auth
def networkOnclick():
	event = Event("mockEvent",network.fig.canvas);
	event.xdata = float(request.args['startX'])
	event.ydata = float(request.args['startY'])
	event.button = int(request.args['type'])
	network.on_button(event)
	event.xdata = float(request.args['endX'])
	event.ydata = float(request.args['endY'])
	event.button = int(request.args['type'])
	network.off_button(event)
	return ""




@app.route("/traces")
@requires_auth
def loadTraces():
	return m.fig_to_html(traces.fig)




@app.route("/info")
@requires_auth
def loadInfo():
	return m.fig_to_html(info.fig)




if __name__ == "__main__":

	initialize()
	app.run(host='0.0.0.0')
	app.debug = True




















