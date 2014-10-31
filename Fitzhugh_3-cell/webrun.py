#!/usr/bin/env python


import mpld3 as m
from mpld3 import plugins
from flask import Flask,request
import pylab as pl
import system as sys
import info as nf
import network3N as netw
import torus as tor
import traces as tra
from WebSupport.Plugins.clickPlugin import ClickPlugin
from WebSupport.Plugins.dragPlugin import DragPlugin
from matplotlib.backend_bases import Event


app = Flask(__name__)

i = nf.info()
n = netw.network(info=i)
system = sys.system(info=i, network=n)
traces = tra.traces(system, n, info=i)
torus = tor.torus(system, n, traces, info=i)
n.system = system
system.traces = traces


plugins.connect(torus.fig, ClickPlugin(eventHandlerURL="updatetorus",radioButtonID="torusRadio"))
plugins.connect(system.fig, DragPlugin(eventHandlerURL="updatesystem",radioButtonID="systemRadio"))
plugins.connect(n.fig, DragPlugin(eventHandlerURL="updatenetwork",radioButtonID="networkRadio"))

sweepphasespace = False;



@app.route("/")
def hello():
	f = open('fitzhugh_3cell.html','r')
	HTML = f.read();
	return HTML

@app.route("/system")
def fun1():
	return m.fig_to_html(system.fig)

@app.route("/torus")
def fun2():
	return m.fig_to_html(torus.fig)

@app.route("/network")
def fun3():
	return m.fig_to_html(n.fig)

@app.route("/traces")
def fun4():
	return m.fig_to_html(traces.fig)

@app.route("/info")
def fun5():
	return m.fig_to_html(i.fig)

@app.route("/pylabshow")
def fun6():
	pl.show()
	return ""

@app.route("/updatetorus")
def fun7():
	if request.args['type']=='sweep':
		global sweepphasespace
		if sweepphasespace==False :
			torus.sweep_phase_space()
			sweepphasespace = True
		return
	elif request.args['type']=='trace':
		print request.args['xval'],request.args['yval']
		torus.click_traces(float(request.args['xval']), float(request.args['yval']))
		return ""
	return ""

@app.route("/updatesystem")
def fun8():
    event = Event("mockEvent",system.fig.canvas)
    event.xdata = float(request.args['startX'])
    event.ydata = float(request.args['startY'])
    system.on_button(event)
    event.xdata = float(request.args['endX'])
    event.ydata = float(request.args['endY'])
    event.button = int(request.args['type'])
    system.off_button(event)
    return ""

@app.route("/updatenetwork")
def fun9():
    event = Event("mockEvent",n.fig.canvas);
    event.xdata = float(request.args['startX'])
    event.ydata = float(request.args['startY'])
    event.button = int(request.args['type'])
    n.on_button(event)
    event.xdata = float(request.args['endX'])
    event.ydata = float(request.args['endY'])
    event.button = int(request.args['type'])
    n.off_button(event)
    return ""

@app.route("/reset")
def fun10():
    pl.close('all')
    global i,n,system,traces,torus,sweepphasespace;
    del i,n,system,traces,torus,sweepphasespace;
    i = nf.info()
    n = netw.network(info=i,system=None)
    system = sys.system(info=i, network=n,traces=None)
    traces = tra.traces(system, n, info=i)
    torus = tor.torus(system, n, traces, info=i)
    n.system = system
    system.traces = traces
    plugins.connect(torus.fig, ClickPlugin(eventHandlerURL="updatetorus",radioButtonID="torusRadio"))
    plugins.connect(system.fig, DragPlugin(eventHandlerURL="updatesystem",radioButtonID="systemRadio"))
    plugins.connect(n.fig, DragPlugin(eventHandlerURL="updatenetwork",radioButtonID="networkRadio"))
    sweepphasespace = False;
    return ""




if __name__ == "__main__":
	#app.run()
    app.run(host='0.0.0.0')
    app.debug = True