#!/usr/bin/env python
from matplotlib.backend_bases import Event

import mpld3 as m
from mpld3 import plugins
from flask import Flask,request
import pylab as pl
import info as nf
import network as netw
import system as sys
import torus as tor
import traces as tra
from WebSupport.Plugins.clickPlugin import ClickPlugin
from WebSupport.Plugins.dragPlugin import DragPlugin

app = Flask(__name__)

pos_info = '+0+600'
pos_tra = '+300+600'
pos_net = '+300+0'
pos_sys = '+0+0'
pos_torus = '+800+0'
i = nf.info(position=pos_info)
n = netw.network(info=i, position=pos_net)
s = sys.system(info=i, position=pos_sys, network=n)
n.system = s
t = tra.traces(s, n, info=i, position=pos_tra)
tr = tor.torus(s, n, t, info=i, position=pos_torus)

plugins.connect(tr.fig, ClickPlugin(tr.fig))
plugins.connect(s.fig,DragPlugin(s.fig))
sweepphasespace=False;



@app.route("/")
def hello():
    f = open('fitzhugh_3cell.html','r')
    HTML = f.read();
    return HTML

@app.route("/system")
def fun1():
    return m.fig_to_html(s.fig)

@app.route("/torus")
def fun2():
    return m.fig_to_html(tr.fig)

@app.route("/network")
def fun3():
    return m.fig_to_html(n.fig)

@app.route("/traces")
def fun4():
    return m.fig_to_html(t.fig)

@app.route("/info")
def fun5():
    return m.fig_to_html(i.fig)

@app.route("/pylabshow")
def fun6():
    pl.show()
    return

@app.route("/updatetorus")
def fun7():
    if request.args['type']=='sweep':
        global sweepphasespace
        if sweepphasespace==False :
            tr.sweep_phase_space()
            sweepphasespace = True
        return
    elif request.args['type']=='trace':
        print request.args['xval'],request.args['yval']
        tr.click_traces(float(request.args['xval']),float(request.args['yval']))
        return
    return

@app.route("/updatesystem")
def fun8():
    event = Event("mockEvent",s.fig.canvas);
    event.xdata = float(request.args['startX'])
    event.ydata = float(request.args['startY'])
    s.on_button(event)
    event.xdata = float(request.args['endX'])
    event.ydata = float(request.args['endY'])
    event.button = int(request.args['type'])
    s.off_button(event)

if __name__ == "__main__":
    app.run()

