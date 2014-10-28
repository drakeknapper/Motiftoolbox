#! /usr/bin/python


import matplotlib.pyplot as plt, mpld3 as m
from flask import Flask
import info as nf
import network as netw
import system as sys
import torus as tor
import traces as tra

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

@app.route("/")
def hello():

    x1=m.fig_to_html(n.fig)
    x2=m.fig_to_html(s.fig)
    x3=m.fig_to_html(i.fig)
    x4=m.fig_to_html(t.fig)
    x5=m.fig_to_html(tr.fig)
    x=x1+"\n"+x2+"\n"+x3+"\n"+x4+"\n"+x5+"\n"
    #print "HTML output : \n",x;
    #return "Hello World"
    return x;


@app.route("/test")
def hello1():
    x = '<h1 style="color:red">The Main Page that opens other plots</h1> \n' \
        '<div id="systemframe" style="position:relative; top: 20px; left: 0px;">' \
        '<iframe id="system" style="border-style: none; border-color: red; border-width: 1px; height:500px; width:500px" src="/system"></iframe>' \
        '</div' \
        '<div id="networkframe" style="position:relative; top: 20px; right: 0px;">' \
        '<iframe id="network" style="border-style: none; border-color: red; border-width: 1px; height:500px; width:500px;" src="/network"></iframe>'\
        '</div'
    return x;

@app.route("/system")
def fun1():
    x2=m.fig_to_html(s.fig)
    return x2;

@app.route("/network")
def fun2():
    x1=m.fig_to_html(n.fig)
    return x1;



if __name__ == "__main__":
    app.run()
    #hello();pl.show()
