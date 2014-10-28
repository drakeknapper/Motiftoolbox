#!/usr/bin/env python



#import Tkinter
import system as sys
import network as netw
import traces as tra
import info as nf
import torus as tor
import pylab as pl

#root = Tkinter.Tk()
#screen_width = root.winfo_screenwidth()
#screen_height = root.winfo_screenheight() 

pos_info = '+0+600'
pos_tra = '+300+600'
pos_net = '+300+0'
pos_sys = '+0+0'
pos_torus = '+800+0'

info = nf.info(position=pos_info)
system = sys.system(info=info, position=pos_sys)
network = netw.network(info=info, position=pos_net, g_inh=0.0055)
traces = tra.traces(system, network, info=info, position=pos_tra)
torus = tor.torus(system, network, traces, info=info, position=pos_torus)

system.traces = traces
network.traces = traces

if pl.get_backend() == 'TkAgg':

	system.fig.tight_layout()
	traces.fig.tight_layout()
	torus.fig.tight_layout()

pl.show()
