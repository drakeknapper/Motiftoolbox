#!/usr/bin/env python

import system as sys
import traces as tra
import info as nf
import pylab as pl


pos_info = '+0+600'
pos_tra = '+300+600'
pos_sys = '+0+0'


i = nf.info(position=pos_info)
s = sys.system(info=i, position=pos_sys)
t = tra.traces(s, info=i, position=pos_tra)

if pl.get_backend() == 'TkAgg':
	s.fig.tight_layout()
	t.fig.tight_layout()

pl.show()
