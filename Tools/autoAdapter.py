#! /usr/bin/python

orbit_0 = """
#include "auto_f2c.h"
#include "math.h"
int func(integer ndim, const doublereal *u, const integer *icp, 
	const doublereal *par, integer ijac, doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
"""

orbit_1 = """
	return 0;
}
int stpnt(integer ndim, doublereal t, doublereal *u, doublereal *par)
{
"""

orbit_2 = """
	return 0;
}
int pvls(integer ndim, const doublereal *u, doublereal *par) { return 0; } 
int bcnd(integer ndim, const doublereal *par, const integer *icp, integer nbc, const doublereal *u0, const doublereal *u1, integer ijac, doublereal *fb, doublereal *dbc)
{ return 0; } 
int icnd(integer ndim, const doublereal *par, const integer *icp, integer nint, const doublereal *u, const doublereal *uold,
		  const doublereal *udot, const doublereal *upold, integer ijac, doublereal *fi, doublereal *dint)
{ return 0; } 
int fopt(integer ndim, const doublereal *u, const integer *icp, const doublereal *par, integer ijac, doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{ return 0; }
"""


adjointICND = """
	return 0;
}
int pvls(integer ndim, const doublereal *u, doublereal *par) { return 0; } 
int bcnd(integer ndim, const doublereal *par, const integer *icp, integer nbc, const doublereal *u0, const doublereal *u1, integer ijac, doublereal *fb, doublereal *dbc)
{ return 0; } 
int icnd(integer ndim, const doublereal *par, const integer *icp, integer nint, const doublereal *u, const doublereal *uold,
		  const doublereal *udot, const doublereal *upold, integer ijac, doublereal *fi, doublereal *dint)
{ 
	integer i, ndim2;
	doublereal ui;

	ndim2 = ndim/2;

	for(i=0; i<ndim2; i++)
	{
		ui = u[i+ndim2];
		fi[i] = ui*ui-1.;
	}
	return 0;
} 
int fopt(integer ndim, const doublereal *u, const integer *icp, const doublereal *par, integer ijac, doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{ return 0; }
"""

orbitSnippets = [orbit_0, orbit_1, orbit_2]

def createAutoCode(filename, diffEqs, params, mode='orbit'):

	f = open(filename, 'w+')

	f.write(orbitSnippets[0])
	f.write(diffEqs)
	f.write(orbitSnippets[1])
	f.write(params)

	if mode == 'orbit':
		f.write(orbitSnippets[2])
	
	elif mode == 'adjoint':
		f.write(adjointICND)

	f.close()



def writeConstantsFile(filename, NDIM, NTST, mode='orbit'):

	if mode == 'orbit':
		constants = """
		dat='orbit'
		NDIM=%i, IPS=2, IRS=0, ILP=0
		ICP=[10, 'PERIOD']
		NTST=%i, NCOL=4, IAD=3, ISP=2, ISW=1, IPLT=4, NBC=0, NINT=0
		NMX=50, NPR=50, MXBF=10, IID=2, ITMX=8, ITNW=7, NWTN=3, JAC=0
		EPSL=1e-07, EPSU=1e-07, EPSS =0.0001
		DS =0.1, DSMIN=0.00001, DSMAX=25.0, IADS=1
		NPAR=5, THL={'PERIOD': 0.0}, THU={}""" % (NDIM, NTST)
	
	elif mode == 'adjoint':
		constants = """
		dat='orbit'
		NDIM=%i, IPS=2, IRS=0, ILP=0
		ICP=[10, 'PERIOD']
		NTST=%i, NCOL=4, IAD=3, ISP=2, ISW=1, IPLT=4, NBC=0, NINT=%i
		NMX=50, NPR=50, MXBF=10, IID=2, ITMX=8, ITNW=7, NWTN=3, JAC=0
		EPSL=1e-07, EPSU=1e-07, EPSS =0.0001
		DS =0.1, DSMIN=0.00001, DSMAX=25.0, IADS=1
		NPAR=5, THL={'PERIOD': 0.0}, THU={}""" % (2*NDIM, NTST, NDIM) # 2NDIM equations, NDIM integral conditions (for each linear adjoint variable)

	f = open(filename, 'w+')

	f.write(constants)
	
	f.close()

###











if __name__ == "__main__":

	createAutoCode('test.c', 'hello;', 'hello;')
