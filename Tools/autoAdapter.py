#! /usr/bin/python

snippet_0 = """
#include "auto_f2c.h"
#include "math.h"
int func(integer ndim, const doublereal *u, const integer *icp, 
	const doublereal *par, integer ijac, doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
"""

snippet_1 = """
	return 0;
}
int stpnt(integer ndim, doublereal t, doublereal *u, doublereal *par)
{
"""

snippet_2 = """
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

snippets = [snippet_0, snippet_1, snippet_2]

def createAutoCode(filename, diffEqs, params):

	f = open(filename, 'w+')

	f.write(snippets[0])
	f.write(diffEqs)
	f.write(snippets[1])
	f.write(params)
	f.write(snippets[2])

	f.close()



def writeConstantsFile(filename, NDIM, NTST):

	f = open(filename, 'w+')

	f.write("""
	dat='orbit'
	NDIM=%i, IPS=2, IRS=0, ILP=0
	ICP=[10, 'PERIOD']
	NTST=%i, NCOL=4, IAD=3, ISP=2, ISW=1, IPLT=4, NBC=0, NINT=0
	NMX=50, NPR=50, MXBF=10, IID=2, ITMX=8, ITNW=7, NWTN=3, JAC=0
	EPSL=1e-07, EPSU=1e-07, EPSS =0.0001
	DS =0.1, DSMIN=0.00001, DSMAX=25.0, IADS=1
	NPAR=5, THL={'PERIOD': 0.0}, THU={}""" % (NDIM, NTST))
	
	f.close()

###











if __name__ == "__main__":

	createAutoCode('test.c', 'hello;', 'hello;')
