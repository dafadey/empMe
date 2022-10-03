#!/usr/bin/python3.9

import pyEmpMe
import math
import time
import numpy

V=1.2

def srcTE(t,z) :
	#print("TE")
	#print(string(z))
	T=13000
	arg = z-V*t #-T*2
	return 0 #math.exp((arg/T)**2)*math.sin(arg/4000*math.pi)

def srcTM(t,z) :
	#print("TM "+str(z) +" "+ str(t))
	T=13000
	arg = z-V*t-168000+T*2
	return 5e-6*math.exp(-(arg/T)**2)*math.sin(2*arg/4000*math.pi)*math.exp(-((t-150000)/50000)**6)

hh = pyEmpMe.hydro2dHandler(single=1)

hh.sourceTE=srcTE
hh.sourceTM=srcTM

hh.setup(velocity=V, Nx=1024, Lx=80000, Nz=2048, Lz=168000, mediaDepth=3000)

count=0
while hh.get_t() < 300000:
	hh.step()
	if count % 100 == 0:
		print("t="+str(hh.get_t()))
	count+=1
	
print('bye')
