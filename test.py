#!/usr/bin/python3

import pyEmpMe
import math
import time
import numpy
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.fft import rfft, irfft

V=1.1

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
	return 5e-6*math.exp(-(arg/T)**2)*math.sin(2*arg/4000*math.pi)#*math.exp(-((t-150000)/70000)**6)

prefix=""
parser=argparse.ArgumentParser("empMe")
parser.add_argument("--nodraw", help="do not draw debug fields", action='store_true')
parser.add_argument("--noshow", help="skip showing plots after run", action='store_true')
parser.add_argument("--doframes", help="make frames for movie", action='store_true')
parser.add_argument("--frames-interval", help="interval between frames")
parser.add_argument("--prefix", help="prefix for filenames")
parser.add_argument("--amplitude", help="incident amplitude")
parser.add_argument("--aperture", help="incident aperture")
parser.add_argument("--tmax", help="incident pulse duration")
parser.add_argument("--nu0", help="nu0 * (1 + nuTratio * T)")
parser.add_argument("--nuTratio", help="nu0 * (1 + nuTratio * T)")
parser.add_argument("--n0", help="media density")
args=parser.parse_args()

Nx=1024
Lx=80000
Nz=2048
Lz=168000

tmax=1130000
DRAW=0 if args.nodraw else 1
SHOW=not(args.noshow)
FRAMES=args.doframes
skipSteps=33 if args.frames_interval == None else int(args.frames_interval)
prefix='' if args.prefix == None else args.prefix
tmax = tmax if args.tmax == None else float(args.tmax)

GPUlist=pyEmpMe.GPUtempsCoolFirst()
print(GPUlist)
#GPUlist[[temperature, devId]]
device1 = GPUlist[0][1]
device2 = GPUlist[1][1]

hh = pyEmpMe.hydro2dHandler(single=1, JHEAT="JJ", device=device1, DRAW=DRAW)
hh_lin = pyEmpMe.hydro2dHandler(single=1, JHEAT="JJ", linear=1, device=device2, DRAW=DRAW)

hh.sourceTE=None #srcTE
hh.sourceTM=None #srcTM

hh_lin.sourceTE=None #srcTE
hh_lin.sourceTM=None #srcTM

cell="tooth_sq.svg" #"asymmetric_sawtooth_3000.svg"
mediaDepth=13000#50
toothDepth=1300#0#100
toothWidth=40300#36723/2#26960#4300#13000
draw_interval=100
mediaNu=1e-5
cell_scaley=1/toothWidth*2200

srcNosc=11
period=4500
srcT=period*srcNosc

hh.setup(srcAmp=[5e-6], srcT=[srcT], srcNosc=[srcNosc], switchOnDelay=0, velocity=V, Nx=Nx, Lx=Lx, Nz=Nz, Lz=Lz, mediaNu=mediaNu, mediaDepth=mediaDepth, toothDepth=toothDepth, cellFilename=cell, cell_scaley=cell_scaley, toothWidth=toothWidth, draw_interval=draw_interval)

hh_lin.setup(srcAmp=[5e-6], srcT=[srcT], srcNosc=[srcNosc], switchOnDelay=0, velocity=V, Nx=Nx, Lx=Lx, Nz=Nz, Lz=Lz, mediaNu=mediaNu, mediaDepth=mediaDepth, toothDepth=toothDepth, cellFilename=cell, cell_scaley=cell_scaley, toothWidth=toothWidth, draw_interval=draw_interval)



hh.addToDrawList('ez')
hh.addToDrawList('ex')
hh.addToDrawList('by')
hh.addToDrawList('T')

count=0
frames_count=0
while hh.get_t() < tmax:
	t=hh.get_t()
	hh.stepAsync()
	hh_lin.stepAsync()
	hh.syncStep();
	hh_lin.syncStep();
	
	if count % draw_interval == 0:
		print("t="+str(hh.get_t()))
		if FRAMES:
			By=numpy.transpose(hh.get_2darray_By())
			By_lin=numpy.transpose(hh_lin.get_2darray_By())
			By_gen = By - By_lin
			Ex=numpy.transpose(hh.get_2darray_Ex())
			Ez=numpy.transpose(hh.get_2darray_Ez())
			Te=numpy.transpose(hh.get_2darray_Te())
			#Jz=numpy.transpose(hh.get_2darray_Jz())
			#Jx=numpy.transpose(hh.get_2darray_Jx())
			#absmax=max(numpy.max(Ex), abs(numpy.min(Ex)))
			#Ex[0,0]=absmax
			#Ex[1,0]=-absmax
			fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(11,15))
			fig.subplots_adjust(hspace=0)
			Nz0=int(Nz*0.0)
			Nz1=int(Nz*1.0)
			imEx = ax1.imshow(Ex[0:int(Nx*0.6),Nz0:Nz1], extent=[Lz*Nz0/Nz,Lz*Nz1/Nz,Lx*0.6,0], aspect=1, cmap=plt.cm.seismic)
			imEz = ax2.imshow(Ez[0:int(Nx*0.6),Nz0:Nz1], extent=[Lz*Nz0/Nz,Lz*Nz1/Nz,Lx*0.6,0], aspect=1, cmap=plt.cm.seismic)
			imBygen = ax3.imshow(By_gen[0:int(Nx*0.6),Nz0:Nz1], extent=[Lz*Nz0/Nz,Lz*Nz1/Nz,Lx*0.6,0], aspect=1, cmap=plt.cm.seismic)
			Nz0=int(Nz*0.0)
			Nz1=int(Nz*1.0)
			imT = ax4.imshow(Te[int(Nx*0.35):int(Nx*0.65),Nz0:Nz1], extent=[Lz*Nz0/Nz,Lz*Nz1/Nz,Lx*0.65,Lx*0.35], aspect=1, cmap=plt.cm.hot)
			
			
			
			imEx.set_clim(-3e-6,3e-6)
			imEz.set_clim(-3e-6,3e-6)
			
			maxBygen = max(numpy.amax(By_gen[int(Nx*0.45):int(Nx*0.55),int(Nz*0.1):int(Nz*0.9)]), abs(numpy.amin(By_gen[int(Nx*0.45):int(Nx*0.55),int(Nz*0.1):int(Nz*0.9)])))
			print("maxBygen="+str(maxBygen))
			imBygen.set_clim(-4e-10, 4e-10)
			#imT.set_clim(0,1.5e-3)
			#Jabsmax = max(abs(numpy.max(Jz)),abs(numpy.min(Jz)))
			#imJ.set_clim(-Jabsmax, Jabsmax)
			plt.setp(ax1.get_xticklabels(), visible=False)
			plt.setp(ax1.get_yticklabels(), visible=False)
			plt.setp(ax2.get_xticklabels(), visible=False)
			plt.setp(ax2.get_yticklabels(), visible=False)
			plt.setp(ax3.get_xticklabels(), visible=False)
			plt.setp(ax3.get_yticklabels(), visible=False)
			plt.setp(ax4.get_xticklabels(), visible=False)
			plt.setp(ax4.get_yticklabels(), visible=False)
			ax1.set(title='$E_x$', xlabel='$z$', ylabel='$x$')
			ax2.set(title='$E_z$', xlabel='$z$', ylabel='$x$')
			ax3.set(title='${B_y}_{gen}$', xlabel='$z$', ylabel='$x$')
			ax4.set(title='$T_e$', xlabel='$z$', ylabel='$x$')
			plt.suptitle("x"+str(skipSteps)+" t="+str(int(t)))
			plt.savefig(prefix+'ex'+str(frames_count).zfill(5)+'.png', bbox_inches='tight')
			plt.close(fig)
			frames_count+=1


	count+=1
	
print('bye')
