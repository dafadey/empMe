#!/usr/bin/python3.9

import pyEmpMe_static
import math
import time
import numpy
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.fft import rfft, irfft

V=1.2

Lz=168000
Lx=80000
Nz=2048
Nx=1024

period=4500

tmax=120000
toff=110000

random.seed(13)

def srcTE(t,z) :
	return 0

def srcTM(t,z) :
	W=111000
	#return 1.5e-5*math.exp(-((z-Lz/2)/W*2)**8) * math.sin(t*math.pi*2/period)*(0.5+0.5*math.cos((t-tmax/2)/tmax*math.pi*2))
	if t < toff:
		return 1.3e-5*math.exp(-((z-Lz/2)/W*2)**8) * math.sin(t*math.pi*2/period)
	else:
		return 0

hh = pyEmpMe_static.hydro2dHandler(single=1, JHEAT=1) #JHEAT=1 : dT/dt~j*j / 0: dT/dt~E*E / 2: dT/dt~j*E

hh.sourceTE=srcTE
hh.sourceTM=srcTM

#hh.setup(velocity=V, Nx=Nx, Lx=Lx, Nz=Nz, Lz=Lz, mediaDepth=13000, mediaNu=0.0001, NUTratio=500000, srcTfactor=300000.0)

hh.setup(velocity=V, Nx=Nx, Lx=Lx, Nz=Nz, Lz=Lz, mediaDepth=13000, mediaNu=0.0001, mediaN0=0.00005, NUTratio=3e4, srcTfactor=100000.0, PMLx=10000, PMLstrength=0.1, diffusion=1)

#hh.setup(velocity=V, Nx=Nx, Lx=Lx, Nz=Nz, Lz=Lz, mediaDepth=13000, mediaN0=0.00005, mediaNu=0.0001, NUTratio=50000, srcTfactor=3000.0)

surf = numpy.zeros(Nz, dtype=numpy.float32)

#level=0
#tilt=0
#for i in range(0,Nz):
#	if i % 5 == 0:
#		tilt = random.randrange(-63,64)
#	level += tilt
#	surf[i] = level

for i in range(0,Nz):
	surf[i] = random.randrange(0,1300)

surf_spec = rfft(surf)
for i in range(0,33):
	surf_spec[i] *= 0
for i in range(177,int(Nz/2)):
	surf_spec[i] *= 0
surf = irfft(surf_spec)

hh.set_surface(surf)

#surf -= numpy.min(surf)
#surf /= numpy.max(surf)
#surf *= 1e-4
#hh.set_surface_Te(surf)

count=0
count_out=0
skipSteps=33

def getOutputCount():
	res=0
	t=0
	count=0
	dt=hh.get_dt()
	while t < tmax:
		if count % skipSteps == 0:
			res += 1
		count += 1
		t += dt
	return res

maxspec = 8*2*math.pi/period
print("max spec is "+str(maxspec))
specN = int(maxspec/(2*math.pi/Lz))
print("sepctral components is "+str(specN))

Te_2d_spec = numpy.zeros((getOutputCount(),specN), dtype=numpy.float32)
Te_2d = numpy.zeros((getOutputCount(),Nz), dtype=numpy.float32)

while hh.get_t() < tmax:
	hh.step()
	if count % skipSteps == 0:
		print("t="+str(hh.get_t()))
	
	if count % skipSteps == 0 :
		therm=hh.get_therm()
		therm_spec = rfft(therm)
		Te_2d[count_out]=therm
		therm_specW=numpy.zeros(specN)
		for i in range(0,specN):
			therm_specW[i] = (therm_spec[i].real**2+therm_spec[i].imag**2)**0.5
		therm_specW /= numpy.max(therm_specW)
		for i in range(0,specN):
			therm_specW[i] = therm_specW[i]**0.5
		Te_2d_spec[count_out] = therm_specW
		
		Ex=numpy.transpose(hh.get_2darray_Ex())
		Te=numpy.transpose(hh.get_2darray_Te())
		#absmax=max(numpy.max(Ex), abs(numpy.min(Ex)))
		#Ex[0,0]=absmax
		#Ex[1,0]=-absmax
		fig, (ax1, ax2) = plt.subplots(2,1,figsize=(11,7))
		fig.subplots_adjust(hspace=0)
		im = ax1.imshow(Ex[0:int(Nx*0.6)], extent=[0,Lz,Lx*0.6,0], aspect=1, cmap=plt.cm.seismic)
		imT = ax2.imshow(Te[int(Nx*0.45):int(Nx*0.55),int(Nz*0.3):int(Nz*0.7)], extent=[Lz*0.3,Lz*0.7,Lx*0.55,Lx*0.45], aspect=1, cmap=plt.cm.hot)
		rect = patches.Rectangle((Lz*0.3, Lx*0.44), Lz*0.4, Lx*0.1, edgecolor='black', facecolor='none')
		ax1.add_patch(rect)
		im.set_clim(-1e-5,1e-5)
		imT.set_clim(0,1e-4)
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.setp(ax1.get_yticklabels(), visible=False)
		plt.setp(ax2.get_xticklabels(), visible=False)
		plt.setp(ax2.get_yticklabels(), visible=False)
		ax1.set(title='$E_x$', xlabel='$z$', ylabel='$x$')
		ax2.set(title='$T_e$', xlabel='$z$', ylabel='$x$')
		plt.savefig('ex'+str(count).zfill(5)+'.png', bbox_inches='tight')
		plt.close(fig)
		count_out+=1
	count+=1

Ex2d=numpy.transpose(hh.get_2darray_Ex())

absmax=max(numpy.max(Ex2d), abs(numpy.min(Ex2d)))

def filter(v):
	f=math.pow(abs(v),0.5)
	if v > 0:
		return f
	else:
		return -f

#for i in range(0,Nz):
#	for j in range(0,Nx):
#		Ex2d[j,i]=filter(Ex2d[j,i])

Ex2d[0,0]=absmax
Ex2d[1,0]=-absmax

fig, ((ax_Te, ax_Te_spec), (ax_Ex, ax_Te_1d)) = plt.subplots(2,2, figsize=(33,33))

xaxis_for_Te=numpy.zeros(Nz)
for i in range(0,Nz):
	xaxis_for_Te[i]=i/Nz*Lz;

crossections=[40000,60000,80000,100000]
pens=['g','b','r','k']
for t_cross, pen in zip(crossections, pens):
	ax_Te_1d.plot(xaxis_for_Te[int(Nz*0.25):int(Nz*0.75)], hh.get_NUTratio()*Te_2d[int(getOutputCount()* t_cross/tmax), int(Nz*0.25):int(Nz*0.75)],pen)

ax_Te_1d.set_yscale('log')
ax_Te_1d.set_ylabel('$NUTratio T(z)$')
ax_Te_1d.set_xlabel('$z$')


for j in range(0,getOutputCount()):
	for i in range(0,Nz):
		Te_2d[j, i]=Te_2d[j, i]**(1/3);

im_Te=ax_Te.imshow(Te_2d,extent=[0,Lz,tmax,0], aspect=Lz/tmax, cmap=plt.cm.seismic)
ax_Te.set_title('$T_e(z,t)$')
ax_Te.set_xlabel('$z$')
ax_Te.set_ylabel('$t$')

for t_cross, pen in zip(crossections, pens):
	ax_Te.add_patch(patches.Polygon([[0, t_cross], [Lz, t_cross]], linewidth=1, edgecolor=pen))

im_Te_spec=ax_Te_spec.imshow(Te_2d_spec,extent=[0,maxspec,tmax,0], aspect=maxspec/tmax, cmap=plt.cm.seismic)
ax_Te_spec.set_title('$~{T_e}(k_z,t)$')
ax_Te_spec.set_xlabel('$k_z$')
ax_Te_spec.set_ylabel('$t$')

for t_cross, pen in zip(crossections, pens):
	ax_Te_spec.add_patch(patches.Polygon([[0, t_cross], [Lz, t_cross]], linewidth=1, edgecolor=pen))

im_Ex=ax_Ex.imshow(Ex2d,extent=[0,Lz,Lx,0], aspect=1, cmap=plt.cm.seismic)
ax_Ex.set_title('$E_x(z,x)$')
ax_Ex.set_xlabel('$z$')
ax_Ex.set_ylabel('$x$')

plt.savefig('therm.png', bbox_inches='tight')
plt.show()
plt.close(fig)
print('bye')
