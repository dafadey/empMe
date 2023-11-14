#!/usr/bin/python3

#how to make a movie
#./test_static.py --amplitude 2e-5 --frames-interval 11 --doframes

import pyEmpMe_static
import math
import time
import numpy
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.fft import rfft, irfft
import argparse

V=1.2

zKoeff=2

Lz=168000*zKoeff
Lx=80000
Nz=int(2048*zKoeff)
Nx=1024
PMLx=10000

mediaNu=0.0001*2
mediaN0=0.00005
NUTratio=3e4/2.0

period=4500

tmax=400000


#tmax=300000
#toff=tmax-31000

amp=1.3e-5
apt=131000

parser=argparse.ArgumentParser("empMe")
parser.add_argument("--dev", help="device to use")
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


DRAW=0 if args.nodraw else 1
SHOW=not(args.noshow)
FRAMES=args.doframes
skipSteps=33 if args.frames_interval == None else int(args.frames_interval)
prefix='' if args.prefix == None else args.prefix
amp = amp if args.amplitude == None else float(args.amplitude)
apt = apt if args.aperture == None else float(args.aperture)
mediaNu = mediaNu if args.nu0 == None else float(args.nu0)
mediaN0 = mediaN0 if args.n0 == None else float(args.n0)
NUTratio = NUTratio if args.nuTratio == None else float(args.nuTratio)
tmax = tmax if args.tmax == None else float(args.tmax)

toff=tmax-130000

if args.dev == None:
	GPUlist=pyEmpMe_static.GPUtempsCoolFirst()
	#GPUlist[[temperature, devId]]
	print(GPUlist)
	device = GPUlist[0][1]
else:
	device = int(args.dev)

skipSteps0=skipSteps;

print("FRAMES="+str(FRAMES))
print("SHOW="+str(SHOW))
print("DRAW="+str(DRAW))
print("amplitude="+str(amp))

def srcTE(t,z) :
	return 0

def srcTM(t,z) :
	W=apt
	if t < toff:
		return amp*math.exp(-((z-Lz/2)/W*2)**16) * math.sin(t*math.pi*2/period)
	else:
		return 0

#hh = pyEmpMe_static.hydro2dHandler(single=1, JHEAT="JJ", device=device, DRAW=DRAW) #JHEAT=1 : dT/dt~j*j / 0: dT/dt~E*E / 2: dT/dt~j*E

hh = pyEmpMe_static.hydro2dHandler(single=1, JHEAT="JE", device=device, DRAW=DRAW) #JHEAT=1 : dT/dt~j*j / 0: dT/dt~E*E / 2: dT/dt~j*E

#JHEAT=0 EE
#JHEAT=1 JJ
#JHEAT=2 JE

hh.sourceTE=srcTE
hh.sourceTM=srcTM

#hh.setup(velocity=V, Nx=Nx, Lx=Lx, Nz=Nz, Lz=Lz, mediaDepth=13000, mediaNu=0.0001, mediaN0=0.00005, NUTratio=3e4, srcTfactor=100000.0, PMLx=10000, PMLstrength=0.1, diffusion=5)

#for nu=nu_0*(1+NUTratio*T)
hh.setup(velocity=V, Nx=Nx, Lx=Lx, Nz=Nz, Lz=Lz, mediaDepth=13000, mediaNu=mediaNu, mediaN0=mediaN0, NUTratio=NUTratio, srcTfactor=10000.0*2.0, PMLx=PMLx, PMLstrength=0.1, diffusion=5)

#for nu=nu_0*(1+NUTratio*sqrt(T))
# dT=src(E,j)*SRCTfactor
#hh.setup(velocity=V, Nx=Nx, Lx=Lx, Nz=Nz, Lz=Lz, mediaDepth=13000, mediaNu=0.0001, mediaN0=0.00005, NUTratio=1.5e2, srcTfactor=10000.0, PMLx=PMLx, PMLstrength=0.1, diffusion=5)

hh.set_dt(hh.get_dt()*0.5);

surf = numpy.zeros(Nz, dtype=numpy.float32)
surf_spec = rfft(surf)

for i in range(int(17 * zKoeff), int(277 * zKoeff)):
	arg=random.randrange(0,1e4)/1e4*2*math.pi
	surf_spec[i] = math.sqrt(zKoeff)*3000 * (math.cos(arg)+1j*math.sin(arg))

surf = irfft(surf_spec)


dp=4414.7551
dopt=4500
#for i in range(0,13):
#	surf[int(Lz/2/Lz*Nz)-i]=-233
#	surf[int((Lz/2+13*dp+0.5*dp)/Lz*Nz)+i]=-233


#hh.set_surface(surf)


#for i in range(0,Nz):
#	surf[i] = math.cos(i/Nz*Lz/(dp*0.7)*2*math.pi)+math.cos(i/Nz*Lz/(dp)*2*math.pi)



#for i in range(int(Nz/2-Nz/5),int(Nz/2+Nz/5)):
#	surf[i] = 0

surf -= numpy.min(surf)
surf /= numpy.max(surf)
surf *= 0.3e-5*7

#surf *= 500
#hh.set_surface(surf)

#for i in range(0,Nz):
#	z=i*Lz/Nz
#	W=31000*zKoeff
#	surf[i] *= math.exp(-((z-Lz/2)/W*2)**16)

hh.set_surface_Te(surf)

count=0
count_out=0

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
dk = 2 * math.pi/Lz
spec0 = 2 * math.pi/period * 0.5
spec1 = 2 * math.pi/period * 2
specN = int(maxspec/dk)
print("sepctral components is "+str(specN))

output_count=getOutputCount();

frames_count=0

Te_2d_spec = numpy.zeros((output_count, specN), dtype = numpy.float32)
Te_2d_max_char=numpy.zeros((output_count, 6), dtype = numpy.float32) # t, value, pos, pos0, pos1, By_above_max
#By_above_max is a maximum value abuve plasma surface in a lyayer with thikness = period
Te_2d = numpy.zeros((output_count,Nz), dtype=numpy.float32)
#ez ex by Te
hh.addToDrawList('ez')
hh.addToDrawList('ex')
hh.addToDrawList('by')
hh.addToDrawList('T')

def get_max_char(spec) :
	max_spec = 0
	max_spec_pos = 0
	for ki in range(int(spec0/dk), int(spec1/dk)) :
		k = ki * dk
		v = abs(spec[ki])
		if max_spec < v :
		  max_spec = v
		  max_spec_pos = k

	max_spec_pos0 = 0
	max_spec_pos1 = 0
	for ki in range(int(spec0/dk), int(spec1/dk)) :
		k = ki * dk
		v = abs(spec[ki])
		if v > max_spec * .5 and k < max_spec_pos and max_spec_pos0 == 0:
			max_spec_pos0 = k
		if v < max_spec * .5 and k > max_spec_pos and max_spec_pos1 == 0:
			max_spec_pos1 = k
	
	return [max_spec, max_spec_pos, max_spec_pos0, max_spec_pos1]

mat_mask=numpy.transpose(hh.get_2darray_matMask())
first_pos_mat_exists=0
for xi in range(0, mat_mask.shape[0]):
	if max(mat_mask[xi]) > 0:
		first_pos_mat_exists = xi
		break

print("first_pos_mat_exists="+str(first_pos_mat_exists))

while hh.get_t() < tmax:
	hh.step()
	if count % skipSteps == 0:
		print("t="+str(hh.get_t()))

	if count % skipSteps == 0 :
		therm=hh.get_therm()
		therm_spec = rfft(therm)
		
		Te_2d[count_out]=therm
		#therm_specW=numpy.zeros(specN)
		#for i in range(0,specN):
		#	therm_specW[i] = (therm_spec[i].real**2+therm_spec[i].imag**2)**0.5
		#therm_specW /= numpy.max(therm_specW)
		#for i in range(0,specN):
		#	therm_specW[i] = therm_specW[i]**0.5
		#Te_2d_spec[count_out] = therm_specW
		Te_2d_spec[count_out]=abs(therm_spec[0:specN])
		
		maxNu=hh.getMaxNu()
		
		print("maxNu="+str(maxNu)+", dt * maxNu="+str(maxNu*hh.get_dt()))

		t=hh.get_t()

		By2d=numpy.transpose(hh.get_2darray_By())

		By_above = By2d[first_pos_mat_exists-3-int(period/Lx*Nx) : first_pos_mat_exists-3]
		#print("By shape is "+str(By_above.shape))
		#print("By above max is "+str(numpy.amax(By_above)))

		Te_2d_max_char[count_out] = [t] + get_max_char(therm_spec) + [numpy.amax(By_above)] # list concatenation
				
		if t > 243000 and t < 267000 and FRAMES :
			skipSteps=3
		else:
			skipSteps=skipSteps0
		
		if FRAMES:
			Ex=numpy.transpose(hh.get_2darray_Ex())
			#Ez=numpy.transpose(hh.get_2darray_Ez())
			Te=numpy.transpose(hh.get_2darray_Te())
			#Jz=numpy.transpose(hh.get_2darray_Jz())
			#Jx=numpy.transpose(hh.get_2darray_Jx())
			#absmax=max(numpy.max(Ex), abs(numpy.min(Ex)))
			#Ex[0,0]=absmax
			#Ex[1,0]=-absmax
			fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(11,11))
			fig.subplots_adjust(hspace=0)
			Nz0=int(Nz*0.2)
			Nz1=int(Nz*0.8)
			im = ax1.imshow(Ex[0:int(Nx*0.6),Nz0:Nz1], extent=[Lz*Nz0/Nz,Lz*Nz1/Nz,Lx*0.6,0], aspect=1, cmap=plt.cm.seismic)
			Nz0=int(Nz*0.3)
			Nz1=int(Nz*0.7)
			imT = ax2.imshow(Te[int(Nx*0.45):int(Nx*0.55),Nz0:Nz1], extent=[Lz*Nz0/Nz,Lz*Nz1/Nz,Lx*0.55,Lx*0.45], aspect=1, cmap=plt.cm.hot)
			#imJ=ax3.imshow(Jz[int(Nx*0.45):int(Nx*0.55),int(Nz*0.3):int(Nz*0.7)], extent=[Lz*0.3,Lz*0.7,Lx*0.55,Lx*0.45], aspect=1, cmap=plt.cm.seismic)
			
			#jz=numpy.zeros(Nz)
			te=numpy.zeros(Nz)
			for i in range(0,Nz):
				for j in range(int(Nx*0.475),int(Nx*0.525)):
					#jz[i] += Ez[j,i]*Jz[j,i]
					if Te[j,i] > 0:
						te[i] += Te[j,i]
						break
			
			#jz *= 1/max(abs(numpy.max(jz)),abs(numpy.min(jz)))*numpy.max(te)
			
			#jz *= 1e10#1e11
			
			#ax3.plot(jz,'b')
			ax3.plot(te,'r')
			#ax3.ticklabel_format(axis='y',style='sci')
			ax3.set_yticklabels([])
			ax3.set_ylabel("$T_e$, max="+str('{:0.3e}'.format(max(te))))
			rect = patches.Rectangle((Lz*Nz0/Nz, Lx*0.44), Lz*(Nz1-Nz0)/Nz, Lx*0.1, edgecolor='black', facecolor='none')
			ax1.add_patch(rect)
			im.set_clim(-3e-6,3e-6)
			#imT.set_clim(0,1.5e-3)
			#Jabsmax = max(abs(numpy.max(Jz)),abs(numpy.min(Jz)))
			#imJ.set_clim(-Jabsmax, Jabsmax)
			plt.setp(ax1.get_xticklabels(), visible=False)
			plt.setp(ax1.get_yticklabels(), visible=False)
			plt.setp(ax2.get_xticklabels(), visible=False)
			plt.setp(ax2.get_yticklabels(), visible=False)
			ax1.set(title='$E_x$', xlabel='$z$', ylabel='$x$')
			ax2.set(title='$T_e$', xlabel='$z$', ylabel='$x$')
			plt.suptitle("x"+str(skipSteps)+" t="+str(int(t)))
			plt.savefig(prefix+'ex'+str(frames_count).zfill(5)+'.png', bbox_inches='tight')
			plt.close(fig)
			frames_count+=1
			
		
		count_out=count_out + 1 if count_out + 1 < output_count else count_out
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

fig, ((ax_Te, ax_Te_spec), (ax_Ex, ax_Te_1d)) = plt.subplots(2,2, figsize=(22,17))

xaxis_for_Te=numpy.zeros(Nz)
for i in range(0,Nz):
	xaxis_for_Te[i]=i/Nz*Lz;

crossections=[40000,60000,80000,263000,toff+Lx/2-PMLx,tmax-1000]
dump_crossection_id=4
pens=['g','b','r','k','orange','purple']
mainfilenamepart=str(amp)+('' if args.tmax == None else '_tmax='+str(tmax))

fileout = open(prefix+'thermspec_amp'+mainfilenamepart+'.dat',"w")
crossection_id=0
for t_cross, pen in zip(crossections, pens):
	ax_Te_1d.plot(xaxis_for_Te[int(Nz*0.25):int(Nz*0.75)], hh.get_NUTratio()*Te_2d[int(getOutputCount()* t_cross/tmax), int(Nz*0.25):int(Nz*0.75)],pen)
	fileout.write(str(t_cross)+'\n')
	nk=numpy.size(Te_2d_spec[int(getOutputCount()* t_cross/tmax)])
	fileout.write(str(maxspec/nk)+'\n')
	fileout.write(str(nk)+'\n')
	for i in range(0,nk):
		fileout.write(str(Te_2d_spec[int(output_count * t_cross / tmax), i])+'\n')
	if crossection_id == dump_crossection_id:
		cross_file=open(prefix+'therm_at_'+str(crossections[crossection_id])+'_amp'+mainfilenamepart+'.dat', "w")
		cross_file.write(str(Lz)+'\n')
		cross_file.write(str(Nz)+'\n')
		for i in range(0,Nz):
			cross_file.write(str(Te_2d[int(output_count * t_cross / tmax), i])+'\n')
		cross_file.close()
	crossection_id += 1
fileout.close()

ax_Te_1d.set_yscale('log')
ax_Te_1d.set_ylabel('$NUTratio T(z)$')
ax_Te_1d.set_xlabel('$z$')


for j in range(0,getOutputCount()):
	Te_2d[j]=Te_2d[j]**(1/3)

im_Te=ax_Te.imshow(Te_2d,extent=[0,Lz,tmax,0], aspect=Lz/tmax, cmap=plt.cm.seismic)
ax_Te.set_title('$T_e(z,t)$')
ax_Te.set_xlabel('$z$')
ax_Te.set_ylabel('$t$')

for t_cross, pen in zip(crossections, pens):
	ax_Te.add_patch(patches.Polygon([[0, t_cross], [Lz, t_cross]], linewidth=1, edgecolor=pen))

for j in range(0,getOutputCount()):
	Te_2d_max = numpy.max(Te_2d_spec[j])
	if Te_2d_max == 0 :
		Te_2d_spec[j] = 0
	else:
		Te_2d_spec[j] /= Te_2d_max 
	Te_2d_spec[j] = Te_2d_spec[j]**0.5

im_Te_spec=ax_Te_spec.imshow(Te_2d_spec,extent=[0,maxspec,tmax,0], aspect=maxspec/tmax, cmap=plt.cm.seismic)
ax_Te_spec.set_title('$~{T_e}(k_z,t)$')
ax_Te_spec.set_xlabel('$k_z$')
ax_Te_spec.set_ylabel('$t$')

for t_cross, pen in zip(crossections, pens):
	ax_Te_spec.add_patch(patches.Polygon([[0, t_cross], [maxspec, t_cross]], linewidth=1, edgecolor=pen))

#ax_Te_spec.plot(Te_2d_max_char[:,2], Te_2d_max_char[:,0], color='white', linestyle='dotted')
#ax_Te_spec.plot(Te_2d_max_char[:,3], Te_2d_max_char[:,0], color='white', linestyle='dashed')
#ax_Te_spec.plot(Te_2d_max_char[:,4], Te_2d_max_char[:,0], color='white', linestyle='dashed')

im_Ex=ax_Ex.imshow(Ex2d,extent=[0,Lz,Lx,0], aspect=1, cmap=plt.cm.seismic)
ax_Ex.set_title('$E_x(z,x)$')
ax_Ex.set_xlabel('$z$')
ax_Ex.set_ylabel('$x$')

filename = prefix+'therm_amp'+mainfilenamepart

maxfile = open(filename+".maxchar.txt", 'w')
for chars in Te_2d_max_char:
	for char in chars:
		maxfile.write(str(char)+' ')
	maxfile.write('\n')
maxfile.close();

description = open(filename+'.desc.txt', 'w')
description.write(hh.getDescription())
description.close()

plt.savefig(filename+'.png', bbox_inches='tight')
if SHOW :
	plt.show()
plt.close(fig)
print('bye')
