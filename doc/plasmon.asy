import grph;
import fft;
import patterns;

print("\textbf{Optical plasmon excitation with grating}");
print("Consider maxwell equations in vacuum$^\mathrm v$ and plasma$^\mathrm p$:",20pt);
print("$\displaystyle \frac {\partial \mathbf H^\mathrm v} {\partial t} = -\mathrm{rot} \mathbf E^\mathrm v$; \ \ \ $\displaystyle \frac {\partial \mathbf E^\mathrm v} {\partial t} = \mathrm{rot} \mathbf H^\mathrm v $",30pt);
print("$\displaystyle \frac {\partial \mathbf H^\mathrm p} {\partial t} = -\mathrm{rot} \mathbf E^\mathrm p$; \ \ \ $\displaystyle \frac {\partial \mathbf E^\mathrm p} {\partial t} = \mathrm{rot} \mathbf H^\mathrm p - \mathbf j $; \ \ \ $\displaystyle \frac{\partial \mathbf j}{\partial t}=\omega_p^2 \mathbf E^\mathrm p$");

{
	add("hatch",hatch(2.1mm));
	real zero = cr-2cm;
	filldraw(shift(1mm,zero)*scale(7cm,2cm)*shift(0,-1)*unitsquare,pattern("hatch"));
	draw((2.5cm,zero+1mm)--(4.5cm,zero+1mm),linewidth(0.7),Arrow);
	label("$k_z$",(3.5cm,zero+1mm),N);
	label("$\mathbf E^\mathrm v$, $\mathbf H^\mathrm v$",(3.5cm,zero+1cm),N);
	label("$\mathbf E^\mathrm p$, $\mathbf H^\mathrm p$",(3.5cm,zero-1cm),S, UnFill);
	draw((0.3cm,zero+0.3cm)--(1.7cm,zero+0.3cm), Arrow);
	label("$z$",(1.5cm,zero+0.3cm),N);
	draw((0.3cm,zero+0.3cm)--(0.3cm,zero+1.7cm), Arrow);
	label("$x$",(0.3cm,zero+1.5cm),E);
	filldraw(shift(0.3cm,zero+0.3cm)*scale(2mm)*unitcircle,white);
	fill(shift(0.3cm,zero+0.3cm)*scale(0.3mm)*unitcircle);
	label("$y$",(0.38cm,zero+0.38cm),NE);

	draw((5cm,zero+0.3cm)--(5cm,zero+1.7cm),linewidth(0.7),Arrow);
	label("$k_\mathrm v$",(5cm,zero+(0.3cm+1.7cm)*0.5),E);
	draw((5cm,zero-0.3cm)--(5cm,zero-1.7cm),linewidth(0.7),Arrow);
	label("$k_\mathrm p$",(5cm,zero-(0.3cm+1.7cm)*0.5),E, UnFill);
}
print("",4.5cm);
print("Now let us go to Fourier space, all vacuum fields assumed to be $\mathrm e^{\mathrm i \omega t - i k_z z - k_\mathrm v x}$,");
print("while plasma fields will have form: $\mathrm e^{\mathrm i \omega t - i k_z z + k_\mathrm p x}$, here $k_\mathrm v > 0$, $k_\mathrm p > 0$", 30pt);
print("Ok, $\mathrm{rot \mathbf f}$ is $\displaystyle \left\{\frac{\partial f_z}{\partial y}-\frac{\partial f_y}{\partial z}; \ \frac{\partial f_x}{\partial z}-\frac{\partial f_z}{\partial x}; \ \frac{\partial f_y}{\partial x}-\frac{\partial f_x}{\partial y}\right\}$",30pt);
print("We will do \textit p-polarized solution, so only $E^\mathrm{v,p}_{x,z}$, $H_y^\mathrm{v,p}$ exists");
print("So our initial equations are transformed to");
print("$\mathrm i \omega H^\mathrm v_y = -k_\mathrm v E^\mathrm v_z + \mathrm i k_z E^\mathrm v_x$; \ \ \ $\color{blue}\mathrm i \omega E^\mathrm v_z = -k_\mathrm v H^\mathrm v_y\color{black}$; \ \ \ $\mathrm i \omega E^\mathrm v_x = \mathrm i k_z H_y^\mathrm v$",25pt);
print("$\mathrm i \omega H^\mathrm p_y = +k_\mathrm p E^\mathrm p_z + \mathrm i k_z E^\mathrm p_x$; \ \ \ $\color{blue}\displaystyle \mathrm i \omega E^\mathrm p_z = +k_\mathrm p H^\mathrm p_y - \frac{E^\mathrm p_z \omega^2_p}{\mathrm i \omega}\color{black}$; \ \ \ $\displaystyle \mathrm i \omega E^\mathrm p_x = \mathrm i k_z H_y^\mathrm p - \frac{E^\mathrm p_x \omega^2_p}{\mathrm i \omega}$",25pt);
print("Ok we perfectly know that $\omega^2 = k_\mathrm z^2-k_\mathrm v^2$, $\omega^2 = \omega_p^2 + k_\mathrm z^2-k_\mathrm p^2$, thus we will not use all the equations");
print("As usually, since the surface is at $x=0$ the following is true: $E^\mathrm v_z = E^\mathrm p_z = E_z$, $H^\mathrm v_y = H^\mathrm p_y = H_y$",27pt);
print("So from the \color{blue}second \color{black} equations on the lines above it follows: \ $\displaystyle \frac{1}{1-\frac{\omega_p^2}{\omega^2}}=-\frac{k_\mathrm v}{k_\mathrm p}$, \ \ ok great signs are good ($\omega < \omega_p$)",25pt);
print("So... $\omega^4\omega_p^2+\omega^4 k_z^2-\omega^6=(\omega_p^4-2\omega^2\omega_p^2+\omega^4)(k_z^2-\omega^2)$; and then we obtain",16pt);
print("$\omega^4\omega_p^2+\omega^4 k_z^2-\omega^6 = \omega_p^4 k_z^2- 2 \omega^2 \omega_p^2 k_z^2 + \omega^4 k_z^2 - \omega_p^4 \omega^2 + 2 \omega_p^2 \omega^4 - \omega^6$",16pt);
print("$0 = \omega_p^4 k_z^2- 2 \omega^2 \omega_p^2 k_z^2 - \omega_p^4 \omega^2 + \omega_p^2 \omega^4$",16pt);
print("So for $\omega^2$ we get quadratic equation: $\omega_p^2 k_z^2- \omega^2 (2 k_z^2 + \omega_p^2) + \omega^4$=0",27pt);
print("Solution is: $\displaystyle \omega^2=\frac{2 k_z^2+\omega_p^2\pm\sqrt{4 k_z^4+\omega_p^4}}{2}$",25pt);
print("Of course we need $-$ to get low frequency optical plasmon: $\color {red}\displaystyle \omega^2=\frac{2 k_z^2+\omega_p^2-\sqrt{4 k_z^4+\omega_p^4}}{2}\color{black}$",33pt);
print("For low frequency we can rewrite it as: $\displaystyle \omega^2\approx k_z^2+\frac{\omega_p^2}{2}-\frac{\omega_p^2}{2}\left[1+\frac{4k_z^4}{\omega_p^4}\right] \approx k_z^2\left[1-\frac{2k_z^2}{\omega_p^2}\right]$",25pt);
print("So $\color {blue} {\displaystyle \omega\approx k_z \left(1-\frac{k_z^2}{\omega_p^2}\right)}$",25pt);

real n0=0.001;
real wp = sqrt(n0);
real T=36000/4;
real w0 = 2*pi/T;

real spp_disp_precise(real kz)
{
	return sqrt(0.5*(2*kz^2+wp^2-sqrt(4*kz^4+wp^4)));
}

real spp_disp(real kz)
{
	return kz*(1-kz^2/wp^2);
}

real V=1.2;

real cherenkov_disp(real kz, real dk)
{
	//cos a = 1/V => k = kz / cos a
	return (kz-dk)*V;
}

plot pdisp;
pdisp.labelx="$k_z$";
pdisp.labely="$\omega$";

path wspp_p, wspp, wch, wopt;
wopt = (0,w0);
wopt = wopt--(w0*1.5,w0);
for(real kz=0; kz<w0*1.5; kz=kz+w0/100)
{
	wspp_p = wspp_p--(kz, spp_disp_precise(kz));
	wspp = wspp--(kz, spp_disp(kz));
	wch = wch--(kz, cherenkov_disp(kz,0));
}

pair spp_res;
{
	real[] res = intersect(wspp_p,wopt);
	spp_res=point(wopt,res[1]);
}

pair ch_res;
{
	real[] res = intersect(wch,wopt);
	ch_res=point(wopt,res[1]);
}
real k1=spp_res.x;
real k2=ch_res.x;
write("k1-k2=",k1-k2);

pdisp.addplot(wspp_p,red);
pdisp.addplot(wspp,blue);
pdisp.addplot(wch,0.5*green);
pdisp.addplot(wopt);
pdisp.addlabel("$k_1$",spp_res,S,blue);
pdisp.addlabel("$k_2$",ch_res,N,0.5*green);

print("We have plasma density $n_0="+string(n0)+"$, so $\omega_p="+string(wp,3)+"$, optical period is $T="+string(T,3)+"$,");
print("so $\omega_{opt}=2\pi/T="+string(w0,3)+"$, ratio $\omega_p/\omega$ is $"+string(wp/w0,3)+"$");

pdisp.show((7cm,4cm));

print("Red is for precise plasmon dispersion, blue is for low frequency approximation,");
print("$k_1-k_2="+string(k1-k2,3)+"$, so grating period is $l_g=2\pi/(k_1-k_2)="+string(2*pi/(k1-k2),3)+"$");
print("");
print("Let us plot plasmon fields:");
real kv = sqrt(k1^2 - w0^2);
real kp = sqrt(wp^2 - w0^2 + k1^2);
real H_y(real H0, real z, real x)
{
	if(x>0)
		return (H0*exp(-(0,1)*k1*z - kv * x)).x;
	else
		return (H0*exp(-(0,1)*k1*z + kp * x)).x;
}
real E_z(real H0, real z, real x)
{
	if(x>0)
		return (-H0*kv/(0,1)/w0*exp(-(0,1)*k1*z - kv * x)).x;
	else
		return (H0*kp/((0,1)*w0 + wp^2/(0,1)/w0)*exp(-(0,1)*k1*z + kp * x)).x;
}
real E_x(real H0, real z, real x)
{
	if(x>0)
		return (-H0*k1/w0*exp(-(0,1)*k1*z - kv * x)).x;
	else
		return (H0*k1/(w0 - wp^2/w0)*exp(-(0,1)*k1*z + kp * x)).x;
}
int Nz = 512;
int Nx = 512;
real[][] Hy = new real[Nz][Nx];
real[][] Ez = new real[Nz][Nx];
real[][] Ex = new real[Nz][Nx];
real Lz=2*pi/k1*3;
real dz=Lz/Nz;
real Lx=2*pi/k1*30;
real dx=Lx/Nx;
for(int i = 0; i != Nz; ++i)
{
	real z=i*dz;
	for(int j = 0; j != Nx; ++j)
	{
		real x=(j-Nx/2)*dx;
		Hy[i][j]=H_y(1.0, z, x);
		Ez[i][j]=E_z(1.0, z, x);
		Ex[i][j]=E_x(1.0, z, x);
	}
}
plot pHy;
pHy.setbounds(0,-Lx/2,Lz,Lx/2);
pHy.setlabels("$z$","$x$","$H_y$");
pHy.addplot(Hy);
pHy.show((7cm,4cm));
plot pEz;
pEz.setbounds(0,-Lx/2,Lz,Lx/2);
pEz.setlabels("$z$","$x$","$E_z$");
pEz.addplot(Ez);
pEz.show((7cm,4cm));
plot pEx;
pEx.setbounds(0,-Lx/2,Lz,Lx/2);
pEx.setlabels("$z$","$x$","$E_x$");
pEx.addplot(Ex);
pEx.show((7cm,4cm));

print("$\mathrm {div} \mathbf E$ is:");
real[][] divE = new real[Nz][Nx];
for(int i = 0; i != Nz; ++i)
{
	for(int j = 0; j != Nx; ++j)
		divE[i][j] = 0.0;
}
for(int i = 1; i != Nz-1; ++i)
{
	for(int j = 1; j != Nx-1; ++j)
		divE[i][j] = -(Ex[i][j+1] - Ex[i][j-1])/dx/2 + (Ez[i+1][j] - Ez[i-1][j])/dz/2;
}
// why "-" dEx/dx ???

real divEmax=0;
real divEmax_vacuum=0;

for(int i = 1; i != Nz-1; ++i)
{
	for(int j = 1; j != Nx-1; ++j)
	{
		divEmax = max(divEmax, abs(divE[i][j]));
		divEmax_vacuum = max(divEmax_vacuum, j < Nx/2 + 2 ? divEmax_vacuum : abs(divE[i][j]));
	}
}

for(int i = 1; i != Nz-1; ++i)
{
	for(int j = 1; j != Nx-1; ++j)
		if(abs(j - Nx/2) < 2)
			divE[i][j] = 0.0;
}

write("div E max = ", divEmax);
write("div E max vacuum = ", divEmax_vacuum);
plot pdivE;
pdivE.setbounds(0,-Lx/2,Lz,Lx/2);
pdivE.setlabels("$z$","$x$","$\mathrm{div} \mathbf E$ in vaccum");
pdivE.addplot(divE);
pdivE.show((7cm,4cm));

