real Emax(real w0, real wp, real a)
{
	return sin(a)*cos(a)/(wp^2/w0^2*sin(a)^2+1);
}
real W(real w0, real wp, real a)
{
	return Emax(w0, wp, a)^2;
}

import graph;

real T=72000;
real nosc=16;
real w0=2*pi/(T/nosc);
real n=5e-5; // plasma density
real wp=sqrt(n);

guide gE, gW;
for(real a=7; a<19; a+=0.5)
{
	gE = gE -- (a, Emax(w0, wp, a/180*pi));
	gW = gW -- (a, W(w0, wp, a/180*pi));
}

void addpic(guide g, string name, pair center=(0,0))
{
	picture p;
	draw(p,g);
	size(p, 5cm, 4cm, point(p, SW), point(p, NE));
	xaxis(p, "$\alpha$", BottomTop, LeftTicks(5.0, 1.0));
	yaxis(p, name, LeftRight, RightTicks);
	add(p.fit(), center);
}

addpic(gE, "$E_{max}$", (0,0));
addpic(gW, "$W$", (0,5cm));
