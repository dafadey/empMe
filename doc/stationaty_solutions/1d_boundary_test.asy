import graph;

int n=100;

real[] ez=new real[n];
real[] b=new real[n+1];

real V=1/cos(5/180*pi);

// variables in which you have no dispersion and
// perfect open boundary iff d x = d new_dz
// new_b = b * V / sqrt(V^2+1)
// new_z = z / sqrt(V^2+1)

real dx=1;
real dz=dx*sqrt(V^2-1);

for(int i=0; i!=n; ++i)
{
  real x=dx*i;
  real L=n*dx;
  real _x = (x-L/2)/10;
  ez[i] = abs(_x)<1 ? cos(_x*pi)+1 : 0;
}
for(int i=0; i!=n+1; ++i)
  b[i]=0;


for(int i=0; i!=45/*try any value here to see how waves pass through boundary*/; ++i)
{
  for(int i=0;i!=n;++i)
    ez[i] += -(b[i+1]-b[i])/dx /V * dz;
  b[0] += -(ez[0]+b[0]/V*sqrt(V^2-1))/dx *V/(V^2-1) * dz;
  for(int i=1;i!=n;++i)
    b[i] += -(ez[i]-ez[i-1])/dx *V/(V^2-1) * dz;
  b[n] += -(b[n]/V*sqrt(V^2-1)-ez[n-1])/dx *V/(V^2-1) * dz;
}

guide ge, gb;
for(int i=0; i!=n;++i)
  ge=ge--((i+0.5)*dx,ez[i]);
for(int i=0; i!=n+1;++i)
  gb=gb--(i*dx,b[i]);
picture p;
draw(p,ge,red);
draw(p,gb,blue);
draw(p,(0,-2));
draw(p,((n+1)*dx,2));
size(p,5cm, 4cm, point(p,SW), point(p,NE));
add(p.fit());
