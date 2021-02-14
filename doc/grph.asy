import graph;
import palette;

real cr=0;

real assymfunc(real x)
{
	real f;
	if(x>0) f=1-(x^2);
	else f=1-13*(x^2);
	f=(f>0)?f:0;
	return f;
}

void addplot(real[] data, picture pic, pen p=defaultpen)
{
	guide g;
	int n=data.length;
	for(int i=0;i<n;++i)	g=g--(i/n,data[i]);
	draw(pic,g,p);
}

struct plot
{
	picture pic;
	string labelx;
	string labely;
	bool is2d = false;
	real xmin=0;
	real ymin=0;
	real xmax=1;
	real ymax=1;
	void setbounds(real _xmin, real _ymin, real _xmax, real _ymax)
	{
		xmin=_xmin;
		ymin=_ymin;
		xmax=_xmax;
		ymax=_ymax;
	}
	bounds range;
	string labelz;
	void setlabels(string lx, string ly, string lz="")
	{
			labelx = lx;
			labely = ly;
			labelz = lz;
	}
	pen[] pal;
	void show(pair dim=(5cm,4cm), pair pos=(0,cr))
	{
		xaxis(pic,labelx,BottomTop,LeftTicks,true);
		yaxis(pic,labely,LeftRight,RightTicks,true);
		size(pic,dim.x,dim.y,point(pic,SW),point(pic,NE));
		write(pos);
		add(pic.fit(),(pos.x-min(pic).x+1cm,pos.y-max(pic).y));
		if(is2d)
		{
			picture bar;
			palette(bar,labelz, range, (0,0), (0.5cm,dim.y), Right, pal);
			add(bar.fit(),(pos.x - min(pic).x + max(pic).x - min(bar).x + 1.5cm, pos.y - max(bar).y));
		}
		cr-=dim.y+1.3cm;
	}
	
	void addlabel(string l, pair pos, pair dir=(0,0), pen p=defaultpen)
	{
		dot(pic, pos, p);
		label(pic, l, pos, dir, p);
	}

	void addplot(real[] data, pen p=defaultpen)
	{
		guide g;
		int n=data.length;
		for(int i=0;i<n;++i)	g=g--(i/n,data[i]);
		draw(pic,g,p);
	}

	void addplot(real[][] data, pen[] p=BWRainbow())
	{
		range = image(pic,data,(xmin,ymin),(xmax,ymax),p);
		pal=p;
		is2d=true;
	}

	void addplot(guide g, pen p=defaultpen)
	{
		draw(pic,g,p);
	}

};

void show(picture pic, pair dim=(5cm,4cm), pair pos=(0,cr))
{
	xaxis(pic,"$t$",BottomTop,LeftTicks,true);
	yaxis(pic,"",LeftRight,RightTicks,true);
	size(pic,dim.x,dim.y,point(pic,SW),point(pic,NE));
	write(pos);
	add(pic.fit(),(pos.x-min(pic).x+1cm,pos.y-max(pic).y));
	cr-=5.3cm;
}

void print(string s, real h=14pt)
{
	label(s,(0,cr),E);
	cr-=h;
}
