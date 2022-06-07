import grph;
import patterns;

print("\textbf{Input parameters for moving window $($test\_mGPU binary$)$}");
printp("This is a default version implementing a moving window case. Compiled with default options $($i.e. without STATIC define$)$. Below is a 2D box representing grid for all fields under simulation:");

{
	real PML=1cm;
	real Lx=7cm;
	real Lz=10cm;
	real depth=2.3cm;
	add("hatch",hatch(2.4mm));
	add("chatch",crosshatch(1.6mm));
	real sz_shift=1cm;
	real sz_subshift=-0.3cm;
	
	pair tlc = (2.5cm,cr-0.5cm);
	
	filldraw(shift(tlc)*scale(Lz,-Lx)*unitsquare, 0.75*red+0.85*green+0.95*blue, black);
	//plasma
	filldraw(shift(tlc)*shift(0,-Lx/2)*scale(Lz,-depth)*unitsquare,pattern("chatch"));

	//PML
	filldraw(shift(tlc)*scale(Lz,-PML)*unitsquare,pattern("hatch"));
	//draw(shift(tlc)*((Lz,0)--(Lz+sz_shift,0)));
	draw(shift(tlc)*((Lz,-PML)--(Lz+sz_shift,-PML)));
	draw(shift(tlc)*((Lz+sz_shift+sz_subshift,0)--(Lz+sz_shift+sz_subshift,-PML)), Arrows);
	label("PMLx",shift(tlc)*(Lz+sz_shift+sz_subshift,-PML/2),E);
	//draw(shift(tlc)*((Lz,-Lx)--(Lz+sz_shift,-Lx)));
	draw(shift(tlc)*((Lz,-Lx+PML)--(Lz+sz_shift,-Lx+PML)));
	draw(shift(tlc)*((Lz+sz_shift+sz_subshift,-Lx)--(Lz+sz_shift+sz_subshift,-Lx+PML)), Arrows);
	label("PMLx",shift(tlc)*(Lz+sz_shift+sz_subshift,-Lx+PML/2),E);

	//plasma
	draw(shift(tlc)*((Lz,-Lx/2)--(Lz+sz_shift,-Lx/2)));
	draw(shift(tlc)*((Lz,-Lx/2-depth)--(Lz+sz_shift,-Lx/2-depth)));
	draw(shift(tlc)*((Lz+sz_shift+sz_subshift,-Lx/2)--(Lz+sz_shift+sz_subshift,-Lx/2-depth)), Arrows);
	label(rotate(90)*"mediaDepth",shift(tlc)*(Lz+sz_shift+sz_subshift,-Lx/2-depth/2),E);

	filldraw(shift(tlc)*shift(0,-Lx+PML)*scale(Lz,-PML)*unitsquare,pattern("hatch"));

	//src
	real src = -Lx/2+Lx/30;
	draw(shift(tlc)*((0,src)--(Lz,src)), dashed);
	draw(shift(tlc)*((Lz+1cm,src+0.5cm)--(Lz,src)), EndArrow);
	label("src", shift(tlc)*(Lz+1cm,src+0.5cm), E);

	//axes
	sz_shift=0.25cm;
	draw(shift(tlc)*((0,sz_shift)--(Lz+0.25cm,sz_shift)), EndArrow);
	label("$z$", shift(tlc)*(Lz+0.25cm,sz_shift),N);
	sz_shift=1cm;
	draw(shift(tlc)*((-sz_shift,0)--(-sz_shift,-Lx-0.25cm)), EndArrow);
	label("$x$", shift(tlc)*(-sz_shift,-Lx-0.25cm),W);

	//sizes
	sz_shift=3cm;
	draw(shift(tlc)*((Lz,0)--(Lz+sz_shift,0)));
	draw(shift(tlc)*((Lz,-Lx)--(Lz+sz_shift,-Lx)));
	draw(shift(tlc)*((Lz+sz_shift+sz_subshift,0)--(Lz+sz_shift+sz_subshift,-Lx)),Arrows);
	label(rotate(90)*"Lx, Nx", shift(tlc)*(Lz+sz_shift+sz_subshift,-0.5*Lx), E);

	sz_shift=1cm;
	draw(shift(tlc)*((0,-Lx)--(0,-Lx-sz_shift)));
	draw(shift(tlc)*((Lz,-Lx)--(Lz,-Lx-sz_shift)));
	draw(shift(tlc)*((0,-Lx-sz_shift-sz_subshift)--(Lz,-Lx-sz_shift-sz_subshift)), Arrows);
	label("Lz, Nz", shift(tlc)*(Lz*0.5, -Lx-sz_shift-sz_subshift), N);
	
	
	cr=tlc.y-Lx-sz_shift-14pt;
}
print("So...");

