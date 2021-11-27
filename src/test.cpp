#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "funcall.h"
#include <ctime>
#include <string>
#include <chrono>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
//#include <zlib.h>
#include "simpledraw.h"
#include "executor.h"

#include "getGPUtemp.h"

//#include "nvml.h"
//#include "tclient.h"

#include "helpers.h"

#include "zstream.h"
#include <unistd.h> // sleep, getpid
#include "geo.h"

struct timer
{
private:
  std::string message;
  std::ostream* out;
  std::chrono::high_resolution_clock::time_point t_start;
public:
  timer(const std::string m, std::ostream& s = std::cout) : message(m), out(&s), t_start(std::chrono::high_resolution_clock::now()) {}
  ~timer()
  {
    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    std::stringstream res;
		*out << message << " " << std::chrono::duration<double>(t_end-t_start).count() << " s" << std::endl;
  }
  timer(const timer& t) : message(t.message), out(t.out), t_start(t.t_start) {}
  const timer &operator=(const timer& t)
  {
    if(this != &t)
    {
      this->message = t.message;
      this->out = t.out;
      this->t_start = std::chrono::high_resolution_clock::now();
    }
    return *this;
  }
};

void dumpData(std::string outputFilename, const hydro2dHandler* hH, std::vector<const FL_DBL*> data, std::vector<std::string> header, std::string description)
{
	timer my("dumpData timig is");
	gzostream gzout(outputFilename);
	gzout << "DESCRIPTION START" << std::endl;
	gzout << description;
	gzout << "DESCRIPTION END" << std::endl;
	gzout << hH->Nz << '\t' << hH->Nx << '\t' << data.size() << std::endl;
	size_t i=0;
	gzout << "z\tx\t";
	for(; i != std::min(data.size(), header.size()); i++)
		gzout << (i ? "\t" : "") << header[i];
	if(i != data.size())
		printf("dump::WARNING: data size and header does not match!!!\n");
	for(; i != data.size(); i++)
		gzout << "\tunknown";
	gzout << std::endl;

	for(int i=0;i<hH->Nz;i++)
	{
		for(int j=0;j<hH->Nx;j++)
		{
			gzout << FL_DBL(i)/hH->dz_1  << '\t' << FL_DBL(j)/hH->dx_1;
			for(size_t k = 0; k != data.size(); k++)
				gzout << "\t" << data[k][i * hH->Nx + j];
			gzout << std::endl;
		}
	}	
	gzout.close();
	return;
}


int halt()
{
	char* a=new char[100];
	printf("paused. enter any number to proceed\n");
	scanf("%s",a);
	int i = atoi(a);
	delete[] a;
	return i;
}

bool cool_dev(int dev)
{
	int t=getGPUtemp(dev);
	if(t < 75)
		return true;
	while(getGPUtemp(dev)>75)
		sleep(7);
	printf("temp is %d\n",t);
	return true;
}

void GPUstep(void* arg)
{
	simpleGPUstep((hydro2dHandler*) arg);
}

int main(int argc, char* argv[])
{
 //struct tms time_buf;
	printf("Hi universe\n");
	/*
	if(nvmlInit() != NVML_SUCCESS)
	{
		std::cerr << "Failed to initialize nvml library to control temperature of GPU boards" << std::endl;
		return -1;
	}
	*/
	bool nonlinVSlin(false);
	
	int Nx = NX;
	int Nz = NZ;
	FL_DBL Lx = LX;
	FL_DBL Lz = LZ;

	int dev(0);
	
	bool toothDir(true);
	bool flip(true);
	bool DRAW(false);
	bool doPhononAbsorbtion(false);
	bool doPolarization(false);
	bool extSource(false);
	int jHeat(2);
  
  FL_DBL NUTratio(NUTRATIO);
	FL_DBL diffusion(DIFFUSION);
	FL_DBL mediaN0(N0);
	FL_DBL mediaTe0(TE0);
	FL_DBL mediaNu(MEDIANU);
	FL_DBL landauDamping(LANDAUDAMPING);
	FL_DBL tMAX(TMAX);
	FL_DBL toothDepth(TOOTHDEPTH);
	FL_DBL toothWidth(TOOTHWIDTH);
	FL_DBL mediaDepth(MEDIADEPTH);
	FL_DBL bound_w2(BOUND_W2);
	FL_DBL bound_beta(BOUND_BETA);
	FL_DBL bound_gamma(BOUND_GAMMA);
	FL_DBL phonon_phw(PHONON_PHW);
	FL_DBL phonon_omega(PHONON_OMEGA);
	FL_DBL phonon_beta(PHONON_BETA);
	FL_DBL mx_1(MX_1);
	FL_DBL my_1(MY_1);
	FL_DBL mz_1(MZ_1);
	
	std::vector<FL_DBL> srcAmp;
	std::vector<FL_DBL> srcT;
	std::vector<FL_DBL> srcNosc;
	std::vector<FL_DBL> srcPhase;
	std::vector<FL_DBL> srcAmpTE;
	std::vector<FL_DBL> srcNoscTE;
	std::vector<FL_DBL> srcTTE;
	std::vector<FL_DBL> srcPhaseTE;

	FL_DBL switchOnDelay = SWITCHONDELAY;
	FL_DBL srcX(SRCX);
	FL_DBL PMLx(PMLX);
	FL_DBL PMLstrength(PMLSTRENGTH);
	FL_DBL cell_scalex(1.0);
	FL_DBL cell_scaley(1.0);
	FL_DBL Tavg(0.0);
	FL_DBL srcTfactor(SRCTFACTOR);
	
	int iter(ITER);

	FL_DBL velocity(VELOCITY);
	std::string outputFilename(OUTPUTFILENAME);
	std::string cellFilename("");
	#define PARSE(ARG) if(name == #ARG) { sval >> ARG; continue;}
	#define PARSE2(ARG, VAR) if(name == #ARG) { sval >> VAR; continue;}
	for(int i=1;i<argc;i++)
	{
		std::string inp = std::string(argv[i]);
		std::cout << "---" << inp << '\n';
		size_t pos = inp.find("=");
		if(pos == std::string::npos)
			printf("you specified parameter wrong way use <name>=<value> format. NOTE: no \'-\' and spaces\n");
		else
		{
			std::string name = inp.substr(0,pos);
			std::stringstream sval;
			sval << inp.substr(pos+1,std::string::npos);
			printf("parameter[%d] has name %s and value %s\n",i-1,name.c_str(), sval.str().c_str());
			PARSE(PMLx);
			PARSE(PMLstrength);
			PARSE(Lx);
			PARSE(Lz);
			PARSE(Nx);
			PARSE(Nz);
			PARSE(switchOnDelay);
			PARSE(srcTfactor);
			PARSE(srcX);
			PARSE(srcT);
			PARSE(srcNosc);
			PARSE(srcNoscTE);
			PARSE(srcTTE);
			PARSE(srcAmp);
			PARSE(srcAmpTE);
			PARSE(srcPhaseTE);
			PARSE(srcPhase);
			PARSE(mediaNu);
			PARSE(landauDamping);
			PARSE(mediaN0);
			PARSE(mediaTe0);
			PARSE(toothWidth);
			PARSE(toothDepth);
			PARSE(velocity);
			PARSE(mediaDepth);
			PARSE2(output,outputFilename);
			PARSE(Tavg);
      PARSE(NUTratio);

			PARSE2(cell,cellFilename);
			PARSE(cell_scalex);
			PARSE(cell_scaley);

			PARSE2(toothDirection,toothDir);
			PARSE(flip);
			PARSE(DRAW);
			PARSE(mz_1);
			PARSE(mx_1);
			PARSE(my_1);
			PARSE(diffusion);
			PARSE2(doPhonon, doPhononAbsorbtion);
			PARSE(doPolarization);
			PARSE(extSource);
      PARSE(jHeat);
			PARSE2(Tmax,tMAX);
			PARSE2(device,dev);
			PARSE(bound_w2);
			PARSE(bound_beta);
			PARSE(bound_gamma);
			PARSE2(phOmega, phonon_omega);
			PARSE2(phWidth, phonon_phw);
			PARSE2(phBeta, phonon_beta);
			PARSE2(silentIterations,iter);
			PARSE(nonlinVSlin);
			printf("ERROR: input parameter %s is unknown\n",name.c_str());
		}
	}
	#undef PARSE
	#undef PARSE2

	if(srcAmp.size()==0)
		srcAmp.push_back(SRCAMP);
	if(srcAmpTE.size()==0)
		srcAmpTE.push_back(SRCAMP_TE);

	#define FILLTO(REF,ARR, VAR)	for(size_t i=ARR.size(); i < REF.size(); i++)	ARR.push_back(VAR);
	//TM
	FL_DBL lastNosc = srcNosc.size() ? srcNosc[srcNosc.size() - 1] : SRCNOSC;
	FILLTO(srcAmp, srcNosc, lastNosc);

	FILLTO(srcAmp, srcPhase, 0);
	
	FL_DBL lastSrcT = srcT.size() ? srcT[srcT.size() - 1] : SRCT;
	FILLTO(srcAmp, srcT, lastSrcT);

	//TE
	lastNosc = srcNoscTE.size() ? srcNoscTE[srcNoscTE.size() - 1] : SRCNOSC;
	FILLTO(srcAmpTE, srcNoscTE, lastNosc);

	FILLTO(srcAmpTE, srcPhaseTE, 0);

	lastSrcT = srcTTE.size() ? srcTTE[srcTTE.size() - 1] : SRCT;
	FILLTO(srcAmpTE, srcTTE, lastSrcT);

	#undef FILLTO

	cell* c_ptr(nullptr);
	if(cellFilename != "")
	{
		c_ptr = new cell();
		if(!(c_ptr->load_from_svg(cellFilename.c_str(), toothWidth)))
		{
			std::cout << "Failed to load cel geometry from " << cellFilename << std::endl;
			delete c_ptr;
			c_ptr = nullptr;
		}
		if(c_ptr)
		{
			c_ptr->scalex((double) cell_scalex);
			c_ptr->scaley((double) cell_scaley);
		}
	}
	
	const FL_DBL dx_1 = FL_DBL(Nx)/Lx;
	const FL_DBL dz_1 = FL_DBL(Nz)/Lz;
	const FL_DBL dt = FPT(0.5)*((FPT(1.0)/dx_1<FPT(1.0)/dz_1)?FPT(1.0)/dx_1:FPT(1.0)/dz_1);
	std::stringstream description;
  #ifdef USEDOUBLE
  description << "accuracy: double\n";
  #else
  description << "accuracy: float\n"; 
  #endif
	description << "grid: Lx=" << Lx << ", Lz=" << Lz << ", Nx=" << Nx << ", Nz=" << Nz << ", dx_1=" << dx_1 << " (dx=" << FPT(1.0)/dx_1 << "), dz_1=" << dz_1 << " (dz=" << FPT(1.0)/dz_1 << "), dt=" << dt << ", PMLx=" << PMLx << ", PMLstrength=" << PMLstrength << std::endl;
	description << "run: Tmax=" << tMAX << ", Tavg=" << Tavg << ", output=" << outputFilename << ", silentIterations(iterations before draw/save)=" << iter << ",\ngeometry: toothWith=" << toothWidth<< ", toothDepth=" << toothDepth << ", mediaDepth=" << mediaDepth;
	if(c_ptr)
		description << ", cellFilename=" << cellFilename << ", scalex=" << cell_scalex << ", scaley=" << cell_scaley;
	description << ", flip=" << (flip?"true":"false") << ", toothDirection=" << (toothDir?"true":"false") << ",\nmedia: n0=" << mediaN0 << ", T0=" << mediaTe0 << ", nu=" << mediaNu << ", diffusion=" << diffusion << "," << " mx_1=" << mx_1 << "," << " my_1=" << my_1 << " mz_1=" << mz_1 << ", NUTratio=" << NUTratio << ",";
	if (doPolarization)
	description << "\n       bound_w2=" << bound_w2 << " bound_beta=" << bound_beta << ", bound_gamma=" << bound_gamma << ",";
	if (doPhononAbsorbtion)
	description << "\n       phOmega=" << phonon_omega << ", phWidth=" << phonon_phw << ", phBeta=" << phonon_beta << ",";
	description << "\nvizualization: DRAW=" << (DRAW?"true":"false") << std::endl;
	if (extSource)
		description << "         therm source is external, srcX=" << srcX;
	else
		description << "         source is electromagnetic (use extSource=1 to switch to thermal source)";
  description << ", srcAmp=" << srcAmp << ", srcT=" << srcT <<  ", srcPhase=" << srcPhase << ", srcNosc=" << srcNosc << ", srcAmpTE=" << srcAmpTE << ", srcTTE=" << srcTTE <<  ", srcPhaseTE=" << srcPhaseTE << ", srcNoscTE=" << srcNoscTE << ", switchOnDelay=" << switchOnDelay << ", srcTfactor=" << srcTfactor << std::endl;
	description << "         velocity=" << velocity << std::endl;
	if (jHeat==1)
	  description << "         electrons heated by current";
	else if (jHeat==2)
	  description << "         electrons heated classicaly by (j, E)";
	else if (jHeat==0)
	  description << "         electrons heated by electric field";
	description << std::endl;
	description << (nonlinVSlin ? "   generated signal = non_linear - linear" : "   generated signal = 2 * (non_linear(A) - 2 * non_linear(A/2))");
	description << std::endl;
	std::cout << description.str();
	//CPU
	// alloc
	FL_DBL* exCPU=new FL_DBL[Nx*Nz];
	FL_DBL* ezCPU=new FL_DBL[Nx*Nz];
	FL_DBL* byCPU=new FL_DBL[Nx*Nz];
	// init
	for(int i=0;i<Nz;i++)
	{
		for(int j=0;j<Nx;j++)
		{
			byCPU[i*Nx+j]=0.0;
			exCPU[i*Nx+j]=0.0;
			ezCPU[i*Nx+j]=0.0;
		}
	}
	if(DRAW)
	{
		printf("visual mode is on\n");
		fadey_init(Nx,Nz,9);
		//fadey_draw(byCPU,Nx,Nz,2);
		//char a;
		//scanf("%c\n",&a);
		printf("init is done\n");
	}
	
	printf("data inited on CPU\n");
	
	//GPU declare
	//alloc CPU copy
	FL_DBL* exGPU=new FL_DBL[Nx*Nz];
	FL_DBL* eyGPU=new FL_DBL[Nx*Nz];
	FL_DBL* ezGPU=new FL_DBL[Nx*Nz];
	FL_DBL* jxGPU=new FL_DBL[Nx*Nz];
	FL_DBL* jyGPU=new FL_DBL[Nx*Nz];
	FL_DBL* jzGPU=new FL_DBL[Nx*Nz];
	FL_DBL* bxGPU=new FL_DBL[Nx*Nz];
	FL_DBL* byGPU=new FL_DBL[Nx*Nz];
	FL_DBL* bzGPU=new FL_DBL[Nx*Nz];
	FL_DBL* TeGPU=new FL_DBL[Nx*Nz];
	FL_DBL* neGPU=new FL_DBL[Nx*Nz];
	FL_DBL* nWeakGPU=new FL_DBL[Nx*Nz];
	FL_DBL* byWeakGPU=new FL_DBL[Nx*Nz];
	FL_DBL* eyWeakGPU=new FL_DBL[Nx*Nz];
	FL_DBL* generationBy=new FL_DBL[Nx*Nz];
	FL_DBL* generationEy=new FL_DBL[Nx*Nz];
	
  std::cout << "NUMBER OF AVAILABLE CUDA DEVICES IS " << CUDA_device_count() << '\n';
  
	//init simple
	hydro2dHandler* hH = new hydro2dHandler(dev /*device*/ % CUDA_device_count(), doPhononAbsorbtion/*doPhononAbsorbtion*/, doPolarization/*doPolarization*/, extSource, jHeat);
	hH->Nx = Nx;
	hH->Nz = Nz;
	hH->Lx = Lx;
	hH->Lz = Lz;
	hH->ix0 = 0;
	hH->iz0 = 0;
	hH->dx_1 = dx_1;
	hH->dz_1 = dz_1;
	hH->dt = dt;
	hH->PMLxmax = PMLx;
	hH->PMLstrength = PMLstrength;
	hH->srcVelocity = velocity;
	hH->srcTfactor = srcTfactor;
	hH->srcT = srcT;
	hH->srcTTE = srcTTE;
	hH->srcNosc = srcNosc;
	hH->srcNoscTE = srcNoscTE;
	hH->srcAmp = srcAmp;
	hH->srcAmpTE = srcAmpTE;
	hH->srcPhase = srcPhase;
	hH->srcPhaseTE = srcPhaseTE;
	hH->mediaN0 = mediaN0;
	hH->mediaTe0 = mediaTe0;
	hH->mediaNu = mediaNu;
	hH->diffusion = diffusion;
	hH->NUTratio = NUTratio;
	hH->media_bound_w2 = bound_w2;
	hH->media_bound_gamma = bound_gamma;
	hH->media_bound_beta = bound_beta;
	hH->mz_1 = mz_1;
	hH->mx_1 = mx_1;
	hH->media_phonon_omega = phonon_omega;
	hH->media_phonon_phw = phonon_phw;
	hH->media_phonon_beta = phonon_beta;
	hH->toothDepth = toothDepth;
	hH->toothWidth = toothWidth;
	hH->mediaDepth = mediaDepth;
	hH->a_cell = c_ptr;
	hH->toothDir = toothDir;
	hH->flip = flip;
	hH->switchOnDelay = switchOnDelay;
	hH->srcX = srcX;
	hH->landauDamping = landauDamping;
	
	simpleGPUinit(hH); // do all allocations and so on...
	
	hydro2dHandler* hH_weak = nullptr;
	if(!extSource)
	{
		hH_weak = new hydro2dHandler(*hH, (dev + 1) % CUDA_device_count());
		if(nonlinVSlin)
			hH_weak->linear = true;
		else {
			if(hH->srcAmp.size() != hH_weak->srcAmp.size())
			{
				std::cerr << "hH->srcAmp.size() != hH_weak->srcAmp.size()\n";
				return -1;
			}
			if(hH->srcAmpTE.size() != hH_weak->srcAmpTE.size())
			{
				std::cerr << "hH->srcAmpTE.size() != hH_weak->srcAmpTE.size()\n";
				return -1;
			}

			for(size_t i=0; i<hH->srcAmp.size(); i++)
				hH_weak->srcAmp[i] = FPT(0.5) * hH->srcAmp[i];

			for(size_t i=0; i<hH->srcAmp.size(); i++)
				hH_weak->srcAmpTE[i] = FPT(0.5) * hH->srcAmpTE[i];

		}
	}
	
	std::cout << "main handler is " << (hH->linear ? "linear" : "nonlinear") << ", weak handler is " << (hH_weak->linear ? "linear" : "nonlinear") << '\n';
	
	printf("USING toothDepth=%g, toothWidth=%g\n",hH->toothDepth,hH->toothWidth);

	if(!hH)
	{
		printf("ERROR: initalizing handler\n");
		return -1;
	}
	GPUsetField(hH, ezCPU, exCPU, byCPU);
	if(hH_weak)
		GPUsetField(hH_weak, ezCPU, exCPU, byCPU);

/*
	if(DRAW)
	{
		GPUgetField(hH, ezGPU, exGPU, byGPU);
		dev_d2h(hH, hH->Te, TeGPU, hH->Nx * hH->Nz);
		dev_d2h(hH, hH->n, nWeakGPU, hH->Nx * hH->Nz);
		if(hH_weak)
			dev_d2h(hH_weak, hH_weak->By, byWeakGPU, hH->Nx * hH->Nz);
		
		fadey_draw(byWeakGPU,Nx,Nz,2);
		fadey_draw(ezGPU,Nx,Nz,3);
		fadey_draw(exGPU,Nx,Nz,4);
		fadey_draw(byGPU,Nx,Nz,5);
		fadey_draw(TeGPU,Nx,Nz,6);
		fadey_draw(nWeakGPU,Nx,Nz,7);
		for(int i=0;i<Nx*Nz;i++)
			generationBy[i]=(byGPU[i]-byWeakGPU[i]*FPT(2.0))*FPT(2.0);
		fadey_draw(generationBy,Nx,Nz,8);
	}
*/

	executor e(2);
	//calculate and save generation
	FL_DBL extraction_factor = nonlinVSlin ? FPT(1.0) : FPT(2.0);
  map_timer tim;
	while(hH->t < tMAX)
	{

		if(hH->step % iter == 0)
		{
			bool have_Tserver=true;
			have_Tserver &= cool_dev(hH->device);
			if(hH_weak)
				have_Tserver &= cool_dev(hH_weak->device);
			
			if(!have_Tserver)
			{
				std::cerr << "please check if nvidiaTserever(tserver) is running" << std::endl;
				break;
			}

      dev_d2h(hH, hH->Bx, bxGPU, hH->Nx * hH->Nz);
      dev_d2h(hH, hH->By, byGPU, hH->Nx * hH->Nz);
      dev_d2h(hH, hH->Bz, bzGPU, hH->Nx * hH->Nz);
			if(hH_weak)
				dev_d2h(hH_weak, hH_weak->By, byWeakGPU, hH_weak->Nx * hH_weak->Nz);
      
      dev_d2h(hH, hH->Ex, exGPU, hH->Nx * hH->Nz);
      dev_d2h(hH, hH->Ey, eyGPU, hH->Nx * hH->Nz);
      dev_d2h(hH, hH->Ez, ezGPU, hH->Nx * hH->Nz);
			if(hH_weak)
				dev_d2h(hH_weak, hH_weak->Ey, eyWeakGPU, hH_weak->Nx * hH_weak->Nz);

      dev_d2h(hH, hH->Jx, jxGPU, hH->Nx * hH->Nz);
      dev_d2h(hH, hH->Jy, jyGPU, hH->Nx * hH->Nz);
      dev_d2h(hH, hH->Jz, jzGPU, hH->Nx * hH->Nz);

      dev_d2h(hH, hH->n, neGPU, hH->Nx * hH->Nz);
			if(hH_weak)
				dev_d2h(hH, hH_weak->n, nWeakGPU, hH->Nx * hH->Nz);
      
      dev_d2h(hH, hH->Te, TeGPU, hH->Nx * hH->Nz);


      FL_DBL bnorm = FL_DBL(0);
      FL_DBL bweaknorm = FL_DBL(0);
      FL_DBL bdiffnorm = FL_DBL(0);
      for(int i=0 ; i!=hH->Nx * hH->Nz; i++)
      {
        bnorm += std::abs(byGPU[i]);
        bweaknorm += std::abs(byWeakGPU[i]);
        if(nonlinVSlin)
					bdiffnorm += std::abs(byGPU[i] - byWeakGPU[i]);
        else
					bdiffnorm += FPT(2.0) * std::abs(byGPU[i] - FPT(2.0) * byWeakGPU[i]);
      }

      FL_DBL maxTe = FPT(0.0);
      FL_DBL maxJ = FPT(0.0);
      FL_DBL surf_J_avg = FPT(0.0);
      FL_DBL surf_Te_avg = FPT(0.0);
      FL_DBL surf_nuf_avg = FPT(0.0);
      FL_DBL surf_norm_T = FPT(0.0);
      FL_DBL surf_norm_J = FPT(0.0);
      for(int i=0; i!=hH->Nz; i++)
      {
        for(int j=0; j!=hH->Nx; j++)
        {
          maxTe = std::max(maxTe, TeGPU[i]);
          if(TeGPU[j+i*hH->Nx] > FPT(1.e-5))
          {
            surf_norm_T += FPT(1.0);
            surf_Te_avg += TeGPU[j+i*hH->Nx];
            surf_nuf_avg += FPT(1.0) + hH->NUTratio * SQRTNUT(TeGPU[j+i*hH->Nx]);
            break;
          }
        }
        for(int j=0; j!=hH->Nx; j++)
        {
          FL_DBL J = sqrt(jxGPU[j+i*hH->Nx] * jxGPU[j+i*hH->Nx] + jzGPU[j+i*hH->Nx] * jzGPU[j+i*hH->Nx]);
          maxJ = std::max(maxJ, J);
          if(J > FPT(0.0))
          {
            surf_norm_J += FPT(1.0);
            surf_J_avg += J;
            break;
          }
        }
      }

      if(surf_norm_T == FPT(0.0))
        surf_norm_T = FPT(1.0);
      surf_Te_avg /= surf_norm_T;
      surf_nuf_avg /= surf_norm_T;

      if(surf_norm_J == FPT(0.0))
        surf_norm_J = FPT(1.0);
      surf_J_avg /= surf_norm_J;

      FL_DBL nuf = FPT(1.0) + hH->NUTratio * SQRTNUT(maxTe);
			std::cout << "\niter #=" << hH->step << ", t=" << hH->t << " (target=" << tMAX << "), bnorm=" << bnorm << ", bweaknorm=" << bweaknorm << ", bdiffnorm=" << bdiffnorm << ",\nmaxT=" << maxTe << ", Te0=" << hH->mediaTe0 << ", nu_factor=" << nuf << ", nu=" << hH->mediaNu * nuf << ",\nt_s_avg=" << surf_Te_avg << ", nuf_s_avg=" << surf_nuf_avg << ", nu_s_avg=" << surf_nuf_avg * hH->mediaNu << ", maxJ=" << maxJ << ", surf_J=" << surf_J_avg << "\n" << std::flush;
      
			if(DRAW)
			{

				fadey_draw(bzGPU, Nx, Nz, 0);
				fadey_draw(bxGPU, Nx, Nz, 1);
				fadey_draw(eyGPU, Nx, Nz, 2);
				fadey_draw(ezGPU, Nx, Nz, 3);
				fadey_draw(exGPU, Nx, Nz, 4);
				fadey_draw(byGPU, Nx, Nz, 5);

				fadey_draw(generationBy, Nx, Nz, 8);
				fadey_draw(generationEy, Nx, Nz, 6);
				fadey_draw(neGPU, Nx, Nz, 7);
				
      }
			for(int i(0); i != Nx * Nz; i++) {
				generationBy[i] = (byGPU[i] - byWeakGPU[i] * extraction_factor) * extraction_factor;
				generationEy[i] = (eyGPU[i] - eyWeakGPU[i] * extraction_factor) * extraction_factor;
			}
		}


    /*
    if(hH->step > 4460)
      iter=1;
  */
    //GPUstep((void*)hH);
    //GPUstep((void*)hH_weak);
    
    e.exec(&GPUstep, (void*) hH, 0);
		if(hH_weak)
			e.exec(&GPUstep, (void*) hH_weak, 1);
		e.sync();
				
		if(hH_weak)
		{
			if(hH_weak->t != hH->t || hH_weak->step != hH->step)
			{
				std::cerr << "weak and main are not syncronized, output filename is " << outputFilename << ", step is " << hH->step <<  " vs " << hH_weak->step << '\n';
				return -1;
			}
		}

		if(hH->step % 100 == 0)
		{
	    tim.start("cool device");
      cool_dev(hH->device);
		  if(hH_weak)
		   	cool_dev(hH_weak->device);
      tim.stop("cool device");
		}

	}
	// ok we reached target time.
	// now we want to average generated magnetic field over some period Tavg
	// the trick here is to average taking into accout movig of the whole field picture with velocity V and shifting (iz0 changes during shifting)
	FL_DBL* avg(nullptr);
	if(Tavg > 0.0)
	{
		std::ofstream therm((outputFilename + "therm.dat").c_str());
		std::cout << "averaging" << std::endl;
		avg = new FL_DBL[hH->Nx * hH->Nz * 4];
		for(int i(0); i != Nx * Nz * 4; i++)
			avg[i] = FL_DBL(0.0);
		FL_DBL Norm(0.0);
		while(hH->t < tMAX + Tavg)
		{
			e.exec(&GPUstep, (void*) hH, 0);
			if(hH_weak)
				e.exec(&GPUstep, (void*) hH_weak, 1);
			e.sync();

			if(hH_weak)
			{
				if(hH_weak->t != hH->t || hH_weak->step != hH->step)
				{
					std::cerr << "weak and main are not syncronized, output filename is " << outputFilename << ", step is " << hH->step <<  " vs " << hH_weak->step << '\n';
					return -1;
				}
			}
			
			if(hH->step % iter == 0)
			{
				cool_dev(hH->device);
				if(hH_weak)
					cool_dev(hH_weak->device);
			}
			
			dev_d2h(hH, hH->By, byGPU, hH->Nx * hH->Nz);
			dev_d2h(hH, hH->Ez, ezGPU, hH->Nx * hH->Nz);
			dev_d2h(hH, hH->Ex, exGPU, hH->Nx * hH->Nz);
			if(hH_weak)
				dev_d2h(hH_weak, hH_weak->By, byWeakGPU, hH->Nx * hH->Nz);
			
			if(hH_weak)
			{
				for(int i(0); i != Nx * Nz; i++)
					generationBy[i] = (byGPU[i] - byWeakGPU[i] * extraction_factor) * extraction_factor;
			}
			else
			{
				for(int i(0); i != Nx * Nz; i++)
					generationBy[i] = byGPU[i];
			}

			FL_DBL dz = FL_DBL(1.0) / hH->dz_1;
			#pragma omp parallel for
			for(int i=0; i < Nz; i++)
			{
				FL_DBL z = (FL_DBL(i) - FL_DBL(hH->iz0)) * dz + hH->t * hH->srcVelocity;
				int ii = (int) floor(z * hH->dz_1);
				FL_DBL z0 = FL_DBL(ii) * dz;
				if(ii <= 0 || ii >= Nz - 1)
					continue;
				for(int j(0); j != Nx; j++)
				{
					avg[j + i * Nx] += hH->dz_1 * (generationBy[j + ii * Nx] * (z0 + dz - z) +
																				 generationBy[j + (ii + 1) * Nx] * (z - z0));
					avg[j + i * Nx + Nx * Nz] += hH->dz_1 * (byGPU[j + ii * Nx] * (z0 + dz - z) +
																								   byGPU[j + (ii + 1) * Nx] * (z - z0));
					avg[j + i * Nx + 2 * Nx * Nz] += hH->dz_1 * (ezGPU[j + ii * Nx] * (z0 + dz - z) +
																								       ezGPU[j + (ii + 1) * Nx] * (z - z0));
					avg[j + i * Nx + 3 * Nx * Nz] += hH->dz_1 * (exGPU[j + ii * Nx] * (z0 + dz - z) +
																								       exGPU[j + (ii + 1) * Nx] * (z - z0));
				}
			
			}
			Norm += FL_DBL(1.0);
			
			dev_d2h(hH, hH->Te, TeGPU, hH->Nx * hH->Nz);
			FL_DBL TeMAX(.0);
			for(int i(0); i != hH->Nx * hH->Nz; i++)
				TeMAX = std::max(TeGPU[i], TeMAX);//std::max(TeGPU[i] == FL_DBL(.0) ? FL_DBL(.0) : TeGPU[i] - mediaTe0, TeMAX);
			therm << hH->t << '\t' << TeMAX << '\n';
		}
		for(int i(0); i != Nx * Nz * 4; i++)
			avg[i] = avg[i] / Norm;
		std::cout << "averaging is done, norm was " << Norm << std::endl;
		therm.close();
	}
	e.finish();

	if(hH_weak)
	{
		if(hH_weak->t != hH->t || hH_weak->step != hH->step)
		{
			std::cerr << "weak and main are not in sync after finish, output filename is " << outputFilename << ", step is " << hH->step <<  " vs " << hH_weak->step << '\n';
			return -1;
		}
	}

	//store fields;
	if(hH_weak)
	{
		if(avg)
			dumpData(outputFilename, hH, std::vector<const FL_DBL*>{generationBy, avg, TeGPU, byGPU, &avg[Nx * Nz], ezGPU, &avg[2 * Nx * Nz],  exGPU, &avg[3 * Nx * Nz], neGPU}, std::vector<std::string>{"B_{generated}", "B_{gen_{avg}}", "Te", "B_y", "B_{y_{avg}}", "E_z", "E_{z_{avg}}", "E_x",  "E_{x_{avg}}", "n_e"}, description.str());
		else
			dumpData(outputFilename, hH, std::vector<const FL_DBL*>{generationBy, generationEy, byGPU, ezGPU, exGPU, eyGPU, bzGPU, bxGPU, neGPU},std::vector<std::string>{"B_{generated}", "E_{generated}", "B_y", "E_z", "E_x", "E_y", "B_z", "B_x", "n_e"}, description.str());
	}
	else
		dumpData(outputFilename, hH, std::vector<const FL_DBL*>{byGPU, ezGPU, exGPU, TeGPU, jzGPU, jxGPU, nWeakGPU},std::vector<std::string>{"B_y", "E_z", "E_x", "Te", "j_z", "j_x", "n_e"}, description.str());

	if(DRAW)
		fadey_close();
	
	if(hH_weak)
		delete hH_weak;

	return 0;
}
