#pragma once
#include "defaults.h"

#include <map>
#include <string>
#include <chrono>
#include <iostream>

#ifndef nullptr
	#define nullptr 0
#endif

void cpustep(int N_x, int N_z, FL_DBL* E_z, FL_DBL* E_x, FL_DBL* B_y, FL_DBL dx_1, FL_DBL dz_1, FL_DBL dt);

struct cell;

extern "C"
struct chrono_data_wrapper
{
	chrono_data_wrapper() : t(std::chrono::high_resolution_clock::now()), dt(.0), count(0) {}
	std::chrono::high_resolution_clock::time_point t;
	double dt;
	int count;
};

extern "C"
struct map_timer
{
	map_timer() : data() {}

  map_timer(const map_timer &obj) : data()
  {
    for(auto& it : obj.data)
      data[it.first] = it.second;
  }

	std::map<std::string, chrono_data_wrapper> data;

	void start(std::string s)
	{
		auto it = data.find(s);
		if(it == data.end())
			data[s] = chrono_data_wrapper();
		else
			it->second.t = std::chrono::high_resolution_clock::now();
	}
	
	void stop(std::string s)
	{
		auto it = data.find(s);
		if(it != data.end())
		{
			std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
			it->second.dt += std::chrono::duration<double>(t1 - it->second.t).count();
			it->second.count ++;
		}
	}
	
	~map_timer()
	{
		for(const auto& it : data)
			std::cout << it.first << " is " << it.second.dt / (double) it.second.count << " seconds per step, " << it.second.dt << " second overall\n";
	}
};

extern "C"
struct hydro2dHandler
{
	map_timer tim;
	const int device; // deny moving to another device by switchig this field
	//holds information about 
	int Nx, Nz; // local
	FL_DBL Lx, Lz; // local
	int ix0,iz0; // position in global mesh
	FL_DBL dx_1,dz_1,dt;
	FL_DBL* Ptz;
	FL_DBL* Ptx;
	FL_DBL* Pz;
	FL_DBL* Px;
	FL_DBL* Ez;
	FL_DBL* Ex;
	FL_DBL* Ez_mid;
	FL_DBL* Ex_mid;
	FL_DBL* By;
	FL_DBL* Phtz;
	FL_DBL* Phtx;
	FL_DBL* Phz;
	FL_DBL* Phx;
	FL_DBL* Jz;
	FL_DBL* Jx;
	FL_DBL* n;
	FL_DBL* Te;
	FL_DBL* dTe;
	FL_DBL* nlx;
	FL_DBL* nlz;
	FL_DBL* vx;
	FL_DBL* vz;
	FL_DBL* PML_Byx; // helper fields required for PML
	FL_DBL* PML_Byz; // helper fields required for PML
	//feed arrays
	FL_DBL* feedPtz;
	FL_DBL* feedPtx;
	FL_DBL* feedPz;
	FL_DBL* feedPx;
	FL_DBL* feedEz;
	FL_DBL* feedEx;
	FL_DBL* feedBy;
	FL_DBL* feedPML_Byx;
	FL_DBL* feedPML_Byz;
	FL_DBL* feedPhtz;
	FL_DBL* feedPhtx;
	FL_DBL* feedPhz;
	FL_DBL* feedPhx;
	FL_DBL* feedJz;
	FL_DBL* feedJx;
	FL_DBL* feedn;
	FL_DBL* feedTe;
	FL_DBL* srcx;
	FL_DBL* srct;
	FL_DBL* host_feed;
	FL_DBL* host_srct;
	double t;
	int step;
	int PMLimin; // x nodes occuied by PML layer
	int PMLimax; // x nodes occuied by PML layer
	int SRCi;
	const int shft; // NOTE: shift is counted in blocks
	const int blockSize;
	FL_DBL PMLxmax; // PML depth
	FL_DBL PMLstrength;
	FL_DBL mediaN0;
	FL_DBL mediaTe0; // effective temperatire of electrons
	FL_DBL diffusion; // thermodiffusion
	FL_DBL mediaNu;
	FL_DBL NUTratio; // this is how nu depends on T : nu = mediaNu * (1.0 + NUTratio * Te)
	FL_DBL landauDamping;
	//effective electron massses
	FL_DBL mz_1;
	FL_DBL mx_1;
	//media geometry
	FL_DBL toothDepth;
	FL_DBL toothWidth;
	FL_DBL mediaDepth;
	cell* a_cell;
	//media electrodynamic specific features
	//parameters of lorents electrons describing bounded electrons in media
	FL_DBL media_bound_w2;
	FL_DBL media_bound_gamma;
	FL_DBL media_bound_beta;
	//phonon absorbtion parameters
	FL_DBL media_phonon_omega;
	FL_DBL media_phonon_phw;
	FL_DBL media_phonon_beta;
	//source parameters
	FL_DBL srcTfactor; // factor for thermal source
	FL_DBL srcAmp;
	FL_DBL srcVelocity;
	FL_DBL srcT;
	FL_DBL srcNosc;
	FL_DBL srcX;
	FL_DBL switchOnDelay;
	const bool doPolarization; // never change this after construction! because reallocation will be needed
	const bool doPhononAbsorbtion; // never change this after construction! because reallocation will be needed
	const bool extSource;
	const int JHEAT; // 0-E^2, 1-j^2, 2-j*E
	bool toothDir; // positive/negative
	bool flip; // if flip then tooths are below surface
	hydro2dHandler(int /*dev*/, bool = false /*doPhononAbsorbtion*/, bool = false /*doPolarization*/, bool = false/*extSource*/, int = 0 /*JHEAT*/);
	hydro2dHandler(const hydro2dHandler& obj, int dev = -1);
private:
	hydro2dHandler();
};

extern "C"
void dev_alloc(hydro2dHandler* hH, void** pparr, int sz);

extern "C"
void dev_h2d(hydro2dHandler* hH, FL_DBL* host_arr, FL_DBL* dev_arr, int sz);

extern "C"
void dev_d2h(const hydro2dHandler* hH, const FL_DBL* dev_arr, FL_DBL* host_arr, int sz);

extern "C"
void GPUsetField(hydro2dHandler* hH, FL_DBL* E_z, FL_DBL* E_x, FL_DBL* B_y);

extern "C"
void GPUgetField(hydro2dHandler* hH, FL_DBL* E_z, FL_DBL* E_x, FL_DBL* B_y);

extern "C"
int simpleGPUstep(hydro2dHandler*);

extern "C"
void simpleGPUinit(hydro2dHandler*); // do all allocations and so on...

