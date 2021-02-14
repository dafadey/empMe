#pragma once
#define USEDOUBLE
#ifdef USEDOUBLE
	#define FL_DBL double
	#define FPT(x) x
  #define ON 1.0
  #define OFF 0.0
#else
	#define FL_DBL float
	#define FPT(x) x##f // Floating Point Type
  #define ON 1.0f
  #define OFF 0.0f
#endif

// sample run scenario
//./test_mGPU DRAW=1 Tmax=500000 Lx=40000 Nx=1024 Nz=2048 Lz=168000 toothWidth=6000 mediaDepth=0 toothDepth=3000 diffusion=20 extSource=1 srcX=100 srcT=36000 srcAmp=0.0000005 toothDirection=0


#define PMLX FPT(4000.0) // 100.0
#define PMLSTRENGTH FPT(0.05) // 0.1
#define SHIFTGRANULARITY 1  // NOTE: shift is counted in blocks
#define HYDROZSRC OFF
#define HYDROXSRC OFF

#define MAINSRC	OFF // ON
#define PRESSURESRC	OFF // ON
#define HSRC	OFF // ON

#define IFNONLIN  OFF // ON

#define BOUND_W2 6.0e-6
#define BOUND_BETA 0.000423//((omegaOpt/1.7)*(omegaOpt/1.7)*70.0)
#define BOUND_GAMMA 0.0027//(omegaOpt/1.7*1.1)

#define PHONON_OMEGA 0.03
#define PHONON_PHW 0.9
#define PHONON_BETA 0.9

#define N0 FPT(1e-4) // 0.001
#define TE0 FPT(3e-3) // 0.03
#define MEDIANU FPT(0.0)
#define LANDAUDAMPING FPT(0.0)//FPT(0.03) // model landau damping
#define DIFFUSION FPT(10.0) // electron thermal diffusion
#define NUTRATIO FPT(.0) // 0.0
#define VELOCITY FPT(1.2) //FPT(1.2)
#define SRCT (FPT(36000.0)) //(FPT(18000.0)*FPT(4.0))
#define SRCNOSC (FPT(8.0))
#define SRCAMP FPT(5.0e-6)
#define SWITCHONDELAY FPT(30000.0)
#define SRCX FPT(100.0)
#define SRCTFACTOR FPT(10000.0) // 300000.0 for j^2; 100.0 for E^2


/*
#define LX FPT(40000.0) //FPT(500.0)
#define LZ FPT(168000.0) //FPT(20000.0*8.0) FPT(20000.0*4.0)
#define NZ (2048)//512*10*4;
//512*10*2;//512*10*4
#define NX (1024)//128
*/

#define LX FPT(80000.0) //FPT(500.0)
#define LZ FPT(672000.0) //FPT(20000.0*8.0) FPT(20000.0*4.0)
#define NZ (8192)//512*10*4;
//512*10*2;//512*10*4
#define NX (2048)//128

#define TMAX FPT(50000.0)
#define TOOTHDEPTH FPT(100.0)
#define TOOTHWIDTH FPT(1500.0)
#define MEDIADEPTH FPT(50.0)

#define MX_1 FPT(1.0)
#define MZ_1 FPT(1.0)

#define OUTPUTFILENAME "data.dat"
#define ITER 1000

#define SQRTNUT //sqrt
