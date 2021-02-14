#include <math_constants.h>
#include <cutil_inline.h>
#include "funcall.h"
#include <iostream>
#include "geo.h"
#include <cmath>

//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION

#define BS 16 // do not use any const variables this can break loops unroll in GPU code

extern "C"
hydro2dHandler::hydro2dHandler(int dev, bool _doPhononAbsorbtion, bool _doPolarization, bool _extSource) :
								device(dev), Nx(NX), Nz(NZ), Lx(LX), Lz(LZ), ix0(0), iz0(0), dx_1(0), dz_1(0), dt(0),
								Ptz(nullptr), Ptx(nullptr), Pz(nullptr), Px(nullptr), Ez(nullptr), Ex(nullptr), By(nullptr), Phtz(nullptr), Phtx(nullptr), Phz(nullptr), Phx(nullptr), Jz(nullptr), Jx(nullptr), n(nullptr), Te(nullptr), dTe(nullptr), nlx(nullptr), nlz(nullptr), vx(nullptr), vz(nullptr), PML_Byx(nullptr), PML_Byz(nullptr), feedPtz(nullptr), feedPtx(nullptr), feedPz(nullptr), feedPx(nullptr), feedEz(nullptr), feedEx(nullptr), feedBy(nullptr), feedPML_Byx(nullptr), feedPML_Byz(nullptr), feedPhtz(nullptr), feedPhtx(nullptr), feedPhz(nullptr), feedPhx(nullptr), feedJz(nullptr), feedJx(nullptr), feedn(nullptr), feedTe(nullptr), srcx(nullptr), srct(nullptr), host_feed(nullptr), host_srct(nullptr),
								t(0), step(0), PMLimin(0), PMLimax(0), SRCi(0), shft(SHIFTGRANULARITY), blockSize(BS), PMLxmax(PMLX), PMLstrength(PMLSTRENGTH),
								mediaN0(N0), mediaTe0(TE0),diffusion(DIFFUSION), mediaNu(MEDIANU), landauDamping(LANDAUDAMPING), mz_1(MZ_1), mx_1(MX_1), toothDepth(TOOTHDEPTH), toothWidth(TOOTHWIDTH), mediaDepth(MEDIADEPTH), a_cell(nullptr), media_bound_w2(BOUND_W2), media_bound_gamma(BOUND_GAMMA), media_bound_beta(BOUND_BETA), media_phonon_omega(PHONON_OMEGA), media_phonon_phw(PHONON_PHW), media_phonon_beta(PHONON_BETA), srcTfactor(SRCTFACTOR), srcAmp(SRCAMP), srcVelocity(VELOCITY), srcT(SRCT), srcNosc(SRCNOSC), srcX(SRCX), switchOnDelay(SWITCHONDELAY),
								doPolarization(_doPolarization),  doPhononAbsorbtion(_doPhononAbsorbtion), extSource(_extSource),
								toothDir(false), flip(false) {}

__device__ FL_DBL
PML(FL_DBL x, FL_DBL PMLx, FL_DBL PMLstrength)
{
	return FPT(1.0)-((x<=PMLx)?(PMLx-x)/PMLx*PMLstrength:FPT(0.0));
}

//Polarization step needed for Bismuth
__global__ void
simpleP_Kernel(int Nx, int Nz, FL_DBL* Ez, FL_DBL* Ex, FL_DBL* Ptz, FL_DBL* Ptx, FL_DBL* Pz, FL_DBL* Px, FL_DBL* n, FL_DBL dt, FL_DBL w2_bound /*media_bound_w2*/, FL_DBL gamma_bound /*media_bound_gamma*/, FL_DBL beta_bound /*media_bound_beta*/)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
		const int im1 = i == 0 ? i : i - 1;
		const int jm1 = j == 0 ? j : j - 1;
				
		const FL_DBL n_z = (n[i*Nx+j] + n[im1*Nx+j]) * FPT(0.5);
		//n_z = n_z > FPT(0.0) ? n_z : FPT(0.0);
		const FL_DBL n_x = (n[i*Nx+j] + n[i*Nx+jm1]) * FPT(0.5);
		//n_x = n_x > FPT(0.0) ? n_x : FPT(0.0);
		
		Ptz[i*Nx+j] = n_z > FPT(0.0) ? (Ptz[i*Nx+j]/dt-w2_bound*Pz[i*Nx+j]-0.5*gamma_bound*Ptz[i*Nx+j]+beta_bound*Ez[i*Nx+j])*dt/(1.0+0.5*gamma_bound*dt) : FPT(0.0);
		Ptx[i*Nx+j] = n_x > FPT(0.0) ? (Ptx[i*Nx+j]/dt-w2_bound*Px[i*Nx+j]-0.5*gamma_bound*Ptx[i*Nx+j]+beta_bound*Ex[i*Nx+j])*dt/(1.0+0.5*gamma_bound*dt) : FPT(0.0);
		Pz[i*Nx+j] += Ptz[i*Nx+j]*dt;
		Px[i*Nx+j] += Ptx[i*Nx+j]*dt;
		Ez[i*Nx+j] += -Ptz[i*Nx+j]*dt;
		Ex[i*Nx+j] += -Ptx[i*Nx+j]*dt;
	}
}

__global__ void
simpleE_Kernel(int Nx, int Nz, FL_DBL* Ez, FL_DBL* Ex, FL_DBL* By, FL_DBL* Jz, FL_DBL* Jx, int ix0, int iz0, FL_DBL dx, FL_DBL dx_1, FL_DBL dz, FL_DBL dz_1, FL_DBL dt, FL_DBL* srct, int srci, int PMLimin, int PMLimax, FL_DBL pml, FL_DBL pmlx)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
		Ez[i*Nx+j] = Ez[i*Nx+j] * (( j+ix0<PMLimin || j+ix0>PMLimax ) ? PML(FL_DBL((j+ix0<PMLimin)?j+ix0:Nx-1-j)*dx, pmlx, pml):FPT(1.0)) + (-(((j+1==Nx)?Ez[i*Nx+j]:By[i*Nx+j+1])-By[i*Nx+j]) * dx_1 - Jz[i*Nx+j] + /*source*/(srct && j + ix0 == srci ? srct[i] : FPT(0.0)))*dt;

// no bottom pml:
//			Ez[i*Nx+j]=Ez[i*Nx+j]*(( j+ix0<PMLimin ) ? PML(FL_DBL((j+ix0<PMLimin)?j+ix0:Nx-1-j)*dx, pmlx, pml):FPT(1.0))+(-(((j+1==Nx)?Ez[i*Nx+j]:By[i*Nx+j+1])-By[i*Nx+j]) * dx_1 - Jz[i*Nx+j] + /*source*/(srct && j + ix0 == srci ? srct[i] : FPT(0.0)))*dt;

		Ex[i*Nx+j] += ((((i+1==Nz)?FPT(0.0):By[(i+1)*Nx+j]) - By[i*Nx+j]) * dz_1 - Jx[i*Nx+j]) * dt;
	}
}

__global__ void
simpleB_Kernel(int Nx, int Nz, FL_DBL* Ez, FL_DBL* Ex, FL_DBL* By, int ix0, int iz0, FL_DBL dx, FL_DBL dx_1, FL_DBL dz, FL_DBL dz_1, FL_DBL dt, int PMLimin, int PMLimax, FL_DBL pml, FL_DBL pmlx, FL_DBL* PML_Byz, FL_DBL* PML_Byx)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
//				By[i*Nx+j]+=(-(Ez[i*Nx+j]-((j==0)?FPT(0.0):Ez[i*Nx+j-1]))*dx_1+(Ex[i*Nx+j]-((i==0)?FPT(0.0):Ex[(i-1)*Nx+j]))*dz_1)*dt;

		if(j+ix0<PMLimin)
		{
			PML_Byz[i*Nx+j]=PML_Byz[i*Nx+j]+(Ex[i*Nx+j]-((i==0)?FPT(0.0):Ex[(i-1)*Nx+j]))*dz_1*dt;
			PML_Byx[i*Nx+j]=PML_Byx[i*Nx+j]*PML(FL_DBL(j+ix0)*dx, pmlx, pml)-(Ez[i*Nx+j]-((j+ix0==0)?-By[i*Nx+j]:Ez[i*Nx+j-1]))*dx_1*dt;
			By[i*Nx+j]=PML_Byz[i*Nx+j]+PML_Byx[i*Nx+j];
		}
		else if(j+ix0>PMLimax)
		{
			PML_Byz[i*Nx+j]=PML_Byz[i*Nx+j]+(Ex[i*Nx+j]-((i==0)?FPT(0.0):Ex[(i-1)*Nx+j]))*dz_1*dt;
			PML_Byx[i*Nx+j]=PML_Byx[i*Nx+j]*PML(FL_DBL(Nx-1-j)*dx, pmlx, pml)-(Ez[i*Nx+j]-Ez[i*Nx+j-1])*dx_1*dt;
			By[i*Nx+j]=PML_Byz[i*Nx+j]+PML_Byx[i*Nx+j];
		}
		else
			By[i*Nx+j]+=(-(Ez[i*Nx+j]-Ez[i*Nx+j-1])*dx_1+(Ex[i*Nx+j]-((i==0)?FPT(0.0):Ex[(i-1)*Nx+j]))*dz_1)*dt;
	}
}

__global__ void
simpleV_Kernel(int Nx, int Nz, FL_DBL* vz, FL_DBL* vx, FL_DBL* Jz, FL_DBL* Jx, FL_DBL* n)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
		int im1 = (i==0)?i:i-1;
		int jm1 = (j==0)?j:j-1;
		FL_DBL n_z = (n[i*Nx+j] + n[im1*Nx+j]) * FPT(0.5);
		FL_DBL n_x = (n[i*Nx+j] + n[i*Nx+jm1]) * FPT(0.5);
				
		vz[i*Nx+j] = (n_z > FPT(0.0)) ? Jz[i*Nx+j] / n_z : FPT(0.0);
		vx[i*Nx+j] = (n_x > FPT(0.0)) ? Jx[i*Nx+j] / n_x : FPT(0.0);
	}
}

__global__ void
simpleNl_Kernel(int Nx, int Nz, FL_DBL* nlz, FL_DBL* nlx, FL_DBL* vz, FL_DBL* vx, FL_DBL* Jz, FL_DBL* Jx, FL_DBL* By, FL_DBL dz_1, FL_DBL dx_1, FL_DBL mz_1, FL_DBL mx_1)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
		int ip1=(i==Nz-1)?i:i+1;
		int im1=(i==0)?i:i-1;
		int jp1=(j==Nx-1)?j:j+1;
		int jm1=(j==0)?j:j-1;

		nlz[i*Nx+j]=
					HYDROZSRC*(Jz[i*Nx+j]*(vz[ip1*Nx+j]-vz[im1*Nx+j])*FPT(0.5)*dz_1+FPT(0.25)*(Jx[i*Nx+j]+Jx[i*Nx+jp1]+Jx[im1*Nx+j]+Jx[im1*Nx+jp1])*(vz[i*Nx+jp1]-vz[i*Nx+jm1])*FPT(0.5)*dx_1)
					+MAINSRC*vz[i*Nx+j]*((Jx[i*Nx+jp1]+Jx[im1*Nx+jp1]-Jx[i*Nx+j]-Jx[im1*Nx+j])*FPT(0.5)*dx_1+(Jz[ip1*Nx+j]-Jz[im1*Nx+j])*FPT(0.5)*dz_1)
					+(Jx[i*Nx+j]+Jx[im1*Nx+j]+Jx[i*Nx+jp1]+Jx[im1*Nx+jp1])*(By[i*Nx+j]+By[i*Nx+jp1])*FPT(0.125)*mz_1*HSRC;


		nlx[i*Nx+j]=
					HYDROXSRC*(FPT(0.25)*(Jz[i*Nx+j]+Jz[i*Nx+jm1]+Jz[ip1*Nx+j]+Jz[ip1*Nx+jm1])*(vx[ip1*Nx+j]-vx[im1*Nx+j])*FPT(0.5)*dz_1+Jx[i*Nx+j]*(vx[i*Nx+jp1]-vx[i*Nx+jm1])*FPT(0.5)*dx_1)
					+PRESSURESRC*vx[i*Nx+j]*((Jx[i*Nx+jp1]-Jx[i*Nx+jm1])*dx_1*FPT(0.5)+(Jz[ip1*Nx+j]+Jz[ip1*Nx+jm1]-Jz[i*Nx+j]-Jz[i*Nx+jm1])*dz_1*FPT(0.5))
					-(Jz[i*Nx+j]+Jz[ip1*Nx+j]+Jz[i*Nx+jm1]+Jz[ip1*Nx+jm1])*(By[ip1*Nx+j]+By[i*Nx+j])*FPT(0.125)*mx_1*HSRC;

	}
}

__global__ void
simpleN_Kernel(int Nx, int Nz, FL_DBL* n, FL_DBL* Jx, FL_DBL* Jz, FL_DBL dz_1, FL_DBL dx_1, FL_DBL dt)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
			int ip1=(i==Nz-1)?i:i+1;
			int jp1=(j==Nx-1)?j:j+1;
			int mask=(n[i*Nx+j]==FPT(0.0))?0:1;
			n[i*Nx+j]+=-((Jx[i*Nx+jp1]-Jx[i*Nx+j])*dx_1+(Jz[ip1*Nx+j]-Jz[i*Nx+j])*dz_1)*dt;
			n[i*Nx+j]=(mask==1)?n[i*Nx+j]:FPT(0.0); // grants n0>0 and 0 outside plasma area
	}
}

__device__ inline FL_DBL sqr(FL_DBL x)
{return x*x;}

__global__ void
simple_dTe_Kernel(int Nx, int Nz, FL_DBL D, FL_DBL* dTe, FL_DBL* Te, FL_DBL* n, FL_DBL* Ex, FL_DBL* Ez, FL_DBL dz_1, FL_DBL dx_1, FL_DBL* srct, FL_DBL* srcx, FL_DBL SRCFACTOR)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
			bool vac = n[i*Nx+j] == FPT(0.0);
			int ip1=(i==Nz-1)?i:i+1;
			int jp1=(j==Nx-1)?j:j+1;
			
			FL_DBL src;
			if(!srct)
				src=(!vac) ? sqr((Ex[i*Nx+j]+Ex[i*Nx+jp1])*FPT(0.5))*FPT(.0) + sqr((Ez[i*Nx+j]+Ez[ip1*Nx+j])*FPT(0.5)) : FPT(0.0);	
			else				
				src = (!vac) ? srct[i] * srcx[j] : FPT(0.0);
			
			ip1 = (i==Nz-1 || n[(i+1)*Nx+j] == FPT(0.0)) ? i : i + 1; //von neumann bc
			jp1 = (j==Nx-1 || n[i*Nx+(j+1)] == FPT(0.0)) ? j : j + 1; //von neumann bc
			int im1 = (i==0 || n[(i-1)*Nx+j] == FPT(0.0)) ? i : i - 1; //von neumann bc
			int jm1 = (j==0 || n[i*Nx+(j-1)] == FPT(0.0)) ? j : j - 1; //von neumann bc
			
			//#define SRCFACTOR FPT(100.0)
			dTe[i*Nx+j] = (!vac) ? SRCFACTOR * src + D * ((Te[i*Nx+jp1]+Te[i*Nx+jm1]-FPT(2.0)*Te[i*Nx+j])*dx_1*dx_1+(Te[ip1*Nx+j]+Te[im1*Nx+j]-FPT(2.0)*Te[i*Nx+j])*dz_1*dz_1) : FPT(0.0);
//                                 
//               0                 
//              000                
//             00 00               
//            00   00              
//           00 /0\ 00             
//          00  080  00            
//         00   080   00           
//        00    \8/    00          
//       00      0      00         
//      00       !       00        
//     00                 00       
//    00        /8\        00      
//   00         \9/         00     
//  00                       00    
// 00000000000000000000000000000   
//                                 
	}
}


__global__ void
simple_Te_Kernel(int Nx, int Nz, FL_DBL* dTe, FL_DBL* Te, FL_DBL dt)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	Te[i*Nx+j] += dTe[i*Nx+j] * dt;
}

__global__ void
simpleJ_Kernel(int Nx, int Nz, FL_DBL* Jz, FL_DBL* Jx, FL_DBL* Ez, FL_DBL* Ex, FL_DBL* nlz, FL_DBL* nlx, FL_DBL* n, FL_DBL* Te, FL_DBL dx_1, FL_DBL dz_1, FL_DBL dt, FL_DBL Nu, FL_DBL landauDamping, FL_DBL* phv_z,  FL_DBL* phv_x, FL_DBL* ph_z, FL_DBL* ph_x, FL_DBL mz_1, FL_DBL mx_1, FL_DBL omegaOpt /*media_phonon_omega*/, FL_DBL PHW /*media_phonon_phw*/, FL_DBL OPT_PH /*media_phonon_beta*/, bool do_phonons)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
		int ip1 = (i == Nz - 1) ? i : i + 1;
		int im1 = (i == 0     ) ? i : i - 1;
		int jp1 = (j == Nx - 1) ? j : j + 1;
		int jm1 = (j == 0     ) ? j : j - 1;

//PROBLEMS APPEARS IF ADDING THIS CODE ON gpu# > 0

		FL_DBL n_z = (n[i*Nx+j] + n[im1*Nx+j]) * FPT(0.5);
		n_z = n_z > FPT(0.0) ? n_z : FPT(0.0);
		FL_DBL n_x = (n[i*Nx+j] + n[i*Nx+jm1]) * FPT(0.5);
		n_x = n_x > FPT(0.0) ? n_x : FPT(0.0);

		FL_DBL Nu_z = (n_z > FPT(0.0)) ? Nu : FPT(0.0);
		FL_DBL Nu_x = (n_x > FPT(0.0)) ? Nu : FPT(0.0);
		
		FL_DBL Jz_t = landauDamping*(Jz[i*Nx+jp1]+Jz[i*Nx+jm1]-FPT(2.0)*Jz[i*Nx+j])*dx_1*dx_1-Jz[i*Nx+j] * Nu_z + (mz_1 * n_z * Ez[i*Nx+j] - IFNONLIN*nlz[i*Nx+j] - mz_1 * (n[i*Nx+j] * Te[i*Nx+j] - n[im1*Nx+j] * Te[im1*Nx+j])*dz_1);
		//FL_DBL Jz_t = n_z * Ez[i*Nx+j] - (n[i*Nx+j] * Te[i*Nx+j] - n[im1*Nx+j] * Te[im1*Nx+j]) * dz_1
		
		FL_DBL Jx_t = landauDamping*(Jx[ip1*Nx+j]+Jx[im1*Nx+j]-FPT(2.0)*Jx[i*Nx+j])*dz_1*dz_1 - Jx[i*Nx+j] * Nu_x + (mx_1 * n_x * Ex[i*Nx+j] - IFNONLIN*nlx[i*Nx+j] - mx_1 * (n[i*Nx+j] * Te[i*Nx+j] - n[i*Nx+jm1] * Te[i*Nx+jm1])*dx_1);
		//FL_DBL Jx_t = n_x * Ex[i*Nx+j] - (n[i*Nx+j] * Te[i*Nx+j] - n[i*Nx+jm1] * Te[i*Nx+jm1]) * dx_1;

		if(do_phonons)
		{
			Jz[i*Nx+j] += (Jz_t - omegaOpt*omegaOpt*ph_z[i*Nx+j]) * dt;
			Jx[i*Nx+j] += (Jx_t - omegaOpt*omegaOpt*ph_x[i*Nx+j]) * dt;
		}
		else
		{
			Jz[i*Nx+j] += Jz_t * dt;
			Jx[i*Nx+j] += Jx_t * dt;
		}
		
		//put zeros at plasma surface
		if(n[i*Nx+j]==FPT(0.0)||n[i*Nx+jm1]==FPT(0.0)) Jx[i*Nx+j]=FPT(0.0);
		if(n[i*Nx+j]==FPT(0.0)||n[im1*Nx+j]==FPT(0.0)) Jz[i*Nx+j]=FPT(0.0);

		if(do_phonons)
		{
			phv_z[i*Nx+j]+=(-omegaOpt*omegaOpt*ph_z[i*Nx+j]-PHW*omegaOpt*phv_z[i*Nx+j]+OPT_PH*(Jz_t-omegaOpt*omegaOpt*ph_z[i*Nx+j]))*dt;
			phv_x[i*Nx+j]+=(-omegaOpt*omegaOpt*ph_x[i*Nx+j]-PHW*omegaOpt*phv_x[i*Nx+j]+OPT_PH*(Jx_t-omegaOpt*omegaOpt*ph_x[i*Nx+j]))*dt;

			ph_z[i*Nx+j]+=phv_z[i*Nx+j]*dt;
			ph_x[i*Nx+j]+=phv_x[i*Nx+j]*dt;		
		}
	}
}

__global__ void
shift_Kernel(int Nx, int Nz, FL_DBL* from, FL_DBL* to, int shift, FL_DBL* feed) // NOTE: shift is counted in blocks
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
	{
		int iold=i+shift*blockDim.y;
		if(iold<Nz)
			to[i*Nx+j]=from[iold*Nx+j];
		else
			to[i*Nx+j]=(!feed)?FPT(0.0):feed[(iold-Nz)*Nx+j];
	}
}

__global__ void
copy_Kernel(int Nx, int Nz, FL_DBL* from, FL_DBL* to) // NOTE: shift is counted in blocks
{
		const int j=threadIdx.x+blockIdx.x*blockDim.x;
		const int i=threadIdx.y+blockIdx.y*blockDim.y;
		if((i>=0)&(i<Nz)&(j>=0)&(j<Nx))
			to[i*Nx+j]=from[i*Nx+j];
}


#ifdef BS
	#undef BS
#endif

// in host code use hydro2dHadnler::blockSize field

//-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION
//-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION
//-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION
//-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION
//-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION
//-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION
//-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION-----HOST-SECTION

inline void shift(hydro2dHandler* hH) // shft is counted in blocks
{
	dim3 grid(-floor(-FL_DBL(hH->Nx)/FL_DBL(hH->blockSize)),-floor(-FL_DBL(hH->Nz)/FL_DBL(hH->blockSize)),1);
	//printf("grid=%d %d %d\n", grid.x, grid.y, grid.z);
	dim3 threads(hH->blockSize,hH->blockSize,1);
	//shift big data

	const bool doPhononAbsorbtion = hH->doPhononAbsorbtion;
	const bool doPolarization = hH->doPolarization;

	if(doPolarization)
	{
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Pz, hH->vz, hH->shft, hH->feedEz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Px, hH->vx, hH->shft, hH->feedEx);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ptz, hH->Pz, hH->shft, hH->feedEz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ptx, hH->Px, hH->shft, hH->feedEx);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->Ptz, hH->shft, hH->feedEz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ex, hH->Ptx, hH->shft, hH->feedEx);
		cudaThreadSynchronize();
	}
	else
	{
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->vz, hH->shft, hH->feedEz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ex, hH->vx, hH->shft, hH->feedEx);
		cudaThreadSynchronize();
	}
	if(doPhononAbsorbtion)
	{
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Phtz, hH->Ez, hH->shft, hH->feedJz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Phtx, hH->Ex, hH->shft, hH->feedJx);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Phz, hH->Phtz, hH->shft, hH->feedJz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Phx, hH->Phtx, hH->shft, hH->feedJx);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jz, hH->Phz, hH->shft, hH->feedJz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jx, hH->Phx, hH->shft, hH->feedJx);
		cudaThreadSynchronize();
	}
	else
	{
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jz, hH->Ez, hH->shft, hH->feedJz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jx, hH->Ex, hH->shft, hH->feedJx);
		cudaThreadSynchronize();
	}
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->By, hH->Jz, hH->shft, hH->feedBy);
	cudaThreadSynchronize();
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->n , hH->Jx, hH->shft, hH->feedn);
	cudaThreadSynchronize();

//shift in PML layer

	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->PML_Byz , hH->By, hH->shft, hH->feedPML_Byz);
	cudaThreadSynchronize();
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->PML_Byx , hH->n, hH->shft, hH->feedPML_Byx);
	cudaThreadSynchronize();
//copy
	copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->By, hH->PML_Byz);
	cudaThreadSynchronize();
	copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->n , hH->PML_Byx);
	cudaThreadSynchronize();

//shift T
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Te , hH->By, hH->shft, hH->feedTe);
	cudaThreadSynchronize();

//roll pointers:
	FL_DBL* new_Ez = 0;
	FL_DBL* new_Ex = 0;
	FL_DBL* new_Pz = 0;
	FL_DBL* new_Px = 0;
	FL_DBL* new_Ptz = 0;
	FL_DBL* new_Ptx = 0;
	
	if(doPolarization)
	{
		new_Pz=hH->vz;
		new_Px=hH->vx;
		new_Ptz=hH->Pz;
		new_Ptx=hH->Px;
		new_Ez=hH->Ptz;
		new_Ex=hH->Ptx;
	}
	else
	{
		new_Ez=hH->vz;
		new_Ex=hH->vx;
	}
	
	FL_DBL* new_Jz = 0;
	FL_DBL* new_Jx = 0;
	FL_DBL* new_Phtz = 0;
	FL_DBL* new_Phtx = 0;
	FL_DBL* new_Phz = 0;
	FL_DBL* new_Phx = 0;

	if(doPhononAbsorbtion)
	{
		new_Phtz=hH->Ez;
		new_Phtx=hH->Ex;
		new_Phz=hH->Phtz;
		new_Phx=hH->Phtx;
		new_Jz=hH->Phz;
		new_Jx=hH->Phx;
	}
	else
	{
		new_Jz=hH->Ez;
		new_Jx=hH->Ex;
	}
	FL_DBL* new_By=hH->Jz;
	FL_DBL* new_n =hH->Jx;
	FL_DBL* new_Te =hH->By;
	FL_DBL* new_vx=hH->Te;
	FL_DBL* new_vz=hH->n;

	if(doPolarization)
	{
		hH->Ptz=new_Ptz;
		hH->Ptx=new_Ptx;
		hH->Pz=new_Pz;
		hH->Px=new_Px;
	}
	hH->Ez=new_Ez;
	hH->Ex=new_Ex;
	if(doPhononAbsorbtion)
	{
		hH->Phtz=new_Phtz;
		hH->Phtx=new_Phtx;
		hH->Phz=new_Phz;
		hH->Phx=new_Phx;
	}
	hH->Jz=new_Jz;
	hH->Jx=new_Jx;
	hH->By=new_By;
	hH->n =new_n ;
	hH->vz=new_vz;
	hH->vx=new_vx;
	hH->Te=new_Te;

	hH->iz0 += hH->shft * hH->blockSize;
}

extern "C"
int simpleGPUstep(hydro2dHandler* hH)
{
	if(hH->t == 0)
	{
		std::cout << "Hi this is first step for device " << hH->device << std::endl;
		std::cout << "grid: Lx=" << hH->Lx << ", Lz=" << hH->Lz << ", Nx=" << hH->Nx << ", Nz=" << hH->Nz << ", dx_1=" << hH->dx_1 << " (dx=" << FPT(1.0)/hH->dx_1 << "), dz_1=" << hH->dz_1 << " (dz=" << FPT(1.0)/hH->dz_1 << "), dt=" << hH->dt << ",\ngeometry: toothWith=" << hH->toothWidth<< ", toothDepth=" << hH->toothDepth << ", mediaDepth=" << hH->mediaDepth << ", flip=" << (hH->flip?"true":"false") << ", toothDirection=" << (hH->toothDir?"true":"false") << ",\nmedia: n0=" << hH->mediaN0 << ", T0=" << hH->mediaTe0 << ", diffusion=" << hH->diffusion << "," << " mz_1=" << hH->mz_1 << "," << " mx_1=" << hH->mx_1 << ",";
		if (hH->doPolarization)
		std::cout << "\n       bound_w2=" << hH->media_bound_w2 << " bound_beta=" << hH->media_bound_beta << ", bound_gamma=" << hH->media_bound_gamma << std::endl;
		if (hH->doPhononAbsorbtion)
		std::cout << "\n       phOmega=" << hH->media_phonon_omega << ", phWidth=" << hH->media_phonon_phw << ", phBeta=" << hH->media_phonon_beta << std::endl;
		if (!hH->doPolarization && !hH->doPolarization)
			std::cout << std::endl;
		if (hH->extSource)
			std::cout << "         therm source is external, srcX=" << hH->srcX;
		else
			std::cout << "         source is electromagnetic";
		std::cout << ", srcAmp=" << hH->srcAmp << ", srcT=" << hH->srcT << ", srcNosc=" << hH->srcNosc << ", switchOnDelay=" << hH->switchOnDelay << ", srcTfactor=" << hH->srcTfactor << std::endl;
		std::cout << "         velocity=" << hH->srcVelocity << std::endl;
	}
	
	if(cudaSetDevice(hH->device))
	{
		printf("ERROR Initializing device #%d\n",hH->device);
		return -1;
	}
	dim3 grid(-floor(-FL_DBL(hH->Nx)/FL_DBL(hH->blockSize)),-floor(-FL_DBL(hH->Nz)/FL_DBL(hH->blockSize)),1);
	//printf("grid=%d %d %d\n", grid.x, grid.y, grid.z);
	dim3 threads(hH->blockSize,hH->blockSize,1);

	FL_DBL delay = hH->switchOnDelay;
	FL_DBL srcZ = hH->srcT * hH->srcVelocity;
	for(int i=0; i != hH->Nz; i++)
	{
		double zzd = (double)(i + hH->iz0 - hH->Nz + hH->blockSize + 1) / (double) hH->dz_1 + (double)(srcZ * 0.5) - hH->t * (double) hH->srcVelocity;
		FL_DBL zz = (FL_DBL(zzd));
		if(hH->extSource)
			hH->host_srct[i] = hH->srcAmp * (FPT(0.5) + FPT(0.5) * tanh((hH->t - delay) / delay)) * ((fabs(zz / srcZ * FPT(2.0)) < FPT(1.0)) ? (FPT(0.5) + FPT(0.5) * cos(zz / srcZ * FPT(2.0) * M_PI)) : FPT(0.0));// * (1.0 + sin(zz / srcZ * hH->srcNosc * M_PI * FPT(2.0)));
		else
			hH->host_srct[i] = hH->srcAmp * (FPT(0.5) + FPT(0.5) * tanh((hH->t - delay) / delay)) * ((fabs(zz / srcZ * FPT(2.0)) < FPT(1.0)) ? (FPT(0.5) + FPT(0.5) * cos(zz / srcZ * FPT(2.0) * M_PI)) : FPT(0.0)) *        sin(zz / srcZ * hH->srcNosc * M_PI * FPT(2.0));
	}
	dev_h2d(hH, hH->host_srct, hH->srct, hH->Nz);

	if(hH->doPolarization)
	{
		simpleP_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->Ex, hH->Ptz, hH->Ptx, hH->Pz, hH->Px, hH->n, hH->dt, hH->media_bound_gamma, hH->media_bound_w2, hH->media_bound_beta);
		cudaThreadSynchronize();
	}
	simpleE_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->Ex, hH->By, hH->Jz, hH->Jx, hH->ix0, hH->iz0, FPT(1.0)/hH->dx_1, hH->dx_1, 1.0/hH->dz_1, hH->dz_1, hH->dt, (hH->extSource ? 0 : hH->srct), hH->SRCi, hH->PMLimin, hH->PMLimax, hH->PMLstrength, hH->PMLxmax);
	cudaThreadSynchronize();
	
	simpleB_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->Ex, hH->By, hH->ix0, hH->iz0, 1.0/ hH->dx_1, hH->dx_1, 1.0/ hH->dz_1, hH->dz_1, hH->dt, hH->PMLimin, hH->PMLimax, hH->PMLstrength, hH->PMLxmax, hH->PML_Byz, hH->PML_Byx);
	cudaThreadSynchronize();
	
	simpleV_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->vz, hH->vx, hH->Jz, hH->Jx, hH->n);
	cudaThreadSynchronize();
	
	simpleNl_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->nlz, hH->nlx, hH->vz, hH->vx, hH->Jz, hH->Jx, hH->By, hH->dz_1, hH->dx_1, hH->mz_1,  hH->mx_1);
	cudaThreadSynchronize();
	
	simpleN_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->n, hH->Jx, hH->Jz, hH->dz_1, hH->dx_1, hH->dt);
	cudaThreadSynchronize();
	
	simple_dTe_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->diffusion, hH->dTe, hH->Te, hH->n, hH->Ex, hH->Ez, hH->dz_1, hH->dx_1, hH->extSource ? hH->srct : 0, hH->srcx, hH->srcTfactor);
	cudaThreadSynchronize();
	
	simple_Te_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->dTe, hH->Te, hH->dt);
	cudaThreadSynchronize();
	
	simpleJ_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jz, hH->Jx, hH->Ez, hH->Ex, hH->nlz, hH->nlx, hH->n, hH->Te, hH->dx_1, hH->dz_1, hH->dt, hH->mediaNu, hH->landauDamping, hH->Phtz, hH->Phtx, hH->Phz, hH->Phx, hH->mz_1,  hH->mx_1, hH->media_phonon_omega, hH->media_phonon_phw, hH->media_phonon_beta, hH->doPhononAbsorbtion);
	cudaThreadSynchronize();
	
	hH->t += (double) hH->dt;
	hH->step++;

	if(hH->t * (double) hH->srcVelocity - (double) hH->iz0 / (double) hH->dz_1 > (double) hH->blockSize / (double) hH->dz_1)
	{
		//manage feed
		FL_DBL* feed = hH->host_feed;//new FL_DBL[hH->Nx * (hH->shft * hH->blockSize)];
		
		for(int i=0;i<hH->Nx * (hH->shft * hH->blockSize);i++) feed[i] = FPT(0.0);

		dev_h2d(hH, feed, hH->feedPML_Byz, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedPML_Byx, hH->Nx * (hH->shft * hH->blockSize));
		if(hH->doPolarization)
		{
			dev_h2d(hH, feed, hH->feedPtz, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPtx, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPz, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPx, hH->Nx * (hH->shft * hH->blockSize));
		}
		dev_h2d(hH, feed, hH->feedEz, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedEx, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedBy, hH->Nx * (hH->shft * hH->blockSize));

		if(hH->doPhononAbsorbtion)
		{
			dev_h2d(hH, feed, hH->feedPhtz, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPhtx, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPhz, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPhx, hH->Nx * (hH->shft * hH->blockSize));
		}
		dev_h2d(hH, feed, hH->feedJz, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedJx, hH->Nx * (hH->shft * hH->blockSize));
		// add grid for n
		if(hH->a_cell)
		{
			//std::cout << "feeding from cell, step is " << hH->step << std::endl;
			int ins(0);
			int outs(0);
			for(int j = 0; j < hH->Nx * hH->shft * hH->blockSize; j++)
				feed[j] = FL_DBL(0.0);
			const cell* c = hH->a_cell;
			const int imin = std::floor(c->get_zmin() / c->get_lz()); 
			const int imax = std::ceil(c->get_zmax() / c->get_lz()); 
			const double dz = FL_DBL(1.0) / hH->dz_1;
			const double dx = FL_DBL(1.0) / hH->dx_1;
			for(int i = 0; i<hH->shft * hH->blockSize; i++)
			{
				int ireal = hH->iz0 + i;
				FL_DBL z = FL_DBL(ireal) * dz;
				for(int cell_id=imin; cell_id!=imax; cell_id++)
				{
					FL_DBL z_cell_local = std::fmod(z, c->get_lz());
					z_cell_local += FL_DBL(cell_id) * c->get_lz();
					for(int j=0;j<hH->Nx;j++)
					{
						FL_DBL x = FL_DBL(j) * dx - hH->Lx * FL_DBL(0.5);
						//std::cout << "calling inside(" << z_cell_local << ", " << x << ")" << std::endl;
						bool res = c->inside(z_cell_local, x);
						ins += res ? 1 : 0;
						outs += res ? 0 : 1;
						//std::cout << (res ? "got yes" : "got no") << std::endl;
						feed[i * hH->Nx + j] = res ? hH->mediaN0 : feed[i * hH->Nx + j];
					}
				}
			}
			//std::cout << "done, ins = " << ins << ", outs = " << outs << std::endl;
		}
		else
		{
			for(int i=0;i<hH->shft * hH->blockSize;i++)
			{
				int ireal=hH->iz0+i;
				FL_DBL z = FL_DBL(ireal) / hH->dz_1;
				FL_DBL levelL = hH->Lx/2.0;// - hH->Lx/4.0;
				FL_DBL levelH = levelL + hH->mediaDepth;

				int toothWidthi = hH->toothWidth * hH->dz_1;
				int toothDepthi = hH->toothDepth * hH->dx_1;
				int profile = FL_DBL(((hH->toothDir) ? ireal : 100000000 - ireal) % toothWidthi) / FL_DBL(toothWidthi) * FL_DBL(toothDepthi);

				if(hH->flip)
				{
					for(int j=0;j<hH->Nx;j++)
						feed[i * hH->Nx + j]=((j < levelH * hH->dx_1 + profile) && (j > levelL * hH->dx_1)) ? hH->mediaN0 : FPT(0.0);
				}
				else
				{
					for(int j=0;j<hH->Nx;j++)
						feed[i * hH->Nx + j]=((j < levelH * hH->dx_1 + toothDepthi) && (j > levelL * hH->dx_1 + profile)) ? hH->mediaN0 : FPT(0.0);
				}
			}
		}
		dev_h2d(hH, feed, hH->feedn, hH->Nx * (hH->shft * hH->blockSize));

		for(int i=0;i<hH->shft * hH->blockSize;i++)
		{
			for(int j=0;j<hH->Nx;j++)
				feed[i * hH->Nx + j]=(feed[i * hH->Nx + j] != FPT(0.0))? hH->mediaTe0 : FPT(0.0);
		}
		dev_h2d(hH, feed, hH->feedTe, hH->Nx * (hH->shft * hH->blockSize));
		shift(hH);
	}
	return 0;
}


extern "C"
static void alloc_main_fields(hydro2dHandler* hH)
{
	printf("allocating fields of size %dx%d\n", hH->Nz, hH->Nx);
		//alloc main fields
	if(hH->doPolarization)
	{
		printf("allocating polarization arrays\n");
		dev_alloc(hH, (void**)&(hH->Ptz), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Ptx), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Pz), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Px), hH->Nx * hH->Nz);
	}
	
	dev_alloc(hH, (void**)&(hH->Ez), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Ex), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->By), hH->Nx * hH->Nz);

	if(hH->doPhononAbsorbtion)
	{
		printf("allocating phonon absorbtion arrays\n");
		dev_alloc(hH, (void**)&(hH->Phtz), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Phtx), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Phz), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Phx), hH->Nx * hH->Nz);
	}

	dev_alloc(hH, (void**)&(hH->Jx), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Jz), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->n), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Te), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->dTe), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->vx), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->vz), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->nlx), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->nlz), hH->Nx * hH->Nz);
	
	//alloc additional fields for PML conditions
	dev_alloc(hH, (void**)&(hH->PML_Byx), hH->Nz * hH->Nx);
	dev_alloc(hH, (void**)&(hH->PML_Byz), hH->Nz * hH->Nx);
	
	//alloc feed arrays
	if(hH->doPolarization)
	{
		dev_alloc(hH, (void**)&(hH->feedPtz), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPtx), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPz), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPx), hH->Nx * (hH->shft * hH->blockSize));
	}
	
	dev_alloc(hH, (void**)&(hH->feedEz), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedEx), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedBy), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedPML_Byz), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedPML_Byx), hH->Nx * (hH->shft * hH->blockSize));
	
	if(hH->doPhononAbsorbtion)
	{
		dev_alloc(hH, (void**)&(hH->feedPhtz), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPhtx), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPhz), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPhx), hH->Nx * (hH->shft * hH->blockSize));
	}
	
	dev_alloc(hH, (void**)&(hH->feedJz), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedJx), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedn), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedTe), hH->Nx * (hH->shft * hH->blockSize));
	//srcx array to speed up calculations
	dev_alloc(hH, (void**)&(hH->srcx), hH->Nx);
	dev_alloc(hH, (void**)&(hH->srct), hH->Nz);
	//hH->srct = 0;
	//allocate some host arrays for external source
	hH->host_feed = new FL_DBL[hH->Nx * (hH->shft * hH->blockSize)];
	hH->host_srct = new FL_DBL[hH->Nz];	
	printf("\tdone\n");
}

extern "C"
void simpleGPUinit(hydro2dHandler* hH) // this return handler for operation
{
	if(cudaSetDevice(hH->device))
	{
		printf("ERROR Initializing device #%d\n",hH->device);
		return;
	}
	hH->PMLimin = -floor(- hH->PMLxmax * hH->dx_1);
	hH->PMLimin = ((hH->PMLimin / hH->blockSize)+(((hH->PMLimin % hH->blockSize)==0)?0:1)) * hH->blockSize;
	hH->PMLimax = hH->Nx-hH->PMLimin;
	hH->SRCi = hH->Nx/2-hH->Nx/16;//hH->PMLimin+1;
	alloc_main_fields(hH);
	// clear helpers
	FL_DBL* zero=new FL_DBL[hH->Nx * hH->Nz];
	for(int i=0;i<hH->Nx * hH->Nz;i++) zero[i]=0.0;
	dev_h2d(hH, zero, hH->PML_Byx, hH->Nz * hH->Nx);
	dev_h2d(hH, zero, hH->PML_Byz, hH->Nz * hH->Nx);
	dev_h2d(hH, zero, hH->Jz, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->Jx, hH->Nx * hH->Nz);
	//we do not care of v and nl arrays.
	FL_DBL* n0 = new FL_DBL[hH->Nx * hH->Nz];
	FL_DBL* T0 = new FL_DBL[hH->Nx * hH->Nz];

	//set the source x:
	FL_DBL* srcx = new FL_DBL[hH->Nx];	
	for(int j = 0; j != hH->Nx; j++)
	{
		FL_DBL x = FL_DBL(j) / hH->dx_1;
		srcx[j] = exp(-pow((x - (hH->Lx/2.0 - hH->Lx/4.0*0.0))/hH->srcX,2.0));
	}
	dev_h2d(hH, srcx, hH->srcx, hH->Nx);
	
	for(int i = 0; i != hH->Nz; i++)
	{
		int ireal=hH->iz0+i;
		FL_DBL z=FL_DBL(ireal) / hH->dz_1;
		//FL_DBL levelL=hH->Lx/3.0;
		FL_DBL levelL=hH->Lx/2.0;
		FL_DBL levelH=hH->Lx/2.0;
		int leveli=(z/(hH->Lz/13.0)-floor(z/(hH->Lz/13.0))>0.5)?int(floor(levelH*hH->dx_1)):int(floor(levelL*hH->dx_1));
		for(int j = 0; j != hH->Nx; j++)
		{
			n0[i * hH->Nx + j] = FPT(0.0);//((j>leveli)&&(j<hH->Nx+hH->Nx/8)) ? hH->mediaN0 : FPT(0.0);
			T0[i * hH->Nx + j] = (n0[i * hH->Nx + j] != FPT(0.0)) ? hH->mediaTe0 : FPT(0.0);
		}
	}

	dev_h2d(hH, n0, hH->n, hH->Nx * hH->Nz);
	dev_h2d(hH, T0, hH->Te, hH->Nx * hH->Nz);
	printf("simpleGPUinit: blockSize=%d, x0=%d, z0=%d, dx_1=%g, dz_1=%g, dt=%g, PMLxmax=%g, PMLimin=%d\n", hH->blockSize, hH->ix0, hH->iz0, hH->dx_1, hH->dz_1, hH->dt, hH->PMLxmax, hH->PMLimin);
	delete[] zero;
	delete[] n0;
	delete[] T0;
	delete[] srcx;
	return;
}

//copy constructor:
extern "C"
hydro2dHandler::hydro2dHandler(const hydro2dHandler &obj, int dev) : device((dev == -1) ? obj.device : dev),
																																		 Nx(obj.Nx),
																																		 Nz(obj.Nz),
																																		 Lx(obj.Lx),
																																		 Lz(obj.Lz),
																																		 ix0(obj.ix0),
																																		 iz0(obj.iz0),
																																		 dx_1(obj.dx_1),
																																		 dz_1(obj.dz_1),
																																		 dt(obj.dt),
																																		 t(obj.t),
																																		 step(obj.step),
																																		 PMLimin(obj.PMLimin),
																																		 PMLimax(obj.PMLimax),
																																		 blockSize(obj.blockSize),
																																		 PMLxmax(obj.PMLxmax),
																																		 PMLstrength(obj.PMLstrength),
																																		 srcVelocity(obj.srcVelocity),
																																		 srcT(obj.srcT), srcNosc(obj.srcNosc),
																																		 srcAmp(obj.srcAmp),
																																		 mediaN0(obj.mediaN0),
																																		 mediaTe0(obj.mediaTe0),
																																		 diffusion(obj.diffusion),
																																		 mx_1(obj.mx_1),
																																		 mz_1(obj.mz_1),
																																		 mediaNu(obj.mediaNu),
																																		 media_bound_w2(obj.media_bound_w2),
																																		 media_bound_beta(obj.media_bound_beta),
																																		 media_bound_gamma(obj.media_bound_gamma),
																																		 media_phonon_omega(obj.media_phonon_omega),
																																		 media_phonon_phw(obj.media_phonon_phw),
																																		 media_phonon_beta(obj.media_phonon_beta),
																																		 toothDepth(obj.toothDepth),
																																		 toothWidth(obj.toothWidth),
																																		 mediaDepth(obj.mediaDepth),
																																		 a_cell(obj.a_cell),
																																		 toothDir(obj.toothDir),
																																		 flip(obj.flip),
																																		 SRCi(obj.SRCi),
																																		 shft(obj.shft),
																																		 doPolarization(obj.doPolarization),
																																		 doPhononAbsorbtion(obj.doPhononAbsorbtion),
																																		 extSource(obj.extSource),
																																		 srcX(obj.srcX),
																																		 srcTfactor(obj.srcTfactor),
																																		 switchOnDelay(obj.switchOnDelay)
{
	alloc_main_fields(this);

	FL_DBL* tmp = new FL_DBL[Nx * Nz];
	#define COPY_FIELD(field, sz)\
	dev_d2h(&obj, obj.field, tmp, sz);\
	dev_h2d(this, tmp, this->field, sz);

	size_t whole(Nx * Nz);
	
	//do srcx separatelly
	COPY_FIELD(srcx, Nx);
	
	if(doPolarization)
	{
		COPY_FIELD(Ptx, whole)
		COPY_FIELD(Ptz, whole)
		COPY_FIELD(Px, whole)
		COPY_FIELD(Pz, whole)
	}
	
	COPY_FIELD(By, whole)
	COPY_FIELD(Ex, whole)
	COPY_FIELD(Ez, whole)
	COPY_FIELD(PML_Byx, whole)
	COPY_FIELD(PML_Byz, whole)
	
	if(doPhononAbsorbtion)
	{
		COPY_FIELD(Phtx, whole)
		COPY_FIELD(Phtz, whole)
		COPY_FIELD(Phx, whole)
		COPY_FIELD(Phz, whole)
	}
	
	COPY_FIELD(Jz, whole)
	COPY_FIELD(Jx, whole)
	COPY_FIELD(n, whole)
	COPY_FIELD(Te, whole)

	#undef COPY_FIELD

	delete [] tmp;
}