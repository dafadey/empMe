#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_timer.h>
#include "funcall.h"
#include <iostream>
#include "geo.h"
#include "helpers.h"
#include <cmath>
#include <mutex>

//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION
//-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION-----DEVICE-SECTION

#define BS 16 // do not use any const variables this can break loops unroll in GPU code

extern "C"
hydro2dHandler::hydro2dHandler(int dev, bool _doPhononAbsorbtion, bool _doPolarization, bool _extSource, eHEATTYPE _JHEAT) :
								tim(), device(dev), Nx(NX), Nz(NZ), Lx(LX), surfaceX(LX/2.0), Lz(LZ), ix0(0), iz0(0), dx_1(0), dz_1(0), dt(0),
								Ptz(nullptr), Pty(nullptr), Ptx(nullptr), Pz(nullptr), Py(nullptr), Px(nullptr), Ez(nullptr), Ey(nullptr), Ex(nullptr), Bz(nullptr), By(nullptr), Bx(nullptr), Phtz(nullptr), Phty(nullptr), Phtx(nullptr), Phz(nullptr), Phy(nullptr), Phx(nullptr), Jz(nullptr), Jy(nullptr), Jx(nullptr), mat_mask(nullptr), n(nullptr), Te(nullptr), dTe(nullptr), nlx(nullptr), nly(nullptr), nlz(nullptr), vx(nullptr), vy(nullptr), vz(nullptr), PML_Byx(nullptr), PML_Byz(nullptr), PML_Eyx(nullptr), PML_Eyz(nullptr), feedPtz(nullptr), feedPtx(nullptr), feedPz(nullptr), feedPx(nullptr), feedEz(nullptr), feedEx(nullptr), feedBy(nullptr), feedPML_Byx(nullptr), feedPML_Byz(nullptr), feedPhtz(nullptr), feedPhtx(nullptr), feedPhz(nullptr), feedPhx(nullptr), feedJz(nullptr), feedJx(nullptr), feed_mat_mask(nullptr), feedn(nullptr), feedTe(nullptr), srcx(nullptr), srct(nullptr), host_feed(nullptr), host_srct(nullptr),
								t(0), step(0), PMLimin(0), PMLimax(0), SRCi(0), shft(SHIFTGRANULARITY), blockSize(BS), PMLxmax(PMLX), PMLstrength(PMLSTRENGTH),
								mediaN0(N0), mediaTe0(TE0),diffusion(DIFFUSION), mediaNu(MEDIANU), NUTratio(NUTRATIO), landauDamping(LANDAUDAMPING), Vdiffusion(VDIFFUSION), mx_1(MX_1), my_1(MY_1), mz_1(MZ_1), toothDepth(TOOTHDEPTH), toothWidth(TOOTHWIDTH), mediaDepth(MEDIADEPTH), a_cell(nullptr), media_bound_w2(BOUND_W2), media_bound_gamma(BOUND_GAMMA), media_bound_beta(BOUND_BETA), media_phonon_omega(PHONON_OMEGA), media_phonon_phw(PHONON_PHW), media_phonon_beta(PHONON_BETA), Byext(.0), srcTfactor(SRCTFACTOR), srcVelocity(VELOCITY), srcX(SRCX), switchOnDelay(SWITCHONDELAY),
								doPolarization(_doPolarization),  doPhononAbsorbtion(_doPhononAbsorbtion), extSource(_extSource), JHEAT(_JHEAT),
								toothDir(false), flip(false), linear(false) {}

__device__ FL_DBL
PML(FL_DBL x, FL_DBL PMLx, FL_DBL PMLstrength)
{
	return FPT(1.0)-((x<=PMLx)?(PMLx-x)/PMLx*PMLstrength:FPT(0.0));
}

//Polarization step needed for Bismuth
__global__ void
simpleP_Kernel(int Nx, int Nz, FL_DBL* Ez, FL_DBL* Ex, FL_DBL* Ptz, FL_DBL* Ptx, FL_DBL* Pz, FL_DBL* Px, FL_DBL* n, FL_DBL dt, FL_DBL w2_bound /*media_bound_w2*/, FL_DBL gamma_bound /*media_bound_gamma*/, FL_DBL beta_bound /*media_bound_beta*/, FL_DBL n0, FL_DBL* mat_mask)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
		const int im1 = i == 0 ? i : i - 1;
		const int jm1 = j == 0 ? j : j - 1;
				
		const FL_DBL n_z = n[i*Nx+j] != FPT(0.0) || n[im1*Nx+j] != FPT(0.0) ? (n[i*Nx+j] + n[im1*Nx+j]) * FPT(0.5) + n0 : FPT(0.0);
		//n_z = n_z > FPT(0.0) ? n_z : FPT(0.0);
		const FL_DBL n_x = mat_mask[i*Nx+j] != FPT(0.0) || mat_mask[i*Nx+jm1] != FPT(0.0) ? (n[i*Nx+j] + n[i*Nx+jm1]) * FPT(0.5) + n0 : FPT(0.0);
		//n_x = n_x > FPT(0.0) ? n_x : FPT(0.0);
		
		if(n_x == FPT(0.0) && n_z == FPT(0.0))
			return;
		
		Ptz[i*Nx+j] = n_z > FPT(0.0) ? (Ptz[i*Nx+j]/dt-w2_bound*Pz[i*Nx+j]-0.5*gamma_bound*Ptz[i*Nx+j]+beta_bound*Ez[i*Nx+j])*dt/(1.0+0.5*gamma_bound*dt) : FPT(0.0);
		Ptx[i*Nx+j] = n_x > FPT(0.0) ? (Ptx[i*Nx+j]/dt-w2_bound*Px[i*Nx+j]-0.5*gamma_bound*Ptx[i*Nx+j]+beta_bound*Ex[i*Nx+j])*dt/(1.0+0.5*gamma_bound*dt) : FPT(0.0);
		Pz[i*Nx+j] += Ptz[i*Nx+j]*dt;
		Px[i*Nx+j] += Ptx[i*Nx+j]*dt;
		Ez[i*Nx+j] += -Ptz[i*Nx+j]*dt;
		Ex[i*Nx+j] += -Ptx[i*Nx+j]*dt;
	}
}

__global__ void
simpleE_Kernel(int Nx, int Nz, FL_DBL* Ez, FL_DBL* Ex, FL_DBL* By, FL_DBL* Jz, FL_DBL* Jx, int ix0, int iz0, FL_DBL dx, FL_DBL dx_1, FL_DBL dz, FL_DBL dz_1, FL_DBL dt, FL_DBL* srct, int srci, int PMLimin, int PMLimax, FL_DBL pml, FL_DBL pmlx, FL_DBL* Ez_mid, FL_DBL* Ex_mid)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
		FL_DBL E_prev = Ez[i*Nx+j];
		FL_DBL pml_z_factor = ( j+ix0<PMLimin || j+ix0>PMLimax ) ? PML(FL_DBL((j+ix0<PMLimin)?j+ix0:Nx-1-j)*dx, pmlx, pml):FPT(1.0);
		#ifdef STATIC
		FL_DBL pml_x_factor = ( i+iz0<PMLimin || i+iz0>PMLimax ) ? PML(FL_DBL((i+iz0<PMLimin)?i+iz0:Nz-1-i)*dz, pmlx, pml):FPT(1.0);
		#else
		FL_DBL pml_x_factor = FL_DBL(1.);
		#endif
		Ez[i*Nx+j] = Ez[i*Nx+j] * pml_z_factor + (-(((j+1==Nx)?Ez[i*Nx+j]:By[i*Nx+j+1])-By[i*Nx+j]) * dx_1 - Jz[i*Nx+j] + /*source*/(srct && j + ix0 == srci ? srct[i] * dx_1 : FPT(0.0)))*dt;
    Ez_mid[i*Nx+j] = (Ez[i*Nx+j] + E_prev) * (FL_DBL) 0.5;

		E_prev = Ex[i*Nx+j];
    Ex[i*Nx+j] = Ex[i*Nx+j] * pml_x_factor + ((((i+1==Nz)?-Ex[i*Nx+j]:By[(i+1)*Nx+j]) - By[i*Nx+j]) * dz_1 - Jx[i*Nx+j]) * dt;
    Ex_mid[i*Nx+j] = (Ex[i*Nx+j] + E_prev) * FPT(0.5);
    
  }
}

__global__ void
simpleB_Kernel(int Nx, int Nz, FL_DBL* Ez, FL_DBL* Ex, FL_DBL* By, int ix0, int iz0, FL_DBL dx, FL_DBL dx_1, FL_DBL dz, FL_DBL dz_1, FL_DBL dt, int PMLimin, int PMLimax, FL_DBL pml, FL_DBL pmlx, FL_DBL* PML_Byz, FL_DBL* PML_Byx)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
//	without PML
//	By[i*Nx+j]+=(-(Ez[i*Nx+j]-((j==0)?FPT(0.0):Ez[i*Nx+j-1]))*dx_1+(Ex[i*Nx+j]-(i == 0 ? By[i * Nx + j] : Ex[(i-1)*Nx+j]))*dz_1)*dt;

		#ifdef STATIC
		bool PML_area = j+ix0<PMLimin || j+ix0>PMLimax || i < PMLimin || i > Nz - PMLimin;
		#else
		bool PML_area = j+ix0<PMLimin || j+ix0>PMLimax;
		#endif
		
		if(PML_area)
		{
			#ifdef STATIC
			if(i < PMLimin)
				PML_Byz[i*Nx+j]=PML_Byz[i*Nx+j]*PML(FL_DBL(i)*dx, pmlx, pml);
			else if(i > Nz - PMLimin)
				PML_Byz[i*Nx+j]=PML_Byz[i*Nx+j]*PML(FL_DBL(Nz-1-i)*dx, pmlx, pml);
			#endif

			if(j+ix0<PMLimin)
				PML_Byx[i*Nx+j] = PML_Byx[i*Nx+j]*PML(FL_DBL(j+ix0)*dx, pmlx, pml);
			else if(j+ix0>Nx-PMLimin)
				PML_Byx[i*Nx+j] = PML_Byx[i*Nx+j]*PML(FL_DBL(Nx-1-j)*dx, pmlx, pml);
			
			
			PML_Byx[i*Nx+j]+=-(Ez[i*Nx+j]-(j + ix0 == 0 ? -PML_Byx[i*Nx+j] : Ez[i*Nx+j-1]))*dx_1*dt;
			PML_Byz[i*Nx+j]+=(Ex[i*Nx+j]-(i == 0 ? PML_Byz[i * Nx + j] : Ex[(i-1)*Nx+j]))*dz_1*dt;
			
			By[i*Nx+j]=PML_Byz[i*Nx+j]+PML_Byx[i*Nx+j];

		}
		else
			By[i*Nx+j] += (-(Ez[i*Nx+j]-(j + ix0 == 0 ? -By[i*Nx+j] : Ez[i*Nx+j-1]))*dx_1+(Ex[i*Nx+j]-(i == 0 ? By[i * Nx + j] : Ex[(i-1)*Nx+j]))*dz_1)*dt;
	}
}

//TE TE TE  Bx Bz
__global__ void
simpleB_TE_Kernel(int Nx, int Nz, FL_DBL* Bz, FL_DBL* Bx, FL_DBL* Ey, int ix0, int iz0, FL_DBL dx, FL_DBL dx_1, FL_DBL dz, FL_DBL dz_1, FL_DBL dt, FL_DBL* srct, int srci, int PMLimin, int PMLimax, FL_DBL pml, FL_DBL pmlx)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
		FL_DBL pml_z_factor = ( j+ix0<PMLimin || j+ix0>Nx - PMLimin ) ? PML(FL_DBL((j+ix0<PMLimin)?j+ix0:Nx-1-j)*dx, pmlx, pml):FPT(1.0);
		Bz[i*Nx+j] = Bz[i*Nx+j] * pml_z_factor + ((Ey[i*Nx+j] - (j==0 ? Bz[i*Nx+j] : Ey[i*Nx+j-1])) * dx_1 + /*source*/(srct && j + ix0 == srci ? srct[i] * dx_1 : FPT(0.0))) * dt;

		#ifdef STATIC
		FL_DBL pml_x_factor = ( i+iz0<PMLimin || i+iz0>PMLimax ) ? PML(FL_DBL((i+iz0<PMLimin)?i+iz0:Nz-1-i)*dz, pmlx, pml):FPT(1.0);
		#else
		FL_DBL pml_x_factor = FPT(1.);
		#endif

    Bx[i*Nx+j] = Bx[i*Nx+j] * pml_x_factor + ((i == 0 ? -Bx[i*Nx+j] : Ey[(i-1)*Nx+j]) - Ey[i*Nx+j]) * dz_1 * dt;
    
  }
}

//TE TE TE  Ey
__global__ void
simpleE_TE_Kernel(int Nx, int Nz, FL_DBL* Bz, FL_DBL* Bx, FL_DBL* Ey, FL_DBL* Jy, int ix0, int iz0, FL_DBL dx, FL_DBL dx_1, FL_DBL dz, FL_DBL dz_1, FL_DBL dt, int PMLimin, int PMLimax, FL_DBL pml, FL_DBL pmlx, FL_DBL* PML_Eyz, FL_DBL* PML_Eyx, FL_DBL* Ey_mid)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
		FL_DBL E_prev = Ey[i*Nx+j];
		FL_DBL E_next = E_prev;

		#ifdef STATIC
		bool PML_area = j+ix0<PMLimin || j+ix0>Nx-PMLimin || i < PMLimin || i > Nz - PMLimin;
		#else
		bool PML_area = j+ix0<PMLimin || j+ix0>Nx-PMLimin;
		#endif
		
		if(PML_area)
		{
			#ifdef STATIC
			if(i < PMLimin)
				PML_Eyz[i*Nx+j]=PML_Eyz[i*Nx+j]*PML(FL_DBL(i)*dx, pmlx, pml);
			else if(i > Nz - PMLimin)
				PML_Eyz[i*Nx+j]=PML_Eyz[i*Nx+j]*PML(FL_DBL(Nz-1-i)*dx, pmlx, pml);
			#endif

			if(j+ix0<PMLimin)
				PML_Eyx[i*Nx+j] = PML_Eyx[i*Nx+j]*PML(FL_DBL(j+ix0)*dx, pmlx, pml);
			else if(j+ix0>Nx-PMLimin)
				PML_Eyx[i*Nx+j] = PML_Eyx[i*Nx+j]*PML(FL_DBL(Nx-1-j)*dx, pmlx, pml);

			PML_Eyz[i*Nx+j] += (Bx[i*Nx+j]-(i+1==Nz ? PML_Eyz[i*Nx+j] : Bx[(i+1)*Nx+j]))*dz_1*dt;
			PML_Eyx[i*Nx+j] += ((j+ix0+1==Nx ? -PML_Eyx[i*Nx+j] : Bz[i*Nx+j+1])-Bz[i*Nx+j])*dx_1*dt;

			E_next=PML_Eyz[i*Nx+j]+PML_Eyx[i*Nx+j] - Jy[i*Nx+j]*dt;
		}
		else
			E_next += (((j+ix0+1==Nx ? -Ey[i*Nx+j] : Bz[i*Nx+j+1])-Bz[i*Nx+j])*dx_1-((i+1==Nz ? Ey[i*Nx+j] : Bx[(i+1)*Nx+j])-Bx[i*Nx+j])*dz_1 - Jy[i*Nx+j])*dt;
		
		Ey_mid[i*Nx+j] = FL_DBL(.5)*(E_next+E_prev);
		Ey[i*Nx+j] = E_next;
	}
}


__global__ void
simpleV_Kernel(int Nx, int Nz, FL_DBL* vx, FL_DBL* vy, FL_DBL* vz, FL_DBL* Jx, FL_DBL* Jy, FL_DBL* Jz, FL_DBL* n, FL_DBL n0)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
    if(Jz[i*Nx+j] == FPT(0.0) && Jx[i*Nx+j] == FPT(0.0) && Jy[i*Nx+j] == FPT(0.0))
    {
      vz[i*Nx+j] = FPT(0.0);
      vx[i*Nx+j] = FPT(0.0);
      vy[i*Nx+j] = FPT(0.0);
      return;
    }
    
		int im1 = i == 0 ? i : i - 1;
		int jm1 = j == 0 ? j : j - 1;
		const FL_DBL n_z = (n[i*Nx+j] + n[im1*Nx+j]) * FPT(0.5) + n0;
		const FL_DBL n_x = (n[i*Nx+j] + n[i*Nx+jm1]) * FPT(0.5) + n0;
		const FL_DBL n_y = n[i*Nx+j] + n0;
				
		vz[i*Nx+j] = (n_z > FPT(0.0)) ? Jz[i*Nx+j] / n_z : FPT(0.0);
		vx[i*Nx+j] = (n_x > FPT(0.0)) ? Jx[i*Nx+j] / n_x : FPT(0.0);
		vy[i*Nx+j] = (n_y > FPT(0.0)) ? Jy[i*Nx+j] / n_y : FPT(0.0);
	}
}

__global__ void
simpleNl_Kernel(int Nx, int Nz, FL_DBL* nlx, FL_DBL* nly, FL_DBL* nlz, FL_DBL* vx, FL_DBL* vy, FL_DBL* vz, FL_DBL* Jx, FL_DBL* Jy, FL_DBL* Jz, FL_DBL* Bx, FL_DBL* By, FL_DBL* Bz, FL_DBL dz_1, FL_DBL dx_1, FL_DBL mx_1, FL_DBL my_1, FL_DBL mz_1, FL_DBL Bxext, FL_DBL Byext, FL_DBL Bzext)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
		int ip1 = i == Nz - 1 ? i : i + 1;
		int im1 = i == 0      ? i : i - 1;
		int jp1 = j == Nx - 1 ? j : j + 1;
		int jm1 = j == 0      ? j : j - 1;

		nlz[i*Nx+j]=
					HYDROZSRC*(Jz[i*Nx+j]*(vz[ip1*Nx+j]-vz[im1*Nx+j])*FPT(0.5)*dz_1+FPT(0.25)*(Jx[i*Nx+j]+Jx[i*Nx+jp1]+Jx[im1*Nx+j]+Jx[im1*Nx+jp1])*(vz[i*Nx+jp1]-vz[i*Nx+jm1])*FPT(0.5)*dx_1)
					+MAINSRC*vz[i*Nx+j]*((Jx[i*Nx+jp1]+Jx[im1*Nx+jp1]-Jx[i*Nx+j]-Jx[im1*Nx+j])*FPT(0.5)*dx_1+(Jz[ip1*Nx+j]-Jz[im1*Nx+j])*FPT(0.5)*dz_1)
					+((Jx[i*Nx+j]+Jx[im1*Nx+j]+Jx[i*Nx+jp1]+Jx[im1*Nx+jp1])*(By[i*Nx+j]+By[i*Nx+jp1]+Byext+Byext)*FPT(0.125)-(Jy[im1*Nx+j]+Jy[i*Nx+j])*(Bx[i*Nx+j]+Bxext)*FPT(0.5))*mz_1*HSRC;


		nlx[i*Nx+j]=
					HYDROXSRC*(FPT(0.25)*(Jz[i*Nx+j]+Jz[i*Nx+jm1]+Jz[ip1*Nx+j]+Jz[ip1*Nx+jm1])*(vx[ip1*Nx+j]-vx[im1*Nx+j])*FPT(0.5)*dz_1+Jx[i*Nx+j]*(vx[i*Nx+jp1]-vx[i*Nx+jm1])*FPT(0.5)*dx_1)
					+PRESSURESRC*vx[i*Nx+j]*((Jx[i*Nx+jp1]-Jx[i*Nx+jm1])*dx_1*FPT(0.5)+(Jz[ip1*Nx+j]+Jz[ip1*Nx+jm1]-Jz[i*Nx+j]-Jz[i*Nx+jm1])*dz_1*FPT(0.5))
					+(-(Jz[i*Nx+j]+Jz[ip1*Nx+j]+Jz[i*Nx+jm1]+Jz[ip1*Nx+jm1])*(By[ip1*Nx+j]+By[i*Nx+j]+Byext+Byext)*FPT(0.125)+(Jy[i*Nx+jm1]+Jy[i*Nx+j])*(Bz[i*Nx+j] + Bzext)*FPT(0.5))*mx_1*HSRC;

		nly[i*Nx+j]=
					HYDROYSRC*(FPT(0.5)*(Jz[i*Nx+j]+Jz[ip1*Nx+j])*(vy[ip1*Nx+j]-vy[im1*Nx+j])*FPT(0.5)*dz_1+FPT(0.5)*(Jx[i*Nx+j]+Jx[i*Nx+jp1])*(vy[i*Nx+jp1]-vy[i*Nx+jm1])*FPT(0.5)*dx_1)
					+MAINSRC*vy[i*Nx+j]*((Jx[i*Nx+jp1]-Jx[i*Nx+j])*dx_1+(Jz[ip1*Nx+j]-Jz[i*Nx+j])*dz_1)
					+((Jz[i*Nx+j]+Jz[ip1*Nx+j])*(Bx[ip1*Nx+j]+Bx[i*Nx+j]+Bxext+Bxext) + (Jx[i*Nx+j]+Jx[i*Nx+jp1])*(Bz[i*Nx+jp1]+Bz[i*Nx+j]+Bzext+Bzext))*FPT(0.25)*my_1*HSRC;
	}
}

__global__ void
simpleN_Kernel(int Nx, int Nz, FL_DBL* n, FL_DBL* Jx, FL_DBL* Jz, FL_DBL dz_1, FL_DBL dx_1, FL_DBL dt, FL_DBL* mat_mask)
{
	const int j = threadIdx.x + blockIdx.x * blockDim.x;
	const int i = threadIdx.y + blockIdx.y * blockDim.y;
  
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
      if (mat_mask[i*Nx+j] == FL_DBL(0.0))
      {
        n[i*Nx+j] = FPT(0.0);
        return;
			}
      int ip1 = (i == Nz - 1) ? i : i + 1;
			int jp1 = (j == Nx - 1) ? j : j + 1;
      
			n[i*Nx+j] += -((Jx[i*Nx+jp1] - Jx[i*Nx+j]) * dx_1 +
                     (Jz[ip1*Nx+j] - Jz[i*Nx+j]) * dz_1) * dt;
			
      n[i*Nx+j] = mat_mask[i*Nx+j] == FL_DBL(0.0) ? FPT(0.0) : n[i*Nx+j]; // grants n0>0 and 0 outside plasma area
	}
}

__device__ inline FL_DBL sqr(FL_DBL x)
{return x*x;}

__global__ void
simple_dTe_Kernel(int Nx, int Nz, FL_DBL D, FL_DBL NuTratio, FL_DBL mediaTe0, FL_DBL* dTe, FL_DBL* Te, FL_DBL* n, FL_DBL* Ex, FL_DBL* Ez, FL_DBL* jx, FL_DBL* jz, FL_DBL dz_1, FL_DBL dx_1, FL_DBL* srct, FL_DBL* srcx, FL_DBL SRCFACTOR, FL_DBL FACTOR_DEPENDS_ON_T, FL_DBL* mat_mask)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
      if(mat_mask[i*Nx+j] == FPT(0.0))
      {
        dTe[i*Nx+j] = FPT(0.0);
        return;
      }
      
			const bool vac = mat_mask[i*Nx+j] == FPT(0.0);
			int ip1=(i==Nz-1)?i:i+1;
			int jp1=(j==Nx-1)?j:j+1;
			
			FL_DBL src;
			if(!srct)
				src=(Ex[i*Nx+j] + Ex[i*Nx+jp1]) * (jx[i*Nx+j] + jx[i*Nx+jp1]) * FPT(0.25)
				    + (Ez[i*Nx+j] + Ez[ip1*Nx+j]) * (jz[i*Nx+j] + jz[ip1*Nx+j]) * FPT(0.25);
			else
				src = srct[i] * srcx[j];
			
			ip1 = (i==Nz-1 || mat_mask[(i+1)*Nx+j] == FPT(0.0)) ? i : i + 1; //von neumann bc
			jp1 = (j==Nx-1 || mat_mask[i*Nx+(j+1)] == FPT(0.0)) ? j : j + 1; //von neumann bc
			int im1 = (i==0 || mat_mask[(i-1)*Nx+j] == FPT(0.0)) ? i : i - 1; //von neumann bc
			int jm1 = (j==0 || mat_mask[i*Nx+(j-1)] == FPT(0.0)) ? j : j - 1; //von neumann bc
      
      const FL_DBL ratio = FPT(1.0) + NuTratio * SQRTNUT(Te[i*Nx+j]);
      const FL_DBL d = D / ratio;
			
			const FL_DBL srcTfactor = SRCFACTOR * (FACTOR_DEPENDS_ON_T * ratio + (FPT(1.) - FACTOR_DEPENDS_ON_T));
			
			//#define SRCFACTOR FPT(100.0)

			dTe[i*Nx+j] = srcTfactor * src + d * ((Te[i*Nx+jp1]+Te[i*Nx+jm1]-FPT(2.0)*Te[i*Nx+j])*dx_1*dx_1+(Te[ip1*Nx+j]+Te[im1*Nx+j]-FPT(2.0)*Te[i*Nx+j])*dz_1*dz_1);

/*
			dTe[i*Nx+j] = SRCFACTOR * src + d * ((Te[i*Nx+jp1]+Te[i*Nx+jm1]-FPT(2.0)*Te[i*Nx+j])*dx_1*dx_1+(Te[ip1*Nx+j]+Te[im1*Nx+j]-FPT(2.0)*Te[i*Nx+j])*dz_1*dz_1);
*/

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
simple_Te_Kernel(int Nx, int Nz, FL_DBL* dTe, FL_DBL* Te, FL_DBL dt, FL_DBL mediaTe0)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
  if(dTe[i*Nx+j] == FPT(0.0))
    return;
	Te[i*Nx+j] += dTe[i*Nx+j] * dt;
	Te[i*Nx+j] = mediaTe0 + Te[i*Nx+j] < FPT(0.) ? -mediaTe0 : Te[i*Nx+j]; 
}

__global__ void
simpleJ_Kernel(int Nx, int Nz, FL_DBL* Jx, FL_DBL* Jy, FL_DBL* Jz, FL_DBL* Ex, FL_DBL* Ey, FL_DBL* Ez, FL_DBL* nlx, FL_DBL* nly, FL_DBL* nlz, FL_DBL* vx, FL_DBL* vy, FL_DBL* vz,FL_DBL* n, FL_DBL* Te, FL_DBL NUTratio, FL_DBL mediaTe0, FL_DBL dx_1, FL_DBL dz_1, FL_DBL dt, FL_DBL Nu, FL_DBL landauDamping, FL_DBL Vdiffusion, FL_DBL* phv_z, FL_DBL* phv_x, FL_DBL* ph_z, FL_DBL* ph_x, FL_DBL mx_1, FL_DBL my_1,FL_DBL mz_1, FL_DBL omegaOpt /*media_phonon_omega*/, FL_DBL PHW /*media_phonon_phw*/, FL_DBL OPT_PH /*media_phonon_beta*/, bool do_phonons, FL_DBL n0, FL_DBL* mat_mask, bool hh_linear, bool neumann)
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
		int ip1 = (i == Nz - 1) ? i : i + 1;
		int im1 = (i == 0     ) ? i : i - 1;
		int jp1 = (j == Nx - 1) ? j : j + 1;
		int jm1 = (j == 0     ) ? j : j - 1;
		
		int ip1vd,im1vd,jp1vd,jm1vd; // indices for vdiffusion
		if(neumann) {
			ip1vd = (i==Nz-1 || mat_mask[(i+1)*Nx+j] == FPT(0.0)) ? i : i + 1; //von neumann bc
			jp1vd = (j==Nx-1 || mat_mask[i*Nx+(j+1)] == FPT(0.0)) ? j : j + 1; //von neumann bc
			im1vd = (i==0 || mat_mask[(i-1)*Nx+j] == FPT(0.0)) ? i : i - 1; //von neumann bc
			jm1vd = (j==0 || mat_mask[i*Nx+(j-1)] == FPT(0.0)) ? j : j - 1; //von neumann bc
		} else {
		  ip1vd = (i == Nz - 1) ? i : i + 1;
		  im1vd = (i == 0     ) ? i : i - 1;
		  jp1vd = (j == Nx - 1) ? j : j + 1;
		  jm1vd = (j == 0     ) ? j : j - 1;
		}
//PROBLEMS APPEAR IF ADDING THIS CODE ON gpu# > 0

		const FL_DBL mask_z = mat_mask[i*Nx+j] == mat_mask[im1*Nx+j] ? mat_mask[i*Nx+j] : FPT(0.0);
		const FL_DBL mask_x = mat_mask[i*Nx+j] == mat_mask[i*Nx+jm1] ? mat_mask[i*Nx+j] : FPT(0.0);
		const FL_DBL mask_y = mat_mask[i*Nx+j];

    if(mask_x == FPT(0.0) && mask_z == FPT(0.0) && mask_y == FPT(0.0))
      return;

    FL_DBL n_z = mask_z > FPT(0.0) ? (n[i*Nx+j] + n[im1*Nx+j]) * FPT(0.5) : FPT(0.0);
    const FL_DBL n0_z = mask_z > FPT(0.0) ? n0 : FPT(0.0);
    n_z = (n_z + n0_z) < FPT(0.0) ? -n0_z : n_z;

		FL_DBL n_x = mask_x > FPT(0.0) ? (n[i*Nx+j] + n[i*Nx+jm1]) * FPT(0.5) : FPT(0.0);
		const FL_DBL n0_x = mask_x > FPT(0.0) ? n0 : FPT(0.0);
    n_x = (n_x + n0_x) < FPT(0.0) ? -n0_x : n_x;

		FL_DBL n_y = n[i*Nx+j];
		const FL_DBL n0_y = mask_y > FPT(0.0) ? n0 : FPT(0.0);
    n_y = (n_y + n0_y) < FPT(0.0) ? -n0_y : n_y;

    n_z = MAINSRC == OFF || hh_linear ? FPT(0.0) : n_z;
    n_x = PRESSURESRC == OFF || hh_linear ? FPT(0.0) : n_x;


    const FL_DBL Te_z = FPT(0.5) * (Te[i*Nx+j] + Te[im1*Nx+j]);


 		const FL_DBL Nu_z = (mask_z == FPT(0.0)) ? FPT(0.0) : (hh_linear ? Nu : Nu * (FPT(1.0) + NUTratio * SQRTNUT(Te_z)) );

/*
 		FL_DBL Nu_z = (mask_z == FPT(0.0)) ? FPT(0.0) : (hh_linear ? Nu : Nu * (FPT(1.0) + FPT(5.0) * NUTratio * SQRTNUT(Te_z)) / (FPT(1.0) + SQRTNUT(Te_z / mediaTe0)));
*/

/*
 		FL_DBL Nu_z = (mask_z == FPT(0.0)) ? FPT(0.0) : (hh_linear ? Nu : Nu * (FPT(1.0) + NUTratio * SQRTNUT(Te_z)) / (FPT(1.0) + FPT(5.0) * SQRTNUT(Te_z / mediaTe0)));
*/
    
    const FL_DBL Te_x = FPT(0.5) * (Te[i*Nx+j] + Te[i*Nx+jm1]);

		const FL_DBL Nu_x = (mask_x == FPT(0.0)) ? FPT(0.0) : (hh_linear ? Nu : Nu * (FPT(1.0) + NUTratio * SQRTNUT(Te_x)) );	

    
    /*
		FL_DBL Nu_x = (mask_x == FPT(0.0)) ? FPT(0.0) : (hh_linear ? Nu : Nu * (FPT(1.0) + FPT(5.0) * NUTratio * SQRTNUT(Te_x)) / (FPT(1.0) + SQRTNUT(Te_x / mediaTe0)));
    */
/*    
		FL_DBL Nu_x = (mask_x == FPT(0.0)) ? FPT(0.0) : (hh_linear ? Nu : Nu * (FPT(1.0) + NUTratio * SQRTNUT(Te_x)) / (FPT(1.0) + FPT(5.0) * SQRTNUT(Te_x / mediaTe0)));	
*/


    const FL_DBL Te_y = Te[i*Nx+j];

		const FL_DBL Nu_y = (mask_y == FPT(0.0)) ? FPT(0.0) : (hh_linear ? Nu : Nu * (FPT(1.0) + NUTratio * SQRTNUT(Te_y)) );

    	
		FL_DBL Jz_t = (n_z + n0_z) * Vdiffusion*((vz[i*Nx+jp1vd]+vz[i*Nx+jm1vd]-FPT(2.0)*vz[i*Nx+j])*dx_1*dx_1+(vz[ip1vd*Nx+j]+vz[im1vd*Nx+j]-FPT(2.0)*vz[i*Nx+j])*dz_1*dz_1)+landauDamping*(Jz[i*Nx+jp1]+Jz[i*Nx+jm1]-FPT(2.0)*Jz[i*Nx+j])*dx_1*dx_1/*-Jz[i*Nx+j] * Nu_z*/ 
                  + mz_1 * (n_z + n0_z) * Ez[i*Nx+j] - (hh_linear ? FPT(0.0) : IFNONLIN*nlz[i*Nx+j]) - mz_1 * ((hh_linear ? FPT(0.0) : (n[i*Nx+j] * Te[i*Nx+j] - n[im1*Nx+j] * Te[im1*Nx+j] + n0_z * (Te[i*Nx+j] - Te[im1*Nx+j]))) + mediaTe0 * (n[i*Nx+j] - n[im1*Nx+j])) * dz_1;


/*
		FL_DBL Jz_t = n0_z * Ez[i*Nx+j] - (n[i*Nx+j] * Te[i*Nx+j] - n[im1*Nx+j] * Te[im1*Nx+j] + n0_z * (Te[i*Nx+j] - Te[im1*Nx+j]) + mediaTe0 * (n[i*Nx+j] - n[im1*Nx+j])) * dz_1;
*/

/*
		FL_DBL Jz_t = n0_z * Ez[i*Nx+j] - n0_z * (Te[i*Nx+j] - Te[im1*Nx+j]) * dz_1;
*/
/*
		FL_DBL Jz_t = n0_z * Ez[i*Nx+j];
*/
    //FL_DBL Jz_t = (n_z + n0_z) * Ez[i*Nx+j] - (n[i*Nx+j] * Te[i*Nx+j] - n[im1*Nx+j] * Te[im1*Nx+j] + mediaTe0 * (n[i*Nx+j] - n[im1*Nx+j])) * dz_1


		FL_DBL Jx_t = (n_x + n0_x) * Vdiffusion*((vx[i*Nx+jp1vd]+vx[i*Nx+jm1vd]-FPT(2.0)*vx[i*Nx+j])*dx_1*dx_1+(vx[ip1vd*Nx+j]+vx[im1vd*Nx+j]-FPT(2.0)*vx[i*Nx+j])*dz_1*dz_1)+landauDamping*(Jx[ip1*Nx+j]+Jx[im1*Nx+j]-FPT(2.0)*Jx[i*Nx+j])*dz_1*dz_1/* - Jx[i*Nx+j] * Nu_x*/
                  + mx_1 * (n_x + n0_x) * Ex[i*Nx+j] - (hh_linear ? FPT(0.0) : IFNONLIN*nlx[i*Nx+j]) - mx_1 * ((hh_linear ? FPT(0.0) : (n[i*Nx+j] * Te[i*Nx+j] - n[i*Nx+jm1] * Te[i*Nx+jm1] + n0_x * (Te[i*Nx+j] - Te[i*Nx+jm1]))) + mediaTe0 * (n[i*Nx+j] - n[i*Nx+jm1])) * dx_1;


/*
		FL_DBL Jx_t = n0_x * Ex[i*Nx+j] - (n[i*Nx+j] * Te[i*Nx+j] - n[i*Nx+jm1] * Te[i*Nx+jm1] + n0_x * (Te[i*Nx+j] - Te[i*Nx+jm1]) + mediaTe0 * (n[i*Nx+j] - n[i*Nx+jm1])) * dx_1;
*/

/*
		FL_DBL Jx_t = n0_x * Ex[i*Nx+j] - n0_x * (Te[i*Nx+j] - Te[i*Nx+jm1]) * dx_1;
*/
/*
		FL_DBL Jx_t = n0_x * Ex[i*Nx+j];
*/
    //FL_DBL Jx_t = (n_x + n0_x) * Ex[i*Nx+j] - (n[i*Nx+j] * Te[i*Nx+j] - n[i*Nx+jm1] * Te[i*Nx+jm1] + mediaTe0 * (n[i*Nx+j] - n[i*Nx+jm1])) * dx_1;



		FL_DBL Jy_t = (n_y + n0_y) * Vdiffusion*((vy[i*Nx+jp1vd]+vy[i*Nx+jm1vd]-FPT(2.0)*vy[i*Nx+j])*dx_1*dx_1+(vy[ip1vd*Nx+j]+vy[im1vd*Nx+j]-FPT(2.0)*vy[i*Nx+j])*dz_1*dz_1)+landauDamping*(Jy[ip1*Nx+j]+Jy[im1*Nx+j]-FPT(2.0)*Jy[i*Nx+j])*dz_1*dz_1/* - Jy[i*Nx+j] * Nu_y*/
                  + my_1 * (n_y + n0_y) * Ey[i*Nx+j] - (hh_linear ? FPT(0.0) : IFNONLIN*nly[i*Nx+j]);


		if(do_phonons)
		{
			Jz[i*Nx+j] += (Jz_t - omegaOpt*omegaOpt*ph_z[i*Nx+j]) * dt;
			Jx[i*Nx+j] += (Jx_t - omegaOpt*omegaOpt*ph_x[i*Nx+j]) * dt;
		}
		else
		{
			Jz[i*Nx+j] += Jz_t * dt;
			Jx[i*Nx+j] += Jx_t * dt;
			Jy[i*Nx+j] += Jy_t * dt;
		}
		
		//put zeros at plasma surface
		Jx[i*Nx+j] = mask_x == FPT(0.0) ? FPT(0.0) : Jx[i*Nx+j];
		Jz[i*Nx+j] = mask_z == FPT(0.0) ? FPT(0.0) : Jz[i*Nx+j];
		Jy[i*Nx+j] = mask_y == FPT(0.0) ? FPT(0.0) : Jy[i*Nx+j];

		if(do_phonons)
		{
			//phv_z[i*Nx+j]+=(-omegaOpt*omegaOpt*ph_z[i*Nx+j]-PHW*omegaOpt*phv_z[i*Nx+j]+OPT_PH*(Jz_t-omegaOpt*omegaOpt*ph_z[i*Nx+j]))*dt;
			phv_z[i*Nx+j]+=(-omegaOpt*omegaOpt*ph_z[i*Nx+j]-PHW*omegaOpt*phv_z[i*Nx+j]+OPT_PH*Jz_t)*dt;
			//phv_x[i*Nx+j]+=(-omegaOpt*omegaOpt*ph_x[i*Nx+j]-PHW*omegaOpt*phv_x[i*Nx+j]+OPT_PH*(Jx_t-omegaOpt*omegaOpt*ph_x[i*Nx+j]))*dt;
			phv_x[i*Nx+j]+=(-omegaOpt*omegaOpt*ph_x[i*Nx+j]-PHW*omegaOpt*phv_x[i*Nx+j]+OPT_PH*Jx_t)*dt;

			ph_z[i*Nx+j]+=phv_z[i*Nx+j]*dt;
			ph_x[i*Nx+j]+=phv_x[i*Nx+j]*dt;		
		}
		
		//Jx[i*Nx+j] = Nu_x * dt > FPT(1.) ? FPT(.0) : Jx[i*Nx+j];
		//Jy[i*Nx+j] = Nu_y * dt > FPT(1.) ? FPT(.0) : Jy[i*Nx+j];
		//Jz[i*Nx+j] = Nu_z * dt > FPT(1.) ? FPT(.0) : Jz[i*Nx+j];
		Jx[i*Nx+j] *= exp(-Nu_x * dt);
		Jy[i*Nx+j] *= exp(-Nu_y * dt);
		Jz[i*Nx+j] *= exp(-Nu_z * dt);
		
	}
}

__global__ void
shift_Kernel(int Nx, int Nz, FL_DBL* from, FL_DBL* to, int shift, FL_DBL* feed) // NOTE: shift is counted in blocks
{
	const int j=threadIdx.x+blockIdx.x*blockDim.x;
	const int i=threadIdx.y+blockIdx.y*blockDim.y;
	if(i >= 0 && i < Nz && j >= 0 && j < Nx)
	{
		int iold = i + shift * blockDim.y;
		if(iold < Nz)
			to[i*Nx+j] = from[iold * Nx + j];
		else
			to[i*Nx+j] = !feed ? FPT(0.0) : feed[(iold - Nz) * Nx + j];
	}
}

__global__ void
copy_Kernel(int Nx, int Nz, FL_DBL* from, FL_DBL* to) // NOTE: shift is counted in blocks
{
		const int j=threadIdx.x+blockIdx.x*blockDim.x;
		const int i=threadIdx.y+blockIdx.y*blockDim.y;
    if(i >= 0 && i < Nz && j >= 0 && j < Nx)
			to[i*Nx+j] = from[i*Nx+j];
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
		//TE:
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Py, hH->vy, hH->shft, hH->feedEy);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Pty, hH->Py, hH->shft, hH->feedEy);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ey, hH->Pty, hH->shft, hH->feedEy);
		cudaThreadSynchronize();
	}
	else
	{
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->vz, hH->shft, hH->feedEz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ex, hH->vx, hH->shft, hH->feedEx);
		cudaThreadSynchronize();
		//TE:
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ey, hH->vy, hH->shft, hH->feedEy);
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
		//TE:
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Phty, hH->Ey, hH->shft, hH->feedEy);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Phy, hH->Phty, hH->shft, hH->feedEy);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jy, hH->Phy, hH->shft, hH->feedEy);
		cudaThreadSynchronize();
	}
	else
	{
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jz, hH->Ez, hH->shft, hH->feedJz);
		cudaThreadSynchronize();
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jx, hH->Ex, hH->shft, hH->feedJx);
		cudaThreadSynchronize();
		//TE:
		shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jy, hH->Ey, hH->shft, hH->feedJy);
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

//shift T
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Te , hH->PML_Byz, hH->shft, hH->feedTe);
	cudaThreadSynchronize();

//shift mask
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->mat_mask , hH->PML_Byx, hH->shft, hH->feed_mat_mask);
	cudaThreadSynchronize();

//TE:
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->PML_Eyx, hH->Jy, hH->shft, hH->feedPML_Eyx);
	cudaThreadSynchronize();
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->PML_Eyz, hH->PML_Eyx, hH->shft, hH->feedPML_Eyz);
	cudaThreadSynchronize();
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Bz, hH->PML_Eyz, hH->shft, hH->feedBz);
	cudaThreadSynchronize();
	shift_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Bx, hH->Bz, hH->shft, hH->feedBx);
	cudaThreadSynchronize();
//roll pointers:
	FL_DBL* new_Ez = 0;
	FL_DBL* new_Ex = 0;
	FL_DBL* new_Pz = 0;
	FL_DBL* new_Px = 0;
	FL_DBL* new_Ptz = 0;
	FL_DBL* new_Ptx = 0;
	//TE:
	FL_DBL* new_Ey = 0;
	FL_DBL* new_Py = 0;
	FL_DBL* new_Pty = 0;
	
	if(doPolarization)
	{
		new_Pz=hH->vz;
		new_Px=hH->vx;
		new_Ptz=hH->Pz;
		new_Ptx=hH->Px;
		new_Ez=hH->Ptz;
		new_Ex=hH->Ptx;
		//TE:
		new_Py=hH->vy;
		new_Pty=hH->Py;
		new_Ey=hH->Pty;
	}
	else
	{
		new_Ez=hH->vz;
		new_Ex=hH->vx;
		//TE:
		new_Ey=hH->vy;
	}
	
	FL_DBL* new_Jz = 0;
	FL_DBL* new_Jx = 0;
	FL_DBL* new_Phtz = 0;
	FL_DBL* new_Phtx = 0;
	FL_DBL* new_Phz = 0;
	FL_DBL* new_Phx = 0;
	//TE:
	FL_DBL* new_Jy = 0;
	FL_DBL* new_Phty = 0;
	FL_DBL* new_Phy = 0;

	if(doPhononAbsorbtion)
	{
		new_Phtz=hH->Ez;
		new_Phtx=hH->Ex;
		new_Phz=hH->Phtz;
		new_Phx=hH->Phtx;
		new_Jz=hH->Phz;
		new_Jx=hH->Phx;
		//TE:
		new_Phty=hH->Ey;
		new_Phy=hH->Phty;
		new_Jy=hH->Phy;
	}
	else
	{
		new_Jz=hH->Ez;
		new_Jx=hH->Ex;
		//TE:
		new_Jy=hH->Ey;
	}
	FL_DBL* new_By=hH->Jz;
	FL_DBL* new_n =hH->Jx;
	FL_DBL* new_PML_Byz=hH->By;
	FL_DBL* new_PML_Byx=hH->n;
	FL_DBL* new_Te=hH->PML_Byz;
	FL_DBL* new_mat_mask=hH->PML_Byx;
	FL_DBL* new_vx=hH->Te;
	FL_DBL* new_vz=hH->mat_mask;
	//TE:
	FL_DBL* new_PML_Eyx=hH->Jy;
	FL_DBL* new_PML_Eyz=hH->PML_Eyx;
	FL_DBL* new_Bz=hH->PML_Eyz;
	FL_DBL* new_Bx=hH->Bz;
	FL_DBL* new_vy=hH->Bx;

  // reassign

	if(doPolarization)
	{
		hH->Ptz=new_Ptz;
		hH->Ptx=new_Ptx;
		hH->Pz=new_Pz;
		hH->Px=new_Px;
		//TE:
		hH->Py=new_Py;
		hH->Pty=new_Pty;
	}
	hH->Ez=new_Ez;
	hH->Ex=new_Ex;
	//TE:
	hH->Ey=new_Ey;
	if(doPhononAbsorbtion)
	{
		hH->Phtz=new_Phtz;
		hH->Phtx=new_Phtx;
		hH->Phz=new_Phz;
		hH->Phx=new_Phx;
		//TE:
		hH->Phty=new_Phty;
		hH->Phy=new_Phy;
	}
	hH->Jz=new_Jz;
	hH->Jx=new_Jx;
	hH->By=new_By;
	hH->n =new_n ;
	hH->vz=new_vz;
	hH->vx=new_vx;
	hH->Te=new_Te;
	hH->mat_mask=new_mat_mask;
	hH->PML_Byz=new_PML_Byz;
	hH->PML_Byx=new_PML_Byx;
	//TE:
	hH->Jy=new_Jy;
	hH->PML_Eyx=new_PML_Eyx;
	hH->PML_Eyz=new_PML_Eyz;
	hH->Bz=new_Bz;
	hH->Bx=new_Bx;
	hH->vy=new_vy;

	hH->iz0 += hH->shft * hH->blockSize;
}




/*
inline void roll(hydro2dHandler* hH) // shft is counted in blocks
{
	dim3 grid(-floor(-FL_DBL(hH->Nx)/FL_DBL(hH->blockSize)),-floor(-FL_DBL(hH->Nz)/FL_DBL(hH->blockSize)),1);
	//printf("grid=%d %d %d\n", grid.x, grid.y, grid.z);
	dim3 threads(hH->blockSize,hH->blockSize,1);
	//shift big data

  copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->vz);
  cudaThreadSynchronize();
  copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ex, hH->vx);
  cudaThreadSynchronize();
  copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jz, hH->Ez);
  cudaThreadSynchronize();
  copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jx, hH->Ex);
  cudaThreadSynchronize();
	copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->By, hH->Jz);
	cudaThreadSynchronize();
	copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->n, hH->Jx);
	cudaThreadSynchronize();
	copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->PML_Byz , hH->By);
	cudaThreadSynchronize();
	copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->PML_Byx , hH->n);
	cudaThreadSynchronize();
	copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Te , hH->PML_Byz);
	cudaThreadSynchronize();
	copy_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->mat_mask , hH->PML_Byx);
	cudaThreadSynchronize();

//roll pointers:
  FL_DBL* new_Ez=hH->vz;
  FL_DBL* new_Ex=hH->vx;
	
  FL_DBL* new_Jz=hH->Ez;
  FL_DBL* new_Jx=hH->Ex;

	FL_DBL* new_By=hH->Jz;
	FL_DBL* new_n =hH->Jx;
	FL_DBL* new_PML_Byz=hH->By;
	FL_DBL* new_PML_Byx=hH->n;
	FL_DBL* new_Te=hH->PML_Byz;
	FL_DBL* new_mat_mask=hH->PML_Byx;
	FL_DBL* new_vx=hH->Te;
	FL_DBL* new_vz=hH->mat_mask;

  // reassign

	hH->Ez=new_Ez;
	hH->Ex=new_Ex;
	hH->Jz=new_Jz;
	hH->Jx=new_Jx;
	hH->By=new_By;
	hH->n =new_n ;
	hH->vz=new_vz;
	hH->vx=new_vx;
	hH->Te=new_Te;
	hH->mat_mask=new_mat_mask;
	hH->PML_Byz=new_PML_Byz;
	hH->PML_Byx=new_PML_Byx;

	hH->iz0 += hH->shft * hH->blockSize;
}
*/


extern "C"
std::string hydro2dHandler::get_description() {
	std::stringstream ss;
	ss << "grid: Lx=" << Lx << ", Lz=" << Lz << ", Nx=" << Nx << ", Nz=" << Nz << ", dx_1=" << dx_1 << " (dx=" << FPT(1.0)/dx_1 << "), dz_1=" << dz_1 << " (dz=" << FPT(1.0)/dz_1 << "), dt=" << dt << ",\ngeometry: " << (a_cell ? "from CELL":"from parameters") << " surfaceX=" << surfaceX << ", toothWith=" << toothWidth<< ", toothDepth=" << toothDepth << ", mediaDepth=" << mediaDepth << ", flip=" << (flip?"true":"false") << ", toothDirection=" << (toothDir?"true":"false") << ",\nmedia: n0=" << mediaN0 << ", T0=" << mediaTe0 << ", nu=" << mediaNu << ", diffusion=" << diffusion << ", Vdiffusion=" << Vdiffusion << "," << " mx_1=" << mx_1 << "," << " my_1=" << my_1 << "," << " mz_1=" << mz_1 << ", NUTratio=" << NUTratio << ",";
	if (doPolarization)
	ss << "\n       bound_w2=" << media_bound_w2 << " bound_beta=" << media_bound_beta << ", bound_gamma=" << media_bound_gamma << std::endl;
	if (doPhononAbsorbtion)
	ss << "\n       phOmega=" << media_phonon_omega << ", phWidth=" << media_phonon_phw << ", phBeta=" << media_phonon_beta << std::endl;
	ss << " current diffusion BC is : " << (neumann_current_diffusion ? "neumann\n" : "dirichlet\n");
	if (!doPolarization && !doPolarization)
		ss << std::endl;
	if (extSource)
		ss << "         therm source is external, srcX=" << srcX;
	else
		ss << "         source is electromagnetic";
	ss << ", Bxext=" << Bxext << ", Byext=" << Byext << ", Bzext=" << Bzext << std::endl;
	#ifdef STATIC
	ss << ", srcAperture=" << srcAperture;
	#endif
	ss << ", srcAmp=" << srcAmp << ", srcT=" << srcT;
	#ifdef STATIC
	ss << ", srcTshift=" << srcTshift;
	#endif
	ss <<  ", srcPhase=" << srcPhase << ", srcNosc=" << srcNosc;
	#ifdef STATIC
	ss << ", srcApertureTE=" << srcApertureTE;
	#endif
	ss << ", srcAmpTE=" << srcAmpTE << ", srcTTE=" << srcTTE;
	#ifdef STATIC
	ss << ", srcTshiftTE=" << srcTshiftTE;
	#endif
	ss <<  ", srcPhaseTE=" << srcPhaseTE << ", srcNoscTE=" << srcNoscTE << ", switchOnDelay=" << switchOnDelay << ", srcTfactor=" << srcTfactor << std::endl;
	ss << "         velocity=" << srcVelocity << std::endl;
	ss << "         electrons heated by " << sHeatType(JHEAT) << std::endl;
	ss << "         mode is " << (linear ? "linear" : "nonlinear") << std::endl;
	return ss.str();
}


extern "C"
int simpleGPUstep(hydro2dHandler* hH)
{
	if(cudaSetDevice(hH->device))
	{
		printf("ERROR Initializing device #%d\n",hH->device);
		return -1;
	}
	
	static std::mutex coutMutex;
	if(hH->t == 0)
	{
		coutMutex.lock();
		std::cout << "Hi this is first step for device " << hH->device << std::endl;
		std::cout << hH->get_description();
    std::cout << std::endl << std::flush; 

		coutMutex.unlock();

    if(cudaSetDevice(hH->device))
    {
      printf("ERROR Initializing device #%d\n",hH->device);
      return -1;
    }
  }
	
	dim3 grid(-floor(-FL_DBL(hH->Nx)/FL_DBL(hH->blockSize)),-floor(-FL_DBL(hH->Nz)/FL_DBL(hH->blockSize)),1);
	//printf("grid=%d %d %d\n", grid.x, grid.y, grid.z);
	dim3 threads(hH->blockSize,hH->blockSize,1);

	//std::cout << "grid: " << -floor(-FL_DBL(hH->Nx)/FL_DBL(hH->blockSize)) << ", " << -floor(-FL_DBL(hH->Nz)/FL_DBL(hH->blockSize)) << '\n' << std::flush;

	//std::cout << "threads: " << hH->blockSize << ", " << hH->blockSize << '\n' << std::flush;

	hH->tim.start("prepare TM source");
	FL_DBL delay = hH->switchOnDelay;
	FL_DBL maxSrcT=0;
	for(int id=0; id < hH->srcT.size(); id++)
		maxSrcT = hH->srcT[id] > maxSrcT ? hH->srcT[id] : maxSrcT;
	for(int id=0; id < hH->srcTTE.size(); id++)
		maxSrcT = hH->srcTTE[id] > maxSrcT ? hH->srcTTE[id] : maxSrcT;
	FL_DBL srcZ = maxSrcT * hH->srcVelocity;
	
	//std::ofstream of_debug("debug.txt");
	
	for(int i=0; i != hH->Nz; i++)
	{
		//zzd is calculated by substraction of big numbers that is why everything is converted to double to make sure we are not loosing precision
		//std::cout << "zzd=" << zzd << ", i=" << i << ", Nz=" << hH->Nz << '\n' << std::flush;
		if(hH->TMsource_callback) {
			hH->host_srct[i] = hH->TMsource_callback(hH->t, (double)(i+hH->iz0) / (double) hH->dz_1, hH->pywrapper);
			//of_debug << hH->host_srct[i] << '\n';
		} else {
			double zzd = (double)(i + hH->iz0 - hH->Nz + hH->blockSize + 1) / (double) hH->dz_1 + (double)(srcZ * 0.5) - hH->t * (double) hH->srcVelocity;
			hH->host_srct[i] = 0;
			for(int id=0; id < hH->srcAmp.size(); id++)
			{
				FL_DBL delay_factor = hH->t < delay ? FPT(0.5) - FPT(0.5) * cos(M_PI * hH->t / delay) : FPT(1.0);
				#ifndef STATIC
				FL_DBL srcZlocal = hH->srcT[id] * hH->srcVelocity;
				FL_DBL zz = (FL_DBL(zzd));
				FL_DBL envelope = (FPT(2.) * std::abs(zz) < std::abs(srcZ)) ? (FPT(0.5) + FPT(0.5) * cos(zz / srcZlocal * FPT(2.0) * M_PI)) : FPT(.0);
				if(hH->extSource)
					hH->host_srct[i] += hH->srcAmp[id] * delay_factor * pow(envelope, FPT(2.0));// * FPT(0.5) * (FPT(1.0) + sin(zz / srcZ * hH->srcNosc * M_PI * FPT(4.0)));
				else
					hH->host_srct[i] += hH->srcAmp[id] * delay_factor * envelope * sin(zz / srcZ * hH->srcNosc[id] * M_PI * FPT(2.0) + hH->srcPhase[id]);
				#else // STATIC
				FL_DBL localT = hH->t - hH->srcTshift[id];
				FL_DBL envelopeT = FPT(2.) * std::abs(localT) < std::abs(hH->srcT[id]) ? (FPT(.5) + FPT(.5) * cos(FPT(2.) * localT / hH->srcT[id] * M_PI)) : FPT(.0); 
				FL_DBL localZ = (FL_DBL) i / hH->dz_1 - hH->Lz * 0.5;
				FL_DBL envelopeAperture = exp(-pow(localZ/hH->srcAperture[id]*2.0, 8));//FPT(2.) * std::abs(localZ) < hH->srcAperture[id] ? FPT(.5) + FPT(.5) * cos(localZ / hH->srcAperture[id] * M_PI * FPT(2.)) : FPT(.0);
				hH->host_srct[i] += hH->srcAmp[id] * delay * envelopeT * envelopeAperture * sin((hH->t - hH->srcTshift[id]) / hH->srcT[id] * M_PI * FPT(2.0) * hH->srcNosc[id] + hH->srcPhase[id]);
				#endif
			}
		}
	}
	//of_debug.close();
	
	dev_h2d(hH, hH->host_srct, hH->srct, hH->Nz);
	hH->tim.stop("prepare source");
	if(hH->doPolarization)
	{
		hH->tim.start("doPolarization");
		simpleP_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->Ex, hH->Ptz, hH->Ptx, hH->Pz, hH->Px, hH->n, hH->dt, hH->media_bound_gamma, hH->media_bound_w2, hH->media_bound_beta, hH->mediaN0, hH->mat_mask);
		cudaThreadSynchronize();
		hH->tim.stop("doPolarization");
	}

	hH->tim.start("simpleE_Kernel");


	simpleE_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->Ex, hH->By, hH->Jz, hH->Jx, hH->ix0, hH->iz0, FPT(1.0)/hH->dx_1, hH->dx_1, 1.0/hH->dz_1, hH->dz_1, hH->dt, hH->extSource ? 0 : hH->srct, hH->SRCi, hH->PMLimin, hH->PMLimax, hH->PMLstrength, hH->PMLxmax, hH->Ez_mid, hH->Ex_mid);
	cudaThreadSynchronize();

	hH->tim.stop("simpleE_Kernel");

	// TE E step
	hH->tim.start("simpleE_TE_Kernel");
	simpleE_TE_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Bz, hH->Bx, hH->Ey, hH->Jy, hH->ix0, hH->iz0, FPT(1.)/hH->dx_1, hH->dx_1, FPT(1.)/hH->dz_1, hH->dz_1, hH->dt, hH->PMLimin, hH->PMLimax, hH->PMLstrength, hH->PMLxmax, hH->PML_Eyz, hH->PML_Eyx, hH->Ey_mid);
	hH->tim.stop("simpleE_TE_Kernel");

	hH->tim.start("prepare TE source");
	for(int i=0; i != hH->Nz; i++)
	{
		//zzd is calculated by substraction of big numbers that is why everything is converted to double to make sure we are not loosing precision
		if(hH->TEsource_callback)
			hH->host_srct[i] = hH->TEsource_callback(hH->t, (double)(i+hH->iz0) / (double) hH->dz_1, hH->pywrapper); 
		else {
			double zzd = (double)(i + hH->iz0 - hH->Nz + hH->blockSize + 1) / (double) hH->dz_1 + (double)(srcZ * 0.5) - hH->t * (double) hH->srcVelocity;
			hH->host_srct[i] = 0;
			for(int id=0; id < hH->srcAmpTE.size(); id++)
			{
				FL_DBL delay_factor = hH->t < delay ? FPT(0.5) - FPT(0.5) * cos(M_PI * hH->t / delay) : FPT(1.0);
				#ifndef STATIC
				FL_DBL srcZlocal = hH->srcTTE[id] * hH->srcVelocity;
				FL_DBL zz = (FL_DBL(zzd));
				FL_DBL envelope = (FPT(2.) * std::abs(zz) < std::abs(srcZ)) ? (FPT(0.5) + FPT(0.5) * cos(zz / srcZlocal * FPT(2.0) * M_PI)) : FPT(.0);
				if(hH->extSource)
					hH->host_srct[i] += hH->srcAmpTE[id] * delay_factor * envelope;// * FPT(0.5) * (FPT(1.0) + sin(zz / srcZ * hH->srcNosc * M_PI * FPT(4.0)));
				else
					hH->host_srct[i] += hH->srcAmpTE[id] * delay_factor * envelope * sin(zz / srcZ * hH->srcNoscTE[id] * M_PI * FPT(2.0) + hH->srcPhaseTE[id]);
				#else // STATIC
				FL_DBL localT = hH->t - hH->srcTshiftTE[id];
				FL_DBL envelopeT = FPT(2.) * std::abs(localT) < std::abs(hH->srcTTE[id]) ? (FPT(.5) + FPT(.5) * cos(FPT(2.) * localT / hH->srcTTE[id] * M_PI)) : FPT(.0); 
				FL_DBL localZ = (FL_DBL) i / hH->dz_1 - hH->Lz * 0.5;
				FL_DBL envelopeAperture = exp(-pow(localZ/hH->srcApertureTE[id]*2.0, 8));//FPT(2.) * std::abs(localZ) < hH->srcApertureTE[id] ? FPT(.5) + FPT(.5) * cos(localZ / hH->srcApertureTE[id] * M_PI * FPT(2.)) : FPT(.0);
				hH->host_srct[i] += hH->srcAmpTE[id] * delay * envelopeT * envelopeAperture * sin((hH->t - hH->srcTshiftTE[id]) / hH->srcTTE[id] * M_PI * FPT(2.0) * hH->srcNoscTE[id] + hH->srcPhaseTE[id]);
				#endif
			}
		}
	}
	dev_h2d(hH, hH->host_srct, hH->srct, hH->Nz);
	hH->tim.stop("prepare TE source");
	
	hH->tim.start("simpleB_Kernel");
	simpleB_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Ez, hH->Ex, hH->By, hH->ix0, hH->iz0, 1.0/ hH->dx_1, hH->dx_1, 1.0/ hH->dz_1, hH->dz_1, hH->dt, hH->PMLimin, hH->PMLimax, hH->PMLstrength, hH->PMLxmax, hH->PML_Byz, hH->PML_Byx);
	cudaThreadSynchronize();
	hH->tim.stop("simpleB_Kernel");
	
	// TE B step
	hH->tim.start("simpleB_TE_Kernel");
	simpleB_TE_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Bz, hH->Bx, hH->Ey, hH->ix0, hH->iz0, FPT(1.)/hH->dx_1, hH->dx_1, FPT(1.)/hH->dz_1, hH->dz_1, hH->dt, hH->extSource ? 0 : hH->srct, hH->SRCi, hH->PMLimin, hH->PMLimax, hH->PMLstrength, hH->PMLxmax);
	cudaThreadSynchronize();
	hH->tim.stop("simpleB_TE_Kernel");

	if(!(hH->linear) && IFNONLIN == ON)
	{
		hH->tim.start("simpleV_Kernel");
		simpleV_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->vx, hH->vy, hH->vz, hH->Jx, hH->Jy, hH->Jz, hH->n, hH->mediaN0);
		cudaThreadSynchronize();
		hH->tim.stop("simpleV_Kernel");
		
		hH->tim.start("simpleNl_Kernel");
		simpleNl_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->nlx, hH->nly, hH->nlz, hH->vx, hH->vy, hH->vz, hH->Jx, hH->Jy, hH->Jz, hH->Bx, hH->By, hH->Bz, hH->dz_1, hH->dx_1, hH->mx_1, hH->my_1,  hH->mz_1, hH->Bxext, hH->Byext, hH->Bzext);
		cudaThreadSynchronize();
		hH->tim.stop("simpleNl_Kernel");
	}
	hH->tim.start("simpleN_Kernel");
	simpleN_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->n, hH->Jx, hH->Jz, hH->dz_1, hH->dx_1, hH->dt, hH->mat_mask);
	cudaThreadSynchronize();
	hH->tim.stop("simpleN_Kernel");

	if(!(hH->linear))
	{
		hH->tim.start("simple_dTe_Kernel");
		if(hH->JHEAT == hydro2dHandler::eHEATTYPE::JJ)
			simple_dTe_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->diffusion, hH->NUTratio, hH->mediaTe0, hH->dTe, hH->Te, hH->n, hH->Jx, hH->Jz, hH->Jx, hH->Jz, hH->dz_1, hH->dx_1, hH->extSource ? hH->srct : 0, hH->srcx, hH->srcTfactor * hH->mediaNu / 1.5 / hH->mediaN0, FPT(1.), hH->mat_mask); // note /1.5 is a fenomenological coefficient
		else if(hH->JHEAT == hydro2dHandler::eHEATTYPE::EE)
			simple_dTe_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->diffusion, hH->NUTratio, hH->mediaTe0, hH->dTe, hH->Te, hH->n, hH->Ex_mid, hH->Ez_mid, hH->Ex_mid, hH->Ez_mid, hH->dz_1, hH->dx_1, hH->extSource ? hH->srct : 0, hH->srcx, hH->srcTfactor, FPT(.0), hH->mat_mask);
		else if(hH->JHEAT == hydro2dHandler::eHEATTYPE::JE)
			simple_dTe_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->diffusion, hH->NUTratio, hH->mediaTe0, hH->dTe, hH->Te, hH->n, hH->Ex_mid, hH->Ez_mid, hH->Jx, hH->Jz, hH->dz_1, hH->dx_1, hH->extSource ? hH->srct : 0, hH->srcx, hH->srcTfactor, FPT(.0), hH->mat_mask);
		
		cudaThreadSynchronize();
		hH->tim.stop("simple_dTe_Kernel");
	
		hH->tim.start("simple_Te_Kernel");
		simple_Te_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->dTe, hH->Te, hH->dt, hH->mediaTe0);
		cudaThreadSynchronize();
		hH->tim.stop("simple_Te_Kernel");
	}

	hH->tim.start("simpleJ_Kernel");
	simpleJ_Kernel<<< grid, threads >>>(hH->Nx, hH->Nz, hH->Jx, hH->Jy, hH->Jz, hH->Ex, hH->Ey, hH->Ez, hH->nlx, hH->nly, hH->nlz, hH->vx, hH->vy, hH->vz, hH->n, hH->Te, hH->NUTratio, hH->mediaTe0, hH->dx_1, hH->dz_1, hH->dt, hH->mediaNu, hH->landauDamping, hH->Vdiffusion, hH->Phtz, hH->Phtx, hH->Phz, hH->Phx, hH->mx_1,  hH->my_1, hH->mz_1, hH->media_phonon_omega, hH->media_phonon_phw, hH->media_phonon_beta, hH->doPhononAbsorbtion, hH->mediaN0, hH->mat_mask, hH->linear, hH->neumann_current_diffusion);
	cudaThreadSynchronize();
	hH->tim.stop("simpleJ_Kernel");

	hH->t += (double) hH->dt;
	hH->step++;

//  roll(hH);

	#ifndef STATIC
	if(hH->srcVelocity > 0 && hH->t * (double) hH->srcVelocity - (double) hH->iz0 / (double) hH->dz_1 > (double) hH->blockSize / (double) hH->dz_1)
	{
		hH->tim.start("feed");
		//manage feed
		FL_DBL* feed = hH->host_feed;//new FL_DBL[hH->Nx * (hH->shft * hH->blockSize)];
		
		for(int i=0;i<hH->Nx * (hH->shft * hH->blockSize);i++) feed[i] = FPT(0.0);

		dev_h2d(hH, feed, hH->feedPML_Byz, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedPML_Byx, hH->Nx * (hH->shft * hH->blockSize));
		//TE
		dev_h2d(hH, feed, hH->feedPML_Eyz, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedPML_Eyx, hH->Nx * (hH->shft * hH->blockSize));
		
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

		//TE
		dev_h2d(hH, feed, hH->feedBz, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedBx, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedEy, hH->Nx * (hH->shft * hH->blockSize));


		if(hH->doPhononAbsorbtion)
		{
			dev_h2d(hH, feed, hH->feedPhtz, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPhtx, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPhz, hH->Nx * (hH->shft * hH->blockSize));
			dev_h2d(hH, feed, hH->feedPhx, hH->Nx * (hH->shft * hH->blockSize));
		}
		dev_h2d(hH, feed, hH->feedJz, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedJy, hH->Nx * (hH->shft * hH->blockSize));
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
				double z = (double) ireal * (double) dz;
				for(int cell_id(imin); cell_id != imax; cell_id++)
				{
					double z_cell_local = std::fmod(z, c->get_lz());
					z_cell_local += (double) cell_id * c->get_lz();
					for(int j=0;j<hH->Nx;j++)
					{
						double x = (double) j * (double) dx + (double) dx * 0.5 - (double) hH->surfaceX;
						//std::cout << "calling inside(" << z_cell_local << ", " << x << ")" << std::endl;
						bool res = c->inside(z_cell_local, x);
						ins += res ? 1 : 0;
						outs += res ? 0 : 1;
						//std::cout << (res ? "got yes" : "got no") << std::endl;
						feed[i * hH->Nx + j] = res ? FPT(1.0) : feed[i * hH->Nx + j];
						//feed[i * hH->Nx + j] = res ? hH->mediaN0 : feed[i * hH->Nx + j];
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
				FL_DBL levelL = hH->surfaceX;
				FL_DBL levelH = levelL + hH->mediaDepth;

				int toothWidthi = hH->toothWidth * hH->dz_1;
				int toothDepthi = hH->toothDepth * hH->dx_1;
				int profile = FL_DBL(((hH->toothDir) ? ireal : 100000000 - ireal) % toothWidthi) / FL_DBL(toothWidthi) * FL_DBL(toothDepthi);

				if(hH->flip)
				{
					for(int j=0;j<hH->Nx;j++)
						feed[i * hH->Nx + j]=((j < levelH * hH->dx_1 + profile) && (j > levelL * hH->dx_1)) ? FPT(1.0) : FPT(0.0);
				}
				else
				{
					for(int j=0;j<hH->Nx;j++)
						feed[i * hH->Nx + j]=((j < levelH * hH->dx_1 + toothDepthi) && (j > levelL * hH->dx_1 + profile)) ? FPT(1.0) : FPT(0.0);
				}
			}
		}

		dev_h2d(hH, feed, hH->feed_mat_mask, hH->Nx * (hH->shft * hH->blockSize));

		for(int i=0;i<hH->shft * hH->blockSize;i++)
		{
			for(int j=0;j<hH->Nx;j++)
				feed[i * hH->Nx + j] = FPT(0.0);//(feed[i * hH->Nx + j] != FPT(0.0))? hH->mediaTe0 : FPT(0.0);
		}
		dev_h2d(hH, feed, hH->feedn, hH->Nx * (hH->shft * hH->blockSize));
		dev_h2d(hH, feed, hH->feedTe, hH->Nx * (hH->shft * hH->blockSize));
		shift(hH);
		hH->tim.stop("feed");
	}
	#endif // STATIC
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
		dev_alloc(hH, (void**)&(hH->Pty), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Ptx), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Pz), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Py), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Px), hH->Nx * hH->Nz);
	}
	
	dev_alloc(hH, (void**)&(hH->Ez), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Ey), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Ex), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Ez_mid), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Ey_mid), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Ex_mid), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Bz), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->By), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Bx), hH->Nx * hH->Nz);

	if(hH->doPhononAbsorbtion)
	{
		printf("allocating phonon absorbtion arrays\n");
		dev_alloc(hH, (void**)&(hH->Phtz), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Phty), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Phtx), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Phz), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Phy), hH->Nx * hH->Nz);
		dev_alloc(hH, (void**)&(hH->Phx), hH->Nx * hH->Nz);
	}

	dev_alloc(hH, (void**)&(hH->Jx), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Jy), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Jz), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->mat_mask), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->n), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->Te), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->dTe), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->vx), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->vy), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->vz), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->nlx), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->nly), hH->Nx * hH->Nz);
	dev_alloc(hH, (void**)&(hH->nlz), hH->Nx * hH->Nz);
	
	//alloc additional fields for PML conditions
	dev_alloc(hH, (void**)&(hH->PML_Byx), hH->Nz * hH->Nx);
	dev_alloc(hH, (void**)&(hH->PML_Byz), hH->Nz * hH->Nx);
	dev_alloc(hH, (void**)&(hH->PML_Eyx), hH->Nz * hH->Nx);
	dev_alloc(hH, (void**)&(hH->PML_Eyz), hH->Nz * hH->Nx);
	
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

	//TE
	dev_alloc(hH, (void**)&(hH->feedBz), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedBx), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedEy), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedPML_Eyz), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedPML_Eyx), hH->Nx * (hH->shft * hH->blockSize));

	
	if(hH->doPhononAbsorbtion)
	{
		dev_alloc(hH, (void**)&(hH->feedPhtz), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPhtx), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPhz), hH->Nx * (hH->shft * hH->blockSize));
		dev_alloc(hH, (void**)&(hH->feedPhx), hH->Nx * (hH->shft * hH->blockSize));
	}
	
	dev_alloc(hH, (void**)&(hH->feedJz), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedJy), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedJx), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feed_mat_mask), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedn), hH->Nx * (hH->shft * hH->blockSize));
	dev_alloc(hH, (void**)&(hH->feedTe), hH->Nx * (hH->shft * hH->blockSize));
	//srcx array to speed up calculations
	dev_alloc(hH, (void**)&(hH->srcx), hH->Nx);
	dev_alloc(hH, (void**)&(hH->srct), hH->Nz);
	//hH->srct = 0;
	//allocate some host arrays for external source
	hH->host_feed = new FL_DBL[hH->Nx * (hH->shft * hH->blockSize)];
	std::cout << "allocating host_srct of size " << hH->Nz << '\n' << std::flush;
	hH->host_srct = new FL_DBL[hH->Nz];	
		
	printf("\tdone\n");
}

extern "C"
int CUDA_device_count()
{
  int res;
  if(cudaGetDeviceCount(&res))
  {
    printf("ERROR getting device count\n");
    return -1;
  }
  return res;
}

extern "C"
void simpleGPUinit(hydro2dHandler* hH) // this return handler for operation
{
	if(cudaSetDevice(hH->device))
	{
		printf("ERROR Initializing device #%d\n",hH->device);
		return;
	}
	std::cout << "simpleGPUinit\n";
	hH->PMLimin = -floor(- hH->PMLxmax * hH->dx_1);
	hH->PMLimin = ((hH->PMLimin / hH->blockSize)+(((hH->PMLimin % hH->blockSize)==0)?0:1)) * hH->blockSize;
	hH->PMLimax = hH->Nx-hH->PMLimin;
	hH->SRCi = hH->PMLimin + hH->Nx/32;//floor(hH->surfaceX * hH->dx_1) - hH->Nx/16;
//	hH->SRCi = hH->PMLimin+1; // usefull for visualization
	alloc_main_fields(hH);

	// clear helpers
	FL_DBL* zero=new FL_DBL[hH->Nx * hH->Nz];
	for(int i=0;i<hH->Nx * hH->Nz;i++) zero[i]=0.0;
	dev_h2d(hH, zero, hH->PML_Byx, hH->Nz * hH->Nx);
	dev_h2d(hH, zero, hH->PML_Byz, hH->Nz * hH->Nx);
	dev_h2d(hH, zero, hH->PML_Eyx, hH->Nz * hH->Nx);
	dev_h2d(hH, zero, hH->PML_Eyz, hH->Nz * hH->Nx);
	dev_h2d(hH, zero, hH->Jz, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->Jy, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->Jx, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->Ez, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->Ex, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->Ey, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->Bz, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->Bx, hH->Nx * hH->Nz);
	dev_h2d(hH, zero, hH->By, hH->Nx * hH->Nz);
	//we do not care of v and nl arrays.
	FL_DBL* n0 = new FL_DBL[hH->Nx * hH->Nz];
	FL_DBL* T0 = new FL_DBL[hH->Nx * hH->Nz];
	FL_DBL* mat_mask = new FL_DBL[hH->Nx * hH->Nz];

	//set the source x:
	FL_DBL* srcx = new FL_DBL[hH->Nx];	
	for(int j = 0; j != hH->Nx; j++)
	{
		FL_DBL x = FL_DBL(j) / hH->dx_1;
		srcx[j] = exp(-pow((x - hH->surfaceX) / hH->srcX, 2.0));
	}
	dev_h2d(hH, srcx, hH->srcx, hH->Nx);
	
	for(int i = 0; i != hH->Nz; i++)
	{
		int ireal=hH->iz0+i;
		FL_DBL z=FL_DBL(ireal) / hH->dz_1;
		//FL_DBL levelL=hH->Lx/3.0;
		FL_DBL levelL=hH->surfaceX;
		FL_DBL levelH=hH->surfaceX;
		int leveli=(z/(hH->Lz/13.0)-floor(z/(hH->Lz/13.0))>0.5)?int(floor(levelH*hH->dx_1)):int(floor(levelL*hH->dx_1));
		for(int j = 0; j != hH->Nx; j++)
		{
      mat_mask[i * hH->Nx + j] = FPT(0.0);
			n0[i * hH->Nx + j] = FPT(0.0);//((j>leveli)&&(j<hH->Nx+hH->Nx/8)) ? hH->mediaN0 : FPT(0.0);
			T0[i * hH->Nx + j] = FPT(0.0);//(n0[i * hH->Nx + j] != FPT(0.0)) ? hH->mediaTe0 : FPT(0.0);
		}
	}
//0.01 * sin(double(i)/double(hH->Nz)*hH->Lz*0.001);

	#ifdef STATIC
	std::vector<FL_DBL> buff(hH->Nx * hH->Nz, FPT(.0));
	if(hH->a_cell)
	{
		//std::cout << "feeding from cell, step is " << hH->step << std::endl;
		int ins(0);
		int outs(0);
		const cell* c = hH->a_cell;
		const int imin = std::floor(c->get_zmin() / c->get_lz()); 
		const int imax = std::ceil(c->get_zmax() / c->get_lz()); 
		const double dz = FL_DBL(1.0) / hH->dz_1;
		const double dx = FL_DBL(1.0) / hH->dx_1;
		for(int i = hH->PMLimin+3; i<hH->Nz-hH->PMLimin-3; i++)
		{
			double z = (double) i * (double) dz;
			double z_cell_local = std::fmod(z, c->get_lz());
			for(int j=0;j<hH->Nx;j++)
			{
				double x = (double) j * (double) dx + (double) dx * 0.5 - (double) hH->surfaceX;
				//std::cout << "calling inside(" << z_cell_local << ", " << x << ")" << std::endl;
				bool res = c->inside(z_cell_local, x);
				ins += res ? 1 : 0;
				outs += res ? 0 : 1;
				//std::cout << (res ? "got yes" : "got no") << std::endl;
				mat_mask[i * hH->Nx + j] = res ? FPT(1.0) : mat_mask[i * hH->Nx + j];
				//feed[i * hH->Nx + j] = res ? hH->mediaN0 : feed[i * hH->Nx + j];
			}
		}
		//std::cout << "done, ins = " << ins << ", outs = " << outs << std::endl;
	}
	else
	{
		for(int i = hH->PMLimin+3; i<hH->Nz-hH->PMLimin-3; i++)
		{
			FL_DBL z = FL_DBL(i) / hH->dz_1;
			FL_DBL levelL = hH->surfaceX;
			FL_DBL levelH = levelL + hH->mediaDepth;

			int toothWidthi = hH->toothWidth * hH->dz_1;
			int toothDepthi = hH->toothDepth * hH->dx_1;
			int profile = FL_DBL(((hH->toothDir) ? i : 100000000 - i) % toothWidthi) / FL_DBL(toothWidthi) * FL_DBL(toothDepthi);
			if(hH->flip)
			{
				for(int j=0;j<hH->Nx;j++)
					mat_mask[i * hH->Nx + j]=((j < levelH * hH->dx_1 + profile) && (j > levelL * hH->dx_1)) ? FPT(1.0) : FPT(0.0);
			}
			else
			{
				for(int j=0;j<hH->Nx;j++)
					mat_mask[i * hH->Nx + j]=((j < levelH * hH->dx_1 + toothDepthi) && (j > levelL * hH->dx_1 + profile)) ? FPT(1.0) : FPT(0.0);
			}
		}
	}
	
	double TdistortionAmp = 0.0001;//0.0003;
	for(int i = 0; i != hH->Nz; i++)
	{
		double z = double(i-hH->Nz/2)/double(hH->Nz)*hH->Lz;
		double value = (double)rand()/(double)RAND_MAX;
		for(int j = 0; j != hH->Nx; j++)
		{
      if(mat_mask[i * hH->Nx + j])
				T0[i * hH->Nx + j] = TdistortionAmp * value * exp(-abs((double(j)/hH->dx_1-hH->surfaceX)/300.));//0.0003 * exp(-pow(z/53000.,4)) * pow(cos(z/hH->toothWidth*M_PI),2);
		}
	}
	#endif

	dev_h2d(hH, mat_mask, hH->mat_mask, hH->Nx * hH->Nz);
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
hydro2dHandler::hydro2dHandler(const hydro2dHandler &obj, int dev) : tim(obj.tim),
																																		 device((dev == -1) ? obj.device : dev),
																																		 Nx(obj.Nx),
																																		 Nz(obj.Nz),
																																		 Lx(obj.Lx),
																																		 surfaceX(obj.surfaceX),
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
																																		 srcT(obj.srcT),
																																		 srcTTE(obj.srcTTE),
																																		 srcTshift(obj.srcTshift),
																																		 srcTshiftTE(obj.srcTshiftTE),
																																		 srcNosc(obj.srcNosc),
																																		 srcNoscTE(obj.srcNoscTE),
																																		 srcAperture(obj.srcAperture),
																																		 srcApertureTE(obj.srcApertureTE),
																																		 srcAmp(obj.srcAmp),
																																		 srcAmpTE(obj.srcAmpTE),
																																		 srcPhase(obj.srcPhase),
																																		 srcPhaseTE(obj.srcPhaseTE),
																																		 mediaN0(obj.mediaN0),
																																		 mediaTe0(obj.mediaTe0),
																																		 diffusion(obj.diffusion),
																																		 Vdiffusion(obj.Vdiffusion),
																																		 landauDamping(obj.landauDamping),
																																		 NUTratio(obj.NUTratio),
																																		 mx_1(obj.mx_1),
																																		 my_1(obj.my_1),
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
																																		 JHEAT(obj.JHEAT),
                                                                     srcX(obj.srcX),
                                                                     Bxext(obj.Bxext),
                                                                     Byext(obj.Byext),
                                                                     Bzext(obj.Bzext),
																																		 srcTfactor(obj.srcTfactor),
																																		 switchOnDelay(obj.switchOnDelay),
																																		 linear(obj.linear),
																																		 neumann_current_diffusion(obj.neumann_current_diffusion)
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
		COPY_FIELD(Pty, whole)
		COPY_FIELD(Ptz, whole)
		COPY_FIELD(Px, whole)
		COPY_FIELD(Py, whole)
		COPY_FIELD(Pz, whole)
	}
	
	COPY_FIELD(Bx, whole)
	COPY_FIELD(By, whole)
	COPY_FIELD(Bz, whole)
	COPY_FIELD(Ex, whole)
	COPY_FIELD(Ey, whole)
	COPY_FIELD(Ez, whole)
	COPY_FIELD(PML_Byx, whole)
	COPY_FIELD(PML_Byz, whole)
	COPY_FIELD(PML_Eyx, whole)
	COPY_FIELD(PML_Eyz, whole)
	
	if(doPhononAbsorbtion)
	{
		COPY_FIELD(Phtx, whole)
		COPY_FIELD(Phty, whole)
		COPY_FIELD(Phtz, whole)
		COPY_FIELD(Phx, whole)
		COPY_FIELD(Phy, whole)
		COPY_FIELD(Phz, whole)
	}
	
	COPY_FIELD(Jz, whole)
	COPY_FIELD(Jy, whole)
	COPY_FIELD(Jx, whole)
	COPY_FIELD(n, whole)
	COPY_FIELD(Te, whole)
	COPY_FIELD(mat_mask, whole)

	#undef COPY_FIELD

	delete [] tmp;
}
