#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_timer.h>
#include "funcall.h"

extern "C"
void dev_alloc(hydro2dHandler* hH, void** pparr, int sz)
{
	if(cudaSetDevice(hH->device))
		printf("ERROR Initializing device #%d\n",hH->device);
	else
	{
		checkCudaErrors( cudaMalloc( pparr, sz*sizeof(FL_DBL)));
		cudaThreadSynchronize();
	}
}

extern "C"
void dev_h2d(hydro2dHandler* hH, FL_DBL* host_arr, FL_DBL* dev_arr, int sz)
{
	//printf("\tcopying source of size %d from %p to %p\n", sz, host_arr, dev_arr);

	if(cudaSetDevice(hH->device))
		printf("ERROR Initializing device #%d\n",hH->device);
	else
	{
		checkCudaErrors( cudaMemcpy( dev_arr, host_arr, sz*sizeof(FL_DBL), cudaMemcpyHostToDevice));
		cudaThreadSynchronize();
	}
	//printf("\t\tdone\n");

}

extern "C"
void dev_d2h(const hydro2dHandler* hH, const FL_DBL* dev_arr, FL_DBL* host_arr, int sz)
{
	if(cudaSetDevice(hH->device))
		printf("ERROR Initializing device #%d\n",hH->device);
	else
	{
		checkCudaErrors( cudaMemcpy( host_arr, dev_arr, sz*sizeof(FL_DBL), cudaMemcpyDeviceToHost));
		cudaThreadSynchronize();
	}
}

extern "C"
void GPUsetField(hydro2dHandler* hH, FL_DBL* E_z, FL_DBL* E_x, FL_DBL* B_y)
{
	if(cudaSetDevice(hH->device))
		printf("ERROR Initializing device #%d\n",hH->device);
	else
	{
		dev_h2d(hH, E_x, hH->Ex, hH->Nx * hH->Nz);
		dev_h2d(hH, E_z, hH->Ez, hH->Nx * hH->Nz);
		dev_h2d(hH, B_y, hH->By, hH->Nx * hH->Nz);
	}
}

extern "C"
void GPUgetField(hydro2dHandler* hH, FL_DBL* E_z, FL_DBL* E_x, FL_DBL* B_y)
{
	if(cudaSetDevice(hH->device))
		printf("ERROR Initializing device #%d\n",hH->device);
	else
	{
		dev_d2h(hH, hH->Ex, E_x, hH->Nx * hH->Nz);
		dev_d2h(hH, hH->Ez, E_z, hH->Nx * hH->Nz);
		dev_d2h(hH, hH->By, B_y, hH->Nx * hH->Nz);
	}
}

