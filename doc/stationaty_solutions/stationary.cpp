#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <simpledraw.h>

double Fz(double z, double T)
{
	double shift = T * 0.5 + T * 0.1;
  return abs(z - shift) < T * 0.5 ? (cos((z - shift) / (T * 0.5) * M_PI) + 1.0) : .0;
}

int main(int argc, char* argv[])
{
  int nx = 512;
  int nz = 1024*16;
  int isurf_top=nx/4;
  int isurf_bot=4*nx/5;
  double T = 72000.*0.5;
  double V = 1./cos(45./180.*M_PI);
  double Lz = 4. * T;
  double dz = Lz / (double) nz;
  double dx = dz / sqrt(V * V - 1.);
  double N0 = 0.2e-6;
  double Lx = dx * (double) nx;
  std::cout << "Lx=" << Lx << " x Lz=" << Lz << std::endl;
  /*
   *  V \partial j_x / \partial z = -n E_x + \partial F / \partial x
   *  V \partial j_z / \partial z = -n E_z + \partial F / \partial z
   *  V \partial E_z / \partial z = j_z - \partial B / \partial x
   *  (V^2-1) \partial E_x / \partial z = V  j_x - \partial E_z / \partial x
   *  (V^2-1) \partial B / \partial z = j_x - V \partial E_z / \partial x
   */
  
  // refer to 1d_boundary_test.asy for more details on dz and dx relation
  // and open boundary condition
  
  /*
                        vacuum                                          
________________________________________^___________________interface   
                          |       here  | jx=0    |                     
                          |                       |                     
                          |                       |                     
  Ez0         jz0        Ez+         jz+          |                     
  -->         -->________-->_________-->          |                     
              |           |            |          |                     
              |           |            |          |                     
              |           |            |          |                     
              |           |            |          |                     
  ^           ^ Ex0,       ^ __________|__________|                     
  | jx0       | B0         | jx+       |Ex+,                            
              |                        | B+                             
              |                        |                                
              |                        |                                
              |________________________|                              
              F0                        F+                              
                                                                        
  */
  
  double* jx = new double[nx+1];
  double* Ez = new double[nx];
  double* jz = new double[nx];
  double* Ex = new double[nx+1];
  double* B = new double[nx+1];
  double* F = new double[nx];
  double* n = new double[nx+1];
  
  for(int i=0; i!=nx; i++)
  {
    Ez[i]=.0;
    jz[i]=.0;
    F[i]=.0;
    n[i]=.0;
  }
  for(int i=0; i!=nx+1; i++)
  {
    jx[i]=.0;
    Ex[i]=.0;
    B[i]=.0;
  }
  
  for(int i=isurf_top;i!=isurf_bot;i++)
  {
    F[i] = exp(- (double) (i - isurf_top) * dx / Lx * 16.0);
    n[i] = N0;
  }
  std::vector<double> signal;
  signal.resize(nz);
  std::cout << "sucessfully inited\n";
  std::vector<std::vector<double>> picture;
  picture.resize(4);
  for(int j=0; j!=nz; j++)
  {
    jx[0] += -n[0] * Ex[0] / V * dz;
    for(int i=1; i!=nx+1; i++)
      jx[i] += (-(n[i-1] + n[i]) * 0.5 * Ex[i] + Fz((double) j * dz, T) * (F[i] - F[i-1]) / dx) / V * dz;
    jx[isurf_top]=.0;
    jx[isurf_bot+1]=.0;
    
    for(int i=0; i!=nx; i++)
      Ez[i] += (jz[i] - (B[i+1] - B[i]) / dx) / V * dz;
    
    for(int i=0; i!=nx; i++)
      jz[i] += (-n[i] * Ez[i] 
               + (Fz(((double) j + 1.0) * dz, T) * F[i]
                  - Fz((double) j * dz, T) * F[i]) / dz) / V * dz;
    
    
    //no need to calculate boundary values for Ex
    for(int i=1; i!=nx+1; i++)
      Ex[i] += (V * jx[i] - (Ez[i] - Ez[i-1]) / dx) / (V * V - 1.) * dz;
    
    
    B[0]  += (jx[0] - V * (Ez[0] + B[0] / V * sqrt(V * V - 1.)) / dx)  / (V * V - 1.) * dz;
    for(int i=1; i!=nx; i++)
      B[i]  += (jx[i] - V * (Ez[i] - Ez[i-1]) / dx)  / (V * V - 1.) * dz;
    B[nx]  += (jx[nx] - V * (B[nx] / V * sqrt(V * V - 1.) - Ez[nx-1]) / dx)  / (V * V - 1.) * dz;
    
    for(int i=0; i!=nx; i++)
			picture[0].push_back(Ez[i]);
    for(int i=0; i!=nx; i++)
			picture[2].push_back(jz[i]);
    for(int i=0; i!=nx+1; i++)
			picture[1].push_back(Ex[i]);
    for(int i=0; i!=nx+1; i++)
			picture[3].push_back(jx[i]);
    //store signal
    signal[j] = B[isurf_top];
  }
  std::cout << nz << " steps done\n";
  std::ofstream fl("signal.dat");
  for(const auto& v : signal)
    fl << v << std::endl;
  fl.close();
  
  fadey_init(nx,nz,4);
  fadey_draw(picture[0].data(),nx,nz,0);
  fadey_draw(picture[1].data(),nx+1,nz,1);
  fadey_draw(picture[2].data(),nx,nz,2);
  fadey_draw(picture[3].data(),nx+1,nz,3);
  char a;
  std::cin >> a;
  fadey_close();
  return 0;
}
