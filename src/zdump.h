#include <zlib.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include "defaults.h"

struct data2d
{
	int Nx;
	int Nz;
	FL_DBL dx;
	FL_DBL dz;
};

//zlib
struct dump
{
	#define CHUNK 16386
	dump() : dumpbuff(), file(), zret(0), zstrm()
	{
		zbuff = new char[CHUNK];
	}
	~dump()
	{
		delete[] zbuff;
	}
	std::string dumpbuff;
	std::ofstream file;
	//zlib stuff:
	char* zbuff;	
	int zret;
	z_stream zstrm;
	
	void myflush(std::string& s, bool last = false)
	{
		size_t i(0);
		if(s.size() >= CHUNK)
		{
			for(; i != s.size() / CHUNK; i++)
			{
				zstrm.avail_in = CHUNK;
				zstrm.next_in = (unsigned char*)&s.c_str()[CHUNK *i];
				
				zstrm.avail_out = 0;
				while(zstrm.avail_out == 0)
				{
					zstrm.avail_out = CHUNK;
					zstrm.next_out = (unsigned char*) zbuff;
					deflate(&zstrm, Z_NO_FLUSH);
					if(CHUNK - zstrm.avail_out)
						file.write(zbuff, CHUNK - zstrm.avail_out);
				}
			}
			s.erase(0, CHUNK * (s.size() / CHUNK));
		}
		if(last)
		{
			zstrm.avail_in = s.size();
			zstrm.next_in = (unsigned char*) s.c_str();
			zstrm.avail_out = 0;
			while(zstrm.avail_out == 0)
			{
				zstrm.avail_out = CHUNK;
				zstrm.next_out = (unsigned char*) zbuff;
				deflate(&zstrm, Z_FINISH);
				if(CHUNK - zstrm.avail_out)
					file.write(zbuff, CHUNK - zstrm.avail_out);
			}
			s.erase();
		}
		return;
	}
	
	void dumpData(std::string filename, const data2d* d2d, std::vector<const FL_DBL*> data, std::vector<std::string> header)
	{
		file.open((filename+".gz").c_str());
		file.close();
		file.open((filename+".gz").c_str(), std::ios::out | std::ios::app | std::ios::binary);

		zstrm.zalloc = Z_NULL;
    zstrm.zfree = Z_NULL;
    zstrm.opaque = Z_NULL;
    /*this is a setup for gzip*/
    int windowBits = 15;
		int GZIP_ENCODING = 16;
    zret = deflateInit2(&zstrm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits | GZIP_ENCODING, 8, Z_DEFAULT_STRATEGY);
    
    if (zret != Z_OK)
    {
			printf("dump::WARNING: failed to initialize zlib\n");
      return;
		}
		char message[100];
		//dump header first:
		sprintf(message,"%d\t%d\t%ld", hH->Nz, hH->Nx, data.size());
		dumpbuff += std::string(message) + "\n";
		dumpbuff += "z\tx";
		size_t i=0;
		for(; i != std::min(data.size(), header.size()); i++)
		{
			sprintf(message,"%s", ("\t" + header[i]).c_str());
			dumpbuff += std::string(message);
		}
		if(i != data.size())
			printf("dump::WARNING: data size and header does not match!!!\n");
		for(; i != data.size(); i++)
			dumpbuff += "\tunknown";

		dumpbuff += "\n";
		//dump raw data:
		for(int i = 0; i != d2d->Nz; i++)
		{
			for(int j = 0; j != d2d->Nx; j++)
			{
				sprintf(message, "%g\t%g", FL_DBL(i) * d2d->dz, FL_DBL(j) * d2d->dx);
				dumpbuff += std::string(message);
				for(size_t k = 0; k != data.size(); k++)
				{
					sprintf(message, "\t%g", data[k][i*hH->Nx+j]);
					dumpbuff += std::string(message);
				}
				dumpbuff += "\n";
			}
			myflush(dumpbuff); // flush data periodically to save memory allocated for mybuff
		}
		myflush(dumpbuff,true);
		deflateEnd(&zstrm);
		file.close();
		return;
	}
	#undef CHUNK
};


int main()
{
	hydro2dHandler *hH = new hydro2dHandler;
	hH->Nx=128;
	hH->Nz=128;
	hH->dx_1=(FL_DBL) hH->Nx;
	hH->dz_1=(FL_DBL) hH->Nz;
	FL_DBL* arrg = new FL_DBL[hH->Nx * hH->Nz];
	for(int i=0;i!=hH->Nz;i++)
	{
		FL_DBL z = FL_DBL(i - hH->Nz/2)/hH->dz_1;
		for(int j=0;j!=hH->Nx;j++)
		{
			FL_DBL x = FL_DBL(j - hH->Nx/2)/hH->dx_1;
			arrg[i*hH->Nx+j] = x*x + z*z;
		}
	}
	
	FL_DBL* arrh = new FL_DBL[hH->Nx * hH->Nz];
	for(int i=0;i!=hH->Nz;i++)
	{
		FL_DBL z = FL_DBL(i - hH->Nz/2)/hH->dz_1;
		for(int j=0;j!=hH->Nx;j++)
		{
			FL_DBL x = FL_DBL(j - hH->Nx/2)/hH->dx_1;
			arrh[i*hH->Nx+j] = x*x - z*z;
		}
	}	
	
	printf("filled array\n");
	dump d;
	d.dumpData("test.dat",hH,std::vector<const FL_DBL*> {arrg, arrh}, std::vector<std::string> {"paraboloid", "hyperboloid"});
	return 0;
}
