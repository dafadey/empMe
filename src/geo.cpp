#include "svg.h"
#include "geometry.h"
#include "geo.h"
#include <limits>
#include <iostream>
#include <vector>

cell::~cell()
{
	if(!has_parent)
		delete geo;
}

bool cell::inside(double x, double y) const
{
	XY p((double)x, (double)y);
//	return geo->is_inside(p);
	return geo->is_inside_simple(p);
}

bool cell::load_from_svg(const char* filename, double _lz)
{
	lz = _lz;
	std::vector<point> bb_nodes;
	std::vector<point> geo_nodes;
  int read_res = svg::import(filename, bb_nodes, 12, 0);
  if(read_res==-1)
  {
    std::cerr << "no bondary box in input file " << filename << std::endl;
		return false;
	}
  else if(read_res==0)
    std::cout << "bounadary box is read" << std::endl;
  else if(read_res==1)
  {
    std::cerr << "path instead of contour is read. giving up" << std::endl;
		return false;
	}
  read_res = svg::import(filename, geo_nodes, 12, 1);
  if(read_res==-1)
  {
    std::cerr << "no geometry in input file " << filename << std::endl;
		return false;
	}
  else if(read_res==0)
    std::cout << "geometry is read" << std::endl;
  else if(read_res==1)
  {
    std::cerr << "path instead of contour is read. giving up" << std::endl;
		return false;
	}
	if(bb_nodes.size() != 4)
	{
    std::cerr << "failed to read bounding box from input svg file. did not you forgot to convert bb rectangle to path object?" << std::endl;
		return false;
	}
	double bbzmin =  std::numeric_limits<double>::max();
	double bbzmax = -std::numeric_limits<double>::max();
	double bbxmin =  std::numeric_limits<double>::max();
	for(const auto& n : bb_nodes)
	{
		bbzmin = std::min(bbzmin, (double) n.x);
		bbzmax = std::max(bbzmax, (double) n.x);
		bbxmin = std::min(bbxmin, (double) n.y);
		std::cout << "\tbb:  " << n.x << ", " << n.y << std::endl;
	}
	double real_lz = bbzmax - bbzmin;
	double f = lz / real_lz;
	geo = new POLY;
	zmin =  std::numeric_limits<double>::max();
	zmax = -std::numeric_limits<double>::max();
	for(const auto& n : geo_nodes)
	{
		std::cout << "\tgeo: " << n.x << ", " << n.y << std::endl;
		double z = ((double) n.x - bbzmin) * f;
		double x = ((double) n.y - bbxmin) * f;
		zmin = std::min(zmin, z);
		zmax = std::max(zmax, z);
		geo->points.push_back(newXY((double)z, (double)x));
	}
	geo->fill_edges();
	geo->fill_tree();
	std::cout << "loaded geometry with zmin = " << zmin
						<< ", zmax = " << zmax
						<< ", number of nodes is " << geo->points.size()
						<< std::endl;
	return true;
}

void cell::scalex(double r)
{
	for(XY* ppt : geo->points)
		ppt->x *= r;
	zmin *= r;
	zmax *= r;
	lz *= r;
}

void cell::scaley(double r)
{
	for(XY* ppt : geo->points)
		ppt->y *= r;
}
