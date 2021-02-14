// this is to provide functionality for creating plasma geometry

// cell is an object which has bounding box and plasma geometry
// bounding box may differ from plasma geometry
// cell has inside(double z, double x) method which is enough to
// deside how to populate plasma in 2d space
// cells create periodical grid along z. each bounding box located after previous
// sharing vertical (along x) boundary
// inside methiod accounts for periodical conditions.
// cell has scale member
// cell has a method to load contours from svg
// cell should have methods to construct primitive geometries
#include "defaults.h"
#include <vector>

struct POLY;

extern "C"
struct cell
{
	cell() : geo(nullptr), lz(1), zmin(0), zmax(0), has_parent(false) {}
	cell(const cell& c) : geo(c.geo), lz(1), zmin(0), zmax(0), has_parent(true) {}
	~cell();
	inline FL_DBL get_lz() const
	{
		return lz;
	}
	inline FL_DBL get_zmin() const
	{
		return zmin;
	}
	inline FL_DBL get_zmax() const
	{
		return zmax;
	}
	void scalex(double /*ratio*/);
	void scaley(double /*ratio*/);
	bool inside(FL_DBL, FL_DBL) const;
	bool load_from_svg(const char*, FL_DBL=1/*lz*/);

private:
	POLY* geo;
	FL_DBL lz;
	FL_DBL zmin;
	FL_DBL zmax;
	bool has_parent;
};
